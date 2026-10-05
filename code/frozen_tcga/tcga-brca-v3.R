#!/usr/bin/env Rscript
# V3: valid adaptive weights + FDR-targeted lambda.
#
# Fixes two root causes found overnight:
#  (1) SKAT-based adaptive weights are INVALID for continuous expression data.
#      SKATBinary expects genotype counts {0,1,2}; on standardized expression it
#      either returns p=1 via SNP-QC ("ALL SNPs have either high missing rates
#      or no-variation") or crashes ("object 'K' not found"), silently falling
#      back to p=1. Result: strongest biological groups got the LARGEST penalty
#      (ESR1's group: 4.91 vs median 0.22) and true drivers were never selected.
#      -> replaced by RMS of univariate logistic Wald z-statistics per group:
#         w_g = sqrt(|g|) / (1 + RMS_z_g), normalized to mean 1. Valid for
#         continuous X / binary y, preserves the thesis's intent.
#  (2) The max(originals - controls) lambda rule is too liberal on real data
#      (implied FDR 0.65 at chosen lambda in a smoke test; 302/1500 selected).
#      -> knockoff+-style per-permutation lambda: max discoveries subject to
#         (1 + #controls)/max(1,#originals) <= 0.10.
#
# Compares: V3-adaptive (new valid weights) vs V3-standard (thesis correlation-
# rank weights), both with the FDR-targeted lambda.
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(Matrix)})
set.seed(20260909)
FDR_TARGET <- 0.10; N_PERM <- 100
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"
dir.create(OUT, showWarnings = FALSE)
logf <- file.path(OUT, "run_log_v3.txt")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

PAM50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
 "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","FOXC1",
 "GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67",
 "MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1",
 "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")
ONCOTYPE_DX <- c("MKI67","AURKA","BIRC5","CCNB1","MYBL2","ERBB2","GRB7","ESR1","PGR",
 "BCL2","SCUBE2","MMP11","CTSV","GSTM1","CD68","BAG1","ACTB","GAPDH","RPLP0","GUSB","TFRC")
MAMMAPRINT_55 <- c("ADGRG6","AKAP2","ALDH4A1","AP2B1","BBC3","CCN4","CCNE2","CDC42BPA",
 "CDCA7","CENPA","CMC2","COL4A2","DCK","DIAPH3","DTL","ECI2","ECT2","ESM1","EXT1",
 "FGF18","FLT1","GMPS","GNAZ","GPR180","GSTM3","IGFBP5","LPCAT1","MCM6","MELK","MMP9",
 "MS4A7","MSANTD3","MTDH","NDC80","NMU","NUSAP1","ORC6","OXCT1","PITRM1","PLAAT1",
 "PRC1","QSOX2","RAB6B","RFC4","RTN4RL1","RUNDC1","SCUBE2","SERF1A","SLC2A3","STK32B",
 "TGFB3","TMEM74B","TSPYL5","UCHL5","ZNF385B")
SIGS <- list(PAM50 = PAM50, OncotypeDX = ONCOTYPE_DX, MammaPrint = MAMMAPRINT_55)
POS_CTRL <- c("ESR1", "PGR", "FOXA1", "GATA3")

define_groups <- function(X) {
  cm <- cor(X, method = "spearman"); dm <- as.dist(1 - abs(cm))
  hc <- hclust(dm, method = "ward.D2"); k <- max(1, round(ncol(X) / 8))
  list(groups = cutree(hc, k = k))
}
# Valid adaptive weights: RMS of univariate logistic Wald z per group.
valid_adaptive_pen <- function(X, y, groups) {
  X <- as.matrix(X)
  z <- apply(X, 2, function(xj) {
    f <- suppressWarnings(glm(y ~ xj, family = binomial()))
    s <- suppressWarnings(summary(f)$coefficients)
    if (nrow(s) < 2 || is.na(s[2, 3])) return(0)
    unname(s[2, 3])
  })
  ug <- sort(unique(groups))
  rms <- sapply(ug, function(g) {
    zj <- z[names(groups)[groups == g]]
    sqrt(mean(zj^2, na.rm = TRUE))
  })
  gs <- sapply(ug, function(g) sum(groups == g))
  rp <- sqrt(gs) / (1 + rms); rp <- rp / mean(rp)
  setNames(rp, ug)[as.character(ug)]
}
# Thesis "standard" weights: correlation-rank x sqrt(size).
standard_pen <- function(X, groups) {
  cm <- cor(as.matrix(X), method = "spearman")
  ug <- sort(unique(groups))
  mc <- sapply(ug, function(g) {
    m <- names(groups)[groups == g]
    if (length(m) < 2) return(0)
    mean(abs(cm[m, m][upper.tri(cm[m, m])]))
  })
  gs <- sapply(ug, function(g) sum(groups == g))
  rp <- rank(-mc, ties.method = "min") * sqrt(gs); rp <- rp / mean(rp)
  setNames(rp, ug)[as.character(ug)]
}
perm_select_v3 <- function(X, y, n_perm, fgroups, gpen, q, cutoff, tag) {
  X <- as.matrix(X); p <- ncol(X); fn <- colnames(X)
  kg <- fgroups + max(fgroups); fp <- c(gpen, gpen)
  hits <- c(); no_med <- c(); nc_med <- c(); nperm_used <- 0
  for (k in seq_len(n_perm)) {
    Xk <- X[sample(nrow(X)), ]
    fit <- suppressWarnings(grpreg(cbind(X, Xk), y, group = c(fgroups, kg),
      penalty = "grLasso", family = "binomial", group.multiplier = fp,
      nlambda = 50, lambda.min = 0.05))
    b <- as.matrix(fit$beta[-1, , drop = FALSE])
    n_o <- colSums(abs(b[1:p, , drop = FALSE]) > 0)
    n_c <- colSums(abs(b[(p + 1):(2 * p), , drop = FALSE]) > 0)
    fdr_hat <- (1 + n_c) / pmax(1, n_o)
    ok <- which(fdr_hat <= q & n_o > 0)
    if (!length(ok)) next
    nperm_used <- nperm_used + 1
    j <- ok[which.max(n_o[ok])]
    no_med <- c(no_med, n_o[j]); nc_med <- c(nc_med, n_c[j])
    si <- which(abs(b[1:p, j]) > 0)
    if (length(si)) hits <- c(hits, rownames(b)[si])
  }
  fc <- setNames(rep(0, p), fn); ct <- table(hits); fc[names(ct)] <- ct
  fr <- fc / n_perm
  logmsg("[%s] perms voting: %d/%d; median at chosen lambda: n_orig=%.0f n_ctrl=%.1f",
    tag, nperm_used, n_perm, median(no_med), median(nc_med))
  list(sel = names(fr)[fr >= cutoff], freq = fr)
}
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }

## ---- data (identical to v1) ----
logmsg("loading RNA...")
rna <- suppressMessages(read_tsv(file.path(DATA, "data_mrna_seq_v2_rsem.txt"),
  comment = "#", show_col_types = FALSE))
genes <- rna$Hugo_Symbol
rna <- as.matrix(rna[, -(1:2)]); rownames(rna) <- genes
clin_s <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_sample.txt"),
  comment = "#", show_col_types = FALSE))
clin_p <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_patient.txt"),
  comment = "#", show_col_types = FALSE))
is_primary <- clin_s$SAMPLE_TYPE == "Primary"
rna <- rna[, colnames(rna) %in% clin_s$SAMPLE_ID[is_primary]]
pat <- substr(colnames(rna), 1, 12)
ord <- order(pat, !(substr(colnames(rna), 14, 16) == "01A"))
keep <- !duplicated(pat[ord]); rna <- rna[, ord[keep]]; pat <- pat[ord[keep]]
er_p <- setNames(clin_p[[grep("ER_STATUS_BY_IHC", names(clin_p), value = TRUE)[1]]],
  clin_p$PATIENT_ID)
er <- toupper(trimws(er_p[pat]))
y_er <- ifelse(er == "POSITIVE", 1L, ifelse(er == "NEGATIVE", 0L, NA_integer_))
keep_s <- !is.na(y_er)
X_all <- log2(t(rna[, keep_s]) + 1); y_all <- y_er[keep_s]
X_all <- X_all[, colSums(is.na(X_all)) == 0]
v0 <- apply(X_all, 2, var)
keep_g <- tapply(seq_len(ncol(X_all)), colnames(X_all), function(ii) ii[which.max(v0[ii])])
X_all <- X_all[, unlist(keep_g)]
X_all <- X_all[, apply(X_all, 2, var) > 1e-6]
v <- apply(X_all, 2, var)
top_genes <- names(sort(v, decreasing = TRUE))[1:min(1500, length(v))]
X <- t(scale(t(X_all[, top_genes]))); X[is.na(X)] <- 0
universe <- colnames(X)
logmsg("analysis matrix: %d x %d; ER+ rate %.2f", nrow(X), ncol(X), mean(y_all))

y <- y_all; n <- length(y)
te <- unlist(lapply(split(seq_len(n), y), function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
tr <- setdiff(seq_len(n), te)
Xtr <- X[tr, ]; ytr <- y[tr]; Xte <- X[te, ]; yte <- y[te]

grp <- define_groups(Xtr)
pen_ad <- valid_adaptive_pen(Xtr, ytr, grp$groups)
pen_st <- standard_pen(Xtr, grp$groups)
for (g in POS_CTRL) {
  gid <- as.character(grp$groups[g])
  logmsg("control %s group %s: adaptive_w=%.3f standard_w=%.3f (medians %.3f / %.3f)",
    g, gid, pen_ad[gid], pen_st[gid], median(pen_ad), median(pen_st))
}

results <- list()
for (cfg in list(list(tag = "adaptive", pen = pen_ad), list(tag = "standard", pen = pen_st))) {
  t0 <- Sys.time()
  sel <- perm_select_v3(Xtr, ytr, N_PERM, grp$groups, cfg$pen, FDR_TARGET, 0.5, cfg$tag)
  auc_te <- NA_real_
  if (length(sel$sel) >= 2) {
    fit <- cv.glmnet(Xtr[, sel$sel, drop = FALSE], ytr, family = "binomial")
    pr <- as.numeric(predict(fit, Xte[, sel$sel, drop = FALSE], s = "lambda.min", type = "response"))
    auc_te <- as.numeric(auc(roc(yte, pr, quiet = TRUE)))
  }
  fr_ord <- sort(sel$freq, decreasing = TRUE)
  pc_rank <- sapply(POS_CTRL, function(g)
    if (g %in% names(fr_ord)) which(names(fr_ord) == g) else NA_integer_)
  ov <- bind_rows(lapply(names(SIGS), function(sg_) {
    sg <- intersect(SIGS[[sg_]], universe); sl <- intersect(sel$sel, universe)
    qq <- length(intersect(sl, sg))
    data.frame(signature = sg_, overlap = qq, sig_in_universe = length(sg),
      selected = length(sl), enrich = round((qq / max(1, length(sl))) /
        (length(sg) / length(universe)), 2),
      p = signif(phyper(qq - 1, length(sg), length(universe) - length(sg),
        length(sl), lower.tail = FALSE), 3))
  }))
  stab <- lapply(1:3, function(b) {
    bi <- unlist(lapply(split(seq_along(ytr), ytr),
      function(ii) sample(ii, floor(0.8 * length(ii)))))
    perm_select_v3(Xtr[bi, ], ytr[bi], 50, grp$groups, cfg$pen, FDR_TARGET, 0.5,
      paste0(cfg$tag, "-stab", b))$sel })
  pw <- combn(3, 2, function(ii) jaccard(stab[[ii[1]]], stab[[ii[2]]]))
  logmsg("[%s] selected %d genes (%.1f min); test AUC %.3f; Jaccard %.3f",
    cfg$tag, length(sel$sel), as.numeric(difftime(Sys.time(), t0, units = "mins")),
    auc_te, mean(pw))
  logmsg("[%s] control ranks: %s", cfg$tag,
    paste(sprintf("%s:%s", POS_CTRL, pc_rank), collapse = " "))
  print(as.data.frame(ov))
  write.csv(data.frame(gene = sel$sel, freq = round(as.numeric(sel$freq[sel$sel]), 3)),
    file.path(OUT, sprintf("selected_ER_IHC_v3_%s.csv", cfg$tag)), row.names = FALSE)
  results[[cfg$tag]] <- list(selected = sel$sel, freq = sel$freq, auc = auc_te,
    overlap = ov, jaccard = mean(pw), pc_rank = pc_rank)
}
saveRDS(results, file.path(OUT, "er_ihc_v3_results.rds"))
logmsg("V3 DONE")

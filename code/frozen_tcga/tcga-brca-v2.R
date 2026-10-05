#!/usr/bin/env Rscript
# V2: FDR-targeted lambda calibration for the permutation-assisted group lasso.
# Same data/preprocessing/groups/SKAT penalties as tcga-brca-real.R, but the
# per-permutation lambda is chosen knockoff+-style: among lambdas, take the one
# with the most original discoveries subject to
#     (1 + #controls selected) / max(1, #originals selected) <= FDR_TARGET
# instead of maximizing the raw (originals - controls) difference.
# Hypothesis: the raw-difference rule is too liberal on real data with strong
# correlation structure, selecting far too many features (302/1500 in v1).
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(SKAT); library(Matrix)})
set.seed(20260909)
FDR_TARGET <- 0.10
N_PERM <- 100
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"
dir.create(OUT, showWarnings = FALSE)
logf <- file.path(OUT, "run_log_v2.txt")
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
skat_pen <- function(X, y, groups) {
  obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
  ug <- sort(unique(groups)); Xm <- as.matrix(X)
  pv <- setNames(numeric(length(ug)), ug); gs <- setNames(numeric(length(ug)), ug)
  for (g in ug) {
    Z <- Xm[, names(groups)[groups == g], drop = FALSE]; gs[as.character(g)] <- ncol(Z)
    pv[as.character(g)] <- tryCatch(suppressWarnings(
      SKATBinary(Z, obj, kernel = "Linear.Weighted")$p.value), error = function(e) 1)
  }
  nlp <- -log(pv + 1e-6); nlp[nlp < 0] <- 0
  rp <- sqrt(gs) * (1 / (nlp + 1)); rp <- rp / mean(rp)
  rp[as.character(sort(unique(groups)))]
}

# One permutation: fit the grpreg path, evaluate every lambda under both the
# v1 rule (max n_orig - n_ctrl) and the v2 rule (max n_orig s.t. FDR_hat <= q).
# Returns the selected original features under the v2 rule (NULL if none qualifies).
perm_one <- function(X, y, fgroups, gpen, q) {
  X <- as.matrix(X); p <- ncol(X); fn <- colnames(X)
  kg <- fgroups + max(fgroups); fp <- c(gpen, gpen)
  Xk <- X[sample(nrow(X)), ]
  fit <- suppressWarnings(grpreg(cbind(X, Xk), y, group = c(fgroups, kg),
    penalty = "grLasso", family = "binomial", group.multiplier = fp,
    nlambda = 50, lambda.min = 0.05))
  b <- as.matrix(fit$beta[-1, , drop = FALSE])
  n_o <- colSums(abs(b[1:p, , drop = FALSE]) > 0)
  n_c <- colSums(abs(b[(p + 1):(2 * p), , drop = FALSE]) > 0)
  fdr_hat <- (1 + n_c) / pmax(1, n_o)
  ok <- which(fdr_hat <= q & n_o > 0)
  v1_pick <- which.max(n_o - n_c)
  if (!length(ok)) return(list(sel = NULL,
    v1 = list(n_o = n_o[v1_pick], n_c = n_c[v1_pick], fdr = fdr_hat[v1_pick])))
  j <- ok[which.max(n_o[ok])]
  si <- which(abs(b[1:p, j]) > 0)
  list(sel = rownames(b)[si],
    v1 = list(n_o = n_o[v1_pick], n_c = n_c[v1_pick], fdr = fdr_hat[v1_pick]),
    v2 = list(n_o = n_o[j], n_c = n_c[j], fdr = fdr_hat[j]))
}

perm_select_v2 <- function(X, y, n_perm, fgroups, gpen, q, cutoff) {
  p <- ncol(X); fn <- colnames(X); hits <- c(); v1f <- c(); v2f <- c()
  for (k in seq_len(n_perm)) {
    r <- perm_one(X, y, fgroups, gpen, q)
    v1f <- c(v1f, r$v1$fdr); v2f <- c(v2f, r$v2$fdr)
    if (!is.null(r$sel)) hits <- c(hits, r$sel)
  }
  fc <- setNames(rep(0, p), fn); ct <- table(hits); fc[names(ct)] <- ct
  fr <- fc / n_perm
  list(sel = names(fr)[fr >= cutoff], freq = fr,
    v1_fdr_med = median(v1f, na.rm = TRUE), v2_fdr_med = median(v2f, na.rm = TRUE))
}
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }

## ---- data (identical preprocessing to v1) ----
logmsg("loading RNA...")
rna <- suppressMessages(read_tsv(file.path(DATA, "data_mrna_seq_v2_rsem.txt"),
  comment = "#", show_col_types = FALSE))
genes <- rna$Hugo_Symbol
rna <- as.matrix(rna[, -(1:2)]); rownames(rna) <- genes
logmsg("RNA: %d genes x %d samples", nrow(rna), ncol(rna))
clin_s <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_sample.txt"),
  comment = "#", show_col_types = FALSE))
clin_p <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_patient.txt"),
  comment = "#", show_col_types = FALSE))
is_primary <- clin_s$SAMPLE_TYPE == "Primary"
samp <- clin_s$SAMPLE_ID[is_primary]
rna <- rna[, colnames(rna) %in% samp]
pat <- substr(colnames(rna), 1, 12)
ord <- order(pat, !(substr(colnames(rna), 14, 16) == "01A"))
keep <- !duplicated(pat[ord]); rna <- rna[, ord[keep]]; pat <- pat[ord[keep]]
er_col_p <- grep("ER_STATUS_BY_IHC", names(clin_p), value = TRUE)[1]
er_p <- setNames(clin_p[[er_col_p]], clin_p$PATIENT_ID)
er <- er_p[pat]; er <- toupper(trimws(er))
y_er <- ifelse(er == "POSITIVE", 1L, ifelse(er == "NEGATIVE", 0L, NA_integer_))
logmsg("primary-tumor samples, one per patient: %d", length(pat))
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
logmsg("analysis matrix: %d samples x %d genes; ER+ rate %.2f", nrow(X), ncol(X), mean(y_all))
logmsg("positive controls in universe: %s",
  paste(sprintf("%s=%s", POS_CTRL, POS_CTRL %in% universe), collapse = " "))

y <- y_all
n <- length(y)
te <- unlist(lapply(split(seq_len(n), y), function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
tr <- setdiff(seq_len(n), te)
Xtr <- X[tr, ]; ytr <- y[tr]; Xte <- X[te, ]; yte <- y[te]
logmsg("train %d (pos %d) | test %d (pos %d)", length(tr), sum(ytr), length(te), sum(yte))

t0 <- Sys.time()
grp <- define_groups(Xtr); pen <- skat_pen(Xtr, ytr, grp$groups)
sel <- perm_select_v2(Xtr, ytr, N_PERM, grp$groups, pen, FDR_TARGET, 0.5)
logmsg("V2 selected %d genes in %.1f min (FDR target %.2f)",
  length(sel$sel), as.numeric(difftime(Sys.time(), t0, units = "mins")), FDR_TARGET)
logmsg("median implied FDR at chosen lambda: v1-rule %.3f | v2-rule %.3f",
  sel$v1_fdr_med, sel$v2_fdr_med)

auc_te <- NA_real_
if (length(sel$sel) >= 2) {
  fit <- cv.glmnet(Xtr[, sel$sel, drop = FALSE], ytr, family = "binomial")
  pr <- as.numeric(predict(fit, Xte[, sel$sel, drop = FALSE], s = "lambda.min", type = "response"))
  auc_te <- as.numeric(auc(roc(yte, pr, quiet = TRUE)))
}
logmsg("held-out test AUC: %.3f", auc_te)

# positive-control ranks
fr_ord <- sort(sel$freq, decreasing = TRUE)
for (g in POS_CTRL) {
  r <- if (g %in% names(fr_ord)) which(names(fr_ord) == g) else NA_integer_
  logmsg("control %s: rank %s / %d, freq %.2f, selected %s", g,
    ifelse(is.na(r), "NA", r), length(fr_ord),
    ifelse(is.na(r), 0, unname(fr_ord[g])), g %in% sel$sel)
}

ov <- bind_rows(lapply(names(SIGS), function(sg_) {
  sg <- intersect(SIGS[[sg_]], universe); sl <- intersect(sel$sel, universe)
  q <- length(intersect(sl, sg))
  data.frame(signature = sg_, overlap = q, sig_in_universe = length(sg),
    selected = length(sl), universe = length(universe),
    enrich = round((q / max(1, length(sl))) / (length(sg) / length(universe)), 2),
    p = signif(phyper(q - 1, length(sg), length(universe) - length(sg),
      length(sl), lower.tail = FALSE), 3))
}))
print(as.data.frame(ov))

# light stability: 3 x 50 perms on 80% subsamples
stab <- lapply(1:3, function(b) {
  bi <- unlist(lapply(split(seq_along(ytr), ytr),
    function(ii) sample(ii, floor(0.8 * length(ii)))))
  perm_select_v2(Xtr[bi, ], ytr[bi], 50, grp$groups, pen, FDR_TARGET, 0.5)$sel })
pw <- combn(3, 2, function(ii) jaccard(stab[[ii[1]]], stab[[ii[2]]]))
logmsg("mean pairwise Jaccard (3x50): %.3f", mean(pw))

# re-threshold curve on the main frequencies (no extra compute)
for (cf in c(0.5, 0.6, 0.7, 0.8, 0.9)) {
  s2 <- names(sel$freq)[sel$freq >= cf]
  pc <- sum(POS_CTRL %in% s2)
  logmsg("cutoff %.1f -> %d genes, positive controls kept %d/%d", cf, length(s2), pc, length(POS_CTRL))
}

write.csv(data.frame(gene = sel$sel, freq = round(as.numeric(sel$freq[sel$sel]), 3)),
  file.path(OUT, "selected_ER_IHC_v2.csv"), row.names = FALSE)
write.csv(ov, file.path(OUT, "overlap_ER_IHC_v2.csv"), row.names = FALSE)
saveRDS(list(selected = sel$sel, freq = sel$freq, auc = auc_te, overlap = ov,
  jaccard = mean(pw), v1_fdr_med = sel$v1_fdr_med, v2_fdr_med = sel$v2_fdr_med,
  config = list(fdr_target = FDR_TARGET, n_perm = N_PERM)),
  file.path(OUT, "er_ihc_v2_results.rds"))
logmsg("V2 DONE")

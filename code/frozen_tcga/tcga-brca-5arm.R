#!/usr/bin/env Rscript
# 5-arm TCGA-BRCA comparison mirroring the thesis simulation design:
#   Baseline LASSO | Early Fusion Std/Adap | Late Fusion Std/Adap
# Modalities: RNA-seq (top 1500 variable genes) + HM450 methylation
#             (top 1500 variable genes, gene-level beta values).
# Endpoint: ER status by IHC. One fixed stratified 70/30 train/test split;
# selection on TRAIN only; held-out AUC on TEST.
# Adaptive weights = V3-validated RMS univariate Wald-z weights (SKAT is invalid
# for continuous data). Lambda via knockoff-inspired FDR-targeted rule (q=0.10).
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
RUN_ID <- if (length(args) >= 1) as.integer(args[1]) else 1L
set.seed(20260909 + RUN_ID)
DO_STAB <- RUN_ID == 1L  # within-run subsample stability only on run 1; across-run Jaccard is the thesis metric
FDR_TARGET <- 0.10; N_PERM <- 80; N_STAB_REP <- 3; N_STAB_PERM <- 30
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"
dir.create(OUT, showWarnings = FALSE)
logf <- file.path(OUT, sprintf("run_log_5arm_run%d.txt", RUN_ID))
sfx <- sprintf("_run%d", RUN_ID)
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
  cutree(hc, k = k)
}
valid_adaptive_pen <- function(X, y, groups) {
  X <- as.matrix(X)
  z <- apply(X, 2, function(xj) {
    f <- suppressWarnings(glm(y ~ xj, family = binomial()))
    s <- suppressWarnings(summary(f)$coefficients)
    if (nrow(s) < 2 || is.na(s[2, 3])) return(0)
    unname(s[2, 3])
  })
  ug <- sort(unique(groups))
  rms <- sapply(ug, function(g) { zj <- z[names(groups)[groups == g]]; sqrt(mean(zj^2, na.rm = TRUE)) })
  gs <- sapply(ug, function(g) sum(groups == g))
  rp <- sqrt(gs) / (1 + rms); rp <- rp / mean(rp)
  setNames(rp, ug)[as.character(ug)]
}
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
# Permutation-assisted group-lasso selection with FDR-targeted lambda.
perm_select <- function(X, y, n_perm, groups, gpen, q, cutoff, tag) {
  X <- as.matrix(X); p <- ncol(X); fn <- colnames(X)
  kg <- groups + max(groups); fp <- c(gpen, gpen)
  hits <- c(); nperm_used <- 0
  for (k in seq_len(n_perm)) {
    Xk <- X[sample(nrow(X)), ]
    fit <- suppressWarnings(grpreg(cbind(X, Xk), y, group = c(groups, kg),
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
    si <- which(abs(b[1:p, j]) > 0)
    if (length(si)) hits <- c(hits, rownames(b)[si])
  }
  fc <- setNames(rep(0, p), fn); ct <- table(hits); fc[names(ct)] <- ct
  fr <- fc / n_perm
  logmsg("[%s] perms voting: %d/%d", tag, nperm_used, n_perm)
  list(sel = names(fr)[fr >= cutoff], freq = fr)
}
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }
strip_prefix <- function(x) sub("^METH_", "", x)

enrichment <- function(sel, universe) {
  sl <- unique(strip_prefix(intersect(sel, universe)))
  bind_rows(lapply(names(SIGS), function(sg_) {
    sg <- intersect(SIGS[[sg_]], universe)
    qq <- length(intersect(sl, sg))
    data.frame(signature = sg_, overlap = qq, sig_in_universe = length(sg),
      selected = length(sl),
      enrich = round((qq / max(1, length(sl))) / (length(sg) / length(universe)), 2),
      p = signif(phyper(qq - 1, length(sg), length(universe) - length(sg),
        length(sl), lower.tail = FALSE), 3))
  }))
}

## ---------------- data ----------------
logmsg("loading RNA...")
rna <- suppressMessages(read_tsv(file.path(DATA, "data_mrna_seq_v2_rsem.txt"),
  comment = "#", show_col_types = FALSE))
genes <- rna$Hugo_Symbol
rna <- as.matrix(rna[, -(1:2)]); rownames(rna) <- genes
clin_s <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_sample.txt"),
  comment = "#", show_col_types = FALSE))
clin_p <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_patient.txt"),
  comment = "#", show_col_types = FALSE))
rna <- rna[, colnames(rna) %in% clin_s$SAMPLE_ID[clin_s$SAMPLE_TYPE == "Primary"]]
pat <- substr(colnames(rna), 1, 12)
ord <- order(pat, !(substr(colnames(rna), 14, 16) == "01A"))
keep <- !duplicated(pat[ord]); rna <- rna[, ord[keep]]; pat <- pat[ord[keep]]
er_p <- setNames(clin_p[[grep("ER_STATUS_BY_IHC", names(clin_p), value = TRUE)[1]]],
  clin_p$PATIENT_ID)
er <- toupper(trimws(er_p[pat]))
y_er <- ifelse(er == "POSITIVE", 1L, ifelse(er == "NEGATIVE", 0L, NA_integer_))
keep_s <- !is.na(y_er)
rna <- rna[, keep_s]; pat_rna <- pat[keep_s]; y_rna <- y_er[keep_s]
logmsg("RNA: %d primary samples with ER labels", length(pat_rna))

logmsg("loading methylation...")
m <- suppressMessages(read_tsv(file.path(DATA, "data_methylation_hm450.txt"),
  show_col_types = FALSE))
mgenes <- m$Hugo_Symbol
mmat <- as.matrix(m[, -(1:2)]); rownames(mmat) <- mgenes
samp <- colnames(mmat)
is_primary_m <- substr(samp, 14, 16) == "01"
mmat <- mmat[, is_primary_m]; samp <- samp[is_primary_m]
pat_m <- substr(samp, 1, 12)
ord_m <- order(pat_m); mmat <- mmat[, ord_m]; pat_m <- pat_m[ord_m]
dup_m <- duplicated(pat_m); mmat <- mmat[, !dup_m]; pat_m <- pat_m[!dup_m]
logmsg("METH: %d genes x %d primary samples", nrow(mmat), ncol(mmat))

common <- intersect(pat_rna, pat_m)
logmsg("overlap: %d patients with RNA + METH + ER label", length(common))
ri <- match(common, pat_rna); mi <- match(common, pat_m)
y <- y_rna[ri]
Xr_all <- log2(t(rna[, ri]) + 1)
Xm_all <- t(mmat[, mi])
logmsg("ER+ rate in overlap set: %.3f", mean(y))

# per-modality QC: drop all-NA genes, median-impute, dedup symbols, variance filter
prep_mod <- function(X, top_n, modname) {
  na_frac <- colMeans(is.na(X))
  logmsg("[%s] genes with >10%% missing: %d", modname, sum(na_frac > 0.10))
  X <- X[, na_frac <= 0.10]
  for (j in seq_len(ncol(X))) { xj <- X[, j]; xj[is.na(xj)] <- median(xj, na.rm = TRUE); X[, j] <- xj }
  v0 <- apply(X, 2, var)
  keep_g <- tapply(seq_len(ncol(X)), colnames(X), function(ii) ii[which.max(v0[ii])])
  X <- X[, unlist(keep_g)]
  X <- X[, apply(X, 2, var) > 1e-6]
  v <- apply(X, 2, var)
  tg <- names(sort(v, decreasing = TRUE))[1:min(top_n, length(v))]
  Xs <- t(scale(t(X[, tg]))); Xs[is.na(Xs)] <- 0
  logmsg("[%s] final: %d x %d", modname, nrow(Xs), ncol(Xs))
  Xs
}
Xr <- prep_mod(Xr_all, 1500, "RNA")
Xm <- prep_mod(Xm_all, 1500, "METH")
colnames(Xm) <- paste0("METH_", colnames(Xm))
X_ef <- cbind(Xr, Xm)
universe <- unique(strip_prefix(colnames(X_ef)))
logmsg("universe genes: %d; total features: %d", length(universe), ncol(X_ef))

## ---------------- train/test split ----------------
n <- length(y)
te <- unlist(lapply(split(seq_len(n), y), function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
tr <- setdiff(seq_len(n), te)
logmsg("train %d (ER+ %.2f) | test %d (ER+ %.2f)", length(tr), mean(y[tr]), length(te), mean(y[te]))

refit_auc <- function(Xtr, ytr, Xte, yte, sel) {
  if (length(sel) < 2) return(list(auc = NA_real_, probs = NULL))
  fit <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial")
  pr <- as.numeric(predict(fit, Xte[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  list(auc = as.numeric(auc(roc(yte, pr, quiet = TRUE))), probs = pr)
}
# per-modality train CV AUC (for late-fusion ensemble weights)
mod_cv_auc <- function(Xtr, ytr, sel) {
  if (length(sel) < 2) return(NA_real_)
  fit <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial", nfolds = 5)
  pr <- as.numeric(predict(fit, Xtr[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  as.numeric(auc(roc(ytr, pr, quiet = TRUE)))
}

stab_jaccard <- function(Xtr, ytr, sel_fun, tag) {
  reps <- lapply(seq_len(N_STAB_REP), function(b) {
    bi <- unlist(lapply(split(seq_along(ytr), ytr),
      function(ii) sample(ii, floor(0.8 * length(ii)))))
    sel_fun(Xtr[bi, , drop = FALSE], ytr[bi], paste0(tag, "-stab", b))
  })
  mean(combn(N_STAB_REP, 2, function(ii) jaccard(reps[[ii[1]]], reps[[ii[2]]])))
}

results <- list()
row <- list()

## ---------------- Arm 1: Baseline LASSO (concatenated) ----------------
t0 <- Sys.time()
Xtr_ef <- X_ef[tr, ]; ytr <- y[tr]; Xte_ef <- X_ef[te, ]; yte <- y[te]
fit_l <- cv.glmnet(Xtr_ef, ytr, family = "binomial")
sel_l <- colnames(Xtr_ef)[which(as.numeric(coef(fit_l, s = "lambda.min"))[-1] != 0)]
pr_l <- as.numeric(predict(fit_l, Xte_ef, s = "lambda.min", type = "response"))
auc_l <- as.numeric(auc(roc(yte, pr_l, quiet = TRUE)))
jac_l <- if (DO_STAB) stab_jaccard(Xtr_ef, ytr, function(Xb, yb, tg) {
  fb <- cv.glmnet(Xb, yb, family = "binomial")
  colnames(Xb)[which(as.numeric(coef(fb, s = "lambda.min"))[-1] != 0)]
}, "lasso") else NA_real_
enr_l <- enrichment(sel_l, universe)
ctrl_l <- sapply(POS_CTRL, function(g) g %in% strip_prefix(sel_l))
logmsg("[LASSO] sel=%d auc=%.3f jacc=%.3f (%.1f min)", length(sel_l), auc_l, jac_l,
  as.numeric(difftime(Sys.time(), t0, units = "mins")))
results$lasso <- list(selected = sel_l, auc = auc_l, jaccard = jac_l, enrich = enr_l)
row$lasso <- c(arm = "Baseline LASSO", n_sel = length(sel_l), auc = auc_l, jaccard = jac_l)

## ---------------- Arms 2-3: Early fusion (Std / Adap) ----------------
grp_ef <- define_groups(Xtr_ef)
pen_ef_ad <- valid_adaptive_pen(Xtr_ef, ytr, grp_ef)
pen_ef_st <- standard_pen(Xtr_ef, grp_ef)
for (cfg in list(list(tag = "EF-Std", pen = pen_ef_st), list(tag = "EF-Adap", pen = pen_ef_ad))) {
  t0 <- Sys.time()
  sel <- perm_select(Xtr_ef, ytr, N_PERM, grp_ef, cfg$pen, FDR_TARGET, 0.5, cfg$tag)
  rf <- refit_auc(Xtr_ef, ytr, Xte_ef, yte, sel$sel)
  jac <- if (DO_STAB) stab_jaccard(Xtr_ef, ytr,
    function(Xb, yb, tg) perm_select(Xb, yb, N_STAB_PERM, define_groups(Xb),
      if (cfg$tag == "EF-Adap") valid_adaptive_pen(Xb, yb, define_groups(Xb)) else standard_pen(Xb, define_groups(Xb)),
      FDR_TARGET, 0.5, tg)$sel, cfg$tag) else NA_real_
  enr <- enrichment(sel$sel, universe)
  ctrl <- sapply(POS_CTRL, function(g) g %in% strip_prefix(sel$sel))
  logmsg("[%s] sel=%d auc=%.3f jacc=%.3f (%.1f min)", cfg$tag, length(sel$sel), rf$auc, jac,
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
  logmsg("[%s] controls: %s", cfg$tag, paste(sprintf("%s:%d", POS_CTRL, as.integer(ctrl)), collapse = " "))
  write.csv(data.frame(feature = sel$sel, freq = round(as.numeric(sel$freq[sel$sel]), 3)),
    file.path(OUT, sprintf("selected_5arm_%s%s.csv", gsub("-", "_", cfg$tag), sfx)), row.names = FALSE)
  results[[cfg$tag]] <- list(selected = sel$sel, freq = sel$freq, auc = rf$auc,
    jaccard = jac, enrich = enr, controls = ctrl)
  row[[cfg$tag]] <- c(arm = cfg$tag, n_sel = length(sel$sel), auc = rf$auc, jaccard = jac)
  saveRDS(results, file.path(OUT, paste0("five_arm_results", sfx, ".rds")))
}

## ---------------- Arms 4-5: Late fusion (Std / Adap) ----------------
mods <- list(RNA = list(Xtr = Xr[tr, , drop = FALSE], Xte = Xr[te, , drop = FALSE]),
             METH = list(Xtr = Xm[tr, , drop = FALSE], Xte = Xm[te, , drop = FALSE]))
for (cfg in list(list(tag = "LF-Std", adap = FALSE), list(tag = "LF-Adap", adap = TRUE))) {
  t0 <- Sys.time()
  sel_m <- list(); probs_m <- list(); auc_m <- c()
  for (mn in names(mods)) {
    Xtrm <- mods[[mn]]$Xtr
    gm <- define_groups(Xtrm)
    penm <- if (cfg$adap) valid_adaptive_pen(Xtrm, ytr, gm) else standard_pen(Xtrm, gm)
    sm <- perm_select(Xtrm, ytr, N_PERM, gm, penm, FDR_TARGET, 0.5, paste0(cfg$tag, "-", mn))
    sel_m[[mn]] <- sm$sel
    rfm <- refit_auc(Xtrm, ytr, mods[[mn]]$Xte, yte, sm$sel)
    probs_m[[mn]] <- rfm$probs; auc_m[mn] <- mod_cv_auc(Xtrm, ytr, sm$sel)
    logmsg("[%s-%s] sel=%d cv_auc=%.3f", cfg$tag, mn, length(sm$sel), auc_m[mn])
  }
  sel_union <- unique(unlist(sel_m))
  # AUC-weighted ensemble (thesis late-fusion prediction rule)
  w <- pmax(0, auc_m - 0.5); w[is.na(w)] <- 0
  if (sum(w) == 0) {
    auc_ens <- NA_real_
  } else {
    w <- w / sum(w)
    P <- sapply(names(mods), function(mn) {
      pr <- probs_m[[mn]]
      if (is.null(pr)) rep(mean(ytr), length(yte)) else pr
    })
    p_ens <- as.numeric(P %*% w)
    auc_ens <- as.numeric(auc(roc(yte, p_ens, quiet = TRUE)))
  }
  jac <- if (DO_STAB) stab_jaccard(Xtr_ef, ytr, function(Xb, yb, tg) {
    su <- c()
    for (mn in names(mods)) {
      Xbm <- Xb[, colnames(mods[[mn]]$Xtr), drop = FALSE]
      gm <- define_groups(Xbm)
      penm <- if (cfg$adap) valid_adaptive_pen(Xbm, yb, gm) else standard_pen(Xbm, gm)
      su <- c(su, perm_select(Xbm, yb, N_STAB_PERM, gm, penm, FDR_TARGET, 0.5, tg)$sel)
    }
    unique(su)
  }, cfg$tag) else NA_real_
  enr <- enrichment(sel_union, universe)
  ctrl <- sapply(POS_CTRL, function(g) g %in% strip_prefix(sel_union))
  logmsg("[%s] union_sel=%d ens_auc=%.3f jacc=%.3f (%.1f min)", cfg$tag, length(sel_union),
    auc_ens, jac, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  logmsg("[%s] controls: %s", cfg$tag, paste(sprintf("%s:%d", POS_CTRL, as.integer(ctrl)), collapse = " "))
  results[[cfg$tag]] <- list(selected = sel_union, per_mod = sel_m, auc = auc_ens,
    jaccard = jac, enrich = enr, controls = ctrl, mod_auc = auc_m)
  row[[cfg$tag]] <- c(arm = cfg$tag, n_sel = length(sel_union), auc = auc_ens, jaccard = jac)
  saveRDS(results, file.path(OUT, paste0("five_arm_results", sfx, ".rds")))
}

## ---------------- summary table ----------------
tab <- bind_rows(lapply(names(row), function(a) {
  r <- results[[a]]
  er <- r$enrich
  gp <- function(sg) { x <- er[er$signature == sg, ]; c(ov = x$overlap, p = x$p) }
  data.frame(
    arm = if (a == "lasso") "Baseline LASSO" else a,
    n_selected = as.integer(row[[a]]["n_sel"]),
    n_RNA = sum(!grepl("^METH_", r$selected)),
    n_METH = sum(grepl("^METH_", r$selected)),
    heldout_AUC = round(as.numeric(row[[a]]["auc"]), 3),
    stability_Jaccard = round(as.numeric(row[[a]]["jaccard"]), 3),
    PAM50_overlap = gp("PAM50")["ov"], PAM50_p = gp("PAM50")["p"],
    Oncotype_overlap = gp("OncotypeDX")["ov"], Oncotype_p = gp("OncotypeDX")["p"],
    MammaPrint_overlap = gp("MammaPrint")["ov"],
    pos_controls = paste(POS_CTRL[as.logical(r$controls)], collapse = ","))
}))
print(tab)
write.csv(tab, file.path(OUT, paste0("five_arm_table", sfx, ".csv")), row.names = FALSE)
saveRDS(results, file.path(OUT, paste0("five_arm_results", sfx, ".rds")))
logmsg("5-ARM DONE")

#!/usr/bin/env Rscript
# Shared library for the TCGA-BRCA adaptive-model stability study (2026-09-22).
# Data protocol is IDENTICAL to tcga-brca-5arm.R so results stay comparable:
# same QC, variance filter, scaling, universe, train/test scheme (new seeds).
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(Matrix)})

FDR_TARGET <- 0.10
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"

PAM50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
 "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","FOXC1",
 "GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67",
 "MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1",
 "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")
ONCOTYPE_DX <- c("MKI67","AURKA","BIRC5","CCNB1","MYBL2","ERBB2","GRB7","ESR1","PGR",
 "BCL2","SCUBE2","MMP9","MMP11","CTSV","GSTM1","CD68","BAG1","ACTB","GAPDH","RPLP0","GUSB","TFRC")
SIGS <- list(PAM50 = PAM50, OncotypeDX = ONCOTYPE_DX)
POS_CTRL <- c("ESR1", "PGR", "FOXA1", "GATA3")

jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }
strip_prefix <- function(x) sub("^METH_", "", x)

## ---------------- data (verbatim protocol from tcga-brca-5arm.R) ----------------
load_data <- function() {
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

  common <- intersect(pat_rna, pat_m)
  ri <- match(common, pat_rna); mi <- match(common, pat_m)
  y <- y_rna[ri]
  Xr_all <- log2(t(rna[, ri]) + 1)
  Xm_all <- t(mmat[, mi])

  prep_mod <- function(X, top_n) {
    na_frac <- colMeans(is.na(X))
    X <- X[, na_frac <= 0.10]
    for (j in seq_len(ncol(X))) { xj <- X[, j]; xj[is.na(xj)] <- median(xj, na.rm = TRUE); X[, j] <- xj }
    v0 <- apply(X, 2, var)
    keep_g <- tapply(seq_len(ncol(X)), colnames(X), function(ii) ii[which.max(v0[ii])])
    X <- X[, unlist(keep_g)]
    X <- X[, apply(X, 2, var) > 1e-6]
    v <- apply(X, 2, var)
    tg <- names(sort(v, decreasing = TRUE))[1:min(top_n, length(v))]
    Xs <- t(scale(t(X[, tg]))); Xs[is.na(Xs)] <- 0
    Xs
  }
  Xr <- prep_mod(Xr_all, 1500)
  Xm <- prep_mod(Xm_all, 1500)
  colnames(Xm) <- paste0("METH_", colnames(Xm))
  X_ef <- cbind(Xr, Xm)
  universe <- unique(strip_prefix(colnames(X_ef)))
  list(Xr = Xr, Xm = Xm, X_ef = X_ef, y = y, universe = universe,
       n = length(y))
}

## Stratified 70/30 split; seed controls the run.
make_split <- function(y, seed) {
  set.seed(seed)
  te <- unlist(lapply(split(seq_len(length(y)), y),
    function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
  list(tr = setdiff(seq_len(length(y)), te), te = te)
}

## ---------------- grouping ----------------
# Per-run grouping (current protocol): single hclust on train.
define_groups <- function(X) {
  cm <- cor(X, method = "spearman"); dm <- as.dist(1 - abs(cm))
  hc <- hclust(dm, method = "ward.D2"); k <- max(1, round(ncol(X) / 8))
  cutree(hc, k = k)
}

# S1: consensus groups. Bootstrap-averaged |Spearman| on the FULL cohort
# (unsupervised: no labels used), one hclust, fixed for every run/perm.
# Group structure is a property of the feature space, not of the split.
consensus_groups <- function(X, B = 10, seed = 1) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  cm_sum <- matrix(0, p, p)
  for (b in seq_len(B)) {
    bi <- sample(n, n, replace = TRUE)
    cm_sum <- cm_sum + abs(suppressWarnings(cor(X[bi, , drop = FALSE], method = "spearman")))
  }
  cm_sum[is.na(cm_sum)] <- 0
  dm <- as.dist(1 - cm_sum / B)
  k <- max(1, round(p / 8))
  cutree(hclust(dm, method = "ward.D2"), k = k)
}

## ---------------- adaptive weights ----------------
# Current protocol: w_g = sqrt(|g|)/(1+RMS_z_g), z = univariate logistic Wald.
valid_adaptive_pen <- function(X, y, groups) {
  X <- as.matrix(X)
  z <- apply(X, 2, function(xj) {
    f <- suppressWarnings(glm(y ~ xj, family = binomial()))
    s <- suppressWarnings(summary(f)$coefficients)
    if (nrow(s) < 2 || is.na(s[2, 3])) return(0)
    unname(s[2, 3])
  })
  group_pen_from_z(z, groups)
}
group_pen_from_z <- function(z, groups) {
  ug <- sort(unique(groups))
  rms <- sapply(ug, function(g) { zj <- z[names(groups)[groups == g]]; sqrt(mean(zj^2, na.rm = TRUE)) })
  gs <- sapply(ug, function(g) sum(groups == g))
  rp <- sqrt(gs) / (1 + rms); rp <- rp / mean(rp)
  setNames(rp, ug)[as.character(ug)]
}
# S2: shrink noisy univariate z's (soft-threshold at |z|=1) before RMS
# aggregation, so weak-signal jitter stops shaking the group penalties.
shrunk_adaptive_pen <- function(X, y, groups, thr = 1.0) {
  X <- as.matrix(X)
  z <- apply(X, 2, function(xj) {
    f <- suppressWarnings(glm(y ~ xj, family = binomial()))
    s <- suppressWarnings(summary(f)$coefficients)
    if (nrow(s) < 2 || is.na(s[2, 3])) return(0)
    unname(s[2, 3])
  })
  zs <- sign(z) * pmax(0, abs(z) - thr)
  group_pen_from_z(zs, groups)
}

## ---------------- selection ----------------
# Permutation-assisted group-lasso selection with FDR-targeted lambda
# (verbatim from tcga-brca-5arm.R).
perm_select <- function(X, y, n_perm, groups, gpen, q, cutoff, tag, logmsg = NULL) {
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
  if (!is.null(logmsg)) logmsg("[%s] perms voting: %d/%d", tag, nperm_used, n_perm)
  list(sel = names(fr)[fr >= cutoff], freq = fr)
}

# S3: stability selection (complementary pairs, Shah & Samworth style).
# B subsamples; within each, perm_select with fewer perms; keep features
# selected in >= pi_thr of subsamples. Groups + penalties fixed (S1+S2).
stab_select <- function(X, y, groups, gpen, B = 30, inner_perm = 15,
                        inner_cut = 0.4, pi_thr = 0.6, seed = 1, tag = "S3") {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X); fn <- colnames(X)
  votes <- setNames(rep(0, p), fn)
  half <- floor(n / 2)
  for (b in seq_len(B / 2)) {
    prm <- sample(n)
    for (s in list(prm[1:half], prm[(half + 1):n])) {
      sel <- perm_select(X[s, , drop = FALSE], y[s], inner_perm, groups, gpen,
                         FDR_TARGET, inner_cut, paste0(tag, "-sub"))$sel
      votes[sel] <- votes[sel] + 1
    }
  }
  list(sel = names(votes)[votes / B >= pi_thr], freq = votes / B)
}

## ---------------- evaluation ----------------
refit_auc <- function(Xtr, ytr, Xte, yte, sel) {
  if (length(sel) < 2) return(NA_real_)
  fit <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial")
  pr <- as.numeric(predict(fit, Xte[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  as.numeric(auc(roc(yte, pr, quiet = TRUE)))
}
mod_cv_auc <- function(Xtr, ytr, sel) {
  if (length(sel) < 2) return(NA_real_)
  fit <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial", nfolds = 5)
  pr <- as.numeric(predict(fit, Xtr[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  as.numeric(auc(roc(ytr, pr, quiet = TRUE)))
}
enrichment <- function(sel, universe) {
  sl <- unique(strip_prefix(intersect(sel, universe)))
  bind_rows(lapply(names(SIGS), function(sg_) {
    sg <- intersect(SIGS[[sg_]], universe)
    qq <- length(intersect(sl, sg))
    data.frame(signature = sg_, overlap = qq, sig_in_universe = length(sg),
      selected = length(sl),
      p = signif(phyper(qq - 1, length(sg), length(universe) - length(sg),
        length(sl), lower.tail = FALSE), 3))
  }))
}
# Cheap "is the endpoint too easy?" baselines: ESR1-only and 4-control-gene
# logistic models, fit on train, AUC on test.
simple_baselines <- function(Xr, y, tr, te) {
  df_tr <- data.frame(y = y[tr], Xr[tr, POS_CTRL, drop = FALSE])
  df_te <- data.frame(Xr[te, POS_CTRL, drop = FALSE])
  auc1 <- tryCatch({
    f1 <- suppressWarnings(glm(y ~ ESR1, data = df_tr, family = binomial()))
    pr <- as.numeric(predict(f1, newdata = df_te, type = "response"))
    as.numeric(auc(roc(y[te], pr, quiet = TRUE)))
  }, error = function(e) NA_real_)
  auc4 <- tryCatch({
    f4 <- suppressWarnings(glm(y ~ ESR1 + PGR + FOXA1 + GATA3, data = df_tr, family = binomial()))
    pr <- as.numeric(predict(f4, newdata = df_te, type = "response"))
    as.numeric(auc(roc(y[te], pr, quiet = TRUE)))
  }, error = function(e) NA_real_)
  c(esr1_only = auc1, four_gene = auc4)
}
pairwise_jaccard <- function(sel_list) {
  cmb <- combn(length(sel_list), 2)
  mean(apply(cmb, 2, function(ii) jaccard(sel_list[[ii[1]]], sel_list[[ii[2]]])))
}

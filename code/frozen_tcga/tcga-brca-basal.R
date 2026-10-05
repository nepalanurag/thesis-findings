#!/usr/bin/env Rscript
# Secondary endpoint: PAM50 Basal-like vs rest (genefu calls) on TCGA-BRCA.
# Same preprocessing + selection engine as tcga-brca-real.R.
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(SKAT); library(Matrix); library(genefu)})
set.seed(20260909)
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"
dir.create(OUT, showWarnings = FALSE)
logf <- file.path(OUT, "run_log_basal.txt")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

define_groups <- function(feature_df) {
  cm <- cor(feature_df, method = "spearman"); dm <- as.dist(1 - abs(cm))
  hc <- hclust(dm, method = "ward.D2"); k <- max(1, round(ncol(feature_df) / 8))
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
perm_grplasso <- function(X, y, n_perm, fgroups, gpen, cutoff) {
  X <- as.matrix(X); p <- ncol(X); fn <- colnames(X)
  kg <- fgroups + max(fgroups); fp <- c(gpen, gpen); hits <- c()
  for (k in seq_len(n_perm)) {
    Xk <- X[sample(nrow(X)), ]
    fit <- suppressWarnings(grpreg(cbind(X, Xk), y, group = c(fgroups, kg),
      penalty = "grLasso", family = "binomial", group.multiplier = fp,
      nlambda = 50, lambda.min = 0.05))
    b <- as.matrix(fit$beta[-1, ])
    d <- colSums(abs(b[1:p, ]) > 0) - colSums(abs(b[(p + 1):(2 * p), ]) > 0)
    oi <- if (all(d <= 0)) NULL else which.max(d)
    if (!is.null(oi)) { si <- which(abs(b[1:p, oi]) > 0)
      if (length(si)) hits <- c(hits, rownames(b)[si]) }
  }
  fc <- setNames(rep(0, p), fn); ct <- table(hits); fc[names(ct)] <- ct
  fr <- fc / n_perm
  list(sel = names(fr)[fr >= cutoff], freq = fr)
}
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }

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

logmsg("loading RNA...")
rna <- suppressMessages(read_tsv(file.path(DATA, "data_mrna_seq_v2_rsem.txt"),
  comment = "#", show_col_types = FALSE))
genes <- rna$Hugo_Symbol
rna <- as.matrix(rna[, -(1:2)]); rownames(rna) <- genes
is_primary <- substr(colnames(rna), 14, 15) == "01"
rna <- rna[, is_primary]
pat <- substr(colnames(rna), 1, 12)
ord <- order(pat, !(substr(colnames(rna), 14, 16) == "01A"))
keep <- !duplicated(pat[ord]); rna <- rna[, ord[keep]]
X_all <- log2(t(rna) + 1)
X_all <- X_all[, colSums(is.na(X_all)) == 0]
v0 <- apply(X_all, 2, var)
keep_g <- tapply(seq_len(ncol(X_all)), colnames(X_all), function(ii) ii[which.max(v0[ii])])
X_all <- X_all[, unlist(keep_g)]
X_all <- X_all[, apply(X_all, 2, var) > 1e-6]
v <- apply(X_all, 2, var)
top_genes <- names(sort(v, decreasing = TRUE))[1:1500]
X <- t(scale(t(X_all[, top_genes]))); X[is.na(X)] <- 0
universe <- colnames(X)

logmsg("computing PAM50 (genefu)...")
data(pam50.robust, package = "genefu")
s <- molecular.subtyping(sbt.model = "pam50", data = t(X_all),
  annot = data.frame(probe = colnames(X_all), EntrezGene.ID = NA,
    Gene.Symbol = colnames(X_all)), do.mapping = FALSE)
pam50 <- setNames(as.character(s$subtype), rownames(X_all))
print(table(pam50))
y <- as.integer(pam50[rownames(X)] == "Basal")
logmsg("Basal: %d / %d", sum(y), length(y))

n <- length(y)
te <- unlist(lapply(split(seq_len(n), y), function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
tr <- setdiff(seq_len(n), te)
Xtr <- X[tr, ]; ytr <- y[tr]; Xte <- X[te, ]; yte <- y[te]
logmsg("train %d (basal %d) | test %d (basal %d)", length(tr), sum(ytr), length(te), sum(yte))
t0 <- Sys.time()
grp <- define_groups(Xtr); pen <- skat_pen(Xtr, ytr, grp$groups)
sel <- perm_grplasso(Xtr, ytr, 100, grp$groups, pen, 0.5)
logmsg("selected %d genes in %.1f min", length(sel$sel),
  as.numeric(difftime(Sys.time(), t0, units = "mins")))
auc_te <- NA_real_
if (length(sel$sel) >= 2) {
  fit <- cv.glmnet(Xtr[, sel$sel, drop = FALSE], ytr, family = "binomial")
  pr <- as.numeric(predict(fit, Xte[, sel$sel, drop = FALSE], s = "lambda.min", type = "response"))
  auc_te <- as.numeric(auc(roc(yte, pr, quiet = TRUE)))
}
logmsg("held-out test AUC: %.3f", auc_te)
ov <- bind_rows(lapply(names(SIGS), function(sg_) {
  sg <- intersect(SIGS[[sg_]], universe); sl <- intersect(sel$sel, universe)
  q <- length(intersect(sl, sg))
  data.frame(signature = sg_, overlap = q, sig_in_universe = length(sg),
    selected = length(sl), universe = length(universe),
    enrich = round((q / length(sl)) / (length(sg) / length(universe)), 2),
    p = signif(phyper(q - 1, length(sg), length(universe) - length(sg),
      length(sl), lower.tail = FALSE), 3))
}))
print(as.data.frame(ov))
stab <- replicate(5, { bi <- unlist(lapply(split(seq_along(ytr), ytr),
  function(ii) sample(ii, floor(0.8 * length(ii)))))
  perm_grplasso(Xtr[bi, ], ytr[bi], 100, grp$groups, pen, 0.5)$sel })
pw <- combn(5, 2, function(ii) jaccard(stab[[ii[1]]], stab[[ii[2]]]))
logmsg("mean pairwise Jaccard: %.3f", mean(pw))
write.csv(data.frame(gene = sel$sel, freq = round(as.numeric(sel$freq[sel$sel]), 3)),
  file.path(OUT, "selected_PAM50_Basal_vs_rest.csv"), row.names = FALSE)
write.csv(ov, file.path(OUT, "overlap_PAM50_Basal_vs_rest.csv"), row.names = FALSE)
saveRDS(list(selected = sel$sel, freq = sel$freq, auc = auc_te, overlap = ov,
  jaccard = mean(pw), pam50_table = table(pam50)), file.path(OUT, "basal_results.rds"))
logmsg("DONE")

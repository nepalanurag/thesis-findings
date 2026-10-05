#!/usr/bin/env Rscript
# Basal-vs-rest endpoint with V3 machinery (valid adaptive weights +
# FDR-targeted lambda). Sanity check: the fixed method should recover basal
# markers (KRT5/14/17, FOXC1...). NOTE the circularity caveat: PAM50 calls are
# computed from this same RNA-seq data, so this is a consistency check, not an
# independent validation like the ER-IHC endpoint.
suppressPackageStartupMessages({library(tidyverse); library(glmnet); library(grpreg);
  library(pROC); library(Matrix); library(genefu)})
set.seed(20260909)
FDR_TARGET <- 0.10; N_PERM <- 100
DATA <- "tcga_data"; OUT <- "tcga_brca_real_results"
dir.create(OUT, showWarnings = FALSE)
logf <- file.path(OUT, "run_log_basal_v3.txt")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

PAM50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
 "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","FOXC1",
 "GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67",
 "MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1",
 "RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")
BASAL_MARKERS <- c("KRT5", "KRT14", "KRT17", "FOXC1", "EGFR")

define_groups <- function(X) {
  cm <- cor(X, method = "spearman"); dm <- as.dist(1 - abs(cm))
  hc <- hclust(dm, method = "ward.D2"); k <- max(1, round(ncol(X) / 8))
  list(groups = cutree(hc, k = k))
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
  rms <- sapply(ug, function(g) {
    zj <- z[names(groups)[groups == g]]
    sqrt(mean(zj^2, na.rm = TRUE))
  })
  gs <- sapply(ug, function(g) sum(groups == g))
  rp <- sqrt(gs) / (1 + rms); rp <- rp / mean(rp)
  setNames(rp, ug)[as.character(ug)]
}
perm_select_v3 <- function(X, y, n_perm, fgroups, gpen, q, cutoff, tag) {
  X <- as.matrix(X); p <- ncol(X); fn <- colnames(X)
  kg <- fgroups + max(fgroups); fp <- c(gpen, gpen)
  hits <- c(); nperm_used <- 0
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
    si <- which(abs(b[1:p, j]) > 0)
    if (length(si)) hits <- c(hits, rownames(b)[si])
  }
  fc <- setNames(rep(0, p), fn); ct <- table(hits); fc[names(ct)] <- ct
  fr <- fc / n_perm
  logmsg("[%s] perms voting: %d/%d", tag, nperm_used, n_perm)
  list(sel = names(fr)[fr >= cutoff], freq = fr)
}
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }

logmsg("loading RNA...")
rna <- suppressMessages(read_tsv(file.path(DATA, "data_mrna_seq_v2_rsem.txt"),
  comment = "#", show_col_types = FALSE))
genes <- rna$Hugo_Symbol
rna <- as.matrix(rna[, -(1:2)]); rownames(rna) <- genes
clin_s <- suppressMessages(read_tsv(file.path(DATA, "data_clinical_sample.txt"),
  comment = "#", show_col_types = FALSE))
is_primary <- clin_s$SAMPLE_TYPE == "Primary"
rna <- rna[, colnames(rna) %in% clin_s$SAMPLE_ID[is_primary]]
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
top_genes <- names(sort(v, decreasing = TRUE))[1:min(1500, length(v))]
X <- t(scale(t(X_all[, top_genes]))); X[is.na(X)] <- 0
universe <- colnames(X)

logmsg("computing PAM50 (genefu)...")
data(pam50.robust, package = "genefu")
s <- molecular.subtyping(sbt.model = "pam50", data = X_all,
  annot = data.frame(probe = colnames(X_all), EntrezGene.ID = NA,
    Gene.Symbol = colnames(X_all)), do.mapping = FALSE)
pam50 <- setNames(as.character(s$subtype), rownames(X_all))
print(table(pam50))
y <- as.integer(pam50[rownames(X)] == "Basal")
logmsg("Basal: %d / %d (%.1f%%)", sum(y), length(y), 100 * mean(y))

n <- length(y)
te <- unlist(lapply(split(seq_len(n), y), function(ii) sample(ii, max(1, round(length(ii) * 0.3)))))
tr <- setdiff(seq_len(n), te)
Xtr <- X[tr, ]; ytr <- y[tr]; Xte <- X[te, ]; yte <- y[te]
logmsg("train %d (basal %d) | test %d (basal %d)", length(tr), sum(ytr), length(te), sum(yte))

t0 <- Sys.time()
grp <- define_groups(Xtr)
pen <- valid_adaptive_pen(Xtr, ytr, grp$groups)
sel <- perm_select_v3(Xtr, ytr, N_PERM, grp$groups, pen, FDR_TARGET, 0.5, "basal")
auc_te <- NA_real_
if (length(sel$sel) >= 2) {
  fit <- cv.glmnet(Xtr[, sel$sel, drop = FALSE], ytr, family = "binomial")
  pr <- as.numeric(predict(fit, Xte[, sel$sel, drop = FALSE], s = "lambda.min", type = "response"))
  auc_te <- as.numeric(auc(roc(yte, pr, quiet = TRUE)))
}
fr_ord <- sort(sel$freq, decreasing = TRUE)
bm_rank <- sapply(BASAL_MARKERS, function(g)
  if (g %in% names(fr_ord)) which(names(fr_ord) == g) else NA_integer_)
sg <- intersect(PAM50, universe); sl <- sel$sel
q <- length(intersect(sl, sg))
enr <- (q / length(sl)) / (length(sg) / length(universe))
pv <- phyper(q - 1, length(sg), length(universe) - length(sg), length(sl), lower.tail = FALSE)
stab <- lapply(1:3, function(b) {
  bi <- unlist(lapply(split(seq_along(ytr), ytr),
    function(ii) sample(ii, floor(0.8 * length(ii)))))
  perm_select_v3(Xtr[bi, ], ytr[bi], 50, grp$groups, pen, FDR_TARGET, 0.5, "basal-stab")$sel })
pw <- combn(3, 2, function(ii) jaccard(stab[[ii[1]]], stab[[ii[2]]]))
logmsg("basal: selected %d genes (%.1f min); test AUC %.3f; Jaccard %.3f",
  length(sel$sel), as.numeric(difftime(Sys.time(), t0, units = "mins")), auc_te, mean(pw))
logmsg("basal marker ranks: %s",
  paste(sprintf("%s:%s", BASAL_MARKERS, bm_rank), collapse = " "))
logmsg("PAM50 overlap: %d/%d, enrich %.2f, p=%.2g", q, length(sg), enr, pv)
write.csv(data.frame(gene = sel$sel, freq = round(as.numeric(sel$freq[sel$sel]), 3)),
  file.path(OUT, "selected_Basal_vs_rest_v3.csv"), row.names = FALSE)
saveRDS(list(selected = sel$sel, freq = sel$freq, auc = auc_te, jaccard = mean(pw),
  bm_rank = bm_rank, pam50_overlap = c(q = q, enrich = unname(enr), p = pv),
  pam50_table = table(pam50)), file.path(OUT, "basal_v3_results.rds"))
logmsg("BASAL V3 DONE")

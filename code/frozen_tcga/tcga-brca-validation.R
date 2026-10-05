#!/usr/bin/env Rscript
# =============================================================================
# TCGA-BRCA Validation of Permutation-Assisted Adaptive Group Lasso
# (Integrative Feature Selection for Multi-Modal Data)
#
# What this does:
#   1. Downloads TCGA-BRCA RNA-seq (STAR counts) + 450K methylation + clinical
#      + published PAM50 subtype calls.
#   2. Runs the thesis method (early-fusion, SKAT-adaptive, permutation-assisted
#      group lasso) to select features associated with Basal-like vs rest.
#   3. Validates the selected signature the way reviewers expect:
#        a) hypergeometric overlap vs published BRCA signatures
#           (PAM50 / Oncotype DX / MammaPrint)
#        b) held-out test AUC (no selection/validation circularity)
#        c) stability across subsamples (Jaccard)
#
# Key corrections vs the draft TCGA-Selection.rmd:
#   - Endpoint is PAM50 subtype (Basal vs rest), NOT vital_status.
#   - Primary-tumor ("-01") samples only; no tumor/normal mixing.
#   - Feature selection on TRAIN split only; AUC on held-out TEST split.
#   - Overlap tests use the number of genes actually tested as background.
#
# Requirements (install once):
#   install.packages(c("tidyverse","glmnet","grpreg","pROC","SKAT","Matrix"))
#   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("TCGAbiolinks","SummarizedExperiment",
#                          "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
#
# Run:  Rscript tcga-brca-validation.R
# Time: data download ~30-60 min first run; analysis ~1-3 h on 4+ cores.
#       Reduce N_PERM / N_TOP_* for a quick pilot (see CONFIG).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(grpreg)
  library(pROC)
  library(SKAT)
  library(Matrix)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

set.seed(20260908)

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
CONFIG <- list(
  n_top_rna    = 1500,   # pre-screen: top variable RNA genes
  n_top_meth   = 1500,   # pre-screen: top variable methylation probes
  n_perm_wt    = 200,    # permutations for adaptive weight p-values
  n_perm_sel   = 100,    # knockoff permutations inside the selector
  sel_cutoff   = 0.5,    # selection frequency cutoff
  test_frac    = 0.30,   # held-out test fraction (stratified)
  n_stab_runs  = 5,      # subsample runs for stability
  out_dir      = "tcga_brca_results"
)
dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(CONFIG$out_dir, "run_log.txt")
logmsg <- function(...) { msg <- sprintf(...); cat(msg, "\n"); cat(msg, "\n", file = logf, append = TRUE) }

# -----------------------------------------------------------------------------
# 1. PUBLISHED BRCA SIGNATURES (HGNC current symbols; see research notes)
# -----------------------------------------------------------------------------
PAM50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1",
  "CDC20","CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T")

ONCOTYPE_DX <- c("MKI67","AURKA","BIRC5","CCNB1","MYBL2",          # proliferation
  "ERBB2","GRB7",                                                  # HER2
  "ESR1","PGR","BCL2","SCUBE2",                                    # estrogen
  "MMP11","CTSV",                                                  # invasion
  "GSTM1","CD68","BAG1",                                           # other
  "ACTB","GAPDH","RPLP0","GUSB","TFRC")                            # reference

# MammaPrint: 70 reporters -> 55 uniquely mappable current HGNC symbols
MAMMAPRINT_55 <- c("ADGRG6","AKAP2","ALDH4A1","AP2B1","BBC3","CCN4","CCNE2",
  "CDC42BPA","CDCA7","CENPA","CMC2","COL4A2","DCK","DIAPH3","DTL","ECI2",
  "ECT2","ESM1","EXT1","FGF18","FLT1","GMPS","GNAZ","GPR180","GSTM3",
  "IGFBP5","LPCAT1","MCM6","MELK","MMP9","MS4A7","MSANTD3","MTDH","NDC80",
  "NMU","NUSAP1","ORC6","OXCT1","PITRM1","PLAAT1","PRC1","QSOX2","RAB6B",
  "RFC4","RTN4RL1","RUNDC1","SCUBE2","SERF1A","SLC2A3","STK32B","TGFB3",
  "TMEM74B","TSPYL5","UCHL5","ZNF385B")

SIGNATURES <- list(PAM50 = PAM50, OncotypeDX = ONCOTYPE_DX, MammaPrint = MAMMAPRINT_55)

# -----------------------------------------------------------------------------
# 2. THESIS METHOD (from Fusion_work.rmd, unchanged logic)
# -----------------------------------------------------------------------------
define_groups <- function(feature_df) {
  cor_matrix  <- cor(feature_df, method = "spearman")
  dist_matrix <- as.dist(1 - abs(cor_matrix))
  hclust_obj  <- hclust(dist_matrix, method = "ward.D2")
  num_groups  <- max(1, round(ncol(feature_df) / 8))
  feature_groups_vec <- cutree(hclust_obj, k = num_groups)
  list(groups = feature_groups_vec, hclust = hclust_obj, cor = cor_matrix)
}

calculate_skat_penalties <- function(X, y, groups) {
  obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
  unique_groups <- sort(unique(groups))
  p_values   <- numeric(length(unique_groups)); names(p_values) <- unique_groups
  group_sizes <- numeric(length(unique_groups))
  X_mat <- as.matrix(X)
  for (g in unique_groups) {
    g_features <- names(groups)[groups == g]
    Z <- X_mat[, g_features, drop = FALSE]
    group_sizes[as.character(g)] <- length(g_features)
    p_values[as.character(g)] <- tryCatch({
      suppressWarnings(SKATBinary(Z, obj, kernel = "Linear.Weighted")$p.value)
    }, error = function(e) 1)
  }
  neg_log_p <- -log(p_values + 1e-6); neg_log_p[neg_log_p < 0] <- 0
  raw_penalty <- sqrt(group_sizes) * (1 / (neg_log_p + 1))
  final_penalty <- raw_penalty / mean(raw_penalty)
  final_penalty[as.character(sort(unique(groups)))]
}

perm_assisted_grplasso_binary <- function(X, y, n_perm, feature_groups, group_penalties, cutoff) {
  X <- as.matrix(X); p <- ncol(X)
  all_feature_names <- colnames(X)
  knockoff_groups <- feature_groups + max(feature_groups)
  full_group_penalties <- c(group_penalties, group_penalties)
  selected_features_all_perms <- c()
  for (k in seq_len(n_perm)) {
    X_knockoff <- X[sample(nrow(X)), ]
    fit <- suppressWarnings(grpreg(cbind(X, X_knockoff), y,
      group = c(feature_groups, knockoff_groups),
      penalty = "grLasso", family = "binomial",
      group.multiplier = full_group_penalties,
      nlambda = 50, lambda.min = 0.05))
    betas <- as.matrix(fit$beta[-1, ])
    n_original <- colSums(abs(betas[1:p, ]) > 0)
    n_knockoff <- colSums(abs(betas[(p + 1):(2 * p), ]) > 0)
    diff <- n_original - n_knockoff
    optimal_idx <- if (all(diff <= 0)) NULL else which.max(diff)
    if (!is.null(optimal_idx)) {
      si <- which(abs(betas[1:p, optimal_idx]) > 0)
      if (length(si) > 0)
        selected_features_all_perms <- c(selected_features_all_perms, rownames(betas)[si])
    }
  }
  freq_counts <- setNames(rep(0, p), all_feature_names)
  if (length(selected_features_all_perms) > 0) {
    ct <- table(selected_features_all_perms)
    freq_counts[names(ct)] <- ct
  }
  feature_frequencies <- freq_counts / n_perm
  list(final_selection = names(feature_frequencies)[feature_frequencies >= cutoff],
       feature_frequencies = feature_frequencies)
}

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (u == 0) return(1)
  length(intersect(a, b)) / u
}

# -----------------------------------------------------------------------------
# 3. DATA: TCGA-BRCA via TCGAbiolinks (GDC harmonized, hg38)
# -----------------------------------------------------------------------------
logmsg("=== STEP 1: querying GDC ===")

# --- 3a. RNA-seq STAR counts (primary tumor only) ---
q_rna <- GDCquery(project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor")
GDCdownload(q_rna, method = "api")
rna_se <- GDCprepare(q_rna, summarizedExperiment = TRUE)
# unstranded counts assay
assay_names <- names(assays(rna_se))
logmsg("RNA assays available: %s", paste(assay_names, collapse = ", "))
counts <- assay(rna_se, "unstranded")
# gene symbols
rna_genes <- as.character(rowData(rna_se)$gene_name)
rownames(counts) <- ifelse(is.na(rna_genes) | rna_genes == "",
  rownames(counts), rna_genes)
# patient barcodes (sample level -> keep one aliquot per patient, prefer vial A)
rna_barcodes_full <- colnames(counts)
rna_patients <- substr(rna_barcodes_full, 1, 12)

# --- 3b. Methylation 450K beta values (primary tumor only) ---
q_meth <- GDCquery(project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  sample.type = "Primary Tumor")
GDCdownload(q_meth, method = "api")
meth_se <- GDCprepare(q_meth, summarizedExperiment = TRUE)
beta <- assay(meth_se, 1)
meth_barcodes_full <- colnames(beta)
meth_patients <- substr(meth_barcodes_full, 1, 12)
logmsg("Methylation probes: %d", nrow(beta))

# --- 3c. PAM50 subtype calls (published marker-paper calls) ---
subtype <- TCGAquery_subtype("BRCA")
# keep: patient barcode + PAM50 call; drop Normal-like (not a tumor intrinsic signal)
subtype$patient <- substr(subtype$patient, 1, 12)
subtype <- subtype[!is.na(subtype$BRCA_Subtype_PAM50) &
                   subtype$BRCA_Subtype_PAM50 != "Normal", ]
logmsg("Subtype table: %d patients with PAM50 calls", nrow(subtype))

# --- 3d. Match patients across RNA x methylation x subtype ---
dedupe_one <- function(patients, full) {
  # one sample per patient: prefer vial A (barcode positions 14-15 == "01A")
  ord <- order(patients, !(substr(full, 14, 16) == "01A"))
  keep <- !duplicated(patients[ord])
  ord[keep]
}
ri <- dedupe_one(rna_patients, rna_barcodes_full)
mi <- dedupe_one(meth_patients, meth_barcodes_full)
rna_patients_u <- rna_patients[ri]; meth_patients_u <- meth_patients[mi]

common <- Reduce(intersect, list(rna_patients_u, meth_patients_u, subtype$patient))
logmsg("Matched patients (RNA x METH x PAM50): %d", length(common))
stopifnot(length(common) >= 50)

counts <- counts[, ri[match(common, rna_patients_u)]]
beta   <- beta[, mi[match(common, meth_patients_u)]]
sub    <- subtype[match(common, subtype$patient), ]
stopifnot(all(substr(colnames(counts), 1, 12) == common))

# --- 3e. Outcome: Basal-like vs rest (binary) ---
y <- as.integer(sub$BRCA_Subtype_PAM50 == "Basal")
logmsg("Outcome: Basal-like = %d / %d (%.1f%%)", sum(y), length(y), 100 * mean(y))
print(table(PAM50 = sub$BRCA_Subtype_PAM50))

# -----------------------------------------------------------------------------
# 4. PREPROCESSING + PRE-SCREENING
# -----------------------------------------------------------------------------
logmsg("=== STEP 2: preprocessing ===")

# RNA: log2(counts+1), drop all-zero / NA genes, top variable
rna_mat <- log2(as.matrix(counts) + 1)
rna_mat <- rna_mat[rowSums(is.na(rna_mat)) == 0, ]
rna_mat <- rna_mat[apply(rna_mat, 1, var) > 1e-6, ]
# collapse duplicated gene symbols (keep highest-variance copy)
if (anyDuplicated(rownames(rna_mat))) {
  v <- apply(rna_mat, 1, var)
  keep <- tapply(seq_len(nrow(rna_mat)), rownames(rna_mat),
                 function(ii) ii[which.max(v[ii])])
  rna_mat <- rna_mat[unlist(keep), ]
}
rna_var <- apply(rna_mat, 1, var)
rna_keep <- names(sort(rna_var, decreasing = TRUE))[seq_len(min(CONFIG$n_top_rna, length(rna_var)))]
X_rna <- t(scale(t(rna_mat[rna_keep, ])))          # samples x genes, z-scored
colnames(X_rna) <- paste0("RNA_", colnames(X_rna))
rna_gene_universe <- sub("^RNA_", "", colnames(X_rna))
logmsg("RNA: %d genes -> top %d variable", nrow(rna_mat), ncol(X_rna))

# Methylation: drop probes with >20% NA, mean-impute rest, top variable
na_frac <- rowMeans(is.na(beta))
beta <- beta[na_frac <= 0.2, ]
for (j in which(rowSums(is.na(beta)) > 0)) {
  beta[j, is.na(beta[j, ])] <- mean(beta[j, ], na.rm = TRUE)
}
meth_var <- apply(beta, 1, var)
meth_keep <- names(sort(meth_var, decreasing = TRUE))[seq_len(min(CONFIG$n_top_meth, length(meth_var)))]
X_meth <- t(scale(t(beta[meth_keep, ])))
colnames(X_meth) <- paste0("METH_", colnames(X_meth))
logmsg("METH: top %d variable probes", ncol(X_meth))

X_early <- cbind(X_rna, X_meth)
logmsg("Early-fusion matrix: %d samples x %d features", nrow(X_early), ncol(X_early))

# sanity: Basal positive control should exist in RNA universe
logmsg("Positive-control check: ESR1/PGR/ERBB2/MKI67 in RNA universe: %s",
  paste(c("ESR1","PGR","ERBB2","MKI67") %in% rna_gene_universe, collapse = ","))

# -----------------------------------------------------------------------------
# 5. TRAIN / TEST SPLIT (stratified; selection NEVER sees the test set)
# -----------------------------------------------------------------------------
set.seed(20260908)
n <- length(y)
test_idx <- unlist(lapply(split(seq_len(n), y), function(ii)
  sample(ii, size = max(1, round(length(ii) * CONFIG$test_frac)))))
train_idx <- setdiff(seq_len(n), test_idx)
logmsg("Train n=%d (Basal %d), Test n=%d (Basal %d)",
  length(train_idx), sum(y[train_idx]), length(test_idx), sum(y[test_idx]))

Xtr <- X_early[train_idx, ]; ytr <- y[train_idx]
Xte <- X_early[test_idx, ];  yte <- y[test_idx]

# -----------------------------------------------------------------------------
# 6. FEATURE SELECTION ON TRAIN (thesis method, early fusion, SKAT-adaptive)
# -----------------------------------------------------------------------------
logmsg("=== STEP 3: feature selection on train ===")
t0 <- Sys.time()
grp <- define_groups(Xtr)
logmsg("Groups: %d", length(unique(grp$groups)))
pen <- calculate_skat_penalties(Xtr, ytr, grp$groups)
sel <- perm_assisted_grplasso_binary(Xtr, ytr, CONFIG$n_perm_sel,
  grp$groups, pen, CONFIG$sel_cutoff)
logmsg("Selection took %.1f min; selected %d features",
  as.numeric(difftime(Sys.time(), t0, units = "mins")), length(sel$final_selection))

sel_genes  <- sub("^RNA_", "", grep("^RNA_", sel$final_selection, value = TRUE))
sel_probes <- sub("^METH_", "", grep("^METH_", sel$final_selection, value = TRUE))
logmsg("Selected: %d RNA genes, %d METH probes", length(sel_genes), length(sel_probes))

write.csv(data.frame(feature = sel$final_selection,
    frequency = round(as.numeric(sel$feature_frequencies[sel$final_selection]), 3)),
  file.path(CONFIG$out_dir, "selected_features.csv"), row.names = FALSE)

# --- map selected methylation probes -> genes (promoter CpGs preferred) ---
probe2gene <- tryCatch({
  suppressPackageStartupMessages(
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann <- ann[rownames(ann) %in% sel_probes, , drop = FALSE]
  # UCSC_RefGene_Name may be ";" separated; keep first gene
  g <- sapply(strsplit(as.character(ann$UCSC_RefGene_Name), ";"),
    function(z) if (length(z) == 0 || all(is.na(z))) NA_character_ else z[1])
  data.frame(probe = rownames(ann), gene = g, row.names = NULL)
}, error = function(e) { logmsg("Probe annotation unavailable: %s", conditionMessage(e)); NULL })
if (!is.null(probe2gene))
  write.csv(probe2gene, file.path(CONFIG$out_dir, "selected_probes_genes.csv"), row.names = FALSE)
sel_meth_genes <- if (!is.null(probe2gene))
  unique(na.omit(probe2gene$gene)) else character(0)

# combined gene-level selected set for signature overlap tests
sel_genes_all <- unique(c(sel_genes, sel_meth_genes))

# -----------------------------------------------------------------------------
# 7. VALIDATION A: hypergeometric overlap vs published signatures
# -----------------------------------------------------------------------------
logmsg("=== STEP 4a: signature overlap tests ===")
hyper_overlap <- function(selected, signature, universe) {
  signature <- intersect(signature, universe)
  selected  <- intersect(selected, universe)
  q <- length(intersect(selected, signature))  # overlap
  m <- length(signature)                       # signature genes in universe
  N <- length(universe); k <- length(selected)
  # P(X >= q) under null
  pval <- phyper(q - 1, m, N - m, k, lower.tail = FALSE)
  enrich <- (q / k) / (m / N)                  # fold enrichment
  data.frame(overlap = q, signature_in_universe = m, selected_in_universe = k,
    universe = N, fold_enrichment = round(enrich, 2), p_value = signif(pval, 3))
}
# background = genes actually tested (RNA universe + mapped METH genes universe)
meth_gene_universe <- if (!is.null(probe2gene)) {
  ann_all <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  unique(na.omit(sapply(strsplit(as.character(
    ann_all$UCSC_RefGene_Name[rownames(ann_all) %in% sub("^METH_","",colnames(X_meth))]),
    ";"), function(z) z[1])))
} else character(0)
universe <- unique(c(rna_gene_universe, meth_gene_universe))
logmsg("Overlap background: %d genes tested", length(universe))

overlap_res <- bind_rows(lapply(names(SIGNATURES), function(s)
  cbind(signature = s, hyper_overlap(sel_genes_all, SIGNATURES[[s]], universe))))
print(overlap_res)
write.csv(overlap_res, file.path(CONFIG$out_dir, "signature_overlap.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. VALIDATION B: held-out TEST AUC (no circularity)
# -----------------------------------------------------------------------------
logmsg("=== STEP 4b: held-out test AUC ===")
heldout_auc <- function(Xtr, ytr, Xte, yte, sel) {
  if (length(sel) < 2) return(list(auc = 0.5, n = length(sel)))
  fit <- cv.glmnet(Xtr[, sel, drop = FALSE], ytr, family = "binomial", alpha = 1)
  pr  <- as.numeric(predict(fit, Xte[, sel, drop = FALSE], s = "lambda.min", type = "response"))
  list(auc = as.numeric(auc(roc(yte, pr, quiet = TRUE))), n = length(sel))
}
res_test <- heldout_auc(Xtr, ytr, Xte, yte, sel$final_selection)
logmsg("Held-out test AUC (selected features, lasso refit): %.3f", res_test$auc)

# baselines on the same test set
bl_rna  <- heldout_auc(Xtr, ytr, Xte, yte,
  colnames(Xtr)[grepl("^RNA_", colnames(Xtr))][seq_len(min(500, sum(grepl("^RNA_", colnames(Xtr)))) )])
logmsg("(sanity) RNA-only top-500-variance baseline AUC: %.3f", bl_rna$auc)

# -----------------------------------------------------------------------------
# 9. VALIDATION C: stability across subsamples (train only)
# -----------------------------------------------------------------------------
logmsg("=== STEP 4c: stability (%d subsamples) ===", CONFIG$n_stab_runs)
stab_lists <- list()
for (b in seq_len(CONFIG$n_stab_runs)) {
  set.seed(1000 + b)
  bi <- unlist(lapply(split(seq_along(ytr), ytr),
    function(ii) sample(ii, size = floor(0.8 * length(ii)))))
  sb <- perm_assisted_grplasso_binary(Xtr[bi, ], ytr[bi], CONFIG$n_perm_sel,
    grp$groups, pen, CONFIG$sel_cutoff)
  stab_lists[[b]] <- sb$final_selection
  logmsg("  subsample %d: %d features", b, length(sb$final_selection))
}
pw <- combn(seq_along(stab_lists), 2, function(ii) jaccard(stab_lists[[ii[1]]], stab_lists[[ii[2]]]))
logmsg("Mean pairwise Jaccard across subsamples: %.3f", mean(pw))

# -----------------------------------------------------------------------------
# 10. SAVE EVERYTHING
# -----------------------------------------------------------------------------
saveRDS(list(config = CONFIG, n_train = length(ytr), n_test = length(yte),
  selected = sel$final_selection, frequencies = sel$feature_frequencies,
  overlap = overlap_res, test_auc = res_test$auc,
  stability_jaccard = mean(pw), stab_lists = stab_lists,
  probe2gene = probe2gene, universe_n = length(universe)),
  file.path(CONFIG$out_dir, "validation_results.rds"))
writeLines(capture.output(sessionInfo()),
  file.path(CONFIG$out_dir, "sessionInfo.txt"))
logmsg("=== DONE. Results in %s/ ===", CONFIG$out_dir)

# Multi-modal feature selection on TCGA-BRCA data with adaptive group lasso
# (inverse-log p-value penalties). Converted from the TCGA-Selection analysis notebook.

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

# --- Data Manipulation & Viz ---
library(tidyverse)
library(pheatmap)
library(knitr)
library(ggplot2)

# --- Modeling ---
library(glmnet)
library(grpreg)
library(pROC)

# --- TCGA Data Access ---
# If missing, install: BiocManager::install(c("curatedTCGAData", "TCGAutils", "MultiAssayExperiment"))
library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)

load_tcga_brca_data <- function() {
  cat("Downloading/Loading TCGA BRCA data... (This may take a moment)\n")
  
  # 1. Fetch Data
  brca <- curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "Methylation*"),
    version = "2.1.1", 
    dry.run = FALSE
  )
  
  # Helper to safe-extract assay to standard matrix
  get_matrix <- function(mae, pattern) {
    found <- grep(pattern, names(assays(mae)), value = TRUE)
    if(length(found) == 0) stop(paste("No assay found for pattern:", pattern))
    
    # Extract raw object
    raw_obj <- assays(mae)[[found[1]]]
    
    # FORCE conversion to standard numeric matrix immediately
    # This fixes issues with DelayedArray/HDF5Matrix
    mat <- as.matrix(raw_obj)
    
    if(!is.numeric(mat)) {
      cat(sprintf("Warning: Assay %s was not numeric. Attempting conversion...\n", pattern))
      mode(mat) <- "numeric"
    }
    return(mat)
  }
  
  cat("Extracting matrices...\n")
  rna <- get_matrix(brca, "RNASeq2GeneNorm")
  meth <- get_matrix(brca, "Methylation")
  
  cat(sprintf("Initial Dimensions -> RNA: %d x %d | Meth: %d x %d\n", 
              nrow(rna), ncol(rna), nrow(meth), ncol(meth)))

  # 2. Robust Sample Matching (Truncated Barcodes)
  # Truncate to first 12 chars (Patient ID only: "TCGA-XX-XXXX")
  # Using 15 chars (Sample) is safer, but 12 guarantees overlap if Vial codes mismatch widely
  clean_barcode <- function(b) substr(b, 1, 12)
  
  rna_pats <- clean_barcode(colnames(rna))
  meth_pats <- clean_barcode(colnames(meth))
  
  common_pats <- intersect(rna_pats, meth_pats)
  
  if(length(common_pats) < 10) {
    stop(sprintf("ERROR: Only %d common patients found. Matching failed.", length(common_pats)))
  }
  
  cat(sprintf("Found %d common patients.\n", length(common_pats)))
  
  # Subset matrices to common patients
  # Note: This logic takes the FIRST sample per patient found
  rna_idx <- match(common_pats, rna_pats)
  meth_idx <- match(common_pats, meth_pats)
  
  rna <- rna[, rna_idx]
  meth <- meth[, meth_idx]
  
  # Rename columns to ensure match
  colnames(rna) <- common_pats
  colnames(meth) <- common_pats
  
  # 3. Extract Outcome
  col_data <- colData(brca)
  col_pats <- clean_barcode(rownames(col_data))
  col_match_idx <- match(common_pats, col_pats)
  
  # Filter out matched metadata
  matched_col_data <- col_data[col_match_idx, ]
  
  if("vital_status" %in% names(matched_col_data)) {
    # Convert 'deceased'/'living' (or 1/0) to numeric binary
    y_raw <- matched_col_data$vital_status
    y <- ifelse(y_raw %in% c("deceased", "1", 1), 1, 0)
    # Handle NA outcomes
    y[is.na(y_raw)] <- NA
  } else {
    stop("Could not find 'vital_status' in clinical data.")
  }
  
  # Final remove NAs from y
  keep_mask <- !is.na(y)
  rna <- rna[, keep_mask]
  meth <- meth[, keep_mask]
  y <- y[keep_mask]
  
  cat(sprintf("Final Analysis Set: %d samples\n", length(y)))

  # 4. Debug & Filter Function
  clean_and_filter <- function(mat, label, n_top=1000) {
    cat(sprintf("Processing %s... Input dim: %d x %d\n", label, nrow(mat), ncol(mat)))
    
    # Check for all-NA rows
    na_counts <- rowSums(is.na(mat))
    if(any(na_counts == ncol(mat))) {
      cat(sprintf("  Removing %d features with 100%% NAs...\n", sum(na_counts == ncol(mat))))
      mat <- mat[na_counts < ncol(mat), ]
    }
    
    # Impute remaining NAs with Row Means
    if(any(is.na(mat))) {
      cat("  Imputing partial NAs with row means...\n")
      row_means <- rowMeans(mat, na.rm=TRUE)
      na_idx <- which(is.na(mat), arr.ind=TRUE)
      mat[na_idx] <- row_means[na_idx[,1]]
    }
    
    # Calculate Variance
    vars <- apply(mat, 1, var)
    
    # Diagnostic print
    cat(sprintf("  Variance Summary: Min=%.4f, Max=%.4f, NAs=%d\n", 
                min(vars, na.rm=TRUE), max(vars, na.rm=TRUE), sum(is.na(vars))))
    
    # Filter
    keep_var <- !is.na(vars) & vars > 1e-6 # Must have some minimal variance
    
    if(sum(keep_var) == 0) {
      stop(sprintf("%s: All features have near-zero variance! Check data scaling.", label))
    }
    
    mat <- mat[keep_var, ]
    vars <- vars[keep_var]
    
    # Top N
    real_n <- min(n_top, nrow(mat))
    top_idx <- order(vars, decreasing=TRUE)[1:real_n]
    
    cat(sprintf("  Kept top %d features.\n", real_n))
    return(mat[top_idx, ])
  }
  
  # Run Cleaning
  X_rna <- t(clean_and_filter(rna, "RNA", 1000))
  X_meth <- t(clean_and_filter(meth, "Methylation", 1000))
  
  colnames(X_rna) <- paste0("RNA_", colnames(X_rna))
  colnames(X_meth) <- paste0("METH_", colnames(X_meth))
  
  return(list(
    X_rna = X_rna,
    X_meth = X_meth,
    y = y
  ))
}

## 3. Methodology: Adaptive P-Value Weighting

### 3.1. Helper: Robust Weight Calculation
#' @title Calculate Weights (Robust Correlation w/ Inverse-Log Logic)
#' @description Computes p-values via permutation and converts them to penalty weights.
#' Uses the 'inverse_log' method to align with concentration inequality theory.
calculate_p_value_weights_robust <- function(X, y, groups, n_perm = 200) {
  
  unique_groups <- levels(groups)
  y_centered <- y - mean(y) 
  X_centered <- scale(X, center = TRUE, scale = FALSE) 
  
  # --- Step 1: Fast Score Calculation (Crossprod) ---
  calc_scores_fast <- function(y_vec, X_mat) {
    # Calculate correlation of ALL features with y at once
    cor_vec <- as.numeric(crossprod(X_mat, y_vec)) 
    
    # Aggregate by group (L2 norm of correlations)
    sapply(unique_groups, function(g) {
      idx <- which(groups == g)
      sqrt(sum(cor_vec[idx]^2))
    })
  }
  
  # --- Step 2: Real & Null Scores ---
  real_scores <- calc_scores_fast(y_centered, X_centered)
  
  null_scores <- matrix(NA, nrow = n_perm, ncol = length(unique_groups))
  for(k in 1:n_perm) {
    null_scores[k, ] <- calc_scores_fast(sample(y_centered), X_centered)
  }
  
  # --- Step 3: Empirical P-values ---
  # Add pseudo-count to avoid p=0
  p_values <- sapply(seq_along(unique_groups), function(i) {
    (sum(null_scores[, i] >= real_scores[i]) + 1) / (n_perm + 1)
  })
  names(p_values) <- unique_groups
  
  # --- Step 4: Inverse-Log Weighting (Concentration Inequality Theory) ---
  # Formula: w = 1 / (-log(p) + 1)
  # Low p (Signal) -> High Surprisal -> Low Weight (Low Penalty)
  surprisal <- -log(p_values)
  weights <- 1 / (surprisal + 1)
  
  return(list(weights = weights, p_values = p_values))
}

### 3.2. Helper: Group Definition (Unsupervised Clustering)
define_feature_groups <- function(X) {
  cat("Defining feature blocks via Hierarchical Clustering...\n")
  # Use Spearman to capture non-linear relationships
  cor_mat <- cor(X, method = "spearman")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hclust_obj <- hclust(dist_mat, method = "ward.D2")
  
  # Cut tree to create ~1 group per 10 features (average block size)
  n_groups <- max(1, round(ncol(X) / 10)) 
  groups <- cutree(hclust_obj, k = n_groups)
  
  return(as.factor(groups))
}

### 3.3. Main Pipeline: Permutation-Assisted CV Group Lasso
run_adaptive_lasso_pipeline <- function(X, y, n_perm = 200) {
  
  # 1. Pre-process
  X <- scale(X) # Z-score standardization is critical for Lasso
  
  # 2. Define Groups (Unsupervised)
  groups <- define_feature_groups(X)
  
  # 3. Calculate Adaptive Weights (Supervised - Permutation)
  cat("Calculating adaptive weights (Permutation)...\n")
  weight_res <- calculate_p_value_weights_robust(X, y, groups, n_perm = n_perm)
  
  # 4. Map Weights to Groups
  group_multipliers <- weight_res$weights[levels(groups)]
  
  # 5. Fit Group Lasso with Cross-Validation
  cat("Fitting CV Group Lasso...\n")
  # cv.grpreg automatically finds the best lambda
  cv_fit <- cv.grpreg(
    X, y, 
    group = groups, 
    penalty = "grLasso", 
    family = "binomial", 
    group.multiplier = group_multipliers,
    seed = 42
  )
  
  # 6. Extract Selected Features (Min MSE)
  coefs <- coef(cv_fit, s = "lambda.min")
  selected <- names(coefs)[coefs != 0]
  selected <- selected[selected != "(Intercept)"]
  
  # 7. Calculate AUC (Internal validation)
  preds <- predict(cv_fit, X, s="lambda.min", type="response")
  roc_obj <- roc(y, as.vector(preds), quiet=TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  return(list(
    selected_features = selected,
    auc = auc_val,
    weights = weight_res$weights,
    model = cv_fit
  ))
}

tcga_data = load_tcga_brca_data()

# Setup Data Variables from the loader list
y <- tcga_data$y
X_rna <- tcga_data$X_rna
X_meth <- tcga_data$X_meth

### 4.1. Early Fusion (Concatenation)
# We merge RNA and Methylation matrices and run the pipeline on the joint set.
X_early <- cbind(X_rna, X_meth)

cat("\n--- Running EARLY FUSION ---\n")
early_res <- run_adaptive_lasso_pipeline(X_early, y, n_perm = 200)

cat("Early Fusion Selected:", length(early_res$selected_features), "features.\n")
cat("Early Fusion AUC:", round(early_res$auc, 3), "\n")


### 4.2. Late Fusion (Ensemble)
# We run the pipeline separately on RNA and Methylation, then combine predictions.

cat("\n--- Running LATE FUSION (Modality A: RNA) ---\n")
late_rna_res <- run_adaptive_lasso_pipeline(X_rna, y, n_perm = 200)

cat("\n--- Running LATE FUSION (Modality B: Methylation) ---\n")
late_meth_res <- run_adaptive_lasso_pipeline(X_meth, y, n_perm = 200)

# Combine Predictions (Average Probability)
pred_rna <- predict(late_rna_res$model, scale(X_rna), s="lambda.min", type="response")
pred_meth <- predict(late_meth_res$model, scale(X_meth), s="lambda.min", type="response")

late_prob <- (pred_rna + pred_meth) / 2
late_auc <- as.numeric(auc(roc(y, as.vector(late_prob), quiet=TRUE)))

cat("\nLate Fusion AUC:", round(late_auc, 3), "\n")


### 4.3. Baseline Lasso (Standard)
# Standard Lasso without adaptive weights or grouping, for comparison.

cat("\n--- Running BASELINE LASSO ---\n")
# alpha=1 is standard Lasso
cv_lasso <- cv.glmnet(scale(X_early), y, family="binomial", alpha=1)

# Calculate AUC
lasso_prob <- predict(cv_lasso, scale(X_early), s="lambda.min", type="response")
lasso_auc <- as.numeric(auc(roc(y, as.vector(lasso_prob), quiet=TRUE)))

# --- FIX: Safe extraction of sparse coefficients ---
lasso_coefs <- coef(cv_lasso, s="lambda.min")

# Convert S4 sparse matrix index to standard integer index
non_zero_indices <- which(lasso_coefs != 0)
lasso_selected <- rownames(lasso_coefs)[non_zero_indices]

# Remove intercept from count
lasso_selected <- lasso_selected[lasso_selected != "(Intercept)"]
lasso_count <- length(lasso_selected)

cat("Baseline Lasso AUC:", round(lasso_auc, 3), "\n")
cat("Baseline Lasso Selected:", lasso_count, "features.\n")

## 5. Results & Visualization

# --- 1. Metric Comparison Table ---
results_df <- data.frame(
  Model = c("Early Fusion", "Late Fusion", "Baseline Lasso"),
  AUC = c(early_res$auc, late_auc, lasso_auc),
  Features_Selected = c(
    length(early_res$selected_features),
    length(late_rna_res$selected_features) + length(late_meth_res$selected_features),
    lasso_count
  )
)

print(kable(results_df, digits=3, caption="TCGA-BRCA Experiment Results"))

# --- 2. Weight Distribution Plot (Interpretation) ---
# Let's visualize the weights assigned in Early Fusion to prove the method worked.

weight_df <- data.frame(
  Block_ID = names(early_res$weights),
  Weight = as.numeric(early_res$weights)
)

# Identify if block is RNA or Methylation (Heuristic: check feature names in that block)
# We regenerate the groups one last time to map names
groups_early <- define_feature_groups(scale(X_early))
feature_names <- colnames(X_early)

block_types <- sapply(names(early_res$weights), function(grp_id) {
  # Get indices of features in this group
  idx <- which(groups_early == grp_id)
  feats <- feature_names[idx]
  
  if(all(grepl("RNA", feats))) return("RNA-Block")
  if(all(grepl("METH", feats))) return("Meth-Block")
  return("Mixed-Block")
})

weight_df$Type <- block_types

# Plot top 40 most significant blocks (Smallest Weight)
# Remember: Small Weight = High Significance
top_blocks <- weight_df %>% 
  arrange(Weight) %>% 
  head(40)

ggplot(top_blocks, aes(x = reorder(Block_ID, -Weight), y = Weight, fill = Type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("RNA-Block" = "steelblue", "Meth-Block" = "orange", "Mixed-Block" = "gray")) +
  labs(
    title = "Top 40 Most Significant Feature Blocks (TCGA-BRCA)",
    subtitle = "Inverse-Log Weights: Lower Bar = More Significant (Less Penalty)",
    y = "Adaptive Penalty Weight",
    x = "Feature Block ID"
  ) +
  theme_minimal()

# --- Data Manipulation & Viz ---
library(tidyverse)
library(pheatmap)
library(knitr)
library(ggplot2)

# --- Modeling ---
library(glmnet)
library(grpreg)
library(pROC)

# --- TCGA Data Access ---
library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)

load_tcga_brca_data <- function() {
  cat("Downloading/Loading TCGA BRCA data...\n")
  
  # 1. Fetch Data
  brca <- curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "Methylation*"),
    version = "2.1.1", 
    dry.run = FALSE
  )
  
  # Helper to safe-extract assay to standard matrix
  get_matrix <- function(mae, pattern) {
    found <- grep(pattern, names(assays(mae)), value = TRUE)
    if(length(found) == 0) stop(paste("No assay found for pattern:", pattern))
    mat <- as.matrix(assays(mae)[[found[1]]])
    if(!is.numeric(mat)) mode(mat) <- "numeric"
    return(mat)
  }
  
  rna <- get_matrix(brca, "RNASeq2GeneNorm")
  meth <- get_matrix(brca, "Methylation")
  
  # 2. Robust Sample Matching
  clean_barcode <- function(b) substr(b, 1, 12)
  rna_pats <- clean_barcode(colnames(rna))
  meth_pats <- clean_barcode(colnames(meth))
  common_pats <- intersect(rna_pats, meth_pats)
  
  cat(sprintf("Found %d common patients.\n", length(common_pats)))
  
  rna <- rna[, match(common_pats, rna_pats)]
  meth <- meth[, match(common_pats, meth_pats)]
  colnames(rna) <- common_pats
  colnames(meth) <- common_pats
  
  # 3. Extract Outcome
  col_data <- colData(brca)
  col_pats <- clean_barcode(rownames(col_data))
  matched_col_data <- col_data[match(common_pats, col_pats), ]
  y_raw <- matched_col_data$vital_status
  y <- ifelse(y_raw %in% c("deceased", "1", 1), 1, 0)
  
  keep_mask <- !is.na(y)
  rna <- rna[, keep_mask]
  meth <- meth[, keep_mask]
  y <- y[keep_mask]
  
  # 4. Filter High Variance Features (Top 1000)
  clean_and_filter <- function(mat, label, n_top=1000) {
    if(any(is.na(mat))) {
      row_means <- rowMeans(mat, na.rm=TRUE)
      na_idx <- which(is.na(mat), arr.ind=TRUE)
      mat[na_idx] <- row_means[na_idx[,1]]
    }
    vars <- apply(mat, 1, var)
    keep_var <- !is.na(vars) & vars > 1e-6
    mat <- mat[keep_var, ]
    vars <- vars[keep_var]
    
    real_n <- min(n_top, nrow(mat))
    top_idx <- order(vars, decreasing=TRUE)[1:real_n]
    return(mat[top_idx, ])
  }
  
  X_rna <- t(clean_and_filter(rna, "RNA", 1000))
  X_meth <- t(clean_and_filter(meth, "Methylation", 1000))
  
  colnames(X_rna) <- paste0("RNA_", colnames(X_rna))
  colnames(X_meth) <- paste0("METH_", colnames(X_meth))
  
  return(list(X_rna = X_rna, X_meth = X_meth, y = y))
}

# Execute Loading
tcga_data <- load_tcga_brca_data()
X_rna <- tcga_data$X_rna
X_meth <- tcga_data$X_meth
y <- tcga_data$y

# --- 1. Define Feature Groups (Hierarchical Clustering) ---
define_feature_groups <- function(X) {
  cor_mat <- cor(X, method = "spearman")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hclust_obj <- hclust(dist_mat, method = "ward.D2")
  # Dynamic group sizing (~10 features per group)
  n_groups <- max(1, round(ncol(X) / 10)) 
  groups <- cutree(hclust_obj, k = n_groups)
  return(as.factor(groups))
}

# --- 2. Calculate Global Adaptive Weights (Inverse-Log) ---
calculate_p_value_weights_robust <- function(X, y, groups, n_perm = 500) {
  unique_groups <- levels(groups)
  y_centered <- y - mean(y) 
  X_centered <- scale(X, center = TRUE, scale = FALSE) 
  
  # Fast Score Calculation (L2 norm of correlations)
  calc_scores_fast <- function(y_vec, X_mat) {
    cor_vec <- as.numeric(crossprod(X_mat, y_vec)) 
    sapply(unique_groups, function(g) {
      sqrt(sum(cor_vec[which(groups == g)]^2))
    })
  }
  
  real_scores <- calc_scores_fast(y_centered, X_centered)
  null_scores <- matrix(NA, nrow = n_perm, ncol = length(unique_groups))
  
  for(k in 1:n_perm) {
    null_scores[k, ] <- calc_scores_fast(sample(y_centered), X_centered)
  }
  
  p_values <- sapply(seq_along(unique_groups), function(i) {
    (sum(null_scores[, i] >= real_scores[i]) + 1) / (n_perm + 1)
  })
  names(p_values) <- unique_groups
  
  # Inverse-Log Weighting: w = 1 / (-log(p) + 1)
  # Low p -> High Surprisal -> Low Weight (Low Penalty)
  weights <- 1 / (-log(p_values) + 1)
  
  return(list(weights = weights))
}

# --- 3. Permutation-Assisted Group Lasso (Knockoff Logic) ---
perm_assisted_grplasso_binary <- function(X, y, n_perm, feature_groups, group_penalties, cutoff) {
  X <- as.matrix(X)
  p <- ncol(X)
  all_feature_names <- colnames(X)
  selected_features_all_perms <- c()
  
  g_int <- as.integer(feature_groups)
  knockoff_groups <- g_int + max(g_int)
  
  # Knockoffs inherit the EXACT SAME penalty as the original group
  full_group_penalties <- c(group_penalties, group_penalties)
  
  # We suppress the internal progress bar for cleaner output
  for (k in 1:n_perm) {
    X_knockoff <- X[sample(nrow(X)), ] # Permute rows to create knockoffs
    
    fit <- suppressWarnings(grpreg(
      cbind(X, X_knockoff), y,
      group = c(g_int, knockoff_groups),
      penalty = "grLasso",
      family = "binomial",
      group.multiplier = full_group_penalties,
      nlambda = 50,
      warn = FALSE
    ))
    
    betas <- as.matrix(fit$beta[-1,]) 
    
    # Selection logic: Find lambda where (Original - Knockoff) count is maximized
    n_original <- colSums(abs(betas[1:p, ]) > 0)
    n_knockoff <- colSums(abs(betas[(p + 1):(2 * p), ]) > 0)
    diff <- n_original - n_knockoff
    
    optimal_idx <- if (all(diff <= 0)) NULL else which.max(diff)

    if (!is.null(optimal_idx)) {
      selected_indices <- which(abs(betas[1:p, optimal_idx]) > 0)
      if (length(selected_indices) > 0) {
        selected_features_all_perms <- c(selected_features_all_perms, rownames(betas)[selected_indices])
      }
    }
  }

  freq_counts <- setNames(rep(0, p), all_feature_names)
  if (length(selected_features_all_perms) > 0) {
    counts_table <- table(selected_features_all_perms)
    freq_counts[names(counts_table)] <- counts_table
  }
  
  feature_frequencies <- freq_counts / n_perm
  final_selection <- names(feature_frequencies)[feature_frequencies >= cutoff]
  
  return(list(final_selection = final_selection))
}

# --- 4. Metric Helpers ---
calculate_jaccard <- function(set1, set2) {
  union_len <- length(union(set1, set2))
  if (union_len == 0) return(1)
  return(length(intersect(set1, set2)) / union_len)
}

calculate_subsequent_jaccard <- function(feature_list) {
  n <- length(feature_list)
  if (n < 2) return(NA)
  scores <- sapply(1:(n-1), function(i) calculate_jaccard(feature_list[[i]], feature_list[[i+1]]))
  return(mean(scores))
}

cat("--- PHASE A: Calculating Global Weights (n_perm=500) ---\n")

# 1. Early Fusion Setup (Joint Matrix)
X_early <- cbind(X_rna, X_meth)
groups_early <- define_feature_groups(X_early)
# Compute weights
w_res_early <- calculate_p_value_weights_robust(X_early, y, groups_early, n_perm = 500)
weights_early <- w_res_early$weights

# 2. Late Fusion Setup (Separate Matrices)
# RNA Leg
groups_rna <- define_feature_groups(X_rna)
weights_rna <- calculate_p_value_weights_robust(X_rna, y, groups_rna, n_perm = 500)$weights

# Methylation Leg
groups_meth <- define_feature_groups(X_meth)
weights_meth <- calculate_p_value_weights_robust(X_meth, y, groups_meth, n_perm = 500)$weights

cat("Global weights successfully calculated.\n")

# Experiment Parameters
N_RUNS <- 5         # Keep at 5 for testing speed, increase to 10+ for final
N_PERM_SELECT <- 50  # Permutations inside the selection algorithm

run_stability_experiment <- function(mode) {
  
  feature_lists <- list()
  auc_scores <- numeric(N_RUNS)
  
  cat(sprintf("\n--- Running Stability Analysis: %s ---\n", mode))
  # Initialize progress bar
  progress_bar <- txtProgressBar(min = 0, max = N_RUNS, style = 3)
  
  for(i in 1:N_RUNS) {
    set.seed(i)
    # 80% Subsample
    train_idx <- sample(seq_len(nrow(X_early)), floor(0.8 * nrow(X_early)))
    y_sub <- y[train_idx]
    
    # --- METHOD LOGIC ---
    if(mode == "Early Fusion") {
      X_sub <- X_early[train_idx, ]
      
      # Use Global Groups & Weights
      res <- perm_assisted_grplasso_binary(
        X = X_sub, y = y_sub,
        n_perm = N_PERM_SELECT,
        feature_groups = groups_early,
        group_penalties = weights_early[levels(groups_early)],
        cutoff = 0.5
      )
      
      sel <- res$final_selection
      feature_lists[[i]] <- sel
      
      # Evaluate AUC (Ridge on selected features)
      if(length(sel) > 1) {
        fit_eval <- cv.glmnet(scale(X_sub[, sel]), y_sub, family="binomial", alpha=0)
        p <- predict(fit_eval, scale(X_sub[, sel]), s="lambda.min", type="response")
        auc_scores[i] <- as.numeric(auc(roc(y_sub, as.vector(p), quiet=TRUE)))
      } else { auc_scores[i] <- 0.5 }
      
    } else if (mode == "Late Fusion") {
      # RNA Leg
      res_rna <- perm_assisted_grplasso_binary(
        X = X_rna[train_idx, ], y = y_sub,
        n_perm = N_PERM_SELECT,
        feature_groups = groups_rna,
        group_penalties = weights_rna[levels(groups_rna)],
        cutoff = 0.5
      )
      
      # Meth Leg
      res_meth <- perm_assisted_grplasso_binary(
        X = X_meth[train_idx, ], y = y_sub,
        n_perm = N_PERM_SELECT,
        feature_groups = groups_meth,
        group_penalties = weights_meth[levels(groups_meth)],
        cutoff = 0.5
      )
      
      # Union Features for Stability
      sel <- unique(c(res_rna$final_selection, res_meth$final_selection))
      feature_lists[[i]] <- sel
      
      # AUC: Ensemble Average
      auc_rna <- 0.5; auc_meth <- 0.5
      
      if(length(res_rna$final_selection)>0) {
         fA <- cv.glmnet(scale(X_rna[train_idx, res_rna$final_selection]), y_sub, family="binomial", alpha=0)
         probs_rna <- predict(fA, scale(X_rna[train_idx, res_rna$final_selection]), s="lambda.min", type="response")
         auc_rna <- as.numeric(auc(roc(y_sub, as.vector(probs_rna), quiet=TRUE)))
      }
      
      if(length(res_meth$final_selection)>0) {
         fB <- cv.glmnet(scale(X_meth[train_idx, res_meth$final_selection]), y_sub, family="binomial", alpha=0)
         # FIXED: Renamed variable from 'pb' to 'probs_meth' to avoid overwriting progress bar
         probs_meth <- predict(fB, scale(X_meth[train_idx, res_meth$final_selection]), s="lambda.min", type="response")
         auc_meth <- as.numeric(auc(roc(y_sub, as.vector(probs_meth), quiet=TRUE)))
      }
      
      auc_scores[i] <- (auc_rna + auc_meth) / 2
      
    } else if (mode == "Baseline Lasso") {
      # Standard Lasso
      fit <- cv.glmnet(scale(X_early[train_idx, ]), y_sub, family="binomial", alpha=1)
      
      coefs <- coef(fit, s="lambda.min")
      sel <- rownames(coefs)[which(coefs != 0)]
      sel <- sel[sel != "(Intercept)"]
      feature_lists[[i]] <- sel
      
      p <- predict(fit, scale(X_early[train_idx, ]), s="lambda.min", type="response")
      auc_scores[i] <- as.numeric(auc(roc(y_sub, as.vector(p), quiet=TRUE)))
    }
    
    setTxtProgressBar(progress_bar, i)
  }
  close(progress_bar)
  
  return(list(
    mode = mode,
    stability = calculate_subsequent_jaccard(feature_lists),
    mean_auc = mean(auc_scores),
    sd_auc = sd(auc_scores),
    avg_features = mean(sapply(feature_lists, length)),
    all_features = feature_lists
  ))
}

# --- Execution ---
res_early <- run_stability_experiment("Early Fusion")

res_late <- run_stability_experiment("Late Fusion")

res_baseline <- run_stability_experiment("Baseline Lasso")

# 1. Summary Table
results_df <- data.frame(
  Method = c("Early Fusion (Adaptive)", "Late Fusion (Adaptive)", "Baseline Lasso"),
  Stability_Jaccard = c(res_early$stability, res_late$stability, res_baseline$stability),
  Mean_AUC = c(res_early$mean_auc, res_late$mean_auc, res_baseline$mean_auc),
  Avg_Features = c(res_early$avg_features, res_late$avg_features, res_baseline$avg_features)
)

print(kable(results_df, digits=3, caption = "Performance & Stability Metrics (TCGA-BRCA)"))

# 2. Visualization
plot_df <- results_df %>%
  pivot_longer(cols = c(Stability_Jaccard, Mean_AUC), names_to = "Metric", values_to = "Value")

ggplot(plot_df, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.9) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("gray40", "forestgreen", "steelblue")) +
  labs(
    title = "Performance vs. Stability Trade-off",
    subtitle = "Comparing Adaptive Fusion Methods against Baseline",
    y = "Score"
  ) +
  theme(legend.position = "bottom")

# 3. Save Results
write.csv(results_df, "TCGA_Final_Results.csv", row.names = FALSE)

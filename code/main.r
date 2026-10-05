
# ARCHIVAL: thesis as written (May 2026). This script contains the SKAT-based
# adaptive weighting issue documented in REPORT.md (SKAT is a rare-variant
# test for genotype counts; on continuous expression it returns p = 1 or
# crashes, so the "adaptive" arms effectively ran with size-only weights).
# Do not use for new work. For the frozen scripts behind the reported
# results, see code/frozen_tcga/ (manifest: code/frozen_tcga/MANIFEST.md).

library(tidyverse)
library(pheatmap)
library(glmnet)
library(grpreg)
library(pROC)
library(SKAT)
library(Matrix)
library(doParallel)
library(foreach)

num_cores <- detectCores(logical = TRUE) - 1
registerDoParallel(cores = num_cores)
# ==============================================================================
# 1. Data Simulation 
# ==============================================================================
simulation_data <- function(EF = c(-1, 1.5, -2, 1.2), 
                            n_ss = 250, 
                            n_x = rep(10, 3), 
                            n_noise = 200, 
                            SD_v = 0.5) {
  
  num_causal <- length(EF) - 1
  x <- matrix(runif(n_ss * num_causal), nrow = n_ss, ncol = num_causal)
  colnames(x) <- paste0("x", 1:num_causal)
  
  v <- list()
  for (i in 1:num_causal) {
    v_i <- matrix(NA, nrow = n_ss, ncol = n_x[i])
    for (j in 1:n_x[i]) {
      noise <- rnorm(n = n_ss, mean = 0, sd = SD_v)
      v_i[, j] <- x[, i] + (0.01 + 0.5 * (j - 1) / (n_x[i] - 1)) * noise
    }
    colnames(v_i) <- paste0("v", i, "_", 1:n_x[i])
    v[[i]] <- v_i
  }
  v_df <- do.call(cbind, v)
  
  w <- matrix(runif(n_ss * n_noise), nrow = n_ss, ncol = n_noise)
  colnames(w) <- paste0("w", 1:n_noise)
  
  if (all(EF[-1] == 0)) {
    prob <- rep(0.5, n_ss)
  } else {
    linear_predictor <- EF[1] + x %*% EF[-1]
    prob <- 1 / (1 + exp(-linear_predictor))
  }
  y <- rbinom(n_ss, 1, prob)
  
  sim_df <- as.data.frame(cbind(v_df, w))
  sim_df$y <- y
  sim_df <- sim_df[, c("y", colnames(v_df), colnames(w))]
  return(sim_df)
}

# ==============================================================================
# 2. Grouping and Penalty Calculations
# ==============================================================================
define_groups <- function(feature_df, linkage_method = "ward.D2", expected_cluster_size = 10) {
  cor_matrix <- cor(feature_df, method = "spearman")
  dist_matrix <- as.dist(1 - abs(cor_matrix))
  hclust_obj <- hclust(dist_matrix, method = linkage_method)
  num_groups <- max(1, round(ncol(feature_df) / expected_cluster_size))
  feature_groups_vec <- cutree(hclust_obj, k = num_groups)
  return(list(groups = feature_groups_vec, hclust = hclust_obj, cor = cor_matrix))
}

calculate_correlation_penalties <- function(feature_df, groups, cor_matrix) {
  unique_groups <- unique(groups)
  group_info <- data.frame(feature = names(groups), group = groups)
  mean_cor_per_group <- sapply(unique_groups, function(g) {
    members <- group_info$feature[group_info$group == g]
    if (length(members) < 2) return(0) 
    sub_cor_matrix <- cor_matrix[members, members]
    mean(abs(sub_cor_matrix[upper.tri(sub_cor_matrix)]))
  })
  group_summary <- data.frame(group = unique_groups, mean_cor = mean_cor_per_group) %>%
    left_join(group_info %>% dplyr::count(group, name = "size"), by = "group") %>%
    mutate(importance_rank = rank(-mean_cor, ties.method = "min"),
           penalty_raw = importance_rank * sqrt(size))
  penalties <- group_summary$penalty_raw / mean(group_summary$penalty_raw)
  names(penalties) <- group_summary$group
  return(penalties[as.character(sort(unique(groups)))])
}

#  w = 1 / (S + 1)
calculate_skat_penalties <- function(X, y, groups, epsilon = 0.1) {
  obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
  unique_groups <- sort(unique(groups))
  p_values <- numeric(length(unique_groups))
  names(p_values) <- unique_groups
  X_mat <- as.matrix(X)
  
  for(g in unique_groups) {
    g_features <- names(groups)[groups == g]
    Z <- X_mat[, g_features, drop=FALSE]
    tryCatch({
      suppressWarnings({
        skat_res <- SKATBinary(Z, obj, kernel = "Linear.Weighted")
      })
      p_values[g] <- skat_res$p.value
    }, error = function(e) { p_values[g] <<- 1 })
  }
  
  p_values <- pmax(pmin(p_values, 1 - 1e-10), 1e-10)
  neg_log_p <- -log(p_values)
  raw_penalty <- 1 / (neg_log_p + epsilon)
  final_penalty <- raw_penalty / mean(raw_penalty)  # normalize
  return(final_penalty)
}

# ==============================================================================
# 3. Parallelized Permutation-Assisted Group Lasso Engine 
# ==============================================================================
perm_assisted_grplasso_binary <- function(X, y, n_perm, feature_groups, group_penalties, cutoff) {
  X <- as.matrix(X)
  p <- ncol(X)
  all_feature_names <- colnames(X)
  knockoff_groups <- feature_groups + max(feature_groups)
  full_group_penalties <- c(group_penalties, group_penalties) 
  
  selected_features_all_perms <- foreach(k = 1:n_perm, .combine = c, .packages = c("grpreg")) %dopar% {
    X_knockoff <- X[sample(nrow(X)), ]
    fit <- grpreg(cbind(X, X_knockoff), y, group = c(feature_groups, knockoff_groups),
                  penalty = "grLasso", family = "binomial", 
                  group.multiplier = full_group_penalties, nlambda = 50, lambda.min = 0.05)
    betas <- as.matrix(fit$beta[-1,]) 
    n_original <- colSums(abs(betas[1:p, ]) > 0)
    n_knockoff <- colSums(abs(betas[(p + 1):(2 * p), ]) > 0)
    diff <- n_original - n_knockoff
    if (all(diff <= 0)) return(NULL)
    optimal_idx <- which.max(diff)
    selected_indices <- which(abs(betas[1:p, optimal_idx]) > 0)
    if (length(selected_indices) > 0) return(rownames(betas)[selected_indices])
    return(NULL)
  }
  
  freq_counts <- setNames(rep(0, p), all_feature_names)
  if (!is.null(selected_features_all_perms) && length(selected_features_all_perms) > 0) {
    counts_table <- table(selected_features_all_perms)
    freq_counts[names(counts_table)] <- counts_table
  }
  
  feature_frequencies <- freq_counts / n_perm
  final_selection <- names(feature_frequencies)[feature_frequencies >= cutoff]
  return(list(final_selection = final_selection, feature_frequencies = feature_frequencies))
}

# ==============================================================================
# 4. Evaluation Helpers
# ==============================================================================
calculate_metrics <- function(selected, true_features) {
  tp <- sum(selected %in% true_features)
  fp <- length(selected) - tp
  fn <- length(true_features) - tp
  
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  sensitivity <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1_score <- ifelse((precision + sensitivity) == 0, 0, 2 * (precision * sensitivity) / (precision + sensitivity))
  fdr <- ifelse((tp + fp) == 0, 0, fp / (tp + fp))
  
  n_v1 <- sum(grepl("^v1_", true_features))
  n_v2 <- sum(grepl("^v2_", true_features))
  n_v3 <- sum(grepl("^v3_", true_features))
  
  power_v1 <- ifelse(n_v1 == 0, 0, sum(grepl("^v1_", selected)) / n_v1)
  power_v2 <- ifelse(n_v2 == 0, 0, sum(grepl("^v2_", selected)) / n_v2)
  power_v3 <- ifelse(n_v3 == 0, 0, sum(grepl("^v3_", selected)) / n_v3)
  
  return(data.frame(F1_Score = f1_score, Sensitivity = sensitivity, FDR = fdr, Num_Selected = length(selected),
                    Power_B1 = power_v1, Power_B2 = power_v2, Power_B3 = power_v3))
}

calculate_pairwise_jaccard <- function(feature_list) {
  n <- length(feature_list)
  if (n < 2) return(NA)
  jaccard_scores <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      set1 <- feature_list[[i]]; set2 <- feature_list[[j]]
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      jaccard_scores <- c(jaccard_scores, ifelse(union == 0, 1, intersection / union))
    }
  }
  return(mean(jaccard_scores))
}

calc_auc <- function(X_mat, sel_vars, y_target) {
  if(length(sel_vars) == 0) return(0.5)
  if(length(sel_vars) == 1) {
    df_single <- data.frame(y = y_target, x = X_mat[, sel_vars])
    fit_glm <- glm(y ~ x, data = df_single, family = "binomial")
    preds <- predict(fit_glm, type = "response")
    return(as.numeric(auc(roc(y_target, as.vector(preds), quiet=TRUE))))
  }
  fit_cv <- cv.glmnet(X_mat[, sel_vars, drop=FALSE], y_target, family="binomial")
  preds <- predict(fit_cv, newx=X_mat[, sel_vars, drop=FALSE], s="lambda.min", type="response")
  return(as.numeric(auc(roc(y_target, as.vector(preds), quiet=TRUE))))
}

get_probs <- function(X_mat, sel_vars, y_target, n) {
  if(length(sel_vars) == 0) return(rep(0.5, n))
  if(length(sel_vars) == 1) {
    df_single <- data.frame(y = y_target, x = X_mat[, sel_vars])
    fit_glm <- glm(y ~ x, data = df_single, family = "binomial")
    return(predict(fit_glm, type = "response"))
  }
  fit <- cv.glmnet(X_mat[, sel_vars, drop=FALSE], y_target, family="binomial")
  return(predict(fit, newx=X_mat[, sel_vars, drop=FALSE], s="lambda.min", type="response"))
}

# ==============================================================================
# 5. Main Simulation Loop
# ==============================================================================
N_FEATURES_PER_BLOCK <- rep(10, 3) 
BASE_INTERCEPT <- -1
BASE_COEFFICIENTS <- c(1.5, -2, 1.2)

param_n_samples <- c(100, 200)
param_effect_multipliers <- c(1.0, 1.5, 2.0)
param_n_noise <- c(200, 300, 400) 

N_RUNS <- 50           
N_PERMUTATIONS <- 300
SELECTION_CUTOFF <- 0.5 

all_experiment_metrics <- list()
run_counter <- 1

cat("--- STARTING CORRECTED 5-MODEL SIMULATION ---\n")

for (n_val in param_n_samples) {
  for (ef_mult in param_effect_multipliers) {
    for (n_noise_val in param_n_noise) {
      
      cat(sprintf("\n*** SCENARIO %d: n=%d, effect=%.1fx, noise=%d ***\n", run_counter, n_val, ef_mult, n_noise_val))
      TRUE_COEFFICIENTS <- c(BASE_INTERCEPT, BASE_COEFFICIENTS * ef_mult)
      run_metrics_list <- list()
      
      run_features <- list(early_std=list(), early_adp=list(), late_std=list(), late_adp=list(), lasso=list())
      
      for (i in 1:N_RUNS) {
        set.seed(i)
        
        full_data <- simulation_data(EF = TRUE_COEFFICIENTS, n_ss = n_val, n_x = N_FEATURES_PER_BLOCK, n_noise = n_noise_val)
        y <- full_data$y
        
        true_signal_features <- colnames(full_data)[grepl("^v[1-3]_", colnames(full_data))]
        modA_features <- c(colnames(full_data)[grepl("^v[1-2]_", colnames(full_data))], 
                           colnames(full_data)[grepl("^w", colnames(full_data))][1:floor(n_noise_val/2)])
        modB_features <- c(colnames(full_data)[grepl("^v3_", colnames(full_data))], 
                           colnames(full_data)[grepl("^w", colnames(full_data))][(floor(n_noise_val/2)+1):n_noise_val])
        
        X_all <- full_data[, c(modA_features, modB_features)]
        X_all_scaled <- scale(as.matrix(X_all))
        X_modA_scaled <- scale(as.matrix(full_data[, modA_features]))
        X_modB_scaled <- scale(as.matrix(full_data[, modB_features]))
        
        grp_struct_early <- define_groups(X_all, expected_cluster_size = 10)
        pen_std_early <- calculate_correlation_penalties(X_all, grp_struct_early$groups, grp_struct_early$cor)
        pen_adp_early <- calculate_skat_penalties(X_all, y, grp_struct_early$groups)
        
        grp_struct_A <- define_groups(full_data[, modA_features], expected_cluster_size = 10)
        pen_std_A <- calculate_correlation_penalties(full_data[, modA_features], grp_struct_A$groups, grp_struct_A$cor)
        pen_adp_A <- calculate_skat_penalties(full_data[, modA_features], y, grp_struct_A$groups)
        
        grp_struct_B <- define_groups(full_data[, modB_features], expected_cluster_size = 10)
        pen_std_B <- calculate_correlation_penalties(full_data[, modB_features], grp_struct_B$groups, grp_struct_B$cor)
        pen_adp_B <- calculate_skat_penalties(full_data[, modB_features], y, grp_struct_B$groups)
        
        res_early_std <- perm_assisted_grplasso_binary(X_all_scaled, y, N_PERMUTATIONS, grp_struct_early$groups, pen_std_early, SELECTION_CUTOFF)
        res_early_adp <- perm_assisted_grplasso_binary(X_all_scaled, y, N_PERMUTATIONS, grp_struct_early$groups, pen_adp_early, SELECTION_CUTOFF)
        
        res_late_std_A <- perm_assisted_grplasso_binary(X_modA_scaled, y, N_PERMUTATIONS, grp_struct_A$groups, pen_std_A, SELECTION_CUTOFF)
        res_late_std_B <- perm_assisted_grplasso_binary(X_modB_scaled, y, N_PERMUTATIONS, grp_struct_B$groups, pen_std_B, SELECTION_CUTOFF)
        sel_late_std <- unique(c(res_late_std_A$final_selection, res_late_std_B$final_selection))
        
        res_late_adp_A <- perm_assisted_grplasso_binary(X_modA_scaled, y, N_PERMUTATIONS, grp_struct_A$groups, pen_adp_A, SELECTION_CUTOFF)
        res_late_adp_B <- perm_assisted_grplasso_binary(X_modB_scaled, y, N_PERMUTATIONS, grp_struct_B$groups, pen_adp_B, SELECTION_CUTOFF)
        sel_late_adp <- unique(c(res_late_adp_A$final_selection, res_late_adp_B$final_selection))
        
        cv_fit_lasso <- cv.glmnet(X_all_scaled, y, family = "binomial", alpha = 1)
        sel_lasso <- rownames(coef(cv_fit_lasso, s = "lambda.min"))[which(coef(cv_fit_lasso, s = "lambda.min") != 0)][-1]
        
        # Store for stability
        run_features$early_std[[i]] <- res_early_std$final_selection
        run_features$early_adp[[i]] <- res_early_adp$final_selection
        run_features$late_std[[i]] <- sel_late_std
        run_features$late_adp[[i]] <- sel_late_adp
        run_features$lasso[[i]] <- sel_lasso
        
        # AUC calculations
        auc_early_std <- calc_auc(X_all_scaled, res_early_std$final_selection, y)
        auc_early_adp <- calc_auc(X_all_scaled, res_early_adp$final_selection, y)
        auc_lasso     <- calc_auc(X_all_scaled, sel_lasso, y)
        
        auc_wA_std <- calc_auc(X_modA_scaled, res_late_std_A$final_selection, y)
        auc_wB_std <- calc_auc(X_modB_scaled, res_late_std_B$final_selection, y)
        wA_std <- auc_wA_std / (auc_wA_std + auc_wB_std + 1e-6); wB_std <- auc_wB_std / (auc_wA_std + auc_wB_std + 1e-6)
        probs_std_A <- get_probs(X_modA_scaled, res_late_std_A$final_selection, y, n_val)
        probs_std_B <- get_probs(X_modB_scaled, res_late_std_B$final_selection, y, n_val)
        auc_late_std <- as.numeric(auc(roc(y, as.vector((wA_std*probs_std_A) + (wB_std*probs_std_B)), quiet=TRUE)))
        
        auc_wA_adp <- calc_auc(X_modA_scaled, res_late_adp_A$final_selection, y)
        auc_wB_adp <- calc_auc(X_modB_scaled, res_late_adp_B$final_selection, y)
        wA_adp <- auc_wA_adp / (auc_wA_adp + auc_wB_adp + 1e-6); wB_adp <- auc_wB_adp / (auc_wA_adp + auc_wB_adp + 1e-6)
        probs_adp_A <- get_probs(X_modA_scaled, res_late_adp_A$final_selection, y, n_val)
        probs_adp_B <- get_probs(X_modB_scaled, res_late_adp_B$final_selection, y, n_val)
        auc_late_adp <- as.numeric(auc(roc(y, as.vector((wA_adp*probs_adp_A) + (wB_adp*probs_adp_B)), quiet=TRUE)))
        
        run_metrics_list[[i]] <- bind_rows(
          data.frame(Model = "Early Fusion (Std)", AUC = auc_early_std, calculate_metrics(res_early_std$final_selection, true_signal_features)),
          data.frame(Model = "Early Fusion (Adap)", AUC = auc_early_adp, calculate_metrics(res_early_adp$final_selection, true_signal_features)),
          data.frame(Model = "Late Fusion (Std)", AUC = auc_late_std, calculate_metrics(sel_late_std, true_signal_features)),
          data.frame(Model = "Late Fusion (Adap)", AUC = auc_late_adp, calculate_metrics(sel_late_adp, true_signal_features)),
          data.frame(Model = "Baseline LASSO", AUC = auc_lasso, calculate_metrics(sel_lasso, true_signal_features))
        )
      } 
      
      # Summary
      summary_metrics <- do.call(rbind, run_metrics_list) %>%
        group_by(Model) %>%
        summarise(across(c(AUC, F1_Score, Sensitivity, FDR, Num_Selected, Power_B1, Power_B2, Power_B3), list(Mean = mean, SD = sd))) %>%
        mutate(n_samples = n_val, effect_multiplier = ef_mult, n_noise = n_noise_val)
      
      stabilities <- c(
        "Early Fusion (Std)" = calculate_pairwise_jaccard(run_features$early_std),
        "Early Fusion (Adap)" = calculate_pairwise_jaccard(run_features$early_adp),
        "Late Fusion (Std)" = calculate_pairwise_jaccard(run_features$late_std),
        "Late Fusion (Adap)" = calculate_pairwise_jaccard(run_features$late_adp),
        "Baseline LASSO" = calculate_pairwise_jaccard(run_features$lasso)
      )
      summary_metrics$Stability <- stabilities[summary_metrics$Model]
      
      all_experiment_metrics[[run_counter]] <- summary_metrics
      run_counter <- run_counter + 1
    }
  }
}

# ==============================================================================
# 6. Final Outputs
# ==============================================================================
final_metrics_df <- do.call(rbind, all_experiment_metrics)
stopImplicitCluster()

print(final_metrics_df)
write_csv(final_metrics_df, "Thesis_Simulation_Corrected_Results.csv")

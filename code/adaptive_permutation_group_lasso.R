# Multi-modal feature selection with adaptive permutation-assisted group lasso.
# Simulation study: correlated features, binary outcome, decoy-based selection.
# Converted from the Fusion_work analysis notebook.

# Data manipulation and visualization
library(tidyverse)
library(pheatmap)
library(knitr)

# Modeling and Evaluation
library(glmnet)
library(grpreg)
library(pROC)

# Adaptive Weights Calculation
library(SKAT) 
library(Matrix)

#' @title Simulate Data with Correlated Features and a Binary Outcome
simulation_data <- function(EF = c(-1, 1.5, -2, 1.2, -1.8, 1.0), 
                              n_ss = 250, 
                              n_x = rep(10, 5), 
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
  
  linear_predictor <- EF[1] + x %*% EF[-1]
  prob <- 1 / (1 + exp(-linear_predictor))
  y <- rbinom(n_ss, 1, prob)
  
  sim_df <- as.data.frame(cbind(v_df, w))
  sim_df$y <- y
  sim_df <- sim_df[, c("y", colnames(v_df), colnames(w))]
  
  return(sim_df)
}

#' @title Define Groups via Clustering
define_groups <- function(feature_df) {
  cor_matrix <- cor(feature_df, method = "spearman")
  dist_matrix <- as.dist(1 - abs(cor_matrix))
  hclust_obj <- hclust(dist_matrix, method = "ward.D2")
  
  num_groups <- max(1, round(ncol(feature_df) / 8))
  feature_groups_vec <- cutree(hclust_obj, k = num_groups)
  
  return(list(groups = feature_groups_vec, hclust = hclust_obj, cor = cor_matrix))
}

#' @title Calculate Penalties based on Correlation (Standard)
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
    mutate(
      importance_rank = rank(-mean_cor, ties.method = "min"),
      penalty_raw = importance_rank * sqrt(size)
    )
  
  penalties <- group_summary$penalty_raw / mean(group_summary$penalty_raw)
  names(penalties) <- group_summary$group
  
  return(penalties[as.character(sort(unique(groups)))])
}

#' @title Calculate Penalties based on SKAT P-values (Adaptive)
#' @description Calculates adaptive penalties: 1 / (-log(p) + 1)
calculate_skat_penalties <- function(X, y, groups) {
  
  # 1. Fit Null Model
  obj <- SKAT_Null_Model(y ~ 1, out_type = "D")
  
  unique_groups <- sort(unique(groups))
  p_values <- numeric(length(unique_groups))
  names(p_values) <- unique_groups
  group_sizes <- numeric(length(unique_groups))
  
  X_mat <- as.matrix(X)
  
  for(g in unique_groups) {
    g_features <- names(groups)[groups == g]
    Z <- X_mat[, g_features, drop=FALSE]
    group_sizes[g] <- length(g_features)
    
    tryCatch({
      suppressWarnings({
        skat_res <- SKATBinary(Z, obj, kernel = "Linear.Weighted")
      })
      p_values[g] <- skat_res$p.value
    }, error = function(e) {
      p_values[g] <<- 1 
    })
  }
  
  # 2. Calculate Penalties
  epsilon <- 1e-6
  neg_log_p <- -log(p_values + epsilon)
  
  # Prevent negatives if p > 1 (rare edge case)
  neg_log_p[neg_log_p < 0] <- 0
  
  # --- MODIFIED FORMULA ---
  # Penalty ~ sqrt(size) * [ 1 / (-log(p) + 1) ]
  raw_penalty <- sqrt(group_sizes) * (1 / (neg_log_p + 1))
  
  # 3. NORMALIZE
  final_penalty <- raw_penalty / mean(raw_penalty)
  print(final_penalty)
  
  return(final_penalty)
}

#' @title Permutation-Assisted Group Lasso
perm_assisted_grplasso_binary <- function(X, y, n_perm, feature_groups, group_penalties, cutoff) {
  X <- as.matrix(X)
  p <- ncol(X)
  all_feature_names <- colnames(X)
  selected_features_all_perms <- c()
  
  knockoff_groups <- feature_groups + max(feature_groups)
  full_group_penalties <- c(group_penalties, group_penalties) 
  
  for (k in 1:n_perm) {
    X_knockoff <- X[sample(nrow(X)), ]
    
    fit <- grpreg(
      cbind(X, X_knockoff), y,
      group = c(feature_groups, knockoff_groups),
      penalty = "grLasso",
      family = "binomial",
      group.multiplier = full_group_penalties,
      nlambda = 50,
      lambda.min = 0.05
    )
    
    betas <- as.matrix(fit$beta[-1,]) 
    
    n_original <- colSums(abs(betas[1:p, ]) > 0)
    n_knockoff <- colSums(abs(betas[(p + 1):(2 * p), ]) > 0)
    diff <- n_original - n_knockoff
    
    if (all(diff <= 0)) {
        optimal_idx <- NULL
    } else {
        optimal_idx <- which.max(diff)
    }

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
  
  return(list(final_selection = final_selection, feature_frequencies = feature_frequencies))
}

calculate_metrics <- function(selected, true) {
  if (length(selected) == 0) {
    tp <- 0; fp <- 0
  } else {
    tp <- sum(selected %in% true)
    fp <- length(selected) - tp
  }
  fn <- length(true) - tp
  
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  sensitivity <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1_score <- ifelse((precision + sensitivity) == 0, 0, 
                     2 * (precision * sensitivity) / (precision + sensitivity))
  fdr <- ifelse((tp + fp) == 0, 0, fp / (tp + fp))
  
  return(data.frame(F1_Score = f1_score, Sensitivity = sensitivity, FDR = fdr, Num_Selected = length(selected)))
}

calculate_subsequent_jaccard <- function(feature_list) {
  n <- length(feature_list)
  if (n < 2) return(NA)
  jaccard_scores <- numeric(n - 1)
  for (i in 1:(n - 1)) {
    set1 <- feature_list[[i]]; set2 <- feature_list[[i+1]]
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    jaccard_scores[i] <- ifelse(union == 0, 1, intersection / union)
  }
  return(mean(jaccard_scores))
}

# --- MAIN SIMULATION LOOP ---

N_FEATURES_PER_BLOCK <- rep(10, 3) 
BASE_INTERCEPT <- -1
BASE_COEFFICIENTS <- c(1.5, -2, 1.2)

# Parameters
param_n_samples <- c(100)
param_effect_multipliers <- c(1)
param_n_noise <- c(200) 

N_RUNS <- 5           
N_PERMUTATIONS <- 25  
SELECTION_CUTOFF <- 0.5 

all_experiment_metrics <- list()
all_experiment_stability <- list()
all_experiment_causal_power <- list()
run_counter <- 1

cat("--- STARTING 5-METHOD SIMULATION ---\n")

for (n_val in param_n_samples) {
  for (ef_mult in param_effect_multipliers) {
    for (n_noise_val in param_n_noise) {
      
      cat(sprintf("\n*** SCENARIO %d: n=%d, effect=%.1fx, noise=%d ***\n",
                  run_counter, n_val, ef_mult, n_noise_val))
      
      TRUE_COEFFICIENTS <- c(BASE_INTERCEPT, BASE_COEFFICIENTS * ef_mult)
      
      run_metrics_list <- list()
      run_features <- list(early_std=list(), early_adp=list(), late_std=list(), late_adp=list(), lasso=list())
      
      for (i in 1:N_RUNS) {
        cat(sprintf("  Run %d / %d\n", i, N_RUNS))
        set.seed(i)
        
        # 1. Data Generation
        full_data <- simulation_data(EF = TRUE_COEFFICIENTS, n_ss = n_val,
                                     n_x = N_FEATURES_PER_BLOCK, n_noise = n_noise_val)
        y <- full_data$y
        
        # 2. Define Features
        true_signal_features <- colnames(full_data)[grepl("^v[1-3]_", colnames(full_data))]
        modA_features <- c(colnames(full_data)[grepl("^v[1-2]_", colnames(full_data))],
                           colnames(full_data)[grepl("^w", colnames(full_data))][1:floor(n_noise_val/2)])
        modB_features <- c(colnames(full_data)[grepl("^v3_", colnames(full_data))],
                           colnames(full_data)[grepl("^w", colnames(full_data))][(floor(n_noise_val/2)+1):n_noise_val])
        
        X_all <- full_data[, c(modA_features, modB_features)]
        X_all_scaled <- scale(as.matrix(X_all))
        X_modA_scaled <- scale(as.matrix(full_data[, modA_features]))
        X_modB_scaled <- scale(as.matrix(full_data[, modB_features]))

        # --- PRE-CALCULATION OF GROUPS AND WEIGHTS ---
        
        # Early Fusion
        grp_struct_early <- define_groups(X_all)
        pen_std_early <- calculate_correlation_penalties(X_all, grp_struct_early$groups, grp_struct_early$cor)
        pen_adp_early <- calculate_skat_penalties(X_all, y, grp_struct_early$groups)
        
        # Late Fusion A
        grp_struct_A <- define_groups(full_data[, modA_features])
        pen_std_A <- calculate_correlation_penalties(full_data[, modA_features], grp_struct_A$groups, grp_struct_A$cor)
        pen_adp_A <- calculate_skat_penalties(full_data[, modA_features], y, grp_struct_A$groups)
        
        # Late Fusion B
        grp_struct_B <- define_groups(full_data[, modB_features])
        pen_std_B <- calculate_correlation_penalties(full_data[, modB_features], grp_struct_B$groups, grp_struct_B$cor)
        pen_adp_B <- calculate_skat_penalties(full_data[, modB_features], y, grp_struct_B$groups)

        # --- RUN MODELS ---
        
        # 1. Early Fusion (Standard)
        res_early_std <- perm_assisted_grplasso_binary(
          X_all_scaled, y, N_PERMUTATIONS, grp_struct_early$groups, pen_std_early, SELECTION_CUTOFF
        )
        
        # 2. Early Fusion (Adaptive SKAT)
        res_early_adp <- perm_assisted_grplasso_binary(
          X_all_scaled, y, N_PERMUTATIONS, grp_struct_early$groups, pen_adp_early, SELECTION_CUTOFF
        )

        # 3. Late Fusion (Standard)
        res_late_std_A <- perm_assisted_grplasso_binary(X_modA_scaled, y, N_PERMUTATIONS, grp_struct_A$groups, pen_std_A, SELECTION_CUTOFF)
        res_late_std_B <- perm_assisted_grplasso_binary(X_modB_scaled, y, N_PERMUTATIONS, grp_struct_B$groups, pen_std_B, SELECTION_CUTOFF)
        sel_late_std <- unique(c(res_late_std_A$final_selection, res_late_std_B$final_selection))

        # 4. Late Fusion (Adaptive SKAT)
        res_late_adp_A <- perm_assisted_grplasso_binary(X_modA_scaled, y, N_PERMUTATIONS, grp_struct_A$groups, pen_adp_A, SELECTION_CUTOFF)
        res_late_adp_B <- perm_assisted_grplasso_binary(X_modB_scaled, y, N_PERMUTATIONS, grp_struct_B$groups, pen_adp_B, SELECTION_CUTOFF)
        sel_late_adp <- unique(c(res_late_adp_A$final_selection, res_late_adp_B$final_selection))

        # 5. Baseline LASSO
        cv_fit_lasso <- cv.glmnet(X_all_scaled, y, family = "binomial", alpha = 1)
        lasso_coefs <- coef(cv_fit_lasso, s = "lambda.min")
        sel_lasso <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
        sel_lasso <- sel_lasso[sel_lasso != "(Intercept)"]

        # Store for Stability Calc
        run_features$early_std[[i]] <- res_early_std$final_selection
        run_features$early_adp[[i]] <- res_early_adp$final_selection
        run_features$late_std[[i]] <- sel_late_std
        run_features$late_adp[[i]] <- sel_late_adp
        run_features$lasso[[i]] <- sel_lasso

        # --- AUC Helper ---
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
        
        # Calculate AUCs
        auc_early_std <- calc_auc(X_all_scaled, res_early_std$final_selection, y)
        auc_early_adp <- calc_auc(X_all_scaled, res_early_adp$final_selection, y)
        auc_lasso     <- calc_auc(X_all_scaled, sel_lasso, y)
        
        # Late Fusion Probabilities
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

        probs_std_A <- get_probs(X_modA_scaled, res_late_std_A$final_selection, y, n_val)
        probs_std_B <- get_probs(X_modB_scaled, res_late_std_B$final_selection, y, n_val)
        auc_late_std <- as.numeric(auc(roc(y, as.vector((probs_std_A + probs_std_B)/2), quiet=TRUE)))
        
        probs_adp_A <- get_probs(X_modA_scaled, res_late_adp_A$final_selection, y, n_val)
        probs_adp_B <- get_probs(X_modB_scaled, res_late_adp_B$final_selection, y, n_val)
        auc_late_adp <- as.numeric(auc(roc(y, as.vector((probs_adp_A + probs_adp_B)/2), quiet=TRUE)))

        # Metrics
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
        summarise(across(c(AUC, F1_Score, Sensitivity, FDR, Num_Selected), list(Mean = mean, SD = sd))) %>%
        mutate(n_samples = n_val, effect_multiplier = ef_mult, n_noise = n_noise_val)
      all_experiment_metrics[[run_counter]] <- summary_metrics
      
      run_counter <- run_counter + 1
    }
  }
}

final_metrics_df <- do.call(rbind, all_experiment_metrics)
print(final_metrics_df)

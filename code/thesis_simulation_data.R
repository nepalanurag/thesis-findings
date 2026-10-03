# Thesis simulation: data generation with causal blocks, correlated features
# and noise features for the group lasso experiments. Converted from the
# Thesis_Sim analysis notebook.

library(tidyverse)
library(reda)

create_simulation_data <- function(n_ss = 100, n_gene = 600, effect_sizes = c(0, -2, 3, -4),
                                   n_x = rep(10, 6), sd_v = 0.3) {

  # --- 1. Generate base causal variables (x1, x2, x3...) ---
  # These are uniformly distributed and will be used to generate correlated features.
  # The number of causal variables is one less than the length of effect_sizes.
  num_causal_vars <- length(effect_sizes) - 1
  num_base_vars <- length(n_x)
  base_vars <- matrix(runif(n_ss * num_base_vars), nrow = n_ss)
  colnames(base_vars) <- paste0("x", 1:num_base_vars)

  # --- 2. Generate correlated features (v1_1, v1_2...) ---
  # These features are created by adding scaled noise to the base variables.
  # Features derived from x1, x2, x3... are considered "causal blocks".
  correlated_features <- list()
  for (i in 1:num_base_vars) {
    v_i <- matrix(NA, nrow = n_ss, ncol = n_x[i])
    for (j in 1:n_x[i]) {
      noise <- rnorm(n = n_ss, mean = 0, sd = sd_v)
      v_i[, j] <- base_vars[, i] + (0.01 + 0.5 * (j - 1) / (n_x[i] - 1)) * noise
    }
    colnames(v_i) <- paste0("v", i, "_", 1:n_x[i])
    correlated_features[[i]] <- v_i
  }
  correlated_features <- do.call(cbind, correlated_features)

  # --- 3. Generate uncorrelated noise features (w1, w2...) ---
  num_noise_features <- n_gene - ncol(correlated_features)
  if (num_noise_features < 0) {
    stop("Total genes (n_gene) is smaller than the number of correlated features.")
  }
  noise_features <- matrix(runif(n_ss * num_noise_features), nrow = n_ss)
  colnames(noise_features) <- paste0("w", 1:num_noise_features)

  # --- 4. Generate Survival Outcome ---
  # We bind an intercept column to the causal base variables.
  causal_matrix <- cbind(rep(1, n_ss), base_vars[, 1:num_causal_vars])
  sim_outcome <- simEventData(z = causal_matrix, zCoef = effect_sizes, recurrent = FALSE)

  # --- 5. Combine into a final data frame ---
  final_data <- as.data.frame(cbind(sim_outcome$time, sim_outcome$event, correlated_features, noise_features))
  colnames(final_data)[1:2] <- c("time", "event")

  return(final_data)
}

sample_size <- 200
effect_coefficients <- c(0, -2, 3, -4) # Intercept, plus effects for x1, x2, x3

# This block will contain the causal features that are truly related to the outcome.
cat("Generating Block 1...\n")
block1 <- create_simulation_data(
  n_ss = sample_size,
  n_gene = 800,
  effect_sizes = effect_coefficients
)

# We need to separate the outcome from the features for the fusion step.
outcome_data <- block1[, c("time", "event")]
block1_features <- block1[, -c(1, 2)]

# Add a prefix to column names to identify their origin
colnames(block1_features) <- paste0("B1_", colnames(block1_features))

# This block will be pure noise, with no actual relationship to the outcome.
# We achieve this by setting its effect sizes to zero.
cat("Generating pure noise...\n")
block2 <- create_simulation_data(
  n_ss = sample_size,
  n_gene = 1200,
  effect_sizes = c(0, 0, 0, 0) # No effects = no relationship to the outcome
)
block2_features <- block2[, -c(1, 2)]

# Add a prefix to column names
colnames(block2_features) <- paste0("B2_", colnames(block2_features))

# Combine the feature blocks side-by-side using cbind().
cat("Applying early fusion to combine the two blocks...\n")
early_fusion_features <- cbind(block1_features, block2_features)

# Finally, combine the fused features with the single outcome data.
# This `early_fusion_data` is the final dataset you would use for model training.
early_fusion_data <- cbind(outcome_data, early_fusion_features)

cat("\n--- Simulation Complete ---\n")
cat("Dimensions of the final early fusion dataset:", dim(early_fusion_data), "\n")
cat("Total features:", ncol(early_fusion_data) - 2, "\n")
print(head(early_fusion_data[, c(1:5, 801:804)]))

library(survival)
library(grpreg)

# --- Strategy: Permutation-Assisted Group Lasso on Early Fusion Data ---
# This method uses permutations to find a robust penalty parameter (lambda)
# by learning the level of penalization required to remove all features
# when no true signal exists.

cat("\n--- Applying Permutation-Assisted Group Lasso ---\n")

# The model needs to know which features belong to which block.
# Group 1: Features from Block 1 (signal)
# Group 2: Features from Block 2 (noise)
group_index <- c(
  rep(1, ncol(block1_features)),
  rep(2, ncol(block2_features))
)

# The `grpsurv` function requires a `Surv` object for the response.
surv_obj <- Surv(outcome_data$time, outcome_data$event)
X_matrix <- as.matrix(early_fusion_features)

# We will shuffle the outcome variable repeatedly to break its true relationship
# with the features. For each shuffled version, we find the smallest lambda
# that shrinks all coefficients to zero. This gives us a distribution of
# lambdas under the null hypothesis (no signal).

n_permutations <- 100
permuted_lambda_max <- numeric(n_permutations)

for (i in 1:n_permutations) {
  # Permute the survival times to break the relationship with features
  permuted_time <- sample(outcome_data$time)
  permuted_surv_obj <- Surv(permuted_time, outcome_data$event)
  # Fit a group lasso path on the permuted data.
  fit_perm <- grpsurv(
    X = X_matrix,
    y = permuted_surv_obj,
    group = group_index,
    penalty = "grLasso",
    nlambda = 3
  )
  
  # Store the largest lambda value (the first in the sequence)
  permuted_lambda_max[i] <- fit_perm$lambda[1]
}

# A robust choice for lambda is one that is unlikely to be selected by chance.
# The median of the permuted lambda.max values is a common and effective choice.
lambda_threshold <- median(permuted_lambda_max)

cat("\nCalculated lambda threshold from permutations:", round(lambda_threshold, 4), "\n")

# Now, we fit the group lasso model on the original, un-permuted data,
# using the single, data-driven lambda value we just determined.
cat("Fitting final group lasso model on original data using the permutation-derived lambda...\n")
final_perm_model <- grpsurv(
  X = X_matrix,
  y = surv_obj,
  group = group_index,
  penalty = "grLasso",
  lambda = lambda_threshold
)

cat("\n--- Permutation-Assisted Model Results ---\n")
selected_groups <- predict(final_perm_model, type = "groups")
cat("Groups selected by the model:", selected_groups, "\n")

selected_coeffs <- coef(final_perm_model)
non_zero_coeffs <- selected_coeffs
print(head(non_zero_coeffs, 10))

# Load necessary libraries for survival analysis and Lasso
library(survival)
library(glmnet)

# For reproducibility of random sampling
set.seed(42)

cat("--- Starting Permutation-Assisted Lasso Feature Selection ---\n\n")

# --- Step 1: Prepare Data for Modeling ---
# Separate the predictor matrix (X) from the survival outcome (y)
X <- as.matrix(early_fusion_data[, -c(1, 2)])
y <- Surv(early_fusion_data$time, early_fusion_data$event)

# --- Step 2: Determine the Optimal Penalty (Lambda) via Permutation ---
# The goal is to find a lambda that is strict enough to penalize all
# coefficients to zero when there is no true signal.

# Define parameters for the permutation process
n_permutations <- 100 # Increase for more stable results, e.g., to 1000
lambda_quantile <- 0.90 # Use a high quantile for a stringent penalty

# Vector to store the optimal lambda from each permutation
permuted_lambda_mins <- numeric(n_permutations)

cat(paste0("Running ", n_permutations, " permutations to find the null lambda distribution...\n"))

# Loop to run cross-validation on permuted data
for (i in 1:n_permutations) {
  # Shuffle the outcome to break the feature-outcome relationship
  y_permuted <- y[sample(nrow(y)), ]
  
  # Run cross-validated Lasso for the Cox Proportional Hazards model
  cv_fit_permuted <- cv.glmnet(X, y_permuted, family = "cox", alpha = 1)
  
  # Store the lambda that resulted in the minimum cross-validated error
  permuted_lambda_mins[i] <- cv_fit_permuted$lambda.min
  
  # Print progress
  if (i %% 10 == 0) cat(paste("  Completed permutation", i, "of", n_permutations, "\n"))
}

# Select the optimal lambda from the upper tail of the null distribution
lambda_optimal <- quantile(permuted_lambda_mins, probs = lambda_quantile)

cat("\n--- Permutation Complete ---\n")
cat("Optimal lambda selected from", lambda_quantile * 100, "th percentile:", round(lambda_optimal, 4), "\n\n")

# --- Step 3: Fit Final Model and Select Features ---
# Now, fit the Lasso model on the ORIGINAL data and apply the penalty we found.
cat("Fitting final Lasso model on original data with the selected lambda...\n")
final_model <- glmnet(X, y, family = "cox", alpha = 1)

# Extract the coefficients at our specific, permutation-assisted lambda
final_coeffs <- coef(final_model, s = lambda_optimal)

# Identify features with non-zero coefficients
selected_features_df <- as.data.frame(as.matrix(final_coeffs))
selected_features_df <- selected_features_df[selected_features_df$s1 != 0, , drop = FALSE]
selected_feature_names <- rownames(selected_features_df)

# --- Step 4: Evaluate Selection Performance ---
cat("\n--- Feature Selection Results ---\n")
cat("Total features selected:", length(selected_feature_names), "\n")
print(selected_features_df)

# Check how many of the selected features are true positives.
# In your simulation, true causal features are `B1_v1_...`, `B1_v2_...`, and `B1_v3_...`
true_causal_prefixes <- c("B1_v1_", "B1_v2_", "B1_v3_")

# Find which selected features are true positives vs. false positives
selected_true_positives <- selected_feature_names[grepl(paste(true_causal_prefixes, collapse="|"), selected_feature_names)]
selected_false_positives <- setdiff(selected_feature_names, selected_true_positives)

cat("\n--- Selection Evaluation ---\n")
cat("Number of correctly selected causal features (True Positives):", length(selected_true_positives), "\n")
print(selected_true_positives)

cat("\nNumber of incorrectly selected noise features (False Positives):", length(selected_false_positives), "\n")
print(selected_false_positives)

# Calculate the recall: What percentage of the actual causal features did we find?
total_causal_features <- 3 * 10 # 3 causal groups, each with 10 features
recall <- length(selected_true_positives) / total_causal_features

cat(sprintf("\nRecall Score: %.2f%% of true causal features were identified (%d / %d).\n",
            100 * recall,
            length(selected_true_positives),
            total_causal_features))

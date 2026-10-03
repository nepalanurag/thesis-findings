# Simulation code

`main.r` is the simulation engine behind the thesis findings. It is one self-contained R script (349 lines).

What it does, in order:

1. **Data simulation** (`simulation_data`): builds a multi-group feature set with a few causal features, correlated decoy features within each group, and noise features, then draws a binary outcome from a logistic model on the causal ones.
2. **Grouping and penalties** (`define_groups`, `calculate_correlation_penalties`, `calculate_skat_penalties`): hierarchical clustering of features into groups and adaptive penalty weights.
3. **Permutation-assisted group lasso engine** (`perm_assisted_grplasso_binary`): runs group lasso (via `grpreg`) against row-permuted decoy copies of each feature in parallel (`doParallel`/`foreach`). A feature is kept only if it beats its own decoys repeatedly.
4. **Evaluation helpers**: precision/recall/F1 against the known causal features, pairwise Jaccard across repeats, AUC on held-out data.
5. **Main simulation loop**: sweeps the pilot scenarios and writes the final outputs.

Packages used: tidyverse, glmnet, grpreg, pROC, SKAT, Matrix, doParallel, foreach. Needs R with those packages installed; the loop is the slow part and uses all but one CPU core.

This is the code that produced the simulation numbers in the thesis README (Early Fusion F1 0.96 vs LASSO F1 0.13 in the pilot scenario).

## Analysis scripts

These came from the analysis notebooks written alongside the thesis. Each is a plain R script with the notebook's code chunks in order:

- `adaptive_permutation_group_lasso.R` — the core simulation study: multi-modal feature selection with adaptive permutation-assisted group lasso on simulated data with correlated features and a binary outcome (grpreg, decoy-based selection).
- `tcga_brca_multimodal_selection.R` — the real-data counterpart: multi-modal feature selection on TCGA-BRCA (RNA-seq + methylation) using adaptive group lasso with inverse-log p-value penalties. Downloads the data via curatedTCGAData.
- `thesis_simulation_data.R` — the simulation data generator used by the experiments: causal blocks, correlated features, and noise features.

#!/usr/bin/env Rscript
# Aggregate the stability study into a summary table + morning report.
suppressPackageStartupMessages(library(tidyverse))
OUT <- "tcga_brca_real_results"

old <- suppressMessages(read_csv(file.path(OUT, "five_arm_summary_table.csv"), show_col_types = FALSE))
ef <- readRDS(file.path(OUT, "stability_ef_results.rds"))
lf <- readRDS(file.path(OUT, "stability_lf_results.rds"))
POS_CTRL <- c("ESR1", "PGR", "FOXA1", "GATA3")

jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }
pairwise_jaccard <- function(sel_list) {
  cmb <- combn(length(sel_list), 2)
  mean(apply(cmb, 2, function(ii) jaccard(sel_list[[ii[1]]], sel_list[[ii[2]]])))
}
msd <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x))
msd3 <- function(x) sprintf("%.3f (%.3f)", mean(x), sd(x))
sig_runs <- function(runs, vtag, sig) {
  sum(sapply(runs, function(r) {
    e <- r[[vtag]]$enrich; e$p[e$signature == sig] < 0.05
  }))
}
ctrl_str <- function(runs, vtag) {
  paste(sapply(POS_CTRL, function(g)
    sprintf("%s %d/5", g, sum(sapply(runs, function(r) r[[vtag]]$controls[[g]])))),
    collapse = "; ")
}

runs_ef <- ef[grep("^run", names(ef))]
runs_lf <- lf[grep("^run", names(lf))]

row_for <- function(label, runs, vtag) {
  nsel <- sapply(runs, function(r) length(r[[vtag]]$sel))
  aucs <- sapply(runs, function(r) r[[vtag]]$auc)
  data.frame(
    variant = label,
    n_selected = msd(nsel),
    heldout_AUC = msd3(aucs),
    stability_Jaccard = round(pairwise_jaccard(lapply(runs, function(r) r[[vtag]]$sel)), 3),
    PAM50_sig_runs = sprintf("%d/5", sig_runs(runs, vtag, "PAM50")),
    Oncotype_sig_runs = sprintf("%d/5", sig_runs(runs, vtag, "OncotypeDX")),
    controls_x_of_5 = ctrl_str(runs, vtag),
    stringsAsFactors = FALSE)
}

tab <- bind_rows(
  data.frame(variant = "EF-Adap (current protocol)", n_selected = old$n_selected[old$arm == "EF-Adap"],
             heldout_AUC = old$heldout_AUC[old$arm == "EF-Adap"],
             stability_Jaccard = old$stability_Jaccard[old$arm == "EF-Adap"],
             PAM50_sig_runs = old$PAM50_sig_runs[old$arm == "EF-Adap"],
             Oncotype_sig_runs = old$Oncotype_sig_runs[old$arm == "EF-Adap"],
             controls_x_of_5 = old$controls_x_of_5[old$arm == "EF-Adap"], stringsAsFactors = FALSE),
  row_for("EF-Adap S1 (consensus groups)", runs_ef, "S1"),
  row_for("EF-Adap S2 (+ shrunk weights)", runs_ef, "S2"),
  row_for("EF-Adap S3 (+ stability selection)", runs_ef, "S3"),
  data.frame(variant = "LF-Adap (current protocol)", n_selected = old$n_selected[old$arm == "LF-Adap"],
             heldout_AUC = old$heldout_AUC[old$arm == "LF-Adap"],
             stability_Jaccard = old$stability_Jaccard[old$arm == "LF-Adap"],
             PAM50_sig_runs = old$PAM50_sig_runs[old$arm == "LF-Adap"],
             Oncotype_sig_runs = old$Oncotype_sig_runs[old$arm == "LF-Adap"],
             controls_x_of_5 = old$controls_x_of_5[old$arm == "LF-Adap"], stringsAsFactors = FALSE),
  row_for("LF-Adap S1 (consensus groups)", runs_lf, "S1"),
  row_for("LF-Adap S2 (+ shrunk weights)", runs_lf, "S2"),
  row_for("LF-Adap S3 (+ stability selection)", runs_lf, "S3"))
write.csv(tab, file.path(OUT, "stability_summary_table.csv"), row.names = FALSE)
print(tab)

# Baselines: ESR1-only / four-gene AUC across the 10 fresh splits.
base_ef <- sapply(runs_ef, function(r) r$baseline)
base_lf <- sapply(runs_lf, function(r) r$baseline)
base_all <- cbind(base_ef, base_lf)
cat(sprintf("ESR1-only AUC: %.3f (%.3f), n=10 splits\n",
  mean(base_all["esr1_only", ]), sd(base_all["esr1_only", ])))
cat(sprintf("Four-gene AUC: %.3f (%.3f), n=10 splits\n",
  mean(base_all["four_gene", ]), sd(base_all["four_gene", ])))

dj <- ef$diagnostic$jaccard_within_split
cat(sprintf("DIAG within-split Jaccard: %.3f\n", dj))

# ---- morning report ----
best_ef <- tab$variant[which.max(tab$stability_Jaccard[1:4])]
best_lf <- tab$variant[4 + which.max(tab$stability_Jaccard[5:8])]
rep_lines <- c(
  "# TCGA-BRCA adaptive-model stability study — morning report (2026-09-23)",
  "",
  "## Question",
  "How to increase across-run selection stability (Jaccard) of the adaptive",
  "fusion models without losing the biology signal (signature enrichment,",
  "positive-control recovery) or prediction (held-out AUC).",
  "",
  "## Diagnostic: where does the instability come from?",
  sprintf("- Within ONE fixed train split, current EF-Adap repeated 3x with different",
  ""),
  sprintf("  seeds gives pairwise Jaccard = %.3f (rep sizes %s).", dj,
          paste(ef$diagnostic$sizes, collapse = ", ")),
  if (dj > 0.5) paste("- Verdict: algorithmic jitter is modest; most across-run instability",
    "(Jaccard 0.23) is genuine sampling variation across splits.") else
    paste("- Verdict: algorithmic jitter is substantial; fixing grouping/weights",
    "should recover stability."),
  "",
  "## What was tried (identical data protocol to the 5-arm study)",
  "- S1 consensus groups: bootstrap-averaged |Spearman| on the full cohort",
  "  (unsupervised, no labels), one hclust, fixed for all runs and perms.",
  "  Removes per-run grouping jitter.",
  "- S2 = S1 + shrunk univariate z weights (soft-threshold |z|=1): weak noisy",
  "  signals stop shaking the group penalties.",
  "- S3 = S2 + stability selection: 30 complementary-pair subsamples, 15 inner",
  "  perms each (inner vote cutoff 0.4), keep features in >=60% of subsamples.",
  "- 5 fresh train/test splits per variant; ESR1-only and 4-gene baselines per split.",
  "",
  "## Results",
  "",
  "```",
  paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  sprintf("- ESR1-only baseline AUC: %.3f (%.3f); four-gene: %.3f (%.3f) over 10 splits.",
          mean(base_all["esr1_only", ]), sd(base_all["esr1_only", ]),
          mean(base_all["four_gene", ]), sd(base_all["four_gene", ])),
  "- Read: if adaptive AUC ≈ ESR1-only AUC, the paper's claim is selection",
  "  quality (stability + biology), not prediction. Frame accordingly.",
  "",
  "## Recommendation",
  sprintf("- Early fusion: adopt **%s**.", best_ef),
  sprintf("- Late fusion: adopt **%s**.", best_lf),
  "- S3 costs ~5-6x compute per run; S1/S2 are free. If S3's Jaccard gain over",
  "  S2 is small, ship S2 (simpler story: fixed consensus groups + shrunk weights).",
  "- Next: fold the winner into the Tier-1 rerun (leakage fix, 20 runs,",
  "  stability-selection comparator arm), then METABRIC external validation.",
  "",
  "## Files",
  "- tcga_brca_real_results/stability_summary_table.csv",
  "- tcga_brca_real_results/stability_ef_results.rds / stability_lf_results.rds",
  "- stability_lib.R, stability_ef.R, stability_lf.R, aggregate_stability.R")
writeLines(rep_lines, "STABILITY_REPORT.md")
cat("wrote STABILITY_REPORT.md\n")

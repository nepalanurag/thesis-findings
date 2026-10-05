#!/usr/bin/env Rscript
# Aggregate the 5 independent runs into a thesis-style summary table + plot:
# per arm: mean +/- SD of held-out AUC and #selected, mean pairwise Jaccard
# across runs (thesis "Stability"), signature-enrichment summaries.
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
OUT <- "tcga_brca_real_results"
ARMS <- c(lasso = "Baseline LASSO", `EF-Std` = "EF-Std", `EF-Adap` = "EF-Adap",
          `LF-Std` = "LF-Std", `LF-Adap` = "LF-Adap")
POS_CTRL <- c("ESR1", "PGR", "FOXA1", "GATA3")
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(1); length(intersect(a, b)) / u }
msd <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  sprintf("%.3f (%.3f)", mean(x), if (length(x) > 1) sd(x) else 0)
}
msd0 <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  sprintf("%.0f (%.0f)", mean(x), if (length(x) > 1) sd(x) else 0)
}

runs <- lapply(1:5, function(r) readRDS(file.path(OUT, sprintf("five_arm_results_run%d.rds", r))))
names(runs) <- paste0("run", 1:5)

tab <- bind_rows(lapply(names(ARMS), function(a) {
  aucs <- sapply(runs, function(rr) rr[[a]]$auc)
  ns <- sapply(runs, function(rr) length(rr[[a]]$selected))
  sets <- lapply(runs, function(rr) rr[[a]]$selected)
  pw <- combn(5, 2, function(ii) jaccard(sets[[ii[1]]], sets[[ii[2]]]))
  enr_pam <- sapply(runs, function(rr) { e <- rr[[a]]$enrich; e$p[e$signature == "PAM50"] })
  enr_onco <- sapply(runs, function(rr) { e <- rr[[a]]$enrich; e$p[e$signature == "OncotypeDX"] })
  ov_pam <- sapply(runs, function(rr) { e <- rr[[a]]$enrich; e$overlap[e$signature == "PAM50"] })
  ov_onco <- sapply(runs, function(rr) { e <- rr[[a]]$enrich; e$overlap[e$signature == "OncotypeDX"] })
  ctr <- sapply(POS_CTRL, function(g) sum(sapply(runs, function(rr)
    g %in% sub("^METH_", "", rr[[a]]$selected))))
  data.frame(
    arm = ARMS[[a]],
    n_selected = msd0(ns),
    heldout_AUC = msd(aucs),
    stability_Jaccard = sprintf("%.3f", mean(pw)),
    PAM50_overlap_med = median(ov_pam),
    PAM50_p_med = signif(median(enr_pam), 3),
    PAM50_sig_runs = sprintf("%d/5", sum(enr_pam < 0.05)),
    Oncotype_overlap_med = median(ov_onco),
    Oncotype_p_med = signif(median(enr_onco), 3),
    Oncotype_sig_runs = sprintf("%d/5", sum(enr_onco < 0.05)),
    controls_x_of_5 = paste(sprintf("%s %d/5", POS_CTRL, ctr), collapse = "; "),
    stringsAsFactors = FALSE)
}))
print(tab)
write.csv(tab, file.path(OUT, "five_arm_summary_table.csv"), row.names = FALSE)

# long-format numbers for the plot
plot_df <- bind_rows(lapply(names(ARMS), function(a) {
  data.frame(arm = ARMS[[a]],
    auc = sapply(runs, function(rr) rr[[a]]$auc),
    n = sapply(runs, function(rr) length(rr[[a]]$selected)))
}))
summ_df <- plot_df %>% group_by(arm) %>%
  summarise(auc_m = mean(auc, na.rm = TRUE), auc_s = sd(auc, na.rm = TRUE),
    n_m = mean(n), n_s = sd(n), .groups = "drop") %>%
  mutate(arm = factor(arm, levels = unname(ARMS)))
summ_df$auc_s[is.nan(summ_df$auc_s)] <- 0

p1 <- ggplot(summ_df, aes(x = arm, y = auc_m, fill = arm)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_errorbar(aes(ymin = auc_m - auc_s, ymax = auc_m + auc_s), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", auc_m)), vjust = -0.6, size = 3.2) +
  labs(title = "Held-out AUC (mean \u00b1 SD, 5 runs)", x = NULL, y = "AUC") +
  ylim(0, 1.05) + theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
p2 <- ggplot(summ_df, aes(x = arm, y = n_m, fill = arm)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_errorbar(aes(ymin = pmax(0, n_m - n_s), ymax = n_m + n_s), width = 0.2) +
  geom_text(aes(label = sprintf("%.0f", n_m)), vjust = -0.6, size = 3.2) +
  labs(title = "# selected features (mean \u00b1 SD, 5 runs)", x = NULL, y = "count") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
fig <- p1 + p2 + patchwork::plot_annotation(
  title = "TCGA-BRCA 5-arm comparison: ER status (RNA + methylation, n=737)",
  subtitle = "Independent train/test splits per run; selection on train only; AUC on held-out test\nEF-Std: 4/5 runs selected 0 features (AUC from the single non-empty run); its Jaccard is inflated by empty-set agreement")
ggsave(file.path(OUT, "five_arm_plot.png"), fig, width = 10, height = 4.6, dpi = 150)
ggsave(file.path(OUT, "five_arm_plot.pdf"), fig, width = 10, height = 4.6)
cat("AGGREGATION DONE\n")

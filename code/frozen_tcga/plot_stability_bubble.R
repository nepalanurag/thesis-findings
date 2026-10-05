#!/usr/bin/env Rscript
# Thesis-style "Prediction vs Reproducibility" bubble plot for the stability study.
# X = selection stability (mean pairwise Jaccard across 5 fresh splits)
# Y = held-out AUC (mean +/- SD)
# Bubble size = mean # selected features
suppressPackageStartupMessages(library(ggplot2))

OUT <- "tcga_brca_real_results"

d <- data.frame(
  arm    = c("EF-Adap (current protocol)", "EF-Adap S1 (consensus groups)",
             "EF-Adap S2 (+ shrunk weights)", "EF-Adap S3 (+ stability selection)",
             "LF-Adap (current protocol)", "LF-Adap S1 (consensus groups)",
             "LF-Adap S2 (+ shrunk weights)", "LF-Adap S3 (+ stability selection)"),
  short  = c("EF current", "EF S1", "EF S2", "EF S3",
             "LF current", "LF S1", "LF S2", "LF S3"),
  family = c(rep("EF-Adap", 4), rep("LF-Adap", 4)),
  jaccard = c(0.229, 0.482, 0.483, 0.680, 0.327, 0.535, 0.557, 0.588),
  auc     = c(0.930, 0.940, 0.940, 0.948, 0.924, 0.946, 0.945, 0.944),
  auc_sd  = c(0.017, 0.008, 0.008, 0.006, 0.015, 0.012, 0.013, 0.016),
  nsel    = c(236.0, 206.6, 206.2, 122.8, 498.0, 437.2, 427.8, 206.8),
  stringsAsFactors = FALSE
)
d$stage <- factor(c("current", "S1", "S2", "S3", "current", "S1", "S2", "S3"),
                  levels = c("current", "S1", "S2", "S3"))

fam_col <- c("EF-Adap" = "#1B9E77", "LF-Adap" = "#377EB8")

# manual label anchors (x, y, hjust, vjust)
lab <- data.frame(
  short = d$short,
  lx = c(0.229, 0.452, 0.505, 0.680, 0.327, 0.535, 0.562, 0.620),
  ly = c(0.9155, 0.9518, 0.9272, 0.9572, 0.9100, 0.9605, 0.9338, 0.9318),
  hj = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0),
  vj = c(1, 0, 1, 0, 1, 0, 0, 0.5)
)
d <- merge(d, lab, by = "short")
d$bold <- d$short == "EF S3"

cap1 <- "Pale bubbles = current protocol; solid = stabilized (S1 consensus groups -> S2 shrunk weights -> S3 stability selection)."
cap2 <- "EF-Adap S3 adopted: Jaccard 0.68, 123 features, AUC 0.948, ESR1/PGR/FOXA1/GATA3 recovered 5/5."
cap3 <- "Caution: LF-Adap S3 gains Jaccard but loses the positive controls in most splits - stability alone must not pick the winner."
caption_txt <- paste(c(paste(strwrap(cap1, width = 120), collapse = "\n"),
                       paste(strwrap(cap2, width = 120), collapse = "\n"),
                       paste(strwrap(cap3, width = 120), collapse = "\n")),
                     collapse = "\n")

p <- ggplot(d, aes(x = jaccard, y = auc)) +
  # progression trajectories: current -> S1 -> S2 -> S3 per family
  geom_path(data = d[order(d$family, d$stage), ],
            aes(group = family, color = family),
            arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
            linewidth = 0.5, linetype = "dashed", alpha = 0.55, show.legend = FALSE) +
  geom_errorbar(aes(ymin = auc - auc_sd, ymax = auc + auc_sd),
                width = 0.012, color = "grey35", linewidth = 0.5) +
  geom_point(aes(size = nsel, fill = family, alpha = stage),
             shape = 21, color = "grey20", stroke = 0.6) +
  geom_text(aes(x = lx, y = ly, label = short, hjust = hj, vjust = vj,
                fontface = ifelse(bold, "bold", "plain")),
            size = 3.6, color = "grey15", show.legend = FALSE) +
  # highlight the adopted spec
  annotate("text", x = 0.680, y = 0.9520, label = "adopted", size = 3.2,
           color = "#1B9E77", fontface = "italic") +
  scale_fill_manual(values = fam_col, name = "Model") +
  scale_color_manual(values = fam_col, guide = "none") +
  scale_alpha_manual(values = c(current = 0.30, S1 = 0.65, S2 = 0.85, S3 = 1.0),
                     guide = "none") +
  scale_size_area(name = "mean # selected\nfeatures", max_size = 24,
                  breaks = c(150, 300, 450)) +
  scale_x_continuous("Selection stability (mean pairwise Jaccard across 5 splits)",
                     limits = c(0.15, 0.76), breaks = seq(0.2, 0.7, 0.1)) +
  scale_y_continuous("Predictive power (held-out AUC, mean +/- SD)",
                     limits = c(0.908, 0.963), breaks = seq(0.91, 0.96, 0.01)) +
  labs(
    title = "Prediction vs Reproducibility: stabilizing the adaptive fusion models",
    subtitle = "TCGA-BRCA, ER-IHC endpoint (RNA + methylation, n = 737); 5 fresh train/test splits per arm",
    caption = caption_txt
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.caption = element_text(hjust = 0, size = 9, color = "grey30"),
        legend.position = "right",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "stability_bubble.png"), p, width = 10, height = 6.5, dpi = 150)
ggsave(file.path(OUT, "stability_bubble.pdf"), p, width = 10, height = 6.5)
cat("wrote", file.path(OUT, "stability_bubble.png"), "and .pdf\n")

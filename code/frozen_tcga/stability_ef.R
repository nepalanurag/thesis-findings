#!/usr/bin/env Rscript
# Stability study: EARLY-FUSION adaptive arm (+ diagnostic). RESUMABLE:
# skips the diagnostic and any (run, variant) already in the .rds, so a
# killed job can be relaunched without losing completed work.
source("stability_lib.R")
OUT <- "tcga_brca_real_results"
rds_path <- file.path(OUT, "stability_ef_results.rds")
logf <- file.path(OUT, "run_log_stability_ef.txt")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
save_results <- function() saveRDS(results, rds_path)

if (file.exists(rds_path)) {
  results <- readRDS(rds_path)
  logmsg("resumed: already have %s", paste(names(results), collapse = ", "))
} else results <- list()

logmsg("loading data...")
D <- load_data()
Xr <- D$Xr; Xm <- D$Xm; X_ef <- D$X_ef; y <- D$y; universe <- D$universe
logmsg("data: n=%d, features=%d", D$n, ncol(X_ef))

if (!("grp_cons_done" %in% names(results))) {
  logmsg("computing consensus groups (full cohort)...")
  t0 <- Sys.time()
  grp_ef_cons <- consensus_groups(X_ef, B = 10, seed = 7)
  logmsg("consensus groups done: %d groups (%.1f min)", length(unique(grp_ef_cons)),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
  results$grp_ef_cons <- grp_ef_cons
  results$grp_cons_done <- TRUE
  save_results()
} else {
  grp_ef_cons <- results$grp_ef_cons
  logmsg("using saved consensus groups (%d groups)", length(unique(grp_ef_cons)))
}

N_RUNS <- 5

## ---------------- DIAGNOSTIC: within-split repeats ----------------
if (!("diagnostic" %in% names(results))) {
  logmsg("=== DIAGNOSTIC: within-split algorithmic jitter ===")
  sp0 <- make_split(y, 20260909 + 1)
  Xtr0 <- X_ef[sp0$tr, ]; ytr0 <- y[sp0$tr]
  grp0 <- define_groups(Xtr0)
  diag_sels <- lapply(1:3, function(r) {
    set.seed(900 + r)
    pen <- valid_adaptive_pen(Xtr0, ytr0, grp0)
    perm_select(Xtr0, ytr0, 80, grp0, pen, FDR_TARGET, 0.5,
                sprintf("DIAG-rep%d", r), logmsg)$sel
  })
  dj <- pairwise_jaccard(diag_sels)
  logmsg("[DIAG] within-split pairwise Jaccard (3 reps): %.3f", dj)
  results$diagnostic <- list(jaccard_within_split = dj, sizes = sapply(diag_sels, length))
  save_results()
} else logmsg("diagnostic already done, skipping")

## ---------------- S1/S2/S3 across 5 fresh runs ----------------
for (i in seq_len(N_RUNS)) {
  rname <- sprintf("run%d", i)
  have <- names(results[[rname]])
  if (all(c("S1", "S2", "S3") %in% have)) { logmsg("[run %d] complete, skipping", i); next }
  seed <- 20260923 + i
  sp <- make_split(y, seed)
  tr <- sp$tr; te <- sp$te
  Xtr <- X_ef[tr, ]; ytr <- y[tr]; Xte <- X_ef[te, ]; yte <- y[te]
  if (!("baseline" %in% have)) {
    results[[rname]]$baseline <- simple_baselines(Xr, y, tr, te)
    b <- results[[rname]]$baseline
    logmsg("[run %d] ESR1-only AUC=%.3f four-gene AUC=%.3f", i, b[["esr1_only"]], b[["four_gene"]])
    save_results()
  }

  if (!("S1" %in% have)) {
    t0 <- Sys.time(); set.seed(seed)
    pen_s1 <- valid_adaptive_pen(Xtr, ytr, grp_ef_cons)
    s1 <- perm_select(Xtr, ytr, 80, grp_ef_cons, pen_s1, FDR_TARGET, 0.5, sprintf("S1-run%d", i), logmsg)
    auc_s1 <- refit_auc(Xtr, ytr, Xte, yte, s1$sel)
    enr_s1 <- enrichment(s1$sel, universe)
    ctrl_s1 <- sapply(POS_CTRL, function(g) g %in% strip_prefix(s1$sel))
    logmsg("[S1 run %d] sel=%d auc=%.3f (%.1f min)", i, length(s1$sel), auc_s1,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
    results[[rname]]$S1 <- list(sel = s1$sel, auc = auc_s1, enrich = enr_s1, controls = ctrl_s1)
    save_results()
  } else logmsg("[run %d] S1 already done, skipping", i)

  if (!("S2" %in% have)) {
    t0 <- Sys.time(); set.seed(seed)
    pen_s2 <- shrunk_adaptive_pen(Xtr, ytr, grp_ef_cons)
    s2 <- perm_select(Xtr, ytr, 80, grp_ef_cons, pen_s2, FDR_TARGET, 0.5, sprintf("S2-run%d", i), logmsg)
    auc_s2 <- refit_auc(Xtr, ytr, Xte, yte, s2$sel)
    enr_s2 <- enrichment(s2$sel, universe)
    ctrl_s2 <- sapply(POS_CTRL, function(g) g %in% strip_prefix(s2$sel))
    logmsg("[S2 run %d] sel=%d auc=%.3f (%.1f min)", i, length(s2$sel), auc_s2,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
    results[[rname]]$S2 <- list(sel = s2$sel, auc = auc_s2, enrich = enr_s2, controls = ctrl_s2)
    save_results()
  } else logmsg("[run %d] S2 already done, skipping", i)

  if (!("S3" %in% have)) {
    t0 <- Sys.time()
    pen_s2b <- shrunk_adaptive_pen(Xtr, ytr, grp_ef_cons)
    s3 <- stab_select(Xtr, ytr, grp_ef_cons, pen_s2b, B = 30, inner_perm = 15,
                      inner_cut = 0.4, pi_thr = 0.6, seed = seed, tag = sprintf("S3-run%d", i))
    auc_s3 <- refit_auc(Xtr, ytr, Xte, yte, s3$sel)
    enr_s3 <- enrichment(s3$sel, universe)
    ctrl_s3 <- sapply(POS_CTRL, function(g) g %in% strip_prefix(s3$sel))
    logmsg("[S3 run %d] sel=%d auc=%.3f (%.1f min)", i, length(s3$sel), auc_s3,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
    results[[rname]]$S3 <- list(sel = s3$sel, auc = auc_s3, enrich = enr_s3, controls = ctrl_s3)
    save_results()
    write.csv(data.frame(feature = s3$sel, freq = round(as.numeric(s3$freq[s3$sel]), 3)),
      file.path(OUT, sprintf("selected_stabEF_S3_run%d.csv", i)), row.names = FALSE)
  } else logmsg("[run %d] S3 already done, skipping", i)
}
logmsg("STABILITY EF DONE")

#!/usr/bin/env Rscript
# Stability study: LATE-FUSION adaptive arm. RESUMABLE (see stability_ef.R).
source("stability_lib.R")
OUT <- "tcga_brca_real_results"
rds_path <- file.path(OUT, "stability_lf_results.rds")
logf <- file.path(OUT, "run_log_stability_lf.txt")
logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, file = logf, append = TRUE) }
save_results <- function() saveRDS(results, rds_path)

if (file.exists(rds_path)) {
  results <- readRDS(rds_path)
  logmsg("resumed: already have %s", paste(names(results), collapse = ", "))
} else results <- list()

logmsg("loading data...")
D <- load_data()
Xr <- D$Xr; Xm <- D$Xm; y <- D$y; universe <- D$universe
logmsg("data: n=%d", D$n)

if (!("grp_cons_done" %in% names(results))) {
  logmsg("computing per-modality consensus groups (full cohort)...")
  t0 <- Sys.time()
  grp_rna_cons <- consensus_groups(Xr, B = 10, seed = 11)
  grp_met_cons <- consensus_groups(Xm, B = 10, seed = 12)
  logmsg("consensus groups done: RNA %d / METH %d groups (%.1f min)",
    length(unique(grp_rna_cons)), length(unique(grp_met_cons)),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
  results$grp_rna_cons <- grp_rna_cons
  results$grp_met_cons <- grp_met_cons
  results$grp_cons_done <- TRUE
  save_results()
} else {
  grp_rna_cons <- results$grp_rna_cons
  grp_met_cons <- results$grp_met_cons
  logmsg("using saved consensus groups")
}

N_RUNS <- 5
mods <- list(RNA = "Xr", METH = "Xm")

run_variant <- function(vtag, i, seed, tr, te, ytr, yte) {
  sel_m <- list(); probs_m <- list(); auc_m <- c()
  for (mn in names(mods)) {
    X <- get(mods[[mn]])
    Xtrm <- X[tr, , drop = FALSE]; Xtem <- X[te, , drop = FALSE]
    gm <- if (mn == "RNA") grp_rna_cons else grp_met_cons
    set.seed(seed)
    penm <- if (vtag == "S1") valid_adaptive_pen(Xtrm, ytr, gm)
            else shrunk_adaptive_pen(Xtrm, ytr, gm)
    if (vtag == "S3") {
      sm <- stab_select(Xtrm, ytr, gm, penm, B = 30, inner_perm = 15,
                        inner_cut = 0.4, pi_thr = 0.6, seed = seed,
                        tag = sprintf("%s-run%d-%s", vtag, i, mn))
      sel <- sm$sel
    } else {
      sm <- perm_select(Xtrm, ytr, 80, gm, penm, FDR_TARGET, 0.5,
                        sprintf("%s-run%d-%s", vtag, i, mn), logmsg)
      sel <- sm$sel
    }
    sel_m[[mn]] <- sel
    if (length(sel) >= 2) {
      fit <- cv.glmnet(Xtrm[, sel, drop = FALSE], ytr, family = "binomial")
      probs_m[[mn]] <- as.numeric(predict(fit, Xtem[, sel, drop = FALSE],
        s = "lambda.min", type = "response"))
    } else probs_m[[mn]] <- NULL
    auc_m[mn] <- mod_cv_auc(Xtrm, ytr, sel)
    logmsg("[%s run %d %s] sel=%d cv_auc=%.3f", vtag, i, mn, length(sel), auc_m[mn])
  }
  sel_union <- unique(unlist(sel_m))
  w <- pmax(0, auc_m - 0.5); w[is.na(w)] <- 0
  auc_ens <- NA_real_
  if (sum(w) > 0) {
    w <- w / sum(w)
    P <- sapply(names(mods), function(mn) {
      pr <- probs_m[[mn]]
      if (is.null(pr)) rep(mean(ytr), length(yte)) else pr
    })
    auc_ens <- as.numeric(auc(roc(yte, as.numeric(P %*% w), quiet = TRUE)))
  }
  enr <- enrichment(sel_union, universe)
  ctrl <- sapply(POS_CTRL, function(g) g %in% strip_prefix(sel_union))
  logmsg("[%s run %d] union_sel=%d ens_auc=%.3f", vtag, i, length(sel_union), auc_ens)
  list(sel = sel_union, per_mod = sel_m, auc = auc_ens, enrich = enr,
       controls = ctrl, mod_auc = auc_m)
}

for (i in seq_len(N_RUNS)) {
  rname <- sprintf("run%d", i)
  have <- names(results[[rname]])
  if (all(c("S1", "S2", "S3") %in% have)) { logmsg("[run %d] complete, skipping", i); next }
  seed <- 20260933 + i
  sp <- make_split(y, seed)
  tr <- sp$tr; te <- sp$te
  ytr <- y[tr]; yte <- y[te]
  if (!("baseline" %in% have)) {
    results[[rname]]$baseline <- simple_baselines(Xr, y, tr, te)
    b <- results[[rname]]$baseline
    logmsg("[run %d] ESR1-only AUC=%.3f four-gene AUC=%.3f", i, b[["esr1_only"]], b[["four_gene"]])
    save_results()
  }
  for (vtag in c("S1", "S2", "S3")) {
    if (vtag %in% names(results[[rname]])) { logmsg("[run %d] %s already done, skipping", i, vtag); next }
    t0 <- Sys.time()
    results[[rname]][[vtag]] <- run_variant(vtag, i, seed, tr, te, ytr, yte)
    logmsg("[%s run %d] total %.1f min", vtag, i, as.numeric(difftime(Sys.time(), t0, units = "mins")))
    save_results()
  }
}
logmsg("STABILITY LF DONE")

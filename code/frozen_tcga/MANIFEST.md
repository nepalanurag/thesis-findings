# Frozen TCGA-BRCA validation scripts — manifest

Scripts run for the real-data validation behind `REPORT.md`, frozen as of
September 2026. Files here are the exact scripts as run (no edits), copied
into the repo for reproducibility. Per-run `.rds` outputs and the raw TCGA
data files are not committed (see "What is not committed" below); the
summary tables and console logs in `logs/` capture each run's results.

How to run: see `REPRODUCE.md` at the repo root.

## REPORT.md table/section -> script

### "As-written vs fixed (RNA only, ER-IHC)"

| REPORT row | Script | Seed | Console log |
|---|---|---|---|
| V1 (thesis as written): 302 genes, AUC 0.974 | `tcga-brca-real.R` | `set.seed(20260908)` | `logs/real_out.txt` |
| V3 (fixed weights + FDR lambda): 224 genes, AUC 0.954 | `tcga-brca-v3.R` (adaptive arm) | `set.seed(20260909)` | `logs/v3_out.txt` |
| V3 standard weights: 75 genes, AUC 0.926 | `tcga-brca-v3.R` (standard arm, same run) | `set.seed(20260909)` | `logs/v3_out.txt` |

`tcga-brca-v2.R` is an intermediate version of the same RNA-only experiment
(FDR-targeted lambda with the older weights); it is kept for provenance but
the REPORT table rows above come from `tcga-brca-real.R` and `tcga-brca-v3.R`.
The 0.9-frequency-cutoff note (152 genes, controls all kept) comes from
re-thresholding the V3 adaptive arm's per-permutation frequencies.

### "Five-arm benchmark (RNA + methylation, ER-IHC, 5 splits)"

| Piece | Script | Seed scheme | Logs / outputs |
|---|---|---|---|
| Per-run analysis, all 5 arms (LASSO, EF-Std, EF-Adap, LF-Std, LF-Adap), one fresh 70/30 split each | `tcga-brca-5arm.R`, run once per run id: `Rscript tcga-brca-5arm.R 1` ... `Rscript tcga-brca-5arm.R 5` | `set.seed(20260909 + RUN_ID)` | `logs/5arm_out.txt` (run 1), `logs/5arm_runs2345_out.txt` (runs 2-5), `logs/5arm_runs345_out.txt` (runs 3-5) |
| Aggregation: means, SDs, pairwise Jaccard, enrichment tallies | `aggregate-5arm.R` | deterministic given the per-run `.rds` files | `logs/five_arm_summary_table.csv` |

`five_arm_summary_table.csv` in `logs/` is the exact table the REPORT
five-arm benchmark rows are built from. Run 1 also carries the within-run
subsample stability diagnostic (`DO_STAB`).

### "Stability study (S1 consensus groups, S2 shrunk weights, S3 stability selection)"

| Piece | Script | Seed scheme | Logs / outputs |
|---|---|---|---|
| Shared library (data loading, consensus groups, selection engine, S1-S3 variants) | `stability_lib.R` | `set.seed(seed)` at each caller entry | — |
| Early-fusion runs: EF-Adap current, S1, S2, S3, 5 fresh splits each | `stability_ef.R` | consensus groups `seed = 7`; run seeds via `set.seed(900 + r)`; per-variant `set.seed(seed)` | `logs/run_log_stability_ef.txt` |
| Late-fusion runs: LF-Adap current, S1, S2, S3, 5 fresh splits each | `stability_lf.R` | same scheme as `stability_ef.R` | `logs/run_log_stability_lf.txt` |
| Aggregation: summary table and morning report | `aggregate_stability.R` | deterministic given the `.rds` files | `logs/stability_summary_table.csv`, `logs/STABILITY_REPORT.md` |

`stability_summary_table.csv` is the exact table the REPORT stability-study
rows are built from. Both scripts are resumable: already-computed
(run, variant) results in the `.rds` are skipped on re-launch.
`plot_stability_bubble.R` builds the stability bubble plot in
`logs/` notes from the summary tables.

### Calibration note (ESR1-only AUC 0.943, four-gene 0.939)

The per-split ESR1-only and four-gene baselines are computed inside
`stability_ef.R` / `stability_lf.R` (`esr1_only`, `four_gene`) and averaged
in `aggregate_stability.R` over the 10 fresh splits (5 EF + 5 LF). See
`logs/STABILITY_REPORT.md`.

### "Secondary endpoint" (PAM50 Basal-like vs rest: 164 genes, AUC 0.994, Jaccard 0.73)

| Script | Seed | Console log |
|---|---|---|
| `tcga-brca-basal-v3.R` (V3 machinery: RMS-z weights + FDR lambda) | `set.seed(20260909)` | `logs/basal_out.txt` |

`tcga-brca-basal.R` is the original basal script (as-written machinery);
the REPORT numbers come from `tcga-brca-basal-v3.R`. `logs/MORNING_REPORT.md`
holds the accompanying notes.

### "Simulation cross-check" (rebuilt thesis pilot)

`fusion_sim_extracted.R` rebuilds the thesis pilot simulation in a fresh
R environment (Early Fusion Std F1 0.96 / sensitivity 0.93 vs LASSO F1
0.13 / sensitivity 0.10). Seeds set per run inside the script.

### To be confirmed

- `tcga-brca-validation.R`: the comprehensive validation script (RNA +
  methylation, ER-IHC and basal endpoints, subsample stability) with its own
  seeds (`set.seed(20260908)` at run start, `set.seed(1000 + b)` in the
  stability bootstrap). Which exact REPORT numbers, if any, trace back to
  this script rather than the v3/5arm/stability scripts above still needs my
  confirmation; it is committed so the full provenance is on record.

## Data

Raw TCGA files (RNA-seq RSEM, HM450 methylation, clinical) are not in this
repo. `dl_rna.py` and `dl_clinical.py` are the exact cBioPortal download
scripts used; expected sizes are ~184 MB (RNA), ~226 MB (methylation),
~2 MB (clinical). See `REPRODUCE.md` for the full setup.

## What is not committed

- `tcga_data/` (~410 MB raw TCGA files) — re-downloadable via `dl_rna.py`
  / `dl_clinical.py` from cBioPortal (see `REPRODUCE.md`).
- `tcga_brca_real_results/*.rds` (per-run result objects) — regenerated by
  running the scripts; the summary CSVs in `logs/` capture the reported
  numbers.

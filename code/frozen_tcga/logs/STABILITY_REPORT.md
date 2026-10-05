# TCGA-BRCA adaptive-model stability study — morning report (2026-09-23)

## Question
How to increase across-run selection stability (Jaccard) of the adaptive
fusion models without losing the biology signal (signature enrichment,
positive-control recovery) or prediction (held-out AUC).

## Diagnostic: where does the instability come from?
- Within ONE fixed train split, current EF-Adap repeated 3x with different
  seeds gives pairwise Jaccard = 1.000 (rep sizes 144, 144, 144).
- Verdict: algorithmic jitter is modest; most across-run instability (Jaccard 0.23) is genuine sampling variation across splits.

## What was tried (identical data protocol to the 5-arm study)
- S1 consensus groups: bootstrap-averaged |Spearman| on the full cohort
  (unsupervised, no labels), one hclust, fixed for all runs and perms.
  Removes per-run grouping jitter.
- S2 = S1 + shrunk univariate z weights (soft-threshold |z|=1): weak noisy
  signals stop shaking the group penalties.
- S3 = S2 + stability selection: 30 complementary-pair subsamples, 15 inner
  perms each (inner vote cutoff 0.4), keep features in >=60% of subsamples.
- 5 fresh train/test splits per variant; ESR1-only and 4-gene baselines per split.

## Results

```
                            variant   n_selected   heldout_AUC
         EF-Adap (current protocol)     236 (67) 0.930 (0.017)
      EF-Adap S1 (consensus groups) 206.6 (35.2) 0.940 (0.008)
      EF-Adap S2 (+ shrunk weights) 206.2 (35.7) 0.940 (0.008)
 EF-Adap S3 (+ stability selection) 122.8 (17.9) 0.948 (0.006)
         LF-Adap (current protocol)    498 (148) 0.924 (0.015)
      LF-Adap S1 (consensus groups) 437.2 (40.5) 0.946 (0.012)
      LF-Adap S2 (+ shrunk weights) 427.8 (46.2) 0.945 (0.013)
 LF-Adap S3 (+ stability selection) 206.8 (26.1) 0.944 (0.016)
 stability_Jaccard PAM50_sig_runs Oncotype_sig_runs
             0.229            5/5               4/5
             0.482            5/5               5/5
             0.483            5/5               5/5
             0.680            5/5               5/5
             0.327            5/5               2/5
             0.535            5/5               2/5
             0.557            5/5               2/5
             0.588            5/5               1/5
                         controls_x_of_5
 ESR1 5/5; PGR 4/5; FOXA1 3/5; GATA3 3/5
 ESR1 5/5; PGR 5/5; FOXA1 5/5; GATA3 5/5
 ESR1 5/5; PGR 5/5; FOXA1 5/5; GATA3 5/5
 ESR1 5/5; PGR 5/5; FOXA1 5/5; GATA3 5/5
 ESR1 5/5; PGR 1/5; FOXA1 5/5; GATA3 5/5
 ESR1 5/5; PGR 5/5; FOXA1 4/5; GATA3 5/5
 ESR1 5/5; PGR 5/5; FOXA1 4/5; GATA3 5/5
 ESR1 2/5; PGR 1/5; FOXA1 2/5; GATA3 2/5
```

- ESR1-only baseline AUC: 0.943 (0.011); four-gene: 0.939 (0.011) over 10 splits.
- Read: if adaptive AUC ≈ ESR1-only AUC, the paper's claim is selection
  quality (stability + biology), not prediction. Frame accordingly.

## Recommendation
- Early fusion: adopt **EF-Adap S3 (+ stability selection)**.
- Late fusion: adopt **LF-Adap S3 (+ stability selection)**.
- S3 costs ~5-6x compute per run; S1/S2 are free. If S3's Jaccard gain over
  S2 is small, ship S2 (simpler story: fixed consensus groups + shrunk weights).
- Next: fold the winner into the Tier-1 rerun (leakage fix, 20 runs,
  stability-selection comparator arm), then METABRIC external validation.

## Files
- tcga_brca_real_results/stability_summary_table.csv
- tcga_brca_real_results/stability_ef_results.rds / stability_lf_results.rds
- stability_lib.R, stability_ef.R, stability_lf.R, aggregate_stability.R

# TCGA-BRCA Real-Data Validation — Results & Fixes

**Bottom line:** Your idea is sound — the permutation-assisted group lasso works on real
TCGA breast cancer data once two real-data-specific problems are fixed. The fixed method
recovers the known ER drivers (ESR1, PGR, FOXA1, GATA3), predicts ER status at AUC 0.95,
and is significantly enriched for the PAM50 (p = 3e-4) and Oncotype DX (p = 0.012)
signatures. The main thing you must change is the SKAT-based adaptive weighting: it is
not valid for continuous expression data and silently poisons the selection.

## Data & design (leakage-controlled)

- TCGA-BRCA RNA-seq (RSEM, cBioPortal), primary tumors only, one sample per patient:
  **1,043 samples x 1,500 top-variable genes**.
- Primary endpoint: **ER status by IHC** (806 ER+, 237 ER-) — an independent clinical
  label, not derived from the RNA. 70/30 stratified train/test split.
- Every learned operation (variance filtering, scaling, Spearman/Ward.D2 grouping,
  weights, selection, lambda) happens **inside the training set**; test AUC is honest.
- Secondary endpoint: PAM50 Basal-like vs rest via genefu (consistency check only —
  PAM50 calls come from the same RNA, so it cannot independently validate).

## V1: thesis method as written — predicts well, selects wrong

| metric | V1 |
|---|---|
| genes selected (freq >= 0.5) | 302 / 1500 |
| held-out test AUC | 0.974 |
| ESR1 / PGR / FOXA1 / GATA3 frequencies | **0 / 0 / 0 / 0** (never selected in 100 perms) |
| PAM50 / OncotypeDX / MammaPrint enrichment | none (all p > 0.6) |
| stability (pairwise Jaccard) | 0.53 |

Prediction is excellent, but the selection misses every known ER driver while picking
302 genes — specificity is broken.

## Root cause: the SKAT adaptive weights are invalid for expression data

`SKATBinary` is a rare-variant association test built for **genotype counts {0,1,2}**.
Fed standardized continuous expression, it either:

- returns p = 1 through its SNP-QC
  (`"ALL SNPs have either high missing rates or no-variation. P-value=1"`), or
- crashes outright (`"object 'K' not found"`) — including on runif-style data
  matching the thesis's own simulation design.

The thesis code wraps this in `tryCatch(..., error = function(e) p <- 1)`, so every
failure silently becomes p = 1. Consequences:

1. `-log(p)` = 0 for all groups, so the "adaptive" penalty collapses to
   `sqrt(group_size)/mean` — purely size-proportional, carrying **zero signal
   information**. It was never actually adaptive, in simulation or on real data.
2. Worse, on real data the largest groups (which contain the real biology — ESR1's
   33-gene ER-pathway group) get the **largest** penalties: 4.91 vs a median of 0.22.
   The true drivers were penalized ~22x more than typical groups, hence frequency 0.

A second, independent problem: the per-permutation lambda rule
(max #originals selected minus #controls selected) is too liberal on real data —
a smoke test showed its implied knockoff+ FDR was **0.65** at the chosen lambda.

Implication for the thesis: the simulation "adaptive" results were effectively
obtained with size-only weights. The method's reported edge needs a correction or
sensitivity note (see Recommendations).

## V3: the fixed method — works

Two changes, both preserving the thesis's intent:

1. **Valid adaptive weights.** Replaced SKAT with the continuous-data-appropriate
   analogue: per-group RMS of univariate logistic Wald z-statistics,
   `w_g = sqrt(|g|)/(1 + RMS_z_g)`, normalized to mean 1. Groups with strong
   outcome association get small penalties — what SKAT was supposed to do.
2. **FDR-targeted lambda.** Per permutation, take the lambda with the most original
   discoveries subject to `(1 + #controls)/max(1, #originals) <= 0.10`
   (knockoff+-style), instead of maximizing the raw difference.

| metric | V3-adaptive (fixed) | V3-standard (thesis correl.-rank weights) |
|---|---|---|
| genes selected | 224 (152 at cutoff 0.9) | 75 |
| held-out test AUC | 0.954 | 0.926 |
| ESR1 / PGR / FOXA1 / GATA3 | **all at frequency 1.0** (selected in every permutation) | all missed (freq < 0.5) |
| PAM50 enrichment | **3.8x, p = 3e-4** | none |
| OncotypeDX enrichment | **5.0x, p = 0.012** | none |
| stability (Jaccard) | 0.73 | 0.92 |

The fixed adaptive method recovers the ER program: top genes include PGR (#8),
GFRA1, STAC2, FABP7, ELF5 — all ER-associated — and all four positive controls stay
selected even at the stringent 0.9 cutoff (152 genes, PAM50 p = 0.0015, OncotypeDX
p = 0.0038). The thesis's "standard" correlation-rank weights, while very stable,
latch onto a tightly-correlated but ER-irrelevant block (S100 genes) and miss the
biology — outcome-blind weighting is the wrong prior here.

## Secondary endpoint: PAM50 Basal-vs-rest (consistency check)

- 164 genes selected; held-out AUC 0.994; Jaccard 0.73.
- Recovers EGFR (freq 1.0), keratins KRT6A/KRT16, and the luminal program
  (SCGB2A2, PIP, TFF1) as negative predictors. ESR1 is the strongest univariate
  discriminator (AUC 0.959, basal tumors are predominantly ER-).
- KRT5/KRT14 were not selected — correctly so: their univariate AUCs for the
  genefu-derived basal calls are only 0.77/0.66 in this dataset. The method
  follows the data, not the textbook. (KRT17/FOXC1 are absent from the
  expression matrix after QC.)
- PAM50 signature overlap 3/14 (p = 0.19, n.s.) — expected: the labels themselves
  come from these genes, so this endpoint can only check consistency, not
  validate. The ER-IHC endpoint above is the real validation.

## Simulation cross-check

Rebuilt the thesis pilot scenario (n=100, 200 noise features, 5 runs x 25 perms)
in R 4.5.3 and ran it end-to-end: fusion methods recover true features far better
than LASSO (Early Fusion Std F1 0.96/Sens 0.93 vs LASSO F1 0.13/Sens 0.10),
matching the thesis's qualitative pattern. (The full 18-scenario generator and
`Thesis_Simulation_Final_Results.csv` were not among the uploaded files.)

## Recommendations for the thesis

1. **Replace the SKAT-based adaptive weights.** They are invalid for continuous
   data (genotype test misapplied). Use the univariate-signal adaptive weights
   (`tcga-brca-v3.R`), or at minimum the correlation-rank standard weights with
   the caveat that they are outcome-blind.
2. **Replace the max-difference lambda rule** with the FDR-targeted
   (knockoff+-style) rule; report the empirical FDR estimate.
3. **Do not claim formal knockoff FDR guarantees.** Row-permuted controls are a
   heuristic, not model-X knockoffs. Report empirical selection frequencies,
   stability, and signature-enrichment p-values instead.
4. **Add a correction note**: the simulation "adaptive" arms effectively ran with
   size-only weights because `SKATBinary` fails on continuous features (verified
   on runif-style simulated data too). Re-run those arms with the fixed weights
   or relabel them.
5. Keep ER IHC (independent label) as the primary real-data endpoint; report
   PAM50 Basal-vs-rest only as a consistency check with the circularity caveat.

## Files (~/workspace/thesis-validation/)

- `tcga-brca-v3.R` — **the fixed pipeline** (recommended): valid adaptive weights
  + FDR-targeted lambda + leakage-controlled evaluation.
- `tcga-brca-real.R` — V1 pipeline (thesis as-written, leakage-fixed) for comparison.
- `tcga-brca-basal-v3.R` — Basal-vs-rest secondary endpoint (V3 machinery).
- `tcga_brca_real_results/` — selected-gene CSVs (`selected_ER_IHC_v3_adaptive.csv`,
  high-confidence sets), overlap tables, `.rds` result objects, run logs.
- `DIAGNOSIS.md` — overnight debugging log with the SKAT evidence.
- `fusion_sim_extracted.R` — extracted thesis pilot simulation (runs as-is).

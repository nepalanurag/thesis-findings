# Methods note: validating the permutation-assisted group lasso on TCGA-BRCA

Companion to the thesis "Integrative Feature Selection for Multi-Modal Data via Permutation-Assisted Group Lasso" (Anurag Nepal, SFSU, May 2026). This note summarizes the real-data validation design, the correction to the thesis method, and the full results behind the figures on the site.

## Data and design

- TCGA-BRCA, primary tumors only, one sample per patient. RNA-seq (RSEM, via cBioPortal), top 1,500 variable genes; methylation (Illumina 450K) at probe level, probes mapped to genes only for interpretation.
- Primary endpoint: ER status by immunohistochemistry (IHC), 806 ER+, 237 ER-, an independent clinical label. Combined cohort with both modalities: n = 737.
- Leakage control: 70/30 stratified train/test splits; every learned operation (variance filtering, scaling, Spearman/Ward.D2 grouping, weights, selection, lambda) inside the training set; held-out test AUC on the test split. The 5-arm and stability studies use 5 fresh splits per arm.
- Biological yardsticks: positive controls (ESR1, PGR, FOXA1, GATA3, canonical ER-pathway factors) and hypergeometric enrichment against the PAM50 (Parker et al. 2009) and Oncotype DX (Paik et al. 2004) signatures, with the number of genes actually tested as background.
- Secondary endpoint: PAM50 Basal-like vs rest (labels derived from the same RNA), reported as a consistency check only.

## The method

Group lasso (Yuan and Lin 2006): penalize whole correlated groups together so a pathway enters or leaves as a unit. Groups from hierarchical clustering of the feature correlation matrix. This matches how genes are organized and is what makes the selection interpretable.

Permutation assist: augment the design matrix with row-permuted copies ("controls") of each feature, preserving marginal distributions while breaking the outcome link. Per permutation, choose the lambda maximizing original discoveries subject to (1 + #controls)/max(1, #originals) <= 0.10, a knockoff-plus style heuristic (Barber and Candes 2015). Keep features selected in at least half the permutations. Row-permuted controls are a heuristic, not model-X knockoffs; I report empirical frequencies, stability, and enrichment p-values, not formal FDR guarantees.

Early fusion: stack modalities into one matrix, one joint model. Late fusion: select within each modality, combine after.

## The SKAT correction

The thesis computed adaptive weights from SKAT p-values. SKAT (Wu et al. 2011) is a rare-variant association test for genotype counts {0,1,2}. On standardized continuous expression it returns p = 1 through its SNP quality filter ("ALL SNPs have either high missing rates or no-variation") or crashes ("object 'K' not found"); the thesis tryCatch mapped every failure to p = 1. I verified this on runif-style simulated data too. So the thesis's "adaptive" arms effectively ran with size-only weights, and on real data the largest groups, which hold the real biology, got the largest penalties (ESR1 group 4.91 vs median 0.22), driving control frequencies to 0. I replaced SKAT with the continuous-data analogue: per-group RMS of univariate logistic Wald z, w = sqrt(group size)/(1 + RMS z), normalized to mean 1. The thesis also used a max-difference lambda rule whose implied false discovery rate was 0.65 on real data; the FDR-targeted rule above replaces it.

## Results

### As-written vs fixed (RNA only, ER-IHC)

| | Genes selected | Held-out AUC | ESR1/PGR/FOXA1/GATA3 | PAM50 | Oncotype DX |
|---|---|---|---|---|---|
| V1 (thesis as written) | 302 | 0.974 | all freq 0 | none | none |
| V3 (fixed weights + FDR lambda) | 224 | 0.954 | all freq 1.0 | 3.8x, p = 3e-4 | 5.0x, p = 0.012 |
| V3 standard weights | 75 | 0.926 | all freq < 0.5 | none | none |

The standard weights latch onto an ER-irrelevant correlated block (S100 genes). Outcome-blind weighting is the wrong prior here. At the stringent 0.9 frequency cutoff the fixed adaptive set (152 genes) keeps all four controls with PAM50 p = 0.0015 and Oncotype DX p = 0.0038.

### Five-arm benchmark (RNA + methylation, ER-IHC, 5 splits)

| Arm | Selected (SD) | AUC (SD) | Jaccard | PAM50 sig | Oncotype sig | Controls |
|---|---|---|---|---|---|---|
| Baseline LASSO | 30 (12) | 0.931 (0.017) | 0.201 | 2/5 | 4/5 | ESR1 5/5; PGR 0/5; FOXA1 0/5; GATA3 2/5 |
| EF-Std | 23 (51) | 0.899 (0.000) | 0.600 | 0/5 | 0/5 | all missed |
| EF-Adap | 236 (67) | 0.930 (0.017) | 0.229 | 5/5 | 4/5 | ESR1 5/5; PGR 4/5; FOXA1 3/5; GATA3 3/5 |
| LF-Std | 107 (47) | 0.917 (0.018) | 0.363 | 0/5 | 0/5 | nearly all missed |
| LF-Adap | 498 (148) | 0.924 (0.015) | 0.327 | 5/5 | 2/5 | ESR1 5/5; PGR 1/5; FOXA1 5/5; GATA3 5/5 |

### Stability study (S1 consensus groups, S2 shrunk weights, S3 stability selection)

| Variant | Selected (SD) | AUC (SD) | Jaccard | PAM50 | Oncotype | Controls |
|---|---|---|---|---|---|---|
| EF-Adap current | 236 (67) | 0.930 (0.017) | 0.229 | 5/5 | 4/5 | ESR1 5/5; PGR 4/5; FOXA1 3/5; GATA3 3/5 |
| EF-Adap S1 | 206.6 (35.2) | 0.940 (0.008) | 0.482 | 5/5 | 5/5 | all 5/5 |
| EF-Adap S2 | 206.2 (35.7) | 0.940 (0.008) | 0.483 | 5/5 | 5/5 | all 5/5 |
| EF-Adap S3 | 122.8 (17.9) | 0.948 (0.006) | 0.680 | 5/5 | 5/5 | all 5/5 |
| LF-Adap current | 498 (148) | 0.924 (0.015) | 0.327 | 5/5 | 2/5 | ESR1 5/5; PGR 1/5; FOXA1 5/5; GATA3 5/5 |
| LF-Adap S1 | 437.2 (40.5) | 0.946 (0.012) | 0.535 | 5/5 | 2/5 | ESR1 5/5; PGR 5/5; FOXA1 4/5; GATA3 5/5 |
| LF-Adap S2 | 427.8 (46.2) | 0.945 (0.013) | 0.557 | 5/5 | 2/5 | ESR1 5/5; PGR 5/5; FOXA1 4/5; GATA3 5/5 |
| LF-Adap S3 | 206.8 (26.1) | 0.944 (0.016) | 0.588 | 5/5 | 1/5 | ESR1 2/5; PGR 1/5; FOXA1 2/5; GATA3 2/5 |

Within one fixed split the method is perfectly reproducible (Jaccard 1.00 across seeds), so the instability was sampling variation, which is what S1-S3 address. S3 costs about 5-6x compute per run; S1/S2 are free. Adopted: EF-Adap S3. Not adopted: LF-Adap S3, which gains stability but loses the positive controls in most splits.

Calibration: ESR1 alone predicts ER status at AUC 0.943 (SD 0.011); a four-gene baseline at 0.939 (0.011), over 10 splits. The claim here is selection quality (stability plus biology), not prediction.

### Secondary endpoint

PAM50 Basal-like vs rest: 164 genes, held-out AUC 0.994, Jaccard 0.73. Recovers EGFR (freq 1.0), KRT6A/KRT16, and luminal markers (SCGB2A2, PIP, TFF1) as negative predictors. KRT5/KRT14 correctly not selected (univariate AUCs 0.77/0.66 for these calls in this dataset).

### Simulation cross-check

Rebuilt the thesis pilot (n=100, 200 noise features, 5 runs x 25 perms) in a fresh R environment: Early Fusion Std F1 0.96 / sensitivity 0.93 vs LASSO F1 0.13 / sensitivity 0.10, matching the thesis's qualitative pattern.

## References

- Barber, R. F. and Candes, E. J. (2015). Controlling the false discovery rate via knockoffs. Annals of Statistics, 43(5), 2055-2085.
- Meinshausen, N. and Buhlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society B, 72(4), 417-473.
- Paik, S. et al. (2004). A multigene assay to predict recurrence of tamoxifen-treated, node-negative breast cancer. New England Journal of Medicine, 351, 2817-2826.
- Parker, J. S. et al. (2009). Supervised risk predictor of breast cancer based on intrinsic subtypes. Journal of Clinical Oncology, 27(8), 1160-1167.
- The Cancer Genome Atlas Network (2012). Comprehensive molecular portraits of human breast tumours. Nature, 490, 61-70.
- Wu, M. C. et al. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. American Journal of Human Genetics, 89(1), 82-93.
- Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society B, 68(1), 49-67.

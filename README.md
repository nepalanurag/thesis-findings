# Thesis findings: permutation-assisted group lasso, validated on TCGA-BRCA

This repo presents my M.S. thesis, "Integrative Feature Selection for Multi-Modal Data via Permutation-Assisted Group Lasso" (Anurag Nepal, San Francisco State University, May 2026), plus the real-data validation and stability studies I ran after defending.

## What the thesis is about

When each patient is measured on thousands of gene expression and methylation features, which small set of features actually matters? The thesis proposes a group lasso with a decoy mechanism: for every feature, add a row-permuted copy of it as a "control." A feature has to beat its own permuted decoys, repeatedly across permutations, to get selected. Early fusion stacks modalities into one joint model; late fusion selects within each modality. Adaptive penalty weights shrink for groups with strong univariate signal.

## Key results (all real, computed Sep 2026)

- **Simulation:** in the thesis pilot scenario (n=100, 200 noise features), fusion methods recovered true features far better than LASSO (Early Fusion F1 0.96, sensitivity 0.93, vs LASSO F1 0.13, sensitivity 0.10).
- **TCGA-BRCA, ER status by IHC (RNA-seq + methylation, n=737, 5 fresh splits per arm):** the stabilized early-fusion adaptive model (EF-Adap S3) selects 122.8 features on average with held-out AUC 0.948 (SD 0.006), selection stability Jaccard 0.680 (up from 0.229), recovers ESR1/PGR/FOXA1/GATA3 in all 5 splits, and is significantly enriched for the PAM50 and Oncotype DX signatures in all 5 splits.
- **Head-to-head vs baselines:** adaptive arms beat plain LASSO (AUC 0.931, Jaccard 0.201, misses most controls) and standard non-adaptive group lasso (EF-Std: AUC 0.899, Jaccard 0.600, but selects nothing biologically meaningful, all controls missed).
- **Honest failure reported:** late-fusion adaptive S3 gained stability (Jaccard 0.588) but lost the positive controls in most splits. Stability alone must not pick the winner.
- **A real bug found and fixed:** the thesis's SKAT-based adaptive weights were invalid for continuous expression data (SKAT is a rare-variant test for genotype counts; it silently returned p=1). Replaced with univariate Wald-z weights plus a false-discovery-rate targeted lambda rule. See REPORT.md.

## The site

`site/` is a hand-built static website (plain HTML/CSS, Plotly from CDN, no build step): key findings, background, methods, and an interactive TCGA-BRCA results page with the prediction-vs-reproducibility bubble plot and the 5-arm benchmark charts. Open `site/index.html` in a browser, or deploy `site/` to any static host (it is deployed on Vercel as anurag-thesis).

## What is where

- `site/` - the static website source
- `REPORT.md` - the longer methods note: design, pipeline, the SKAT correction, full result tables, references
- The underlying analysis scripts live in my private working directory; the report tables in `REPORT.md` reproduce their outputs exactly.

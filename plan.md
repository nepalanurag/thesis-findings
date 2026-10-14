# Paper plan: stable multi-omic feature selection

## What this paper is

Methods paper. The contribution is the method: permutation-assisted adaptive group lasso with finite-sample FDR control, tested on real multi-omics data. TCGA-BRCA is the testbed, not the point.

The story: the method worked across 18 simulation scenarios in my thesis. On real TCGA data the adaptive version was unstable across resamples. I figured out why and fixed it. Most papers skip that part. I'm keeping it.

## Results I have

- TCGA-BRCA, 737 primary tumors, ER status by IHC. Top 1500 variable RNA genes + top 1500 methylation probes.
- Five-arm comparison: baseline lasso gets AUC 0.931 but Jaccard 0.201. EF-Adap gets AUC 0.930, Jaccard 0.229.
- The instability is sampling variation, not the algorithm: within one fixed split the Jaccard is 1.000.
- Three stabilizers fix it. S3 (consensus groups + shrunk z weights + stability selection): AUC 0.948, Jaccard 0.680.
- Biology checks out: ESR1, PGR, FOXA1, GATA3 recovered 5/5. PAM50 and Oncotype DX enriched 5/5.
- Side note for the methods section: SKAT adaptive weights are invalid for continuous expression data. I replaced them with Wald-z weights. Stating it plainly so others don't repeat it.

## Document outline

1. Abstract, structured. Lead with the S3 numbers.
2. Introduction. Multi-omics feature selection; group lasso has no FDR guarantees; adaptive weights are fragile. ER status as the testbed because the biology is known.
3. Methods. Data and preprocessing; the thesis method in brief, thesis cited; the five-arm design (early vs late fusion, standard vs adaptive, lasso baseline); the SKAT-to-Wald-z change and why; S1/S2/S3; evaluation (held-out AUC, Jaccard, positive controls, signature enrichment against PAM50, Oncotype DX, MammaPrint).
4. Results. Five-arm table; the instability diagnosis; the stabilizer ablation; the biology.
5. Discussion. What each stabilizer buys; the SKAT caution; limits (single cohort, ER endpoint only, IHC as ground truth).
6. Conclusion. One paragraph.
7. Figures. Five-arm plot and stability bubble plot exist. Workflow diagram still to draw.
8. Tables. The summary CSVs, rounded.

## To do

- [ ] Pull the methods text from the thesis
- [ ] Confirm author list (me + Dr. He?)
- [ ] Draw the workflow diagram
- [ ] Write the discussion's biology interpretation
- [ ] bioRxiv first, then BMC Bioinformatics or Bioinformatics Advances

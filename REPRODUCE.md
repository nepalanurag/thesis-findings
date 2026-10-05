# Reproducing the TCGA-BRCA results

The frozen scripts behind every table in `REPORT.md` live in
`code/frozen_tcga/`, with a table-by-table mapping in
`code/frozen_tcga/MANIFEST.md`. This file documents the environment, the
data, and how to rerun everything.

## Environment

R (>= 4.3) with these packages:

- CRAN: `tidyverse`, `glmnet`, `grpreg`, `pROC`, `Matrix`, `patchwork`, `SKAT`
- Bioconductor: `genefu` (used only if available, for the basal endpoint;
  the scripts degrade gracefully without it)

No `renv.lock` is committed. Bioconductor packages are unpinned; if a
rerun produces different groupings or weights on a much newer Bioconductor,
that is the place to look. `curatedTCGAData` (version `2.1.1`) is pinned
only in the superseded notebook-era script
`code/tcga_brca_multimodal_selection.R`, which is now archival and not
part of the frozen validation pipeline.

Install, for example:

```r
install.packages(c("tidyverse", "glmnet", "grpreg", "pROC", "Matrix", "patchwork", "SKAT"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("genefu")  # optional; basal endpoint only
```

## Data

Raw TCGA files are not in this repo (they are ~410 MB). Download them from
cBioPortal into a `tcga_data/` directory alongside the scripts:

```bash
python3 dl_rna.py data_mrna_seq_v2_rsem.txt        # ~184 MB, batched API fetch; allow 30-60 min
python3 dl_clinical.py data_clinical_sample.txt    # ~0.5 MB
# methylation: data_methylation_hm450.txt (~226 MB) from the brca_tcga
#   study on cBioPortal (HM450 beta values, probes x samples)
# clinical patient file: data_clinical_patient.txt (~1.7 MB, ER status by IHC)
```

Expected layout (`DATA <- "tcga_data"` in the scripts):

```
tcga_data/
  data_mrna_seq_v2_rsem.txt        # RNA-seq V2 RSEM, genes x samples
  data_methylation_hm450.txt       # Illumina 450K methylation, probes x samples
  data_clinical_patient.txt        # ER status by IHC (patient-level)
  data_clinical_sample.txt
```

The scripts expect the current working directory to be
`code/frozen_tcga/` so that relative paths to `tcga_data/` and the
output directory `tcga_brca_real_results/` resolve.

## Rerunning (in order)

```bash
cd code/frozen_tcga

# 1. As-written vs fixed, RNA only (REPORT "As-written vs fixed")
Rscript tcga-brca-real.R        # V1: ~10 min for the selection step
Rscript tcga-brca-v3.R          # V3 adaptive + standard arms: ~10 min total

# 2. Five-arm benchmark (REPORT "Five-arm benchmark")
for i in 1 2 3 4 5; do
  Rscript tcga-brca-5arm.R $i   # one fresh 70/30 split per run id; tens of minutes each
done
Rscript aggregate-5arm.R        # builds five_arm_summary_table.csv + plot

# 3. Stability study (REPORT "Stability study")
Rscript stability_ef.R          # early fusion, S1-S3: ~2 h (resumable)
Rscript stability_lf.R          # late fusion, S1-S3: ~2 h (resumable)
Rscript aggregate_stability.R   # builds stability_summary_table.csv + report

# 4. Secondary endpoint (REPORT "Secondary endpoint")
Rscript tcga-brca-basal-v3.R    # ~5-10 min

# 5. Simulation cross-check (REPORT "Simulation cross-check")
Rscript fusion_sim_extracted.R  # no TCGA data needed
```

Both stability scripts are resumable: completed (run, variant) results in
`tcga_brca_real_results/stability_{ef,lf}_results.rds` are skipped if you
relaunch a killed job.

## Seed scheme

All seeds are fixed in the scripts, so a rerun on the same data and
package versions reproduces the runs bit-for-bit:

| Script | Seed scheme |
|---|---|
| `tcga-brca-real.R` | `set.seed(20260908)` |
| `tcga-brca-v3.R` | `set.seed(20260909)` |
| `tcga-brca-5arm.R` | `set.seed(20260909 + RUN_ID)`, RUN_ID = 1..5 on the command line |
| `stability_ef.R`, `stability_lf.R` | consensus groups `seed = 7`; per-run seeds via `set.seed(900 + r)` and `set.seed(seed)` inside `stability_lib.R` |
| `tcga-brca-basal-v3.R` | `set.seed(20260909)` |
| `tcga-brca-validation.R` | `set.seed(20260908)` at start; `set.seed(1000 + b)` in the stability bootstrap |
| `aggregate-*.R` | deterministic given the per-run `.rds` files |

The scripts also write `run_log*.txt` files into `tcga_brca_real_results/`
as they go; `code/frozen_tcga/logs/` holds the logs from the original
runs for comparison.

## What maps to what

See `code/frozen_tcga/MANIFEST.md` for the script -> REPORT table mapping.
If any mapping is still uncertain after reading it, it is marked
"to be confirmed" rather than guessed.

# CAFEH-UFP-DNAm

R code for a targeted analysis of DNA methylation and traffic-related ultrafine particle exposure in the Community Assessment of Freeway Exposure and Health (CAFEH) study. The primary analysis was restricted to 558 candidate CpG sites previously associated with cardiovascular health, as measured by the American Heart Association's Life's Essential 8 score. A secondary analysis examined epigenetic aging.

## Repository contents

| Script | Analysis |
|---|---|
| `candidate_CpG_limma_analysis.R` | Primary outcome: methylation at candidate CpG sites |
| `epigenetic_aging_analysis.R` | Secondary outcome: epigenetic age acceleration and pace of aging |

## Primary analysis: candidate CpG sites

The script `candidate_CpG_limma_analysis.R`:

- Compares DNA methylation (M-values) between high and low particle number concentration (PNC) exposure groups, with low exposure as the reference
- Fits two linear models using `limma` on the same complete-case sample:
  - **Model 1 (adjusted):** PNC + age + sex + education + BMI category + estimated blood cell proportions
  - **Model 2 (unadjusted):** PNC only
- Applies Benjamini-Hochberg false discovery rate (FDR) correction
- Produces results tables, Q-Q plots, volcano plots, an adjusted vs. unadjusted comparison, and a model summary

## Secondary analysis: epigenetic aging

The script `epigenetic_aging_analysis.R`:

- Links epigenetic clock output (keyed by array barcode) to participant data using the EPIC sample legend, and checks that age and sex agree between sources
- Analyzes five outcomes: age acceleration (AA) for the Horvath, Hannum, and PhenoAge clocks, Horvath intrinsic epigenetic age acceleration (IEAA), and DunedinPACE. AA is the residual from regressing each clock on chronological age
- Compares high vs. low TAA-PNC (reference = low) using linear regression in three models:
  - **Model 1:** age + sex
  - **Model 2:** Model 1 + race/ethnicity + smoking (ever vs. never) + BMI
  - **Model 3:** Model 2 + estimated CD8+ T cell, CD4+ T cell, NK cell, B cell, and monocyte proportions
- Handles missing covariates (BMI, smoking) with multiple imputation by chained equations (`mice`, predictive mean matching, 20 imputations), pooled with Rubin's rules
- Applies Benjamini-Hochberg FDR across the five outcomes within each model
- Runs a complete-case analysis and sensitivity analyses: continuous ln(TAA-PNC) per IQR, excluding age acceleration outliers (|z| > 3), excluding statin users, and additional adjustment for statin use or array slide
- Writes a clock quality check, results tables, a forest plot, and session information to `results_epigenetic_aging/`

## Requirements

R with the packages `data.table`, `limma`, `ggplot2`, `ggrepel`, `readxl`, and `mice`. Session information is saved automatically when each script runs.

## Input files (not included)

Primary analysis:

- `M_values_cprop_metadata.csv`: M-values, cell proportions, and covariates
- `CVD_906_CpGs.csv`: curated list of 609 candidate CpG IDs associated with Life's Essential 8 cardiovascular health (column `IlmnID`)

Secondary analysis:

- `CAFEH_epigenetic_clock_results.rds`: epigenetic clock estimates, age acceleration, and cell proportions per array
- `2025_07_30_EPIC_Tufts_SampleLegend_4.xlsx`: links array barcodes (BeadChip and position) to participant IDs
- `dna_dataset.csv`: exposure and covariate data

## Data availability

Individual-level data are not included in this repository. DNA methylation data are available from the NCBI Gene Expression Omnibus (GEO; accession number GSE353206). Other de-identified data may be available from the CAFEH principal investigator upon reasonable request, subject to institutional review board approval and a data use agreement.

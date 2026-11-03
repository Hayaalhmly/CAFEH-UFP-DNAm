# CAFEH-UFP-DNAm

R code for a targeted analysis of DNA methylation and traffic-related ultrafine particle exposure in the Community Assessment of Freeway Exposure and Health (CAFEH) study. Analyses were restricted to 558 candidate CpG sites previously associated with air pollution.

## Analysis overview

The script `candidate_CpG_limma_analysis.R`:

- Compares DNA methylation (M-values) between high and low particle number concentration (PNC) exposure groups, with low exposure as the reference
- Fits two linear models using `limma` on the same complete-case sample:
  - **Model 1 (adjusted):** PNC + age + sex + education + BMI category + estimated blood cell proportions
  - **Model 2 (unadjusted):** PNC only
- Applies Benjamini-Hochberg false discovery rate (FDR) correction
- Produces results tables, Q-Q plots, volcano plots, an adjusted vs. unadjusted comparison, and a model summary

## Requirements

R with the packages `data.table`, `limma`, `ggplot2`, and `ggrepel`. Session information is saved automatically when the script runs.

## Input files (not included)

- `M_values_cprop_metadata.csv`: M-values, cell proportions, and covariates
- `airpollution_CpGs.csv`: list of candidate CpG IDs

## Data availability

Individual-level data are not included in this repository. DNA methylation data have been submitted to the NCBI Gene Expression Omnibus (GEO; accession number forthcoming). Other de-identified data may be available from the CAFEH principal investigator upon reasonable request, subject to institutional review board approval and a data use agreement.

#!/usr/bin/env Rscript

# ============================================================
# TARGETED CARDIOVASCULAR HEALTH (LIFE'S ESSENTIAL 8) CpG ANALYSIS
# High vs Low PNC Exposure
#
# MODEL 1: Adjusted
#   PNC + age + sex + education + BMI + cell proportions
#
# MODEL 2: Unadjusted
#   PNC exposure only
#
# BOTH models use the SAME analytic sample.
#
# Analysis is restricted to CpGs associated with Life's Essential 8 cardiovascular health in:
# CVD_906_CpGs.csv
# ============================================================


# ============================================================
# LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})


# ============================================================
# USER SETTINGS
# ============================================================

in_file <- "M_values_cprop_metadata.csv"
candidate_file <- "CVD_906_CpGs.csv"

# Exposure
expo_var <- "pnc_exposure"

# Covariates
age_var <- "age_years"
sex_var <- "sex"
edu_var <- "edu_cat"
bmi_var <- "BMI_category"

# Cell proportions
cell_candidates <- c(
  "CD8T",
  "CD4T",
  "NK",
  "Bcell",
  "Mono",
  "Gran",
  "Neu",
  "Eos"
)

# Output prefix
out_prefix <- "HighLow_CVD_CpGs"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

stop_if_missing <- function(x) {
  
  if (!file.exists(x)) {
    stop("Missing file: ", x)
  }
  
}


# ------------------------------------------------------------
# Collapse duplicate CpG columns
# ------------------------------------------------------------

collapse_duplicate_cpg_cols <- function(dt) {
  
  cpg_names <- grep(
    "^cg",
    names(dt),
    value = TRUE
  )
  
  if (length(cpg_names) == 0) {
    stop("No CpG columns starting with 'cg' found.")
  }
  
  dup_cpg <- cpg_names[
    duplicated(cpg_names)
  ]
  
  if (length(dup_cpg) == 0) {
    return(dt)
  }
  
  dup_unique <- unique(dup_cpg)
  
  message(
    "Duplicated CpGs found: ",
    length(dup_unique)
  )
  
  for (nm in dup_unique) {
    
    idx <- which(
      names(dt) == nm
    )
    
    dt[
      ,
      (nm) := rowMeans(
        .SD,
        na.rm = TRUE
      ),
      .SDcols = idx
    ]
    
    dt[
      ,
      (idx[-1]) := NULL
    ]
  }
  
  return(dt)
}


# ------------------------------------------------------------
# Q-Q plot
# ------------------------------------------------------------

make_qqplot <- function(
    p,
    out_png,
    title
) {
  
  p <- p[
    is.finite(p) &
      !is.na(p) &
      p > 0 &
      p <= 1
  ]
  
  if (length(p) == 0) {
    warning("No valid P-values for Q-Q plot.")
    return(NULL)
  }
  
  p <- sort(p)
  
  expected <- -log10(
    ppoints(length(p))
  )
  
  observed <- -log10(p)
  
  png(
    out_png,
    width = 1200,
    height = 1200,
    res = 150
  )
  
  plot(
    expected,
    observed,
    pch = 16,
    cex = 0.7,
    xlab = "Expected -log10(P)",
    ylab = "Observed -log10(P)",
    main = title
  )
  
  abline(
    0,
    1,
    lwd = 2
  )
  
  dev.off()
}


# ------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------

make_volcano <- function(
    res,
    model_name,
    out_png
) {
  
  plot_res <- copy(res)
  
  # Avoid problems if a P-value is exactly zero
  plot_res[
    ,
    neglog10P := -log10(
      pmax(P.Value, .Machine$double.xmin)
    )
  ]
  
  # Significance categories
  plot_res[
    ,
    significance := fifelse(
      adj.P.Val < 0.05,
      "FDR < 0.05",
      fifelse(
        P.Value < 0.05,
        "Nominal P < 0.05",
        "Not significant"
      )
    )
  ]
  
  n_nominal <- sum(
    plot_res$P.Value < 0.05,
    na.rm = TRUE
  )
  
  n_fdr05 <- sum(
    plot_res$adj.P.Val < 0.05,
    na.rm = TRUE
  )
  
  n_fdr10 <- sum(
    plot_res$adj.P.Val < 0.10,
    na.rm = TRUE
  )
  
  
  p <- ggplot(
    plot_res,
    aes(
      x = logFC,
      y = neglog10P
    )
  ) +
    
    geom_point(
      aes(
        color = significance
      ),
      size = 2.5,
      alpha = 0.80
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    # Label nominally significant CpGs
    geom_text_repel(
      data = plot_res[
        P.Value < 0.05
      ],
      aes(
        label = CpG
      ),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0
    ) +
    
    scale_color_manual(
      values = c(
        "Not significant" = "grey70",
        "Nominal P < 0.05" = "#E69F00",
        "FDR < 0.05" = "#D55E00"
      )
    ) +
    
    labs(
      title = paste0(
        model_name,
        ": CVD CpGs"
      ),
      
      subtitle = paste0(
        nrow(plot_res),
        " CpGs tested | ",
        n_nominal,
        " P < 0.05 | ",
        n_fdr05,
        " FDR < 0.05 | ",
        n_fdr10,
        " FDR < 0.10"
      ),
      
      x = "M-value difference (High vs Low PNC)",
      
      y = expression(
        -log[10](italic(P))
      ),
      
      color = NULL
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      legend.position = "top",
      
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      
      plot.subtitle = element_text(
        size = 10
      )
    )
  
  
  print(p)
  
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
}


# ============================================================
# START
# ============================================================

message("")
message("==============================================")
message("TARGETED CARDIOVASCULAR HEALTH (LIFE'S ESSENTIAL 8) CpG ANALYSIS")
message("==============================================")


# ============================================================
# CHECK FILES
# ============================================================

stop_if_missing(in_file)
stop_if_missing(candidate_file)


# ============================================================
# READ DATA
# ============================================================

message(
  "Reading: ",
  in_file
)

dt <- fread(
  in_file
)


message(
  "Initial samples: ",
  nrow(dt)
)


# ============================================================
# SAMPLE ID
# ============================================================

if (!("sample_id" %in% names(dt))) {
  
  message(
    "sample_id not found; renaming first column."
  )
  
  setnames(
    dt,
    1,
    "sample_id"
  )
  
}


dt[
  ,
  sample_id := as.character(sample_id)
]


# ============================================================
# COLLAPSE DUPLICATE CpGs
# ============================================================

dt <- collapse_duplicate_cpg_cols(dt)


# ============================================================
# IDENTIFY CpGs IN METHYLATION MATRIX
# ============================================================

all_cpg_cols <- grep(
  "^cg",
  names(dt),
  value = TRUE
)


message(
  "Total CpGs in methylation matrix: ",
  length(all_cpg_cols)
)


# ============================================================
# READ CANDIDATE CpGs
# ============================================================

cand <- fread(
  candidate_file
)


# Assume first column contains CpG IDs
cand_ids <- as.character(
  cand[[1]]
)


# Remove blank/missing values
cand_ids <- cand_ids[
  !is.na(cand_ids) &
    cand_ids != ""
]


# Remove duplicates
cand_ids <- unique(
  cand_ids
)


message(
  "Unique CVD CpGs in candidate file: ",
  length(cand_ids)
)


# ============================================================
# KEEP ONLY CANDIDATE CpGs
# ============================================================

cpg_cols <- intersect(
  cand_ids,
  all_cpg_cols
)


message(
  "CVD CpGs found in methylation matrix: ",
  length(cpg_cols),
  " / ",
  length(cand_ids)
)


if (length(cpg_cols) == 0) {
  
  stop(
    "No CVD CpGs were found in the methylation matrix."
  )
  
}


# ============================================================
# SAVE MISSING CANDIDATE CpGs
# ============================================================

missing_cpgs <- setdiff(
  cand_ids,
  all_cpg_cols
)


message(
  "CVD CpGs missing from matrix: ",
  length(missing_cpgs)
)


if (length(missing_cpgs) > 0) {
  
  fwrite(
    data.table(
      CpG = missing_cpgs
    ),
    paste0(
      out_prefix,
      "_missing_candidate_CpGs.csv"
    )
  )
  
}


# ============================================================
# CHECK COVARIATES
# ============================================================

cell_vars <- intersect(
  cell_candidates,
  names(dt)
)


need_covars <- c(
  expo_var,
  age_var,
  sex_var,
  edu_var,
  bmi_var
)


missing_covars <- setdiff(
  need_covars,
  names(dt)
)


if (length(missing_covars) > 0) {
  
  stop(
    "Missing required covariates: ",
    paste(
      missing_covars,
      collapse = ", "
    )
  )
  
}


if (length(cell_vars) == 0) {
  
  stop(
    "No cell proportion columns found."
  )
  
}


message(
  "Cell covariates: ",
  paste(
    cell_vars,
    collapse = ", "
  )
)


# ============================================================
# BUILD COVARIATE DATA
# ============================================================

covar <- dt[
  ,
  c(
    "sample_id",
    expo_var,
    age_var,
    sex_var,
    edu_var,
    bmi_var,
    cell_vars
  ),
  with = FALSE
]


# ============================================================
# FORMAT VARIABLES
# ============================================================

covar[[expo_var]] <- factor(
  covar[[expo_var]]
)


if (nlevels(covar[[expo_var]]) != 2) {
  
  stop(
    "Exposure must have exactly two levels. Current levels: ",
    paste(
      levels(covar[[expo_var]]),
      collapse = ", "
    )
  )
  
}


# Set LOW as reference
lvl <- levels(
  covar[[expo_var]]
)


if (
  any(
    grepl(
      "^low$",
      lvl,
      ignore.case = TRUE
    )
  )
) {
  
  low_level <- lvl[
    grep(
      "^low$",
      lvl,
      ignore.case = TRUE
    )[1]
  ]
  
  covar[[expo_var]] <- relevel(
    covar[[expo_var]],
    ref = low_level
  )
  
}


message(
  "Reference exposure group: ",
  levels(covar[[expo_var]])[1]
)


message(
  "Comparison exposure group: ",
  levels(covar[[expo_var]])[2]
)


# Factor covariates
covar[[sex_var]] <- factor(
  covar[[sex_var]]
)

covar[[edu_var]] <- factor(
  covar[[edu_var]]
)

covar[[bmi_var]] <- factor(
  covar[[bmi_var]]
)


# ============================================================
# COMPLETE CASE ANALYTIC SAMPLE
#
# IMPORTANT:
# We use the SAME samples for adjusted and unadjusted models.
# ============================================================

complete_model <- complete.cases(
  covar
)


message(
  "Complete samples for adjusted covariates: ",
  sum(complete_model),
  " / ",
  nrow(covar)
)


if (!all(complete_model)) {
  
  message(
    "Removing ",
    sum(!complete_model),
    " samples with missing covariates."
  )
  
  covar <- covar[
    complete_model,
  ]
  
  dt <- dt[
    complete_model,
  ]
  
}


message(
  "Final analytic sample size for BOTH models: ",
  nrow(dt)
)


# ============================================================
# CREATE CANDIDATE CpG M-VALUE MATRIX
#
# CpGs x Samples
# ============================================================

Mmat <- t(
  as.matrix(
    dt[
      ,
      ..cpg_cols
    ]
  )
)


message(
  "CVD CpGs before M-value filtering: ",
  nrow(Mmat)
)


# ============================================================
# REMOVE CpGs WITH NON-FINITE VALUES
# ============================================================

finite_cpg <- apply(
  Mmat,
  1,
  function(x) {
    all(is.finite(x))
  }
)


if (!all(finite_cpg)) {
  
  message(
    "Removing ",
    sum(!finite_cpg),
    " CpGs with non-finite M-values."
  )
  
  Mmat <- Mmat[
    finite_cpg,
    ,
    drop = FALSE
  ]
  
}


message(
  "FINAL CVD CpGs entering BOTH models: ",
  nrow(Mmat)
)


# ============================================================
# MODEL 1
# ADJUSTED MODEL
# ============================================================

message("")
message("==============================================")
message("MODEL 1: ADJUSTED")
message("PNC + age + sex + education + BMI + cells")
message("==============================================")


design_formula_adj <- as.formula(
  paste0(
    "~ ",
    expo_var,
    " + ",
    age_var,
    " + ",
    sex_var,
    " + ",
    edu_var,
    " + ",
    bmi_var,
    " + ",
    paste(
      cell_vars,
      collapse = " + "
    )
  )
)


message(
  "Adjusted formula: ",
  deparse(design_formula_adj)
)


design_adj <- model.matrix(
  design_formula_adj,
  data = covar
)


# ============================================================
# RUN ADJUSTED LIMMA
# ============================================================

fit_adj <- lmFit(
  Mmat,
  design_adj
)


fit_adj <- eBayes(
  fit_adj
)


# Identify exposure coefficient
coef_adj <- grep(
  paste0(
    "^",
    expo_var
  ),
  colnames(design_adj),
  value = TRUE
)


if (length(coef_adj) != 1) {
  
  stop(
    "Could not identify adjusted exposure coefficient."
  )
  
}


message(
  "Adjusted exposure coefficient: ",
  coef_adj
)


# ============================================================
# EXTRACT ADJUSTED RESULTS
# ============================================================

res_adj <- topTable(
  fit_adj,
  coef = coef_adj,
  number = Inf,
  sort.by = "P",
  adjust.method = "BH"
)


res_adj$CpG <- rownames(
  res_adj
)


res_adj <- as.data.table(
  res_adj
)


setcolorder(
  res_adj,
  c(
    "CpG",
    setdiff(
      names(res_adj),
      "CpG"
    )
  )
)


# ============================================================
# ADJUSTED SIGNIFICANCE COUNTS
# ============================================================

adj_nominal <- sum(
  res_adj$P.Value < 0.05,
  na.rm = TRUE
)


adj_fdr05 <- sum(
  res_adj$adj.P.Val < 0.05,
  na.rm = TRUE
)


adj_fdr10 <- sum(
  res_adj$adj.P.Val < 0.10,
  na.rm = TRUE
)


message(
  "Adjusted: Nominal P < 0.05 = ",
  adj_nominal
)


message(
  "Adjusted: FDR < 0.05 = ",
  adj_fdr05
)


message(
  "Adjusted: FDR < 0.10 = ",
  adj_fdr10
)


# ============================================================
# SAVE ADJUSTED RESULTS
# ============================================================

fwrite(
  res_adj,
  paste0(
    out_prefix,
    "_ADJUSTED_ALL.csv"
  )
)


fwrite(
  res_adj[
    P.Value < 0.05
  ],
  paste0(
    out_prefix,
    "_ADJUSTED_nominal_P005.csv"
  )
)


fwrite(
  res_adj[
    adj.P.Val < 0.05
  ],
  paste0(
    out_prefix,
    "_ADJUSTED_FDR005.csv"
  )
)


fwrite(
  res_adj[
    adj.P.Val < 0.10
  ],
  paste0(
    out_prefix,
    "_ADJUSTED_FDR010.csv"
  )
)


# ============================================================
# ADJUSTED Q-Q PLOT
# ============================================================

make_qqplot(
  res_adj$P.Value,
  paste0(
    out_prefix,
    "_ADJUSTED_QQ.png"
  ),
  "Q-Q Plot: Adjusted Candidate CpG Model"
)


# ============================================================
# ADJUSTED VOLCANO
# ============================================================

make_volcano(
  res = res_adj,
  model_name = "Adjusted model",
  out_png = paste0(
    out_prefix,
    "_ADJUSTED_volcano.png"
  )
)


# ============================================================
# MODEL 2
# UNADJUSTED MODEL
# PNC EXPOSURE ONLY
# ============================================================

message("")
message("==============================================")
message("MODEL 2: UNADJUSTED")
message("PNC EXPOSURE ONLY")
message("==============================================")


design_formula_unadj <- as.formula(
  paste0(
    "~ ",
    expo_var
  )
)


message(
  "Unadjusted formula: ",
  deparse(design_formula_unadj)
)


design_unadj <- model.matrix(
  design_formula_unadj,
  data = covar
)


# ============================================================
# RUN UNADJUSTED LIMMA
# ============================================================

fit_unadj <- lmFit(
  Mmat,
  design_unadj
)


fit_unadj <- eBayes(
  fit_unadj
)


# Identify exposure coefficient
coef_unadj <- grep(
  paste0(
    "^",
    expo_var
  ),
  colnames(design_unadj),
  value = TRUE
)


if (length(coef_unadj) != 1) {
  
  stop(
    "Could not identify unadjusted exposure coefficient."
  )
  
}


message(
  "Unadjusted exposure coefficient: ",
  coef_unadj
)


# ============================================================
# EXTRACT UNADJUSTED RESULTS
# ============================================================

res_unadj <- topTable(
  fit_unadj,
  coef = coef_unadj,
  number = Inf,
  sort.by = "P",
  adjust.method = "BH"
)


res_unadj$CpG <- rownames(
  res_unadj
)


res_unadj <- as.data.table(
  res_unadj
)


setcolorder(
  res_unadj,
  c(
    "CpG",
    setdiff(
      names(res_unadj),
      "CpG"
    )
  )
)


# ============================================================
# UNADJUSTED SIGNIFICANCE COUNTS
# ============================================================

unadj_nominal <- sum(
  res_unadj$P.Value < 0.05,
  na.rm = TRUE
)


unadj_fdr05 <- sum(
  res_unadj$adj.P.Val < 0.05,
  na.rm = TRUE
)


unadj_fdr10 <- sum(
  res_unadj$adj.P.Val < 0.10,
  na.rm = TRUE
)


message(
  "Unadjusted: Nominal P < 0.05 = ",
  unadj_nominal
)


message(
  "Unadjusted: FDR < 0.05 = ",
  unadj_fdr05
)


message(
  "Unadjusted: FDR < 0.10 = ",
  unadj_fdr10
)


# ============================================================
# SAVE UNADJUSTED RESULTS
# ============================================================

fwrite(
  res_unadj,
  paste0(
    out_prefix,
    "_UNADJUSTED_ALL.csv"
  )
)


fwrite(
  res_unadj[
    P.Value < 0.05
  ],
  paste0(
    out_prefix,
    "_UNADJUSTED_nominal_P005.csv"
  )
)


fwrite(
  res_unadj[
    adj.P.Val < 0.05
  ],
  paste0(
    out_prefix,
    "_UNADJUSTED_FDR005.csv"
  )
)


fwrite(
  res_unadj[
    adj.P.Val < 0.10
  ],
  paste0(
    out_prefix,
    "_UNADJUSTED_FDR010.csv"
  )
)


# ============================================================
# UNADJUSTED Q-Q PLOT
# ============================================================

make_qqplot(
  res_unadj$P.Value,
  paste0(
    out_prefix,
    "_UNADJUSTED_QQ.png"
  ),
  "Q-Q Plot: Unadjusted Candidate CpG Model"
)


# ============================================================
# UNADJUSTED VOLCANO PLOT
# ============================================================

make_volcano(
  res = res_unadj,
  model_name = "Unadjusted model",
  out_png = paste0(
    out_prefix,
    "_UNADJUSTED_volcano.png"
  )
)


# ============================================================
# COMPARE ADJUSTED VS UNADJUSTED
# ============================================================

message("")
message("==============================================")
message("COMPARING ADJUSTED VS UNADJUSTED")
message("==============================================")


comparison <- merge(
  
  res_adj[
    ,
    .(
      CpG,
      logFC_adjusted = logFC,
      P_adjusted = P.Value,
      FDR_adjusted = adj.P.Val
    )
  ],
  
  res_unadj[
    ,
    .(
      CpG,
      logFC_unadjusted = logFC,
      P_unadjusted = P.Value,
      FDR_unadjusted = adj.P.Val
    )
  ],
  
  by = "CpG",
  all = TRUE
)


# ============================================================
# ADD COMPARISON FLAGS
# ============================================================

comparison[
  ,
  nominal_adjusted :=
    P_adjusted < 0.05
]


comparison[
  ,
  nominal_unadjusted :=
    P_unadjusted < 0.05
]


comparison[
  ,
  FDR05_adjusted :=
    FDR_adjusted < 0.05
]


comparison[
  ,
  FDR05_unadjusted :=
    FDR_unadjusted < 0.05
]


comparison[
  ,
  FDR10_adjusted :=
    FDR_adjusted < 0.10
]


comparison[
  ,
  FDR10_unadjusted :=
    FDR_unadjusted < 0.10
]


# ============================================================
# SAVE COMPARISON TABLE
# ============================================================

fwrite(
  comparison,
  paste0(
    out_prefix,
    "_ADJUSTED_vs_UNADJUSTED.csv"
  )
)


# ============================================================
# OVERLAP COUNTS
# ============================================================

nominal_both <- sum(
  comparison$nominal_adjusted &
    comparison$nominal_unadjusted,
  na.rm = TRUE
)


fdr05_both <- sum(
  comparison$FDR05_adjusted &
    comparison$FDR05_unadjusted,
  na.rm = TRUE
)


fdr10_both <- sum(
  comparison$FDR10_adjusted &
    comparison$FDR10_unadjusted,
  na.rm = TRUE
)


message(
  "Nominal P < 0.05 in BOTH models: ",
  nominal_both
)


message(
  "FDR < 0.05 in BOTH models: ",
  fdr05_both
)


message(
  "FDR < 0.10 in BOTH models: ",
  fdr10_both
)


# ============================================================
# SAVE SUMMARY
# ============================================================

summary_file <- paste0(
  out_prefix,
  "_MODEL_SUMMARY.txt"
)


sink(
  summary_file
)


cat(
  "TARGETED CARDIOVASCULAR HEALTH (LIFE'S ESSENTIAL 8) CpG ANALYSIS\n"
)


cat(
  "========================================\n\n"
)


cat(
  "Candidate file: ",
  candidate_file,
  "\n",
  sep = ""
)


cat(
  "CVD CpGs in candidate list: ",
  length(cand_ids),
  "\n",
  sep = ""
)


cat(
  "CVD CpGs found in methylation matrix: ",
  length(cpg_cols),
  "\n",
  sep = ""
)


cat(
  "CVD CpGs entering models: ",
  nrow(Mmat),
  "\n",
  sep = ""
)


cat(
  "Analytic sample size: ",
  nrow(dt),
  "\n\n",
  sep = ""
)


cat(
  "Exposure reference: ",
  levels(covar[[expo_var]])[1],
  "\n",
  sep = ""
)


cat(
  "Exposure comparison: ",
  levels(covar[[expo_var]])[2],
  "\n\n",
  sep = ""
)


cat(
  "MODEL 1: ADJUSTED\n"
)


cat(
  "PNC + age + sex + education + BMI + cell proportions\n"
)


cat(
  "Nominal P < 0.05: ",
  adj_nominal,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.05: ",
  adj_fdr05,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.10: ",
  adj_fdr10,
  "\n\n",
  sep = ""
)


cat(
  "MODEL 2: UNADJUSTED\n"
)


cat(
  "PNC exposure only\n"
)


cat(
  "Nominal P < 0.05: ",
  unadj_nominal,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.05: ",
  unadj_fdr05,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.10: ",
  unadj_fdr10,
  "\n\n",
  sep = ""
)


cat(
  "OVERLAP\n"
)


cat(
  "Nominal P < 0.05 in both: ",
  nominal_both,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.05 in both: ",
  fdr05_both,
  "\n",
  sep = ""
)


cat(
  "FDR < 0.10 in both: ",
  fdr10_both,
  "\n",
  sep = ""
)


sink()


# ============================================================
# SESSION INFORMATION
# ============================================================

sink(
  paste0(
    out_prefix,
    "_sessionInfo.txt"
  )
)


print(
  sessionInfo()
)


sink()


# ============================================================
# FINAL SCREEN OUTPUT
# ============================================================

message("")
message("==============================================")
message("ALL DONE")
message("==============================================")


message(
  "Samples: ",
  nrow(dt)
)


message(
  "CVD CpGs tested: ",
  nrow(Mmat)
)


message("")
message("ADJUSTED MODEL:")

message(
  "  Nominal P < 0.05: ",
  adj_nominal
)

message(
  "  FDR < 0.05: ",
  adj_fdr05
)

message(
  "  FDR < 0.10: ",
  adj_fdr10
)


message("")
message("UNADJUSTED MODEL:")

message(
  "  Nominal P < 0.05: ",
  unadj_nominal
)

message(
  "  FDR < 0.05: ",
  unadj_fdr05
)

message(
  "  FDR < 0.10: ",
  unadj_fdr10
)


message("")
message("OVERLAP:")

message(
  "  Nominal in both: ",
  nominal_both
)

message(
  "  FDR < 0.05 in both: ",
  fdr05_both
)

message(
  "  FDR < 0.10 in both: ",
  fdr10_both
)


message("==============================================")

# =============================================================================
# Epigenetic aging and traffic-related ultrafine particle exposure in CAFEH
# Secondary outcome analysis
#
# Compares epigenetic age acceleration (Horvath, Hannum, PhenoAge), Horvath
# intrinsic epigenetic age acceleration (IEAA), and DunedinPACE between high
# and low time-activity-adjusted particle number concentration (TAA-PNC)
# groups (reference = low).
#
# Models (linear regression):
#   Model 1: exposure + age + sex
#   Model 2: Model 1 + race/ethnicity + smoking (ever vs never) + BMI
#   Model 3: Model 2 + CD8T, CD4T, NK, B cell and monocyte proportions
# Missing covariates (BMI, smoking) are multiply imputed with mice
# (predictive mean matching, 20 imputations) and pooled with Rubin's rules.
# Benjamini-Hochberg FDR is applied across the five outcomes within each model.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(mice)
  library(ggplot2)
})

set.seed(2026)

# ---- Input files (not included in the repository) ---------------------------
clock_file  <- "CAFEH_epigenetic_clock_results.rds"            # clock output, keyed by array barcode
legend_file <- "2025_07_30_EPIC_Tufts_SampleLegend_4.xlsx"     # barcode -> participant ID
pheno_file  <- "dna_dataset.csv"                               # exposure and covariates
out_dir     <- "results_epigenetic_aging"
dir.create(out_dir, showWarnings = FALSE)

M_IMP <- 20   # number of imputations

# ---- 1. Load and link data ----------------------------------------------------
clocks <- as.data.table(readRDS(clock_file))

legend <- as.data.table(read_excel(legend_file))
legend[, SampleID := paste0(BeadChip, "_", Position)]
legend[, ParticipantID := as.integer(Sample_Name)]

clocks <- merge(clocks, legend[, .(SampleID, ParticipantID)], by = "SampleID", all.x = TRUE)
stopifnot(!anyNA(clocks$ParticipantID))

pheno <- fread(pheno_file, na.strings = c("NA", ""), colClasses = list(character = c("statin_use", "obese")))
d <- merge(clocks, pheno, by = "ParticipantID")

# Linkage checks: all samples matched, and age and sex agree between sources
stopifnot(nrow(d) == 96,
          all(d$Age == d$age_years),
          all(d$Female == as.integer(d$sex == "Female")))

# ---- 2. Variables ------------------------------------------------------------
d[, high       := as.integer(pnc_exposure == "High")]
d[, female     := as.integer(sex == "Female")]
d[, race       := relevel(factor(race), ref = "White")]
d[, ever_smoke := fifelse(is.na(smoking_status), NA_integer_,
                          as.integer(smoking_status != "Non-smoker"))]
d[, BMI        := as.numeric(BMI)]
d[, statin     := fcase(statin_use %in% c("1", "Yes"), 1L,
                        statin_use %in% c("0", "No"),  0L)]
d[, slide      := factor(sub("_.*", "", SampleID))]

iqr_ln_pnc <- IQR(d$ln_pnc)
d[, pnc_iqr := ln_pnc / iqr_ln_pnc]

outcomes <- c("Horvath AA"   = "AA_Horvath",
              "Hannum AA"    = "AA_Hannum",
              "PhenoAge AA"  = "AA_Pheno",
              "Horvath IEAA" = "IEAA_Horvath",
              "DunedinPACE"  = "PACE")
cells <- c("CD8T", "CD4T", "NK", "Bcell", "Mono")

# Outliers: |z| > 3 on any age acceleration measure
z <- d[, lapply(.SD, function(v) abs(as.numeric(scale(v))) > 3),
       .SDcols = c("AA_Horvath", "AA_Hannum", "AA_Pheno", "IEAA_Horvath")]
d[, outlier := as.integer(rowSums(z) > 0)]

cat("n =", nrow(d), "| High =", sum(d$high), "| Low =", sum(1 - d$high), "\n")
cat("IQR of ln(TAA-PNC):", round(iqr_ln_pnc, 3), "\n")
cat("Missing: BMI =", sum(is.na(d$BMI)), "| smoking =", sum(is.na(d$ever_smoke)), "\n")
cat("Outliers:", paste(d[outlier == 1, ParticipantID], collapse = ", "), "\n")

# ---- 3. Clock quality check ---------------------------------------------------
qc <- rbindlist(lapply(c(Horvath = "DNAmAge_Horvath", Hannum = "DNAmAge_Hannum", PhenoAge = "DNAmAge_Pheno"),
  function(k) data.table(r_with_age  = cor(d[[k]], d$Age),
                         mean_offset = mean(d[[k]] - d$Age),
                         median_abs_error = median(abs(d[[k]] - d$Age)))),
  idcol = "clock")
qc <- rbind(qc, data.table(clock = "DunedinPACE", r_with_age = cor(d$PACE, d$Age),
                           mean_offset = mean(d$PACE), median_abs_error = sd(d$PACE)))
fwrite(qc, file.path(out_dir, "clock_qc.csv"))
print(qc)

# ---- 4. Model formulas --------------------------------------------------------
m1 <- c("Age", "female")
m2 <- c(m1, "race", "ever_smoke", "BMI")
m3 <- c(m2, cells)
models <- list(M1 = m1, M2 = m2, M3 = m3)

fml <- function(y, x, covs) reformulate(c(x, covs), response = y)

# ---- 5. Multiple imputation ---------------------------------------------------
imp_vars <- c("Age", "female", "race", "high", "pnc_iqr", "statin", "BMI", "ever_smoke",
              unname(outcomes), cells, "outlier", "slide", "ParticipantID")
imp_dat <- as.data.frame(d[, ..imp_vars])

meth <- make.method(imp_dat)
meth[c("BMI", "ever_smoke")] <- "pmm"
pred <- make.predictorMatrix(imp_dat)
pred[, c("outlier", "slide", "ParticipantID")] <- 0   # carried along, not used as predictors

imp <- mice(imp_dat, m = M_IMP, maxit = 20, method = meth, predictorMatrix = pred,
            printFlag = FALSE, seed = 2026)

# Fit a model in each imputed dataset and pool with Rubin's rules
pool_fit <- function(f, x, keep = NULL) {
  fits <- lapply(seq_len(M_IMP), function(i) {
    dd <- complete(imp, i)
    if (!is.null(keep)) dd <- dd[keep(dd), ]
    lm(f, data = dd)
  })
  s <- summary(pool(as.mira(fits)), conf.int = TRUE)
  r <- s[grepl(paste0("^", x), s$term), ][1, ]
  n <- nrow(if (is.null(keep)) imp_dat else imp_dat[keep(imp_dat), ])
  data.table(est = r$estimate, lo = r$`2.5 %`, hi = r$`97.5 %`, p = r$p.value, n = n)
}

# Complete-case fit
cc_fit <- function(f, x, data = d) {
  m  <- lm(f, data = data)
  ci <- confint(m)
  data.table(est = coef(m)[[x]], lo = ci[x, 1], hi = ci[x, 2],
             p = summary(m)$coefficients[x, 4], n = nobs(m))
}

add_fdr <- function(res) res[, q := p.adjust(p, "BH"), by = model][]

# ---- 6. Main analysis (multiple imputation) ----------------------------------
res_mi <- rbindlist(lapply(names(models), function(mn)
  rbindlist(lapply(names(outcomes), function(o)
    cbind(data.table(model = mn, outcome = o),
          pool_fit(fml(outcomes[[o]], "high", models[[mn]]), "high"))))))
res_mi <- add_fdr(res_mi)
fwrite(res_mi, file.path(out_dir, "results_main_MI.csv"))

# ---- 7. Complete-case analysis -----------------------------------------------
res_cc <- rbindlist(lapply(names(models), function(mn)
  rbindlist(lapply(names(outcomes), function(o)
    cbind(data.table(model = mn, outcome = o),
          cc_fit(fml(outcomes[[o]], "high", models[[mn]]), "high"))))))
res_cc <- add_fdr(res_cc)
fwrite(res_cc, file.path(out_dir, "results_complete_case.csv"))

# ---- 8. Sensitivity analyses (Model 2 covariates, multiple imputation) -------
sens_list <- list(
  "Continuous ln(TAA-PNC), per IQR"      = list(x = "pnc_iqr", covs = m2,                 keep = NULL),
  "Excluding age acceleration outliers"  = list(x = "high",    covs = m2,                 keep = function(dd) dd$outlier == 0),
  "Excluding statin users"               = list(x = "high",    covs = m2,                 keep = function(dd) dd$statin != 1),
  "Additionally adjusted for statin use" = list(x = "high",    covs = c(m2, "statin"),    keep = NULL),
  "Additionally adjusted for array slide"= list(x = "high",    covs = c(m2, "slide"),     keep = NULL)
)
res_sens <- rbindlist(lapply(names(sens_list), function(a) {
  s <- sens_list[[a]]
  rbindlist(lapply(names(outcomes), function(o)
    cbind(data.table(analysis = a, outcome = o),
          pool_fit(fml(outcomes[[o]], s$x, s$covs), s$x, s$keep))))
}))
fwrite(res_sens, file.path(out_dir, "results_sensitivity.csv"))

# ---- 9. Forest plot (Models 2 and 3, multiple imputation) --------------------
plot_dat <- res_mi[model %in% c("M2", "M3")]
plot_dat[, model_lab := factor(fifelse(model == "M2", "Fully adjusted", "Fully adjusted + cell types"),
                               levels = c("Fully adjusted", "Fully adjusted + cell types"))]
plot_dat[, panel := factor(fifelse(outcome == "DunedinPACE",
                                   "Difference in pace of aging",
                                   "Difference in age acceleration (years)"),
                           levels = c("Difference in age acceleration (years)", "Difference in pace of aging"))]
lev <- rev(names(outcomes))
plot_dat[, y := match(outcome, lev) + fifelse(model == "M2", -0.15, 0.15)]

p <- ggplot(plot_dat, aes(x = est, y = y, colour = model_lab, shape = model_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.8) +
  geom_point(size = 2.6) +
  geom_blank(data = data.frame(panel = factor("Difference in pace of aging", levels = levels(plot_dat$panel)),
                               est = 0, y = c(0.4, 1.6), model_lab = plot_dat$model_lab[1])) +
  scale_y_continuous(breaks = seq_along(lev), labels = lev) +
  facet_wrap(~ panel, scales = "free", nrow = 1) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834")) +
  scale_shape_manual(values = c(16, 15)) +
  labs(x = "High vs Low TAA-PNC", y = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(size = 9))

ggsave(file.path(out_dir, "Figure_epigenetic_age_forest.pdf"),  p, width = 7.2, height = 3.4)
ggsave(file.path(out_dir, "Figure_epigenetic_age_forest.tiff"), p, width = 7.2, height = 3.4,
       dpi = 600, compression = "lzw")

# ---- 10. Session information --------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Done. Results written to", out_dir, "\n")

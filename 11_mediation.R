# 11_mediation_WL.R
#
# Does immune exclusion mediate the effect of risk group on overall survival? Cox attenuation check and causal mediation analysis.
# Input: results1/survival_with_drug_sensitivity.rds  |  Output: results1/mediation_exclusion_result.rds

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("mediation", quietly = TRUE)) install.packages("mediation")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(survival)
library(mediation)
library(dplyr)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# LOAD DATA
survival_drug <- readRDS("results1/survival_with_drug_sensitivity.rds")

stopifnot(all(c("OS", "OS.time", "risk_group", "Exclusion", "Subtype_Selected") %in% colnames(survival_drug)))
cat("[CHECK] n =", nrow(survival_drug), "\n")

n_events <- sum(survival_drug$OS == 1, na.rm = TRUE)
cat("[CHECK] OS events:", n_events, "/", sum(!is.na(survival_drug$OS)), "\n")
if (n_events < 50) warning("Fewer than 50 events -- survival/mediation estimates may be unstable")

med_data <- survival_drug %>%
  filter(!is.na(OS), !is.na(OS.time), !is.na(Exclusion), !is.na(Subtype_Selected)) %>%
  mutate(risk_group = factor(risk_group, levels = c("Low", "High")))

# SANITY CHECK: complete cases
cat("[CHECK] Complete cases for mediation:", nrow(med_data), "\n")
if (nrow(med_data) < 500) warning("More missingness than expected -- verify upstream merges")

# CHECKPOINT
saveRDS(med_data, "results1/mediation_analysis_data.rds")
cat("[SAVED] results1/mediation_analysis_data.rds\n\n")

# SEQUENTIAL COX ATTENUATION CHECK (simple, standard first look)
# Logic: if Exclusion mediates risk_group's effect on survival, adding
# Exclusion to the model should shrink risk_group's log(HR) toward zero.
med_data <- readRDS("results1/mediation_analysis_data.rds")

cox_unadjusted      <- coxph(Surv(OS.time, OS) ~ risk_group, data = med_data)
cox_mediator_added  <- coxph(Surv(OS.time, OS) ~ risk_group + Exclusion, data = med_data)

cat("\n=== Model A: risk_group only ===\n")
print(summary(cox_unadjusted)$coefficients)
cat("\n=== Model B: risk_group + Exclusion ===\n")
print(summary(cox_mediator_added)$coefficients)

logHR_A <- coef(cox_unadjusted)["risk_groupHigh"]
logHR_B <- coef(cox_mediator_added)["risk_groupHigh"]
pct_attenuation <- 100 * (logHR_A - logHR_B) / logHR_A

cat("\n[CHECK] log(HR) attenuation after adding Exclusion:", round(pct_attenuation, 1), "%\n")
cat("[CHECK] Interpretation: a LARGE POSITIVE attenuation (e.g. >20%) would\n")
cat("        suggest mediation; a small/negative value (as found here) means\n")
cat("        Exclusion does not explain risk_group's survival effect.\n")

# SUBTYPE-ADJUSTED COX MODEL FOR EXCLUSION (this is the correct
cox_full <- coxph(Surv(OS.time, OS) ~ risk_group + Exclusion + Subtype_Selected,
                  data = med_data)

cat("\n=== Model C: risk_group + Exclusion + subtype (subtype-adjusted) ===\n")
print(summary(cox_full)$coefficients)

cat("\n[VERIFY] Exclusion p-value, subtype-adjusted Cox model:",
    summary(cox_full)$coefficients["Exclusion", "Pr(>|z|)"], "\n")

saveRDS(cox_full, "results1/cox_full_subtype_adjusted.rds")
cox_full <- readRDS("results1/cox_full_subtype_adjusted.rds")   # checkpoint reload

# FORMAL CAUSAL MEDIATION ANALYSIS
# Mediator model: does risk_group predict Exclusion (subtype-adjusted)?
# Outcome model: parametric AFT survival model (mediate() does not support
# coxph directly; survreg is the standard substitute for survival mediation).
# NOTE: survreg's AFT parameterization has the OPPOSITE sign convention from
# Cox HRs -- do not directly compare coefficient signs across Step 2 and 3.
mediator_model <- lm(Exclusion ~ risk_group + Subtype_Selected, data = med_data)
outcome_model  <- survreg(Surv(OS.time, OS) ~ risk_group + Exclusion + Subtype_Selected,
                          data = med_data, dist = "weibull")

cat("\n=== Mediator model: risk_group -> Exclusion (subtype-adjusted) ===\n")
print(summary(mediator_model)$coefficients["risk_groupHigh", ])

cat("\n=== Outcome model: survival ~ risk_group + Exclusion + subtype ===\n")
print(summary(outcome_model)$table)

mediator_p <- summary(mediator_model)$coefficients["risk_groupHigh", "Pr(>|t|)"]
if (mediator_p > 0.001) warning("Mediator model p-value looks weaker than script 13's original Exclusion finding -- verify data consistency")

med_result <- mediate(
  mediator_model, outcome_model,
  treat = "risk_group", mediator = "Exclusion",
  sims = 1000
)

cat("\n=== Formal mediation results ===\n")
print(summary(med_result))
cat("\n[CHECK] ACME = mediated/indirect effect; ADE = direct effect;\n")
cat("        Prop. Mediated = fraction of total effect explained by Exclusion.\n")
cat("        A non-significant ACME with CI spanning zero means NO mediation.\n")

# SAVE
saveRDS(med_result, "results1/mediation_exclusion_result.rds")
cat("\n[SAVED] results1/mediation_exclusion_result.rds\n")

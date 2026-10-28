# 07_tmb_hrd_WL.R
#
# Molecular subtype (TCGA PanCancer Atlas), TMB (GDC MAF) and HRD (Marquard et al. 2015) by risk group, with subtype-adjusted and within-subtype tests.
# Input: results1/full_survival_with_risk.rds  |  Output: results1/risk_subtype_tmb_hrd.rds

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(TCGAbiolinks)
library(dplyr)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

survival_patient <- readRDS("results1/full_survival_with_risk.rds")   # n=541
stopifnot(nrow(survival_patient) == 541)

# PULL TCGA-UCEC MOLECULAR SUBTYPES (Levine et al. 2013 / PanCancerAtlas)
all_subtypes <- PanCancerAtlas_subtypes()

ucec_subtypes <- all_subtypes %>%
  filter(cancer.type == "UCEC") %>%
  mutate(tcga_patient = substr(pan.samplesID, 1, 12)) %>%
  dplyr::select(tcga_patient, Subtype_Selected)

stopifnot(sum(duplicated(ucec_subtypes$tcga_patient)) == 0)

risk_subtype <- merge(survival_patient, ucec_subtypes,
                      by.x = "patient_id", by.y = "tcga_patient", all.x = TRUE)

# SANITY CHECK 1: subtype merge
stopifnot(nrow(risk_subtype) == nrow(survival_patient))
n_subtype <- sum(!is.na(risk_subtype$Subtype_Selected))
cat("[CHECK] Patients with subtype call:", n_subtype, "/", nrow(risk_subtype), "\n")
if (n_subtype < 500) warning("Fewer subtype matches than expected (~532) -- verify barcode format")

# Foundational confound check -- expect highly significant association
subtype_chisq <- chisq.test(table(risk_subtype$Subtype_Selected, risk_subtype$risk_group))
print(subtype_chisq)
if (subtype_chisq$p.value > 0.05) warning("Subtype-risk association not significant -- unexpected, re-verify upstream data")

saveRDS(risk_subtype, "results1/checkpoint_11a_risk_subtype.rds")
cat("[SAVED] results1/checkpoint_11a_risk_subtype.rds\n\n")

# TMB -- VIA GDC ENSEMBLE MAF
# (GDCquery_Maf() is deprecated; use base GDCquery/GDCdownload/GDCprepare)

# Reload checkpoint explicitly (do not rely on in-memory risk_subtype,
# in case this section is re-run independently after an environment clear)
risk_subtype <- readRDS("results1/checkpoint_11a_risk_subtype.rds")

if (!file.exists("results/ucec_maf.rds")) {
  maf_query <- GDCquery(
    project      = "TCGA-UCEC",
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(maf_query, method = "api", directory = "GDCdata_maf")
  ucec_maf <- GDCprepare(maf_query, directory = "GDCdata_maf")
  saveRDS(ucec_maf, "results/ucec_maf.rds")
}

# SANITY CHECK 2: MAF pulled correctly
ucec_maf <- readRDS("results/ucec_maf.rds")
stopifnot(nrow(ucec_maf) > 0)
cat("[CHECK] MAF rows:", nrow(ucec_maf), "| unique samples:", length(unique(ucec_maf$Tumor_Sample_Barcode)), "\n")

nonsyn_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"
)

tmb_data <- ucec_maf %>%
  filter(Variant_Classification %in% nonsyn_classes) %>%
  mutate(tcga_patient = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  group_by(tcga_patient) %>%
  summarise(n_mutations = n(), .groups = "drop") %>%
  mutate(TMB = n_mutations / 38)   # mutations per Mb, 38 Mb standard exome size

risk_subtype_tmb <- merge(risk_subtype, tmb_data, by.x = "patient_id", by.y = "tcga_patient", all.x = TRUE)

# SANITY CHECK 3: TMB merge
stopifnot(nrow(risk_subtype_tmb) == nrow(risk_subtype))
n_tmb <- sum(!is.na(risk_subtype_tmb$TMB))
cat("[CHECK] Patients with TMB:", n_tmb, "\n")
if (n_tmb < 480) warning("Fewer TMB matches than expected (~511) -- verify MAF barcode overlap")

saveRDS(risk_subtype_tmb, "results1/checkpoint_11b_risk_subtype_tmb.rds")
cat("[SAVED] results1/checkpoint_11b_risk_subtype_tmb.rds\n\n")

# HRD -- MARQUARD ET AL. 2015 PER-SAMPLE GENOMIC SCAR SCORES

# Reload checkpoint explicitly
risk_subtype_tmb <- readRDS("results1/checkpoint_11b_risk_subtype_tmb.rds")

hrd_url <- "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs40364-015-0033-4/MediaObjects/40364_2015_33_MOESM8_ESM.txt"
if (!file.exists("data/marquard_hrd_scores.txt")) {
  dir.create("data", showWarnings = FALSE)
  download.file(hrd_url, destfile = "data/marquard_hrd_scores.txt", mode = "wb")
}
hrd_raw <- read.delim("data/marquard_hrd_scores.txt", stringsAsFactors = FALSE)

# SANITY CHECK 4: HRD source table
stopifnot(nrow(hrd_raw) == 5371)   # known size of Marquard et al. pan-cancer table
cat("[CHECK] Marquard HRD table rows:", nrow(hrd_raw), "(expect 5371)\n")

hrd_data <- hrd_raw %>%
  rename(tcga_patient = Tumor, HRD_LOH = HRD.LOH) %>%
  mutate(HRD_score = NtAI + LST + HRD_LOH) %>%
  dplyr::select(tcga_patient, NtAI, LST, HRD_LOH, HRD_score)

risk_full <- merge(risk_subtype_tmb, hrd_data, by.x = "patient_id", by.y = "tcga_patient", all.x = TRUE)

# SANITY CHECK 5: HRD merge
stopifnot(nrow(risk_full) == nrow(risk_subtype_tmb))
n_hrd <- sum(!is.na(risk_full$HRD_score))
cat("[CHECK] Patients with HRD:", n_hrd, "\n")
if (n_hrd < 380) warning("Fewer HRD matches than expected (~410) -- verify UCEC coverage in Marquard table")

# SUBTYPE-ADJUSTED MODELS -- THE KEY ANALYSIS
risk_full$logTMB <- log1p(risk_full$TMB)

tmb_unadj <- wilcox.test(TMB ~ risk_group, data = risk_full)
tmb_adj   <- lm(logTMB ~ risk_group + Subtype_Selected, data = risk_full)

hrd_unadj <- wilcox.test(HRD_score ~ risk_group, data = risk_full)
hrd_adj   <- lm(HRD_score ~ risk_group + Subtype_Selected, data = risk_full)

cat("\n=== TMB: unadjusted Wilcoxon ===\n"); print(tmb_unadj)
cat("\n=== TMB: subtype-adjusted model ===\n"); print(summary(tmb_adj))
cat("\n=== HRD: unadjusted Wilcoxon ===\n"); print(hrd_unadj)
cat("\n=== HRD: subtype-adjusted model ===\n"); print(summary(hrd_adj))

# WITHIN-SUBTYPE STRATIFIED TESTS (this is the source of the
subtypes <- unique(na.omit(risk_full$Subtype_Selected))

stratified_results <- data.frame(
  Subtype = character(), n = integer(),
  TMB_p = numeric(), HRD_p = numeric()
)

for (s in subtypes) {
  sub_data <- risk_full[risk_full$Subtype_Selected == s & !is.na(risk_full$Subtype_Selected), ]
  tmb_p <- tryCatch(wilcox.test(TMB ~ risk_group, data = sub_data)$p.value, error = function(e) NA)
  hrd_p <- tryCatch(wilcox.test(HRD_score ~ risk_group, data = sub_data)$p.value, error = function(e) NA)
  stratified_results <- rbind(stratified_results, data.frame(
    Subtype = s, n = nrow(sub_data), TMB_p = tmb_p, HRD_p = hrd_p
  ))
}

print(stratified_results)
cat("[VERIFY] TMB p range across strata:", range(stratified_results$TMB_p, na.rm = TRUE), "\n")
cat("[VERIFY] HRD p range across strata:", range(stratified_results$HRD_p, na.rm = TRUE), "\n")

saveRDS(stratified_results, "results1/tmb_hrd_stratified_by_subtype.rds")

hrd_interaction <- lm(HRD_score ~ risk_group * Subtype_Selected, data = risk_full)

cat("\n=== HRD: risk_group * Subtype_Selected interaction model ===\n")
print(summary(hrd_interaction))

# Formal model comparison: does adding the interaction term significantly
# improve fit over the main-effects-only model (hrd_adj)? This is the actual
# test of "is this a real interaction or a stratum-level fluke" -- an F-test
# on nested models, not just eyeballing one stratum's p-value.
cat("\n=== HRD: ANOVA -- main-effects model vs. interaction model ===\n")
hrd_interaction_test <- anova(hrd_adj, hrd_interaction)
print(hrd_interaction_test)

hrd_interaction_coefs <- summary(hrd_interaction)$coefficients
cn_low_interaction_row <- grep("risk_groupHigh:Subtype_SelectedUCEC.CN_LOW",
                               rownames(hrd_interaction_coefs), value = TRUE)
if (length(cn_low_interaction_row) == 1) {
  cat("\n[VERIFY] risk_group x CN_LOW interaction term:\n")
  print(hrd_interaction_coefs[cn_low_interaction_row, , drop = FALSE])
} else {
  warning("Could not locate the risk_group x CN_LOW interaction coefficient by name -- check factor level labels in the printed summary above.")
}

saveRDS(list(hrd_interaction_model = hrd_interaction,
             anova_vs_main_effects = hrd_interaction_test),
        "results1/hrd_interaction_model.rds")
cat("[SAVED] results1/hrd_interaction_model.rds\n\n")

# SAVE
saveRDS(risk_full, "results1/risk_subtype_tmb_hrd.rds")
risk_full <- readRDS("results1/risk_subtype_tmb_hrd.rds")   # checkpoint reload
cat("\nSaved: results1/risk_subtype_tmb_hrd.rds  (", nrow(risk_full), "patients )\n")

# 13_cptac_external_validation_WL.R
#
# External replication in CPTAC-UCEC: risk score reconstruction, TME scores and checkpoint genes by risk group, and TCGA vs CPTAC subtype composition.
# Input: ucec_cptac_2020/ucec_cptac_2020/ (cBioPortal download), results1/lasso_coefficients_dataframe.rds, results1/risk_score_cutoff.rds  |  Output: cptac_validation/

library(dplyr)

set.seed(123)

# LOAD CPTAC DATA (assumes manual download + extraction already done)
cptac_dir <- "ucec_cptac_2020/ucec_cptac_2020"   # nested folder from manual extraction
stopifnot(dir.exists(cptac_dir))

cptac_expr_raw      <- read.delim(file.path(cptac_dir, "data_mrna_seq_rsem.txt"), check.names = FALSE)
cptac_clin_patient  <- read.delim(file.path(cptac_dir, "data_clinical_patient.txt"), comment.char = "#", check.names = FALSE)
cptac_clin_sample   <- read.delim(file.path(cptac_dir, "data_clinical_sample.txt"),  comment.char = "#", check.names = FALSE)

# SANITY CHECK: basic shapes
stopifnot(ncol(cptac_expr_raw) > 90, nrow(cptac_clin_sample) > 80)
cat("[CHECK] Expression matrix:", nrow(cptac_expr_raw), "genes x", ncol(cptac_expr_raw) - 2, "patients\n")
cat("[CHECK] Clinical sample rows:", nrow(cptac_clin_sample), "\n")

# SANITY CHECK: confirm no survival data exists (documented scope decision,
survival_cols <- grep("SURVIVAL|STATUS|OS_|VITAL|FOLLOW|DAYS|DEATH|ALIVE|DECEASED",
                      c(colnames(cptac_clin_patient), colnames(cptac_clin_sample)),
                      value = TRUE, ignore.case = TRUE)
cat("[CHECK] Survival-related columns found:", if(length(survival_cols)==0) "NONE (expected)" else paste(survival_cols, collapse=", "), "\n")

dir.create("cptac_validation", showWarnings = FALSE)
saveRDS(cptac_expr_raw, "cptac_validation/cptac_expr_raw.rds")
saveRDS(cptac_clin_patient, "cptac_validation/cptac_clin_patient.rds")
saveRDS(cptac_clin_sample, "cptac_validation/cptac_clin_sample.rds")
cat("[SAVED] cptac_validation/cptac_expr_raw.rds (+ clinical files)\n\n")

cptac_expr_raw <- readRDS("cptac_validation/cptac_expr_raw.rds")
cptac_clin_patient <- readRDS("cptac_validation/cptac_clin_patient.rds")
cptac_clin_sample <- readRDS("cptac_validation/cptac_clin_sample.rds")

sig_genes <- c("TIMP3", "PTGIS", "PPP1R16B", "DIO2", "CTSW", "ZNF683", "GZMK", "IHH", "SPN")
checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2", "CD27", "ICOS", "IDO1")

missing_sig <- setdiff(sig_genes, cptac_expr_raw$Hugo_Symbol)
missing_checkpoint <- setdiff(checkpoint_genes, cptac_expr_raw$Hugo_Symbol)

cat("[CHECK] Signature genes (9 total) missing from CPTAC:",
    if (length(missing_sig) == 0) "NONE" else paste(missing_sig, collapse = ", "), "\n")
cat("[CHECK]", length(sig_genes) - length(missing_sig), "/", length(sig_genes),
    "signature genes available in CPTAC\n")
cat("[CHECK] Missing checkpoint genes:", if(length(missing_checkpoint)==0) "NONE" else paste(missing_checkpoint, collapse=", "), "\n")
if (length(missing_sig) > 2) {
  warning("More than 2 of 9 signature genes missing from CPTAC -- risk score reconstruction may be unreliable; inspect missing_sig before proceeding")
}

# RECOMPUTE RISK SCORE (using whichever signature genes are available
# in CPTAC's RSEM matrix -- count and identity determined above, not assumed)
lasso_coef_df <- readRDS("results1/lasso_coefficients_dataframe.rds")
cutoff <- readRDS("results1/risk_score_cutoff.rds")

available_genes <- intersect(lasso_coef_df$Gene, cptac_expr_raw$Hugo_Symbol)
cat("[CHECK]", length(available_genes), "/", nrow(lasso_coef_df),
    "signature genes available for CPTAC risk-score reconstruction:",
    paste(available_genes, collapse = ", "), "\n")
if (length(available_genes) < 7) {
  warning("Fewer than 7 signature genes available in CPTAC -- risk score reconstruction may be unreliable")
}

coef_subset <- lasso_coef_df %>% filter(Gene %in% available_genes)
coef_subset <- coef_subset[match(available_genes, coef_subset$Gene), ]
coef_vector <- coef_subset$Coefficient
names(coef_vector) <- coef_subset$Gene
cat("[CHECK] Coefficients used (", length(available_genes), "/", nrow(lasso_coef_df),
    " genes available):\n", sep = "")
print(coef_vector)

cptac_expr_mat <- cptac_expr_raw %>%
  filter(Hugo_Symbol %in% available_genes) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
rownames(cptac_expr_mat) <- cptac_expr_mat$Hugo_Symbol
cptac_expr_mat <- as.matrix(cptac_expr_mat %>% dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id))
cptac_expr_mat <- cptac_expr_mat[available_genes, ]
stopifnot(all(rownames(cptac_expr_mat) == names(coef_vector)))

risk_score_cptac <- as.vector(t(cptac_expr_mat) %*% coef_vector)
names(risk_score_cptac) <- colnames(cptac_expr_mat)

cptac_risk_df <- data.frame(patient_id = names(risk_score_cptac), risk_score = risk_score_cptac)

# DEFINE RISK GROUPS -- BOTH APPROACHES, MEDIAN SPLIT AS PRIMARY
cptac_risk_df$risk_group_tcga_cutoff <- factor(
  ifelse(cptac_risk_df$risk_score >= cutoff, "High", "Low"), levels = c("Low", "High"))
cptac_median <- median(cptac_risk_df$risk_score)
cptac_risk_df$risk_group_median_split <- factor(
  ifelse(cptac_risk_df$risk_score >= cptac_median, "High", "Low"), levels = c("Low", "High"))

cat("\n[CHECK] Group sizes, TCGA cutoff (secondary check):\n"); print(table(cptac_risk_df$risk_group_tcga_cutoff))
cat("[CHECK] Group sizes, CPTAC median split (PRIMARY):\n");     print(table(cptac_risk_df$risk_group_median_split))
if (min(table(cptac_risk_df$risk_group_tcga_cutoff)) < 10) {
  cat("[NOTE] TCGA-cutoff split is unbalanced as expected -- median split used as primary throughout.\n")
}

saveRDS(cptac_risk_df, "cptac_validation/cptac_risk_scores.rds")
cat("[SAVED] cptac_validation/cptac_risk_scores.rds\n\n")

# MERGE WITH CPTAC'S PRECOMPUTED TME/SUBTYPE DATA
# (independently-derived EPIC/CIBERSORT/ESTIMATE scores + GENOMICS_SUBTYPE --
# NOT recomputed from our own pipeline, which strengthens this as a genuine
# external check rather than a re-run of our own method on new data)
cptac_risk_df <- readRDS("cptac_validation/cptac_risk_scores.rds")

cptac_merged <- merge(cptac_risk_df, cptac_clin_sample, by.x = "patient_id", by.y = "SAMPLE_ID")
stopifnot(nrow(cptac_merged) > 80)
cat("[CHECK] Merged CPTAC cohort n =", nrow(cptac_merged), "\n")

cptac_merged$risk_group <- cptac_merged$risk_group_median_split
cptac_merged$GENOMICS_SUBTYPE <- factor(cptac_merged$GENOMICS_SUBTYPE)
cat("[CHECK] GENOMICS_SUBTYPE distribution:\n"); print(table(cptac_merged$GENOMICS_SUBTYPE, useNA = "always"))

# TEST TME SCORES BY RISK GROUP
run_risk_group_test_cptac <- function(data, var) {
  formula_unadj <- as.formula(paste0("`", var, "` ~ risk_group"))
  wt <- wilcox.test(formula_unadj, data = data)
  formula_adj <- as.formula(paste0("`", var, "` ~ risk_group + GENOMICS_SUBTYPE"))
  lm_adj <- lm(formula_adj, data = data)
  adj_coef <- summary(lm_adj)$coefficients["risk_groupHigh", ]
  data.frame(
    variable     = var,
    unadj_p      = wt$p.value,
    median_Low   = median(data[[var]][data$risk_group == "Low"],  na.rm = TRUE),
    median_High  = median(data[[var]][data$risk_group == "High"], na.rm = TRUE),
    adj_estimate = adj_coef["Estimate"],
    adj_p        = adj_coef["Pr(>|t|)"]
  )
}

# NOTE: "CIBERSORT_MACROPHAGES _M2" carries a space before "_M2" that looks
tme_vars <- c("EPIC_CAFS", "CIBERSORT_MACROPHAGES _M2", "ESTIMATE_STROMALSCORE", "ESTIMATE_IMMUNESCORE")
missing_tme_vars <- setdiff(tme_vars, colnames(cptac_merged))
if (length(missing_tme_vars) > 0) {
  cat("[CHECK] tme_vars NOT found in cptac_merged:", paste(missing_tme_vars, collapse = ", "), "\n")
  cat("[CHECK] Available columns containing 'MACROPHAGES':",
      paste(grep("MACROPHAGES", colnames(cptac_merged), value = TRUE), collapse = ", "), "\n")
}
stopifnot(all(tme_vars %in% colnames(cptac_merged)))

tme_results <- do.call(rbind, lapply(tme_vars, run_risk_group_test_cptac, data = cptac_merged))
tme_results$unadj_p_fdr <- p.adjust(tme_results$unadj_p, method = "fdr")
tme_results$adj_p_fdr   <- p.adjust(tme_results$adj_p, method = "fdr")

cat("\n=== CPTAC-UCEC external validation: TME scores by risk group ===\n")
print(tme_results, digits = 4)

cptac_cp_expr <- cptac_expr_raw %>%
  filter(Hugo_Symbol %in% checkpoint_genes) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
rownames(cptac_cp_expr) <- cptac_cp_expr$Hugo_Symbol
cptac_cp_expr <- as.data.frame(t(cptac_cp_expr %>% dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)))
cptac_cp_expr$patient_id <- rownames(cptac_cp_expr)

cptac_merged_cp <- merge(cptac_merged, cptac_cp_expr, by = "patient_id")
cat("[CHECK] Checkpoint-merged CPTAC cohort n =", nrow(cptac_merged_cp), "\n")
if (nrow(cptac_merged_cp) < 80) {
  warning("Checkpoint-merged CPTAC n dropped below 80 -- verify the merge, since this is well below the raw cohort size")
}

checkpoint_results_cptac <- do.call(rbind, lapply(checkpoint_genes, run_risk_group_test_cptac, data = cptac_merged_cp))
checkpoint_results_cptac$unadj_p_fdr <- p.adjust(checkpoint_results_cptac$unadj_p, method = "fdr")
checkpoint_results_cptac$adj_p_fdr   <- p.adjust(checkpoint_results_cptac$adj_p, method = "fdr")

cat("\n=== CPTAC-UCEC external validation: checkpoint genes by risk group ===\n")
print(checkpoint_results_cptac, digits = 4)

n_lower_in_high <- sum(checkpoint_results_cptac$median_High < checkpoint_results_cptac$median_Low)
cat("\n[CHECK]", n_lower_in_high, "/ 10 checkpoint genes lower in High-risk (TCGA had 10/10)\n")

# SAVE
saveRDS(cptac_merged_cp, "cptac_validation/cptac_final_validation_data.rds")
cptac_merged_cp <- readRDS("cptac_validation/cptac_final_validation_data.rds")   # checkpoint reload
write.csv(tme_results, "cptac_validation/cptac_tme_validation_results.csv", row.names = FALSE)
write.csv(checkpoint_results_cptac, "cptac_validation/cptac_checkpoint_validation_results.csv", row.names = FALSE)

cat("\n[SAVED] cptac_validation/cptac_final_validation_data.rds\n")
cat("[SAVED] cptac_validation/cptac_tme_validation_results.csv\n")
cat("[SAVED] cptac_validation/cptac_checkpoint_validation_results.csv\n")

# CHECK: DOES CPTAC HAVE A DIFFERENT SUBTYPE MIX THAN TCGA?
# If so, that alone could explain partial replication (e.g. M2 macrophage
# direction flip) without needing "low power" as the only explanation.
library(dplyr)

risk_full <- readRDS("results1/risk_subtype_tmb_hrd.rds")
cptac_merged <- readRDS("cptac_validation/cptac_final_validation_data.rds")

# TCGA subtype distribution
tcga_subtype_pct <- risk_full %>%
  filter(!is.na(Subtype_Selected)) %>%
  count(Subtype_Selected) %>%
  mutate(pct = round(100 * n / sum(n), 1), cohort = "TCGA-UCEC")

# CPTAC subtype distribution (map their labels to the same naming)
subtype_map <- c("CNV_high" = "UCEC.CN_HIGH", "CNV_low" = "UCEC.CN_LOW",
                 "MSI-H" = "UCEC.MSI", "POLE" = "UCEC.POLE")

cptac_subtype_pct <- cptac_merged %>%
  filter(!is.na(GENOMICS_SUBTYPE), GENOMICS_SUBTYPE != "") %>%
  mutate(Subtype_Selected = recode(GENOMICS_SUBTYPE, !!!subtype_map)) %>%
  count(Subtype_Selected) %>%
  mutate(pct = round(100 * n / sum(n), 1), cohort = "CPTAC-UCEC")

cat("[CHECK] TCGA subtype distribution:\n"); print(tcga_subtype_pct)
cat("\n[CHECK] CPTAC subtype distribution:\n"); print(cptac_subtype_pct)

# Formal comparison: does subtype composition differ between cohorts?
comparison_table <- full_join(
  tcga_subtype_pct %>% dplyr::select(Subtype_Selected, TCGA_n = n, TCGA_pct = pct),
  cptac_subtype_pct %>% dplyr::select(Subtype_Selected, CPTAC_n = n, CPTAC_pct = pct),
  by = "Subtype_Selected"
)
cat("\n=== Side-by-side comparison ===\n"); print(comparison_table)

# Chi-square test: is the subtype distribution significantly different
# between the two cohorts?
count_matrix <- comparison_table %>%
  dplyr::select(TCGA_n, CPTAC_n) %>%
  as.matrix()
rownames(count_matrix) <- comparison_table$Subtype_Selected
count_matrix[is.na(count_matrix)] <- 0

chisq_result <- chisq.test(count_matrix)
cat("\n=== Chi-square test: TCGA vs. CPTAC subtype composition ===\n")
print(chisq_result)

fisher.test(count_matrix, simulate.p.value = TRUE, B = 10000)

if (chisq_result$p.value < 0.05) {
  cat("\n[FINDING] Subtype composition DIFFERS significantly between cohorts.\n")
  cat("This is a plausible contributor to partial replication differences\n")
  cat("(e.g. M2 macrophage direction), alongside reduced statistical power.\n")
} else {
  cat("\n[FINDING] No significant difference in subtype composition between cohorts.\n")
  cat("Partial replication differences are more likely attributable to power (n=95) alone.\n")
}

# VERIFICATION: CPTAC external validation
tme_results <- read.csv("cptac_validation/cptac_tme_validation_results.csv")
checkpoint_results_cptac <- read.csv("cptac_validation/cptac_checkpoint_validation_results.csv")

cat("=== TME results ===\n")
print(tme_results, digits = 4)

cat("\n=== Checkpoint gene results ===\n")
print(checkpoint_results_cptac, digits = 4)

# Check 1: direction consistency (all 10 lower in High-risk?)
n_lower <- sum(checkpoint_results_cptac$median_High < checkpoint_results_cptac$median_Low)
cat("\n[VERIFY] Checkpoint genes lower in High-risk:", n_lower, "/ 10\n")

# Check 2: how many survive FDR + subtype adjustment, and what's the p range?
sig_checkpoints <- checkpoint_results_cptac[checkpoint_results_cptac$adj_p_fdr < 0.05, ]
cat("[VERIFY] Checkpoint genes significant after FDR+subtype adjustment:", nrow(sig_checkpoints), "/ 10\n")
cat("[VERIFY] Adjusted+FDR p range among significant ones:",
    range(sig_checkpoints$adj_p_fdr), "\n")

# Check 3: ESTIMATE_IMMUNESCORE
cat("\n[VERIFY] ESTIMATE_IMMUNESCORE row:\n")
print(tme_results[tme_results$variable == "ESTIMATE_IMMUNESCORE", ], digits = 4)

# Check 4: EPIC_CAFS (CAF enrichment)
cat("\n[VERIFY] EPIC_CAFS row:\n")
print(tme_results[tme_results$variable == "EPIC_CAFS", ], digits = 4)

# Check 5: CIBERSORT M2 macrophage direction
cat("\n[VERIFY] CIBERSORT M2 macrophage row:\n")
print(tme_results[grepl("MACROPHAGES", tme_results$variable), ], digits = 4)

# Check 6: cohort n
cptac_final <- readRDS("cptac_validation/cptac_final_validation_data.rds")
cat("\n[VERIFY] CPTAC cohort n =", nrow(cptac_final), "\n")

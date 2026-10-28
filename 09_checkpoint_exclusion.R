# 09_checkpoint_exclusion_WL.R
#
# Immune checkpoint genes and TIDE dysfunction/exclusion components by risk group.
# Input: results1/survival_with_tide_ips.rds, ips_calc/IPS_genes.txt  |  Output: results1/survival_with_checkpoints.rds

library(dplyr)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# LOAD DEPENDENCIES EXPLICITLY
survival_full <- readRDS("results1/survival_with_tide_ips.rds")
expr_vst <- readRDS("results/expr_tumor_vst_genesymbol.rds")
colnames(expr_vst) <- substr(colnames(expr_vst), 1, 12)

stopifnot(nrow(survival_full) == 541)
stopifnot(all(c("Dysfunction", "Exclusion", "CAF", "TAM.M2", "MDSC", "CTL.flag") %in% colnames(survival_full)))

# IMMUNE CHECKPOINT GENE PANEL
ips_genes <- read.table("ips_calc/IPS_genes.txt", header = TRUE, sep = "\t", check.names = FALSE)
checkpoint_panel <- ips_genes %>% dplyr::filter(CLASS == "CP")
checkpoint_genes <- checkpoint_panel$GENE

stopifnot(length(checkpoint_genes) == 10)
cat("[CHECK] Checkpoint panel (n=10):", paste(checkpoint_genes, collapse = ", "), "\n")

missing_cp <- setdiff(checkpoint_genes, rownames(expr_vst))
stopifnot(length(missing_cp) == 0)

common_patients <- intersect(colnames(expr_vst), survival_full$patient_id)
cp_expr <- as.data.frame(t(expr_vst[checkpoint_genes, common_patients]))
cp_expr$patient_id <- rownames(cp_expr)

dup_names <- intersect(checkpoint_genes, colnames(survival_full))
if (length(dup_names) > 0) {
  cat("[CHECK] Dropping pre-existing columns to avoid merge collision:", paste(dup_names, collapse = ", "), "\n")
  survival_full <- survival_full %>% dplyr::select(-all_of(dup_names))
}

survival_cp <- merge(survival_full, cp_expr, by = "patient_id")
stopifnot(nrow(survival_cp) == 541)
stopifnot(!any(grepl("\\.x$|\\.y$", colnames(survival_cp))))   # confirm no silent collision

# CHECKPOINT (downstream of risk_group -- results1/)
saveRDS(survival_cp, "results1/checkpoint_13a_survival_cp.rds")
cat("[SAVED] results1/checkpoint_13a_survival_cp.rds\n\n")

# RISK_GROUP VS. EACH CHECKPOINT GENE (UNADJUSTED + SUBTYPE-ADJUSTED)
survival_cp <- readRDS("results1/checkpoint_13a_survival_cp.rds")

run_risk_group_test <- function(data, var, adjust_for_subtype = TRUE) {
  formula_unadj <- as.formula(paste(var, "~ risk_group"))
  wt <- wilcox.test(formula_unadj, data = data)
  out <- data.frame(
    variable    = var,
    unadj_p     = wt$p.value,
    median_Low  = median(data[[var]][data$risk_group == "Low"],  na.rm = TRUE),
    median_High = median(data[[var]][data$risk_group == "High"], na.rm = TRUE)
  )
  if (adjust_for_subtype) {
    formula_adj <- as.formula(paste(var, "~ risk_group + Subtype_Selected"))
    lm_adj <- lm(formula_adj, data = data)
    adj_coef <- summary(lm_adj)$coefficients["risk_groupHigh", ]
    out$adj_estimate <- adj_coef["Estimate"]
    out$adj_p        <- adj_coef["Pr(>|t|)"]
  }
  out
}

checkpoint_results <- do.call(rbind, lapply(checkpoint_genes, run_risk_group_test, data = survival_cp))
checkpoint_results$unadj_p_fdr <- p.adjust(checkpoint_results$unadj_p, method = "fdr")
checkpoint_results$adj_p_fdr   <- p.adjust(checkpoint_results$adj_p, method = "fdr")

cat("\n=== Checkpoint gene results (FDR-corrected across 10 genes) ===\n")
print(checkpoint_results, digits = 4)

# SANITY CHECK: expect direction consistent with prior run
n_lower_in_high <- sum(checkpoint_results$median_High < checkpoint_results$median_Low)
cat("[CHECK]", n_lower_in_high, "/ 10 checkpoint genes lower in High-risk (expect 10)\n")

# DYSFUNCTION VS. EXCLUSION MECHANISM BREAKDOWN
tide_components <- c("Dysfunction", "Exclusion", "CAF", "TAM.M2", "MDSC")
component_results <- do.call(rbind, lapply(tide_components, run_risk_group_test, data = survival_cp))
component_results$unadj_p_fdr <- p.adjust(component_results$unadj_p, method = "fdr")
component_results$adj_p_fdr   <- p.adjust(component_results$adj_p, method = "fdr")

cat("\n=== TIDE sub-component results ===\n")
print(component_results, digits = 4)

# Primary mechanism split (CTL.flag TRUE = Dysfunction primary; FALSE = Exclusion primary)
ctl_flag_table <- table(survival_cp$risk_group, survival_cp$CTL.flag)
ctl_flag_chisq <- chisq.test(ctl_flag_table)
cat("\n=== Primary mechanism (CTL.flag) by risk group ===\n")
print(ctl_flag_table)
print(ctl_flag_chisq)

# SAVE FINAL OBJECT + RESULTS TABLES
saveRDS(survival_cp, "results1/survival_with_checkpoints.rds")
survival_cp <- readRDS("results1/survival_with_checkpoints.rds")   # checkpoint reload
write.csv(checkpoint_results, "results1/checkpoint_gene_results.csv", row.names = FALSE)
write.csv(component_results, "results1/tide_component_results.csv", row.names = FALSE)

cat("\n[SAVED] results1/survival_with_checkpoints.rds\n")
cat("[SAVED] results1/checkpoint_gene_results.csv\n")
cat("[SAVED] results1/tide_component_results.csv\n")

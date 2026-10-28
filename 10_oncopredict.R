# 10_oncopredict_WL.R
#
# Predicted drug sensitivity (oncoPredict, GDSC2) for endometrial-cancer-relevant drugs by risk group.
# Requires DataFiles.zip from https://osf.io/c6tfx/ unzipped into oncopredict_data/.
# Input: results1/survival_with_checkpoints.rds  |  Output: results1/survival_with_drug_sensitivity.rds, results1/drug_sensitivity_results.csv

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
  BiocManager::install("GenomicFeatures", update = FALSE, ask = FALSE)
}
if (!requireNamespace("oncoPredict", quietly = TRUE)) install.packages("oncoPredict")

library(oncoPredict)
library(dplyr)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# DOWNLOAD + LOAD GDSC2 TRAINING DATA
train_dir <- "oncopredict_data/DataFiles/Training Data/"

if (!dir.exists(train_dir)) {
  if (!file.exists("DataFiles.zip")) {
    tryCatch(
      download.file("https://osf.io/c6tfx/files/osfstorage",
                    destfile = "DataFiles.zip", mode = "wb"),
      error = function(e) NULL
    )
  }
  
  valid_zip <- file.exists("DataFiles.zip") &&
    tryCatch({ unzip("DataFiles.zip", list = TRUE); TRUE }, error = function(e) FALSE)
  
  if (!valid_zip) {
    stop(
      "GDSC2 training data not found or not a valid zip.\n",
      "Manually download 'DataFiles.zip' from https://osf.io/c6tfx/ in your\n",
      "browser, place it in the working directory, and re-run this script."
    )
  }
  
  unzip("DataFiles.zip", exdir = "oncopredict_data")
}

stopifnot(dir.exists(train_dir))

GDSC2_Expr <- readRDS(file.path(train_dir, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
GDSC2_Res  <- readRDS(file.path(train_dir, "GDSC2_Res.rds"))
GDSC2_Res  <- exp(GDSC2_Res)   # stored log-transformed; calcPhenotype expects raw IC50 scale

# SANITY CHECK: training data shape
stopifnot(dim(GDSC2_Expr)[1] > 15000, dim(GDSC2_Res)[2] > 150)
cat("[CHECK] GDSC2_Expr:", paste(dim(GDSC2_Expr), collapse = " x "), "\n")
cat("[CHECK] GDSC2_Res:",  paste(dim(GDSC2_Res),  collapse = " x "), "(cell lines x drugs)\n")

# PREPARE TEST EXPRESSION DATA
expr_vst <- readRDS("results/expr_tumor_vst_genesymbol.rds")
colnames(expr_vst) <- substr(colnames(expr_vst), 1, 12)

survival_cp <- readRDS("results1/survival_with_checkpoints.rds")
common_patients <- intersect(colnames(expr_vst), survival_cp$patient_id)
stopifnot(length(common_patients) == 541)

testExpr <- as.matrix(expr_vst[, common_patients])
cat("[CHECK] testExpr:", paste(dim(testExpr), collapse = " x "), "\n")

if (!file.exists("results/checkpoint_09a_drug_predictions.rds")) {
  dir.create("results/oncopredict_out", showWarnings = FALSE)
  
  calcPhenotype(
    trainingExprData         = GDSC2_Expr,
    trainingPtype             = GDSC2_Res,
    testExprData              = testExpr,
    batchCorrect               = "eb",
    powerTransformPhenotype   = TRUE,
    removeLowVaryingGenes     = 0.2,
    minNumSamples              = 10,
    printOutput                = TRUE,
    removeLowVaringGenesFrom  = "rawData"
  )
  
  # calcPhenotype writes its own output file rather than returning an object
  pred_path <- "calcPhenotype_Output/DrugPredictions.csv"
  stopifnot(file.exists(pred_path))
  
  drug_predictions <- read.csv("calcPhenotype_Output/DrugPredictions.csv", row.names = 1, check.names = FALSE)
  
  # SANITY CHECK: prediction output shape
  stopifnot(nrow(drug_predictions) == 541)
  cat("[CHECK] Drug predictions:", nrow(drug_predictions), "patients x", ncol(drug_predictions), "drugs\n")
  
  drug_predictions$patient_id <- rownames(drug_predictions)
  
  # CHECKPOINT (expression-only -- results/, always save, never rerun)
  saveRDS(drug_predictions, "results/checkpoint_09a_drug_predictions.rds")
  cat("[SAVED] results/checkpoint_09a_drug_predictions.rds\n\n")
} else {
  cat("[CHECK] Drug predictions loaded from cache (calcPhenotype skipped)\n\n")
}

# TEST CLINICALLY RELEVANT, EC-EVIDENCED DRUGS BY RISK GROUP
drug_predictions <- readRDS("results/checkpoint_09a_drug_predictions.rds")
drug_predictions <- as.data.frame(drug_predictions)
drug_predictions$patient_id <- rownames(drug_predictions)

candidate_terms <- c("Paclitaxel", "Cisplatin", "Docetaxel", "Oxaliplatin",
                     "Epirubicin", "Olaparib", "Talazoparib",
                     "Niraparib", "Tamoxifen", "Fulvestrant")
matched_cols <- unique(unlist(lapply(candidate_terms, function(t) grep(t, colnames(drug_predictions), value = TRUE))))

cat("[CHECK] Matched EC-evidenced drug columns (n =", length(matched_cols), "):\n")
print(matched_cols)

stopifnot(length(grep("Carboplatin", colnames(drug_predictions), ignore.case = TRUE)) == 0)
stopifnot(length(grep("Doxorubicin", colnames(drug_predictions), ignore.case = TRUE)) == 0)

survival_drug <- merge(survival_cp, drug_predictions, by = "patient_id")
stopifnot(nrow(survival_drug) == 541)

run_risk_group_test <- function(data, var) {
  formula_unadj <- as.formula(paste0("`", var, "` ~ risk_group"))
  wt <- wilcox.test(formula_unadj, data = data)
  formula_adj <- as.formula(paste0("`", var, "` ~ risk_group + Subtype_Selected"))
  lm_adj <- lm(formula_adj, data = data)
  adj_coef <- summary(lm_adj)$coefficients["risk_groupHigh", ]
  data.frame(
    drug         = var,
    unadj_p      = wt$p.value,
    median_Low   = median(data[[var]][data$risk_group == "Low"],  na.rm = TRUE),
    median_High  = median(data[[var]][data$risk_group == "High"], na.rm = TRUE),
    adj_estimate = adj_coef["Estimate"],
    adj_p        = adj_coef["Pr(>|t|)"]
  )
}

drug_results <- do.call(rbind, lapply(matched_cols, run_risk_group_test, data = survival_drug))
drug_results$unadj_p_fdr <- p.adjust(drug_results$unadj_p, method = "fdr")
drug_results$adj_p_fdr   <- p.adjust(drug_results$adj_p, method = "fdr")

cat("\n=== Drug sensitivity results (lower predicted IC50 = more sensitive) ===\n")
print(drug_results, digits = 4)

saveRDS(survival_drug, "results1/survival_with_drug_sensitivity.rds")
survival_drug <- readRDS("results1/survival_with_drug_sensitivity.rds")   # checkpoint reload
write.csv(drug_results, "results1/drug_sensitivity_results.csv", row.names = FALSE)

cat("\n[SAVED] results1/survival_with_drug_sensitivity.rds\n[SAVED] results1/drug_sensitivity_results.csv\n")

# VERIFICATION: drug sensitivity results
drug_results <- read.csv("results1/drug_sensitivity_results.csv")
print(drug_results, digits = 4)

# Check 1: how many show higher IC50 (lower sensitivity) in High-risk
cat("\n[VERIFY] Show median_High > median_Low (higher IC50 in High-risk):",
    sum(drug_results$median_High > drug_results$median_Low), "/", nrow(drug_results), "\n")

# Check 2: adjusted p range, after FDR correction
cat("[VERIFY] Adjusted+FDR p range:", range(drug_results$adj_p_fdr), "\n")

# Check 3: consistency for the 3 duplicated-batch drugs
dup_drugs <- drug_results[grepl("Oxaliplatin|Docetaxel|Fulvestrant", drug_results$drug), ]
cat("\n[VERIFY] Duplicate-batch drugs (direction + significance):\n")
print(dup_drugs[, c("drug", "median_Low", "median_High", "adj_p", "adj_p_fdr")], digits = 4)

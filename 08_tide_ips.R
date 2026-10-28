# 08_tide_ips_WL.R
#
# TIDE (TIDEpy) and Immunophenoscore (IPS) by risk group, unadjusted and subtype-adjusted.
# Input: results/expr_tumor_vst_genesymbol.rds, results1/full_survival_with_risk.rds  |  Output: results1/survival_with_tide_ips.rds

library(dplyr)
library(reticulate)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# PREPARE EXPRESSION MATRIX (TIDE NEEDS MEDIAN-CENTERED, IPS NEEDS RAW)
expr_vst <- readRDS("results/expr_tumor_vst_genesymbol.rds")
colnames(expr_vst) <- substr(colnames(expr_vst), 1, 12)

survival_patient <- readRDS("results1/full_survival_with_risk.rds")

common_patients <- intersect(colnames(expr_vst), survival_patient$patient_id)
stopifnot(length(common_patients) == 541)
cat("[CHECK] Common patients for expression-based scoring:", length(common_patients), "\n")

expr_cohort <- expr_vst[, common_patients]

# TIDE -- MEDIAN-CENTER (PER TIDE'S DOCUMENTED REQUIREMENT), EXPORT
gene_medians <- apply(expr_cohort, 1, median, na.rm = TRUE)
expr_tide_centered <- sweep(expr_cohort, 1, gene_medians, FUN = "-")

# SANITY CHECK: centering worked
med_check <- summary(apply(expr_tide_centered, 1, median))
stopifnot(abs(med_check["Median"]) < 1e-8)
cat("[CHECK] Gene-wise median after centering (expect ~0):", med_check["Median"], "\n")

tide_export <- data.frame(GeneSymbol = rownames(expr_tide_centered), expr_tide_centered, check.names = FALSE)
write.table(tide_export, file = "results/tide_input_ucec_n541.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# RUN TIDE VIA TIDEPY (PYTHON, THROUGH RETICULATE)

# 3a. Ensure Miniconda is installed
# NOTE: on Windows, install_miniconda() uses pyenv-win, which can crash with
# BrokenProcessPool errors under conda's newer experimental "rattler" solver.
# Forcing CONDA_SOLVER="classic" before install avoids this.
Sys.setenv(CONDA_SOLVER = "classic")

if (!tryCatch({ reticulate::conda_binary(); TRUE }, error = function(e) FALSE)) {
  reticulate::install_miniconda()
}

# 3b. Ensure the tide_env conda environment exists
existing_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
if (!"tide_env" %in% existing_envs) {
  reticulate::conda_create("tide_env", python_version = "3.11")
}

# 3c. Ensure tidepy is installed in tide_env
reticulate::use_condaenv("tide_env", required = TRUE)
tidepy_available <- tryCatch({
  reticulate::py_run_string("import tidepy")
  TRUE
}, error = function(e) FALSE)

if (!tidepy_available) {
  reticulate::conda_install("tide_env", packages = "tidepy", pip = TRUE)
  reticulate::use_condaenv("tide_env", required = TRUE)
}

# SANITY CHECK: environment ready
py_config_check <- reticulate::py_config()
cat("[CHECK] Python in use:", py_config_check$python, "\n")
reticulate::py_run_string("import tidepy; print('tidepy OK:', tidepy.__file__)")

# TIDE scoring is expression-only (independent of the signature/risk_group),
# and the TIDEpy call is the slowest step in this script -- skip it if the
# output already exists so re-running this script after a signature update
# never needlessly recomputes it.
if (!file.exists("results/tide_scores_ucec_n541.csv")) {
  py_run_string("
import pandas as pd
from tidepy.pred import TIDE

df = pd.read_csv('results/tide_input_ucec_n541.txt', sep='\t', index_col=0)
assert df.shape[1] == 541, f'Expected 541 patients, got {df.shape[1]}'

result = TIDE(df, cancer='Other', pretreat=False, ignore_norm=True, vthres=0.)
result.to_csv('results/tide_scores_ucec_n541.csv')
")
}

tide_scores <- read.csv("results/tide_scores_ucec_n541.csv", row.names = 1)
tide_scores$patient_id <- rownames(tide_scores)

# SANITY CHECK: TIDE output shape
stopifnot(nrow(tide_scores) == 541)
cat("[CHECK] TIDE scores:", nrow(tide_scores), "patients x", ncol(tide_scores) - 1, "columns\n")

survival_tide <- merge(survival_patient, tide_scores, by = "patient_id")
stopifnot(sum(!is.na(survival_tide$TIDE)) == 541)

# CHECKPOINT (downstream of risk_group -- results1/)
saveRDS(survival_tide, "results1/checkpoint_12a_survival_tide.rds")
cat("[SAVED] results1/checkpoint_12a_survival_tide.rds\n\n")

# IPS -- UNCENTERED VST INPUT (DIFFERENT NORMALIZATION FROM TIDE)
survival_tide <- readRDS("results1/checkpoint_12a_survival_tide.rds")

dir.create("ips_calc", showWarnings = FALSE)
ips_expr_export <- data.frame(GeneSymbol = rownames(expr_cohort), expr_cohort, check.names = FALSE)
write.table(ips_expr_export, file = "ips_calc/EXPR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

if (!file.exists("ips_calc/IPS.R")) {
  download.file("https://raw.githubusercontent.com/icbi-lab/Immunophenogram/master/IPS.R",
                destfile = "ips_calc/IPS.R", mode = "wb")
}
if (!file.exists("ips_calc/IPS_genes.txt")) {
  download.file("https://raw.githubusercontent.com/icbi-lab/Immunophenogram/master/IPS_genes.txt",
                destfile = "ips_calc/IPS_genes.txt", mode = "wb")
}

ipsmap <- function(x) {
  if (x <= 0) return(0)
  if (x >= 3) return(10)
  return(round(x * 10 / 3, digits = 0))
}

# IPS scoring is also expression-only (independent of the signature) -- skip
# the loop if the output already exists.
if (!file.exists("results/ips_scores_ucec_n541.txt")) {
  gene_expression <- read.table("ips_calc/EXPR.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
  IPSG <- read.table("ips_calc/IPS_genes.txt", header = TRUE, sep = "\t", check.names = FALSE)
  unique_ips_genes <- as.vector(unique(IPSG$NAME))
  sample_names <- names(gene_expression)
  
  missing_genes <- setdiff(as.vector(IPSG$GENE), rownames(gene_expression))
  cat("[CHECK] IPS genes missing from expression matrix:", paste(missing_genes, collapse = ", "), "\n")
  
  IPS <- MHC <- CP <- EC <- SC <- AZ <- numeric(length(sample_names))
  for (i in seq_along(sample_names)) {
    GE  <- gene_expression[[i]]
    Z1  <- (gene_expression[as.vector(IPSG$GENE), i] - mean(GE)) / sd(GE)
    W1  <- IPSG$WEIGHT
    MIG <- WEIGHT <- numeric(length(unique_ips_genes))
    for (k in seq_along(unique_ips_genes)) {
      gen <- unique_ips_genes[k]
      MIG[k]    <- mean(Z1[which(as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(W1[which(as.vector(IPSG$NAME) == gen)])
    }
    WG <- MIG * WEIGHT
    MHC[i] <- mean(WG[1:10]); EC[i] <- mean(WG[21:24]); SC[i] <- mean(WG[25:26]); CP[i] <- mean(WG[11:20])
    AZ[i]  <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }
  
  ips_df <- data.frame(patient_id = sample_names, MHC = MHC, EC = EC, SC = SC, CP = CP, AZ = AZ, IPS = IPS)
  
  # SANITY CHECK: IPS output
  stopifnot(nrow(ips_df) == 541, all(ips_df$IPS >= 0 & ips_df$IPS <= 10))
  cat("[CHECK] IPS scores: n=", nrow(ips_df), " range=[", min(ips_df$IPS), ",", max(ips_df$IPS), "]\n")
  
  write.table(ips_df, file = "results/ips_scores_ucec_n541.txt", row.names = FALSE, quote = FALSE, sep = "\t")
} else {
  ips_df <- read.table("results/ips_scores_ucec_n541.txt", header = TRUE, sep = "\t", check.names = FALSE)
  cat("[CHECK] IPS scores loaded from cache: n=", nrow(ips_df), " range=[", min(ips_df$IPS), ",", max(ips_df$IPS), "]\n")
}

# MERGE, TEST, SAVE
survival_full <- merge(survival_tide, ips_df, by = "patient_id")
stopifnot(nrow(survival_full) == 541)

subtype_lookup <- readRDS("results1/checkpoint_11a_risk_subtype.rds") %>%
  dplyr::select(patient_id, Subtype_Selected)
survival_full <- merge(survival_full, subtype_lookup, by = "patient_id", all.x = TRUE)

cat("\n=== TIDE: unadjusted ===\n");        print(wilcox.test(TIDE ~ risk_group, data = survival_full))
cat("\n=== TIDE: subtype-adjusted ===\n");  print(summary(lm(TIDE ~ risk_group + Subtype_Selected, data = survival_full))$coefficients["risk_groupHigh", ])
cat("\n=== IPS: unadjusted ===\n");         print(wilcox.test(IPS ~ risk_group, data = survival_full))
cat("\n=== IPS: subtype-adjusted ===\n");   print(summary(lm(IPS ~ risk_group + Subtype_Selected, data = survival_full))$coefficients["risk_groupHigh", ])

saveRDS(survival_full, "results1/survival_with_tide_ips.rds")
survival_full <- readRDS("results1/survival_with_tide_ips.rds")   # checkpoint reload
cat("\n[SAVED] results1/survival_with_tide_ips.rds  (", nrow(survival_full), "patients )\n")

# VARIANCE IN TIDE EXPLAINED BY SUBTYPE ALONE
tide_subtype_only <- lm(TIDE ~ Subtype_Selected, data = survival_full)
cat("\n[VERIFY] R-squared, subtype alone, TIDE model:",
    summary(tide_subtype_only)$r.squared, "\n")

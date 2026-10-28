# 03_cox_lasso_WL.R
#
# Prognostic model: univariate Cox screening, LASSO-Cox (10-fold CV), risk score, Kaplan-Meier, time-dependent ROC and subtype-adjusted multivariable Cox.
# Input: results1/ (from 02)  |  Output: results1/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks", update = FALSE, ask = FALSE)
}

library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(timeROC)
library(broom)
library(TCGAbiolinks)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# LOAD PRE-SPLIT COHORT + CANDIDATE STROMAL DEG EXPRESSION
stromal_deg_expr <- readRDS("results1/stromal_deg_expr_matrix.rds")            # 541 x 501
split            <- readRDS("results1/analytic_cohort_train_test_split.rds")   # train/test patient IDs
survival_patient <- readRDS("results/survival_patient.rds")   # raw input
clinical_patient <- readRDS("results/clinical_patient.rds")   # raw input
survival_patient$patient_id <- substr(survival_patient$sample, 1, 12)

survival_patient <- survival_patient[match(rownames(stromal_deg_expr), survival_patient$patient_id), ]
clinical_patient  <- clinical_patient[match(rownames(stromal_deg_expr), clinical_patient$patient_id), ]

stopifnot(all(rownames(stromal_deg_expr) == survival_patient$patient_id))
stopifnot(all(rownames(stromal_deg_expr) == clinical_patient$patient_id))
stopifnot(all(survival_patient$patient_id == clinical_patient$patient_id))

dim(stromal_deg_expr)   # expect 541 x 501

train_expr <- stromal_deg_expr[split$train_patients, ]
test_expr  <- stromal_deg_expr[split$test_patients, ]
dim(train_expr)   # expect 405 x 501
dim(test_expr)    # expect 136 x 501

train_survival <- survival_patient[match(split$train_patients, survival_patient$patient_id), ]
test_survival  <- survival_patient[match(split$test_patients,  survival_patient$patient_id), ]
dim(train_survival)   # expect 405 x 7
dim(test_survival)    # expect 136 x 7

# UNIVARIATE COX SCREENING (TRAINING COHORT ONLY)
surv_object <- Surv(time = train_survival$OS.time, event = train_survival$OS)

uni_cox_results <- lapply(colnames(train_expr), function(gene) {
  cox_model   <- coxph(surv_object ~ train_expr[, gene])
  cox_summary <- summary(cox_model)
  data.frame(
    Gene    = gene,
    HR      = cox_summary$conf.int[, "exp(coef)"],
    lower95 = cox_summary$conf.int[, "lower .95"],
    upper95 = cox_summary$conf.int[, "upper .95"],
    p.value = cox_summary$coefficients[, "Pr(>|z|)"]
  )
})
uni_cox_results <- do.call(rbind, uni_cox_results)
dim(uni_cox_results)   # expect 501 x 5

uni_cox_results$FDR <- p.adjust(uni_cox_results$p.value, method = "fdr")

uni_sig_genes <- uni_cox_results %>% filter(p.value < 0.01)
dim(uni_sig_genes)   # expect 32 x 6

saveRDS(uni_cox_results, "results1/univariate_cox_results_train.rds")
saveRDS(uni_sig_genes,   "results1/univariate_significant_genes_train.rds")

# LASSO-COX (10-FOLD CV)
uni_cox_results <- readRDS("results1/univariate_cox_results_train.rds")
uni_sig_genes <- readRDS("results1/univariate_significant_genes_train.rds")
lasso_genes <- uni_sig_genes$Gene
x <- as.matrix(train_expr[, lasso_genes])
y <- Surv(train_survival$OS.time, train_survival$OS)
dim(x)   # 405 x 32

set.seed(123)
foldid <- sample(rep(1:10, length.out = nrow(x)))
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, foldid = foldid, cox.ties = "breslow")
saveRDS(cvfit, "results1/cvfit.rds")

cvfit$lambda.min
cvfit$lambda.1se

lasso_coef <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
lasso_genes

saveRDS(lasso_genes, "results1/lasso_selected_genes_train.rds")
saveRDS(lasso_coef,  "results1/lasso_coefficients_train.rds")

# FINAL COEFFICIENT TABLE
lasso_coef_df <- as.data.frame(as.matrix(lasso_coef))
lasso_coef_df$Gene <- rownames(lasso_coef_df)
colnames(lasso_coef_df)[1] <- "Coefficient"
lasso_coef_df <- lasso_coef_df %>% filter(Coefficient != 0)

lasso_summary <- lasso_coef_df %>%
  left_join(
    uni_cox_results %>% dplyr::select(Gene, HR, lower95, upper95, p.value),
    by = "Gene"
  ) %>%
  mutate(
    Role   = ifelse(Coefficient > 0, "Risk", "Protective"),
    `95% CI` = sprintf("%.3f-%.3f", lower95, upper95)
  ) %>%
  dplyr::select(Gene, HR, `95% CI`, p.value, Coefficient, Role)

lasso_summary
saveRDS(lasso_coef_df, "results1/lasso_coefficients_dataframe.rds")

# RISK SCORE CALCULATION (TRAIN / TEST / FULL)
selected_genes <- lasso_coef_df$Gene
lasso_coef_df  <- lasso_coef_df[match(selected_genes, lasso_coef_df$Gene), ]
coef_vector    <- lasso_coef_df$Coefficient
names(coef_vector) <- lasso_coef_df$Gene

train_selected <- train_expr[, selected_genes]
test_selected  <- test_expr[, selected_genes]

train_survival$risk_score <- as.vector(as.matrix(train_selected) %*% coef_vector)
test_survival$risk_score  <- as.vector(as.matrix(test_selected)  %*% coef_vector)

cutoff <- median(train_survival$risk_score)
saveRDS(cutoff, "results1/risk_score_cutoff.rds")

train_survival$risk_group <- factor(
  ifelse(train_survival$risk_score >= cutoff, "High", "Low"),
  levels = c("Low", "High")
)
test_survival$risk_group <- factor(
  ifelse(test_survival$risk_score >= cutoff, "High", "Low"),
  levels = c("Low", "High")
)

saveRDS(train_survival, "results1/train_survival_with_risk.rds")
saveRDS(test_survival,  "results1/test_survival_with_risk.rds")

# KAPLAN-MEIER (TRAIN / TEST)
km_train <- survfit(Surv(OS.time, OS) ~ risk_group, data = train_survival)
ggsurvplot(km_train, data = train_survival, pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
           legend.title = "Risk Group", legend.labs = c("Low", "High"),
           xlab = "Time (days)", ylab = "Overall Survival Probability",
           title = "Kaplan-Meier Survival Curve (Training Cohort)",
           palette = c("#00BFC4", "#F8766D"))
saveRDS(km_train, "results1/km_train.rds")

km_test <- survfit(Surv(OS.time, OS) ~ risk_group, data = test_survival)
ggsurvplot(km_test, data = test_survival, pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
           legend.title = "Risk Group", legend.labs = c("Low", "High"),
           xlab = "Time (days)", ylab = "Overall Survival Probability",
           title = "Kaplan-Meier Survival Curve (Testing Cohort)",
           palette = c("#00BFC4", "#F8766D"))
saveRDS(km_test, "results1/km_test.rds")

# UNIVARIATE COX ON RISK SCORE (TRAIN / TEST)
cox_train <- coxph(Surv(OS.time, OS) ~ risk_score, data = train_survival)
summary(cox_train)
saveRDS(cox_train, "results1/cox_train.rds")

cox_test <- coxph(Surv(OS.time, OS) ~ risk_score, data = test_survival)
summary(cox_test)
saveRDS(cox_test, "results1/cox_test.rds")

cox_train_results <- tidy(cox_train, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value)
cox_test_results <- tidy(cox_test, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value)

saveRDS(cox_train_results, "results1/cox_train_results.rds")
saveRDS(cox_test_results,  "results1/cox_test_results.rds")

# TIME-DEPENDENT ROC (TRAIN / TEST)
roc_train <- timeROC(
  T = train_survival$OS.time, delta = train_survival$OS,
  marker = train_survival$risk_score, cause = 1,
  weighting = "marginal", times = c(365, 1095, 1825), iid = TRUE
)
roc_train$AUC
saveRDS(roc_train, "results1/roc_train.rds")

roc_test <- timeROC(
  T = test_survival$OS.time, delta = test_survival$OS,
  marker = test_survival$risk_score, cause = 1,
  weighting = "marginal", times = c(365, 1095, 1825), iid = TRUE
)
roc_test$AUC
saveRDS(roc_test, "results1/roc_test.rds")

# FULL-COHORT RISK SCORE, KM, COX, ROC
full_selected <- stromal_deg_expr[, selected_genes]
dim(full_selected)   # 541 x 9

survival_patient$risk_score <- as.vector(as.matrix(full_selected) %*% coef_vector)
summary(survival_patient$risk_score)

survival_patient$risk_group <- factor(
  ifelse(survival_patient$risk_score >= cutoff, "High", "Low"),
  levels = c("Low", "High")
)
table(survival_patient$risk_group)

saveRDS(survival_patient, "results1/full_survival_with_risk.rds")

km_full <- survfit(Surv(OS.time, OS) ~ risk_group, data = survival_patient)
ggsurvplot(km_full, data = survival_patient, pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
           legend.title = "Risk Group", legend.labs = c("Low", "High"),
           xlab = "Time (days)", ylab = "Overall Survival Probability",
           title = "Kaplan-Meier Survival Curve (Entire TCGA-UCEC Cohort)",
           palette = c("#00BFC4", "#F8766D"))
saveRDS(km_full, "results1/km_full.rds")

cox_full <- coxph(Surv(OS.time, OS) ~ risk_score, data = survival_patient)
summary(cox_full)
saveRDS(cox_full, "results1/cox_full.rds")

roc_full <- timeROC(
  T = survival_patient$OS.time, delta = survival_patient$OS,
  marker = survival_patient$risk_score, cause = 1,
  weighting = "marginal", times = c(365, 1095, 1825), iid = TRUE
)
roc_full$AUC
saveRDS(roc_full, "results1/roc_full.rds")

# MULTIVARIABLE COX
clinical_selected <- clinical_patient[, c("patient_id", "age_at_index.demographic", "figo_stage.diagnoses")]

multi_data_base <- merge(survival_patient, clinical_selected, by = "patient_id")
multi_data_base$age   <- as.numeric(multi_data_base$age_at_index.demographic)
multi_data_base$stage <- multi_data_base$figo_stage.diagnoses

multi_data_base <- multi_data_base %>%
  filter(!is.na(risk_score), !is.na(age), !is.na(stage))

multi_data_base$stage_group <- factor(
  ifelse(
    multi_data_base$stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                 "Stage IIIC1", "Stage IIIC2", "Stage IV", "Stage IVA", "Stage IVB"),
    "Advanced", "Early"
  ),
  levels = c("Early", "Advanced")
)

dim(multi_data_base)
saveRDS(multi_data_base, "results1/multi_data_base_n539.rds")

uni_age   <- coxph(Surv(OS.time, OS) ~ age, data = multi_data_base)
uni_stage <- coxph(Surv(OS.time, OS) ~ stage_group, data = multi_data_base)
uni_risk  <- coxph(Surv(OS.time, OS) ~ risk_score, data = multi_data_base)
print(summary(uni_risk))
saveRDS(uni_risk, "results1/uni_risk_n539.rds")

all_subtypes <- PanCancerAtlas_subtypes()
ucec_subtypes <- all_subtypes %>%
  filter(cancer.type == "UCEC") %>%
  mutate(tcga_patient = substr(pan.samplesID, 1, 12)) %>%
  dplyr::select(tcga_patient, Subtype_Selected)
stopifnot(sum(duplicated(ucec_subtypes$tcga_patient)) == 0)

multi_data_subtype <- merge(
  multi_data_base, ucec_subtypes,
  by.x = "patient_id", by.y = "tcga_patient", all.x = TRUE
)
stopifnot(nrow(multi_data_subtype) == nrow(multi_data_base))

n_subtype <- sum(!is.na(multi_data_subtype$Subtype_Selected))
cat("[CHECK] Patients with subtype call:", n_subtype, "/", nrow(multi_data_subtype), "\n")
if (n_subtype < 500) warning("Fewer subtype matches than expected (~530) -- verify barcode format")

multi_data <- multi_data_subtype %>% filter(!is.na(Subtype_Selected))
dim(multi_data)

saveRDS(multi_data, "results1/multi_data.rds")
multi_data <- readRDS("results1/multi_data.rds")

# Nested 3-model comparison, all on n=530
model1 <- coxph(Surv(OS.time, OS) ~ risk_score, data = multi_data)
model2 <- coxph(Surv(OS.time, OS) ~ risk_score + age + stage_group, data = multi_data)
model3 <- coxph(Surv(OS.time, OS) ~ risk_score + age + stage_group + Subtype_Selected, data = multi_data)

cat("\n=== Model 1: risk_score alone (n=", nrow(multi_data), ") ===\n"); print(summary(model1)$coefficients)
cat("\n=== Model 2: + age + stage (n=", nrow(multi_data), ") ===\n"); print(summary(model2)$coefficients)
cat("\n=== Model 3: + age + stage + subtype (n=", nrow(multi_data), ") ===\n"); print(summary(model3)$coefficients)

lrt_subtype <- anova(model2, model3)
cat("\n=== LRT: does molecular subtype improve Model 2? ===\n"); print(lrt_subtype)

extract_risk_row <- function(model, model_name) {
  tidy_model <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "risk_score")
  data.frame(
    Model = model_name, n = model$n, events = model$nevent,
    HR = tidy_model$estimate, lower95 = tidy_model$conf.low, upper95 = tidy_model$conf.high,
    p.value = tidy_model$p.value, Concordance = summary(model)$concordance[1]
  )
}
subtype_cox_summary <- rbind(
  extract_risk_row(model1, "Model 1: risk_score alone"),
  extract_risk_row(model2, "Model 2: + age + stage"),
  extract_risk_row(model3, "Model 3: + age + stage + subtype")
)
print(subtype_cox_summary, digits = 4)

saveRDS(list(model1 = model1, model2 = model2, model3 = model3, lrt = lrt_subtype),
        "results1/subtype_nested_cox_models.rds")
write.csv(subtype_cox_summary, "results1/subtype_adjusted_cox_summary.csv", row.names = FALSE)

multi_cox <- model3
summary(multi_cox)
saveRDS(multi_cox, "results1/multivariable_cox.rds")

ph_test <- cox.zph(multi_cox)
ph_test

saveRDS(ph_test, "results1/ph_test.rds")

par(mfrow = c(2, 2))
plot(ph_test)
par(mfrow = c(1, 1))

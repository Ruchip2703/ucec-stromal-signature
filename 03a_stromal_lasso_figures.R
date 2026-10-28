# 03a_stromal_lasso_figures_WL.R
#
# Main figures for the signature: stromal score distribution, volcano, heatmaps, forest plots, LASSO path/CV, KM, AUROC, Cox forest plots, nomogram, calibration and DCA.
# Input: results1/ (from 02, 03)  |  Output: plots1/

library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(survival)
library(survminer)
library(glmnet)
library(timeROC)
library(rms)
library(dcurves)
library(broom)

dir.create("plots1", showWarnings = FALSE, recursive = TRUE)

# LOAD ALL REQUIRED INPUTS EXPLICITLY (from results1/)
tum_est_df <- readRDS("results1/tumor_estimate_scores_grouped.rds")  # 545 x 7,
stromal_deg <- readRDS("results1/stromal_deg_full.rds")      # 19,616 x 6
stromal_sig_deg <- readRDS("results1/stromal_sig_deg.rds")   # 501 x 6
expr_tum_df <- readRDS("results1/expr_tumor_vst_genesymbol.rds")    # 19,616 x 545
stromal_deg_expr <- readRDS("results1/stromal_deg_expr_matrix.rds") # 541 x 501

uni_sig_genes <- readRDS("results1/univariate_significant_genes_train.rds")  # 32 x 6
cvfit <- readRDS("results1/cvfit.rds")                                  # LASSO cv.glmnet fit
lasso_coef_df <- readRDS("results1/lasso_coefficients_dataframe.rds")   # 9 x 2
cutoff <- readRDS("results1/risk_score_cutoff.rds")                     # single number

train_survival <- readRDS("results1/train_survival_with_risk.rds")     # 405 x n
test_survival <- readRDS("results1/test_survival_with_risk.rds")       # 136 x n
survival_patient <- readRDS("results1/full_survival_with_risk.rds")    # 541 x n, full cohort

roc_train <- readRDS("results1/roc_train.rds")
roc_test <- readRDS("results1/roc_test.rds")
roc_full <- readRDS("results1/roc_full.rds")

multi_data      <- readRDS("results1/multi_data.rds")            # 530 x 15 (subtype-complete subset) -- Figure 4A/4B/4C
multi_data_base <- readRDS("results1/multi_data_base_n539.rds")  # 539 x 14 (age+stage-complete, no subtype filter) -- Figure 5A/5B/5C
multi_cox <- readRDS("results1/multivariable_cox.rds")  # final multivariable Cox model

n_sig_genes <- nrow(lasso_coef_df)
cat("[CHECK] Final signature gene count:", n_sig_genes, "\n")   # expect 9

# STROMAL SCORE DISTRIBUTION
group_counts <- table(tum_est_df$stromal_group)

tum_est_df$stromal_group_label <- factor(
  tum_est_df$stromal_group,
  levels = c("High", "Low"),
  labels = c(
    paste0("High (n=", group_counts["High"], ")"),
    paste0("Low (n=", group_counts["Low"], ")")
  )
)

p_stromal_dist <- ggplot(tum_est_df, aes(x = stromal_group_label, y = StromalScore, fill = stromal_group_label)) +
  geom_boxplot(width = 0.55, alpha = 0.9, outlier.shape = NA, color = "black", linewidth = 0.7) +
  geom_jitter(width = 0.18, size = 1, alpha = 0.25, color = "black") +
  scale_fill_manual(values = c("salmon", "#6BBFB0")) +
  labs(title = "Distribution of Stromal Scores", x = "", y = "ESTIMATE Stromal Score") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "plain"),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none"
  )

p_stromal_dist

ggsave("plots1/STROMAL SCORE DISTRIBUTION.pdf", p_stromal_dist, width = 6, height = 5, units = "in")
ggsave("plots1/STROMAL SCORE DISTRIBUTION.png", p_stromal_dist, width = 6, height = 5, units = "in", dpi = 300)

# VOLCANO PLOT OF STROMAL DEGs
stromal_deg$significance <- "Not Significant"
stromal_deg$significance[stromal_deg$adj.P.Val < 0.05 & stromal_deg$logFC > 1] <- "Upregulated"
stromal_deg$significance[stromal_deg$adj.P.Val < 0.05 & stromal_deg$logFC < -1] <- "Downregulated"

p_volcano_deg <- ggplot(stromal_deg, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5", "Not Significant" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic(base_size = 13) +
  labs(title = "Stromal-Associated Differentially Expressed Genes", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.9, size = 12, face = "bold"), legend.title = element_blank())

p_volcano_deg

ggsave("plots1/Stromal-Associated Differentially Expressed Genes.pdf", p_volcano_deg, width = 6, height = 5, units = "in")
ggsave("plots1/Stromal-Associated Differentially Expressed Genes.png", p_volcano_deg, width = 6, height = 5, units = "in", dpi = 300)

# TOP STROMAL DEG HEATMAP (top 15 up + top 15 down)
top_up <- stromal_sig_deg %>% filter(logFC > 1) %>% arrange(adj.P.Val) %>% head(15)
top_down <- stromal_sig_deg %>% filter(logFC < -1) %>% arrange(adj.P.Val) %>% head(15)
top_genes <- rbind(top_up, top_down)
dim(top_genes)

heatmap_matrix <- expr_tum_df[rownames(expr_tum_df) %in% rownames(top_genes), ]
heatmap_matrix <- heatmap_matrix[rownames(top_genes), ]
dim(heatmap_matrix)

heatmap_patient <- substr(colnames(heatmap_matrix), 1, 12)
match_idx <- match(heatmap_patient, tum_est_df$patient_id)

heatmap_matrix <- heatmap_matrix[, !is.na(match_idx)]
match_idx <- match(substr(colnames(heatmap_matrix), 1, 12), tum_est_df$patient_id)

sample_order <- order(tum_est_df$stromal_group[match_idx])
heatmap_matrix <- heatmap_matrix[, sample_order]
colnames(heatmap_matrix) <- substr(colnames(heatmap_matrix), 1, 12)

annotation_col <- data.frame(StromalGroup = tum_est_df$stromal_group[match_idx][sample_order])
rownames(annotation_col) <- colnames(heatmap_matrix)

stopifnot(all(colnames(heatmap_matrix) == rownames(annotation_col)))

hm <- pheatmap(
  heatmap_matrix, scale = "row", show_colnames = FALSE, fontsize_row = 9,
  clustering_method = "ward.D2", cluster_cols = FALSE, annotation_col = annotation_col,
  annotation_colors = list(StromalGroup = c(High = "salmon", Low = "#6BBFB0")),
  color = colorRampPalette(c("#3B6FB6", "#F7F7F7", "#D73027"))(100),
  border_color = NA, main = "Top Stromal-Associated Differentially Expressed Genes"
)

pdf("plots1/Top Stromal-Associated Differentially Expressed Genes.pdf", width = 10, height = 7)
grid::grid.draw(hm$gtable)
dev.off()

png("plots1/Top Stromal-Associated Differentially Expressed Genes.png", width = 10, height = 7, units = "in", res = 300)
grid::grid.draw(hm$gtable)
dev.off()

# FOREST PLOT -- TOP PROGNOSTIC STROMAL GENES (post univariate Cox)
top_risk <- uni_sig_genes %>% filter(HR > 1) %>% arrange(p.value) %>% head(7)
top_protective <- uni_sig_genes %>% filter(HR < 1) %>% arrange(p.value) %>% head(8)
forest_df <- rbind(top_risk, top_protective)
forest_df$Gene <- factor(forest_df$Gene, levels = rev(forest_df$Gene))

p_forest_unicox <- ggplot(forest_df, aes(x = HR, y = Gene, color = HR > 1)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.2, linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("#009E73", "#D55E00")) +
  labs(title = "Top Prognostic Stromal-Associated Genes", x = "Hazard Ratio", y = "") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.title.x = element_text(face = "bold"), axis.text = element_text(color = "black"),
        legend.position = "none")

p_forest_unicox

ggsave("plots1/Top Prognostic Stromal-Associated Genes after uni cox.pdf", p_forest_unicox, width = 6, height = 6, units = "in")
ggsave("plots1/Top Prognostic Stromal-Associated Genes after uni cox.png", p_forest_unicox, width = 6, height = 6, units = "in", dpi = 300)

# LASSO COEFFICIENT PATH
plot(cvfit$glmnet.fit, xvar = "lambda", label = FALSE, lwd = 1)
abline(v = -log(cvfit$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = -log(cvfit$lambda.1se), col = "blue", lty = 2, lwd = 2)
legend("bottomleft", legend = c("lambda.min", "lambda.1se"), col = c("red", "blue"), lty = 2, lwd = 2, cex = 0.5)

dev.copy(pdf, "plots1/LASSO Coefficient Path.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/LASSO Coefficient Path.png", width = 6, height = 5, units = "in", res = 300)
dev.off()

# CROSS-VALIDATION CURVE
library(glmnet)

# PDF version
pdf("plots1/Cross Validation Curve lasso.pdf", width = 7, height = 5.5)
par(mar = c(5, 5, 4, 6) + 0.1)  # extra room on the right (6th value) for the legend
plot(cvfit)
abline(v = -log(cvfit$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = -log(cvfit$lambda.1se), col = "blue", lty = 2, lwd = 2)
legend("bottomright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, lwd = 2, cex = 0.7,
       bty = "n", xpd = TRUE)  # xpd=TRUE lets legend sit outside the plot box if needed
dev.off()

# PNG version
png("plots1/Cross Validation Curve lasso.png", width = 7, height = 5.5, units = "in", res = 300)
par(mar = c(5, 5, 4, 6) + 0.1)
plot(cvfit)
abline(v = -log(cvfit$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = -log(cvfit$lambda.1se), col = "blue", lty = 2, lwd = 2)
legend("bottomright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, lwd = 2, cex = 0.7,
       bty = "n", xpd = TRUE)
dev.off()

# FINAL SIGNATURE BARPLOT (coefficients)
p_sig_barplot <- ggplot(lasso_coef_df, aes(x = reorder(Gene, Coefficient), y = Coefficient,
                                           fill = ifelse(Coefficient > 0, "Risk-associated", "Protective"))) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Protective" = "#6BBFB0", "Risk-associated" = "salmon"),
                    name = NULL) +
  labs(title = "Final Stromal-Associated Prognostic Signature", x = "", y = "Coefficient") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title.y = element_blank(), axis.title.x = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom")

p_sig_barplot

ggsave("plots1/lasso gene Prognostic Signature.pdf", p_sig_barplot, width = 8, height = 6, units = "in")
ggsave("plots1/lasso gene Prognostic Signature.png", p_sig_barplot, width = 7, height = 6, units = "in", dpi = 300)

# KM CURVES (TRAIN / TEST / FULL) -- publication-styled versions
km_theme <- theme_classic(base_size = 13)

fit_train <- survfit(Surv(OS.time, OS) ~ risk_group, data = train_survival)
ggsurvplot(
  fit_train, data = train_survival, pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table = TRUE, palette = c("#4DBBD5", "#E64B35"), legend.title = "",
  legend.labs = c("Low Risk", "High Risk"), title = "TRAIN Cohort Survival",
  xlab = "Time (Days)", ylab = "Overall Survival Probability", break.time.by = 1000, xlim = c(0, 6000),
  risk.table.height = 0.3, risk.table.y.text.col = TRUE, risk.table.y.text = TRUE,
  ggtheme = km_theme, font.title = c(15, "bold"), font.x = c(12, "bold"), font.y = c(12, "bold"),
  censor.size = 2.5, lwd = 1.6
)
dev.copy(pdf, "plots1/Train Km.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Train KM.png", width = 6, height = 5, units = "in", res = 300); dev.off()

fit_test <- survfit(Surv(OS.time, OS) ~ risk_group, data = test_survival)
ggsurvplot(
  fit_test, data = test_survival, pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table = TRUE, palette = c("#4DBBD5", "#E64B35"), legend.title = "",
  legend.labs = c("Low Risk", "High Risk"), title = "TEST Cohort Survival",
  xlab = "Time (Days)", ylab = "Overall Survival Probability", break.time.by = 1000, xlim = c(0, 4000),
  risk.table.height = 0.3, risk.table.y.text.col = TRUE, risk.table.y.text = TRUE,
  ggtheme = km_theme, font.title = c(15, "bold"), font.x = c(12, "bold"), font.y = c(12, "bold"),
  censor.size = 2.5, lwd = 1.6
)
dev.copy(pdf, "plots1/Test KM.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Test KM.png", width = 6, height = 5, units = "in", res = 300); dev.off()

fit_full <- survfit(Surv(OS.time, OS) ~ risk_group, data = survival_patient)
ggsurvplot(
  fit_full, data = survival_patient, pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table = TRUE, palette = c("#4DBBD5", "#E64B35"), legend.title = "",
  legend.labs = c("Low Risk", "High Risk"), title = "Full TCGA Cohort Survival",
  xlab = "Time (Days)", ylab = "Overall Survival Probability", break.time.by = 1000, xlim = c(0, 6000),
  risk.table.height = 0.3, risk.table.y.text.col = TRUE, risk.table.y.text = TRUE,
  ggtheme = km_theme, font.title = c(15, "bold"), font.x = c(12, "bold"), font.y = c(12, "bold"),
  censor.size = 2.5, lwd = 1.6
)
dev.copy(pdf, "plots1/Full KM.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Full KM.png", width = 6, height = 5, units = "in", res = 300); dev.off()

# TIME-DEPENDENT AUROC (FULL / TEST / TRAIN)
plot(roc_full, time = 365, col = "#E64B35", lwd = 2, title = FALSE)
plot(roc_full, time = 1095, col = "#4DBBD5", lwd = 2, add = TRUE)
plot(roc_full, time = 1825, col = "#00A020", lwd = 2, add = TRUE)
title(main = "Time-dependent ROC Curve (Full Cohort)")
legend("bottomright", legend = c(
  paste0("1-year AUC = ", round(roc_full$AUC[1], 3)),
  paste0("3-year AUC = ", round(roc_full$AUC[2], 3)),
  paste0("5-year AUC = ", round(roc_full$AUC[3], 3))
), col = c("#E64B35", "#4DBBD5", "#00A020"), lwd = 2, cex = 0.7, bty = "n")
dev.copy(pdf, "plots1/Full AUROC.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Full AUROC.png", width = 6, height = 5, units = "in", res = 300); dev.off()

plot(roc_test, time = 365, col = "#E64B35", lwd = 2, title = FALSE)
plot(roc_test, time = 1095, col = "#4DBBD5", lwd = 2, add = TRUE)
plot(roc_test, time = 1825, col = "#00A020", lwd = 2, add = TRUE)
title(main = "Time-dependent ROC Curve (Test Cohort)")
legend("bottomright", legend = c(
  paste0("1-year AUC = ", round(roc_test$AUC[1], 3)),
  paste0("3-year AUC = ", round(roc_test$AUC[2], 3)),
  paste0("5-year AUC = ", round(roc_test$AUC[3], 3))
), col = c("#E64B35", "#4DBBD5", "#00A020"), lwd = 2, cex = 0.7, bty = "n")
dev.copy(pdf, "plots1/Test AUROC.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Test AUROC.png", width = 6, height = 5, units = "in", res = 300); dev.off()

plot(roc_train, time = 365, col = "#E64B35", lwd = 2, title = FALSE)
plot(roc_train, time = 1095, col = "#4DBBD5", lwd = 2, add = TRUE)
plot(roc_train, time = 1825, col = "#00A020", lwd = 2, add = TRUE)
title(main = "Time-dependent ROC Curve (Train Cohort)")
legend("bottomright", legend = c(
  paste0("1-year AUC = ", round(roc_train$AUC[1], 3)),
  paste0("3-year AUC = ", round(roc_train$AUC[2], 3)),
  paste0("5-year AUC = ", round(roc_train$AUC[3], 3))
), col = c("#E64B35", "#4DBBD5", "#00A020"), lwd = 2, cex = 0.7, bty = "n")
dev.copy(pdf, "plots1/Train AUROC.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Train AUROC.png", width = 6, height = 5, units = "in", res = 300); dev.off()

# RISK SCORE DISTRIBUTION (full cohort, ordered)
risk_order <- survival_patient[order(survival_patient$risk_score), ]
risk_order$Index <- 1:nrow(risk_order)

par(bg = "white")
plot(
  risk_order$Index, risk_order$risk_score,
  col = ifelse(risk_order$risk_group == "High", "salmon", "#6BBFB0"),
  pch = 16, cex = 0.55, xlab = "Patients", ylab = "Risk score",
  main = "Risk Score Distribution", cex.lab = 1.0, cex.main = 1.15,
  font.lab = 1, font.main = 2, las = 1, bty = "l"
)
abline(h = cutoff, lty = 2, lwd = 1.5, col = "grey40")
abline(v = sum(risk_order$risk_group == "Low"), lty = 3, lwd = 1.2, col = "grey55")

dev.copy(pdf, "plots1/Risk Score Distr.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/Risk Score Distr.png", width = 6, height = 5, units = "in", res = 300); dev.off()

# SURVIVAL STATUS DISTRIBUTION
plot(
  risk_order$Index, risk_order$OS.time,
  col = ifelse(risk_order$OS == 1, "salmon", "#6BBFB0"),
  pch = 16, cex = 0.6, xlab = "Patients", ylab = "Survival Time (Days)",
  main = "Survival Status Distribution", cex.lab = 1.0, cex.main = 1.15, font.lab = 1, font.main = 2
)
legend("topright", legend = c("Alive", "Dead"), col = c("#6BBFB0", "salmon"), pch = 16, bty = "n", cex = 0.6)

dev.copy(pdf, "plots1/survival status.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/survival status.png", width = 6, height = 5, units = "in", res = 300); dev.off()

# HEATMAP OF FINAL SIGNATURE
rownames(stromal_deg_expr) <- substr(rownames(stromal_deg_expr), 1, 12)

risk_order2 <- survival_patient %>% arrange(risk_group, risk_score)
final_genes <- lasso_coef_df$Gene

heatmap_matrix2 <- stromal_deg_expr[risk_order2$patient_id, final_genes]
heatmap_matrix2 <- t(heatmap_matrix2)

annotation_col2 <- data.frame(RiskGroup = factor(risk_order2$risk_group, levels = c("Low", "High")))
rownames(annotation_col2) <- risk_order2$patient_id

stopifnot(all(colnames(heatmap_matrix2) == rownames(annotation_col2)))

hm2 <- pheatmap(
  heatmap_matrix2, scale = "row", cluster_cols = FALSE, show_colnames = FALSE,
  fontsize = 9, fontsize_row = 10, annotation_col = annotation_col2,
  annotation_colors = list(RiskGroup = c(Low = "#6BBFB0", High = "salmon")),
  clustering_method = "ward.D2", color = colorRampPalette(c("#3B6FB6", "#F7F7F7", "#D73027"))(100),
  border_color = NA, main = "Final Prognostic Signature"
)

pdf("plots1/heatmap of prognostic signature.pdf", width = 10, height = 7)
grid::grid.draw(hm2$gtable)
dev.off()

png("plots1/heatmap of prognostic signature.png", width = 10, height = 7, units = "in", res = 300)
grid::grid.draw(hm2$gtable)
dev.off()

# UNIVARIATE COX (risk_score, age, stage_group)
variables <- c("risk_score", "age", "stage_group")

uni_clinical <- list()
for (var in variables) {
  formula_uni <- as.formula(paste("Surv(OS.time, OS) ~", var))
  model_uni <- coxph(formula_uni, data = multi_data)
  summary_uni <- summary(model_uni)
  uni_clinical[[var]] <- data.frame(
    Variable = var,
    HR = summary_uni$coefficients[, "exp(coef)"],
    lower95 = summary_uni$conf.int[, "lower .95"],
    upper95 = summary_uni$conf.int[, "upper .95"],
    p.value = summary_uni$coefficients[, "Pr(>|z|)"]
  )
}

uni_clinical_df <- do.call(rbind, uni_clinical)
uni_clinical_df$Variable <- c("Risk score", "Age", "Advanced stage")
uni_clinical_df$Variable <- factor(uni_clinical_df$Variable, levels = rev(uni_clinical_df$Variable))

xlim_lower <- max(0.1, min(uni_clinical_df$lower95) * 0.9)
xlim_upper <- max(uni_clinical_df$upper95) * 1.1

p_unicox_4a <- ggplot(uni_clinical_df, aes(x = HR, y = Variable)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.15, linewidth = 0.5, color = "black") +
  geom_point(size = 2, color = "#E64B35") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  coord_cartesian(xlim = c(xlim_lower, xlim_upper)) +
  theme_classic(base_size = 12) +
  labs(title = "Univariate Cox Regression", x = "Hazard ratio", y = NULL) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
        axis.text = element_text(color = "black", size = 11), axis.title.x = element_text(size = 12))

p_unicox_4a

ggsave("plots1/unicox after lasso.pdf", p_unicox_4a, width = 6, height = 5, units = "in")
ggsave("plots1/unicox after lasso.png", p_unicox_4a, width = 6, height = 5, units = "in", dpi = 300)

# MULTIVARIABLE COX -- AGE + STAGE (n=530, subtype-complete subset)
cox_models <- readRDS("results1/subtype_nested_cox_models.rds")
multi_cox_A <- cox_models$model2   # age + stage only

summary_A <- summary(multi_cox_A)
raw_terms_A <- rownames(summary_A$coefficients)

nice_label <- function(term) {
  if (term == "risk_score") return("Risk Score")
  if (term == "age") return("Age")
  if (term == "stage_groupAdvanced") return("Advanced Stage")
  if (grepl("^Subtype_Selected", term)) {
    return(paste0("Subtype: ", sub("^Subtype_SelectedUCEC\\.", "", term)))
  }
  term
}

df_A <- data.frame(
  Variable = raw_terms_A,
  HR = summary_A$coefficients[, "exp(coef)"],
  lower95 = summary_A$conf.int[, "lower .95"],
  upper95 = summary_A$conf.int[, "upper .95"],
  p.value = summary_A$coefficients[, "Pr(>|z|)"]
)
df_A$Variable <- vapply(raw_terms_A, nice_label, character(1))
df_A$Variable <- factor(df_A$Variable, levels = rev(df_A$Variable))

p_multicox_4A <- ggplot(df_A, aes(x = HR, y = Variable)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.15, linewidth = 0.5, color = "black") +
  geom_point(size = 2, color = "#E64B35") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(title = "Multivariable Cox Regression (age + stage)", x = "Hazard Ratio", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
        axis.text = element_text(color = "black", size = 11), axis.title.x = element_text(size = 12))

p_multicox_4A
ggsave("plots1/multicox age+stage.pdf", p_multicox_4A, width = 7, height = 5.5, units = "in")
ggsave("plots1/multicox age+stage.png", p_multicox_4A, width = 7, height = 5.5, units = "in", dpi = 300)

# MULTIVARIABLE COX -- AGE + STAGE + SUBTYPE (n=530, same subset)
multi_cox_B <- cox_models$model3   # risk_score + age + stage + subtype, n=530
summary_B <- summary(multi_cox_B)
raw_terms_B <- rownames(summary_B$coefficients)

df_B <- data.frame(
  Variable = raw_terms_B,
  HR = summary_B$coefficients[, "exp(coef)"],
  lower95 = summary_B$conf.int[, "lower .95"],
  upper95 = summary_B$conf.int[, "upper .95"],
  p.value = summary_B$coefficients[, "Pr(>|z|)"]
)
df_B$Variable <- vapply(raw_terms_B, nice_label, character(1))
df_B$Variable <- factor(df_B$Variable, levels = rev(df_B$Variable))

p_multicox_4B <- ggplot(df_B, aes(x = HR, y = Variable)) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.15, linewidth = 0.5, color = "black") +
  geom_point(size = 2, color = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  theme_classic(base_size = 13) +
  labs(title = "Multivariable Cox Regression (subtype-adjusted)", x = "Hazard Ratio", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
        axis.text = element_text(color = "black", size = 11), axis.title.x = element_text(size = 12))

p_multicox_4B
ggsave("plots1/multicox subtype-adjusted.pdf", p_multicox_4B, width = 7, height = 5.5, units = "in")
ggsave("plots1/multicox subtype-adjusted.png", p_multicox_4B, width = 7, height = 5.5, units = "in", dpi = 300)

summary_A
summary_B

# RISK SCORE BY CLINICAL STAGE
set.seed(123)

p_stage_4c <- ggplot(multi_data, aes(x = stage_group, y = risk_score, fill = stage_group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, linewidth = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.18, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4, hjust = -1) +
  coord_cartesian(ylim = c(min(multi_data$risk_score), max(multi_data$risk_score) + 0.3)) +
  scale_fill_manual(values = c(Early = "#90EE90", Advanced = "#E76F51")) +
  labs(title = "Risk score by Clinical stage", x = NULL, y = "Risk Score") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
        axis.text = element_text(color = "black", size = 11), axis.title.y = element_text(size = 12),
        legend.position = "none")

p_stage_4c

ggsave("plots1/risk score distr by stage.pdf", p_stage_4c, width = 6, height = 5, units = "in")
ggsave("plots1/risk score distr by stage.png", p_stage_4c, width = 6, height = 5, units = "in", dpi = 300)

# NOMOGRAM
# NOTE: deliberately built on risk_score + age + stage_group only (NOT the
# subtype-adjusted model from Figure 4B). A nomogram is a practical clinical
# tool meant to be usable with routinely available inputs -- molecular
# subtyping isn't always available at the point of care, so keeping the
# simpler 3-variable version here is an intentional design choice, not an
# oversight. The subtype-adjusted result belongs in the Results/Table text
# (Figure 4B) as the evidence that risk_score's prognostic value holds
# beyond subtype; the nomogram is the deployable tool built from what's
# realistically available.
library(Hmisc)

label(multi_data_base$risk_score)  <- "Risk Score"
label(multi_data_base$age)         <- "Age"
label(multi_data_base$stage_group) <- "Stage Group"

dd <- datadist(multi_data_base)
options(datadist = "dd")

cox_rms <- cph(
  Surv(OS.time, OS) ~ risk_score + age + stage_group,
  data = multi_data_base, x = TRUE, y = TRUE, surv = TRUE
)

surv <- Survival(cox_rms)

nom <- nomogram(
  cox_rms,
  fun = list(function(x) surv(365, x), function(x) surv(1095, x), function(x) surv(1825, x)),
  funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"),
  lp = FALSE
)

print(nom)
print(cox_rms)

pdf("plots1/Nomogram1.pdf", width = 11, height = 7)
par(mar = c(5, 5, 4, 2))
plot(nom, xfrac = 0.28, cex.axis = 0.7, cex.var = 1.0, col.grid = grey(0.9), lmgp = 0.3, tcl = -0.3)
title(main = "Nomogram for Predicting Overall Survival", font.main = 1, cex.main = 1)
dev.off()

png("plots1/Nomogram1.png", width = 11, height = 7, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))
plot(nom, xfrac = 0.28, cex.axis = 0.7, cex.var = 1.0, col.grid = grey(0.9), lmgp = 0.3, tcl = -0.3)
title(main = "Nomogram for Predicting Overall Survival", font.main = 1, cex.main = 1)
dev.off()

# CALIBRATION (1-/3-/5-YEAR)
set.seed(123)
cox_rms_1yr <- cph(Surv(OS.time, OS) ~ risk_score + age + stage_group, data = multi_data_base,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
cal1 <- calibrate(cox_rms_1yr, cmethod = "KM", method = "boot", u = 365, m = 70, B = 1000)
cal1
plot(cal1, xlab = "Predicted 1-year Survival", ylab = "Observed 1-year Survival", lwd = 2, subtitles = FALSE)
dev.copy(pdf, "plots1/cal1.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/cal1.png", width = 6, height = 5, units = "in", res = 300); dev.off()

set.seed(123)
cox_rms_3yr <- cph(Surv(OS.time, OS) ~ risk_score + age + stage_group, data = multi_data_base,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095)
cal3 <- calibrate(cox_rms_3yr, cmethod = "KM", method = "boot", u = 1095, m = 70, B = 1000)
cal3
plot(cal3, xlab = "Predicted 3-year Survival", ylab = "Observed 3-year Survival", lwd = 2, subtitles = FALSE)
dev.copy(pdf, "plots1/cal3.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/cal3.png", width = 6, height = 5, units = "in", res = 300); dev.off()

set.seed(123)
cox_rms_5yr <- cph(Surv(OS.time, OS) ~ risk_score + age + stage_group, data = multi_data_base,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = 1825)
cal5 <- calibrate(cox_rms_5yr, cmethod = "KM", method = "boot", u = 1825, m = 70, B = 1000)
cal5
plot(cal5, xlab = "Predicted 5-year Survival", ylab = "Observed 5-year Survival", lwd = 2, subtitles = FALSE)
dev.copy(pdf, "plots1/cal5.pdf", width = 6, height = 5); dev.off()
dev.copy(png, filename = "plots1/cal5.png", width = 6, height = 5, units = "in", res = 300); dev.off()

# DECISION CURVE ANALYSIS
cox_dca <- coxph(Surv(OS.time, OS) ~ risk_score + age + stage_group, data = multi_data_base, x = TRUE)

multi_data_dca <- broom::augment(
  cox_dca, newdata = multi_data_base %>% mutate(OS.time = 1825), type.predict = "expected"
) %>%
  mutate(pred_5yr_event = 1 - exp(-.fitted)) %>%
  dplyr::select(-.fitted, -.se.fit)

summary(multi_data_dca$pred_5yr_event)
range(multi_data_dca$pred_5yr_event)

dca_model <- dca(
  Surv(OS.time, OS) ~ pred_5yr_event, data = multi_data_dca, time = 1825,
  thresholds = seq(0.01, 0.40, 0.01), label = list(pred_5yr_event = "Nomogram")
)

pdf("plots1/DCA.pdf", width = 6, height = 5)
plot(dca_model, smooth = TRUE, linewidth = 1.5)
dev.off()

png("plots1/DCA.png", width = 6, height = 5, units = "in", res = 300)
plot(dca_model, smooth = TRUE, linewidth = 1.5)
dev.off()

# Reset the global datadist option set above (Figure 5A), so it doesn't leak
# into any later interactive session or script run after this one.
options(datadist = NULL)

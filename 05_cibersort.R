# 05_cibersort_WL.R
#
# CIBERSORT (LM22) immune deconvolution and association with risk group and risk score.
# Requires CIBERSORT.R and LM22.txt in the working directory (obtained from the CIBERSORT website; not redistributed here).
# Input: results1/expr_tumor_vst_genesymbol.rds, results1/full_survival_with_risk.rds  |  Output: results1/, plots1/

set.seed(123)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(tidyr)
library(ggsignif)
library(rstatix)
library(coin)     # required by rstatix::wilcox_effsize()
library(corrplot)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)
dir.create("plots1", showWarnings = FALSE, recursive = TRUE)

col_low            <- "#3B6FB6"   # blue   -- low risk
col_high           <- "#D73027"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

col_categorical_8 <- c(
  "#2a78d6",  # 1 blue
  "#eb6834",  # 2 orange
  "#1baf7a",  # 3 aqua
  "#eda100",  # 4 yellow
  "#e87ba4",  # 5 magenta
  "#008300",  # 6 green
  "#4a3aa7",  # 7 violet
  "#e34948"   # 8 red
)

# LOAD EXPRESSION MATRIX + RISK-SCORED SURVIVAL DATA
expr_tumor_ciber <- readRDS("results1/expr_tumor_vst_genesymbol.rds")
survival_patient_ciber <- readRDS("results1/full_survival_with_risk.rds")
dim(expr_tumor_ciber)           # expect 19,616 x 545
dim(survival_patient_ciber)     # expect 541 x 9

# MATCH PATIENT IDs, SUBSET TO COMMON PATIENTS
colnames(expr_tumor_ciber) <- substr(colnames(expr_tumor_ciber), 1, 12)

common_patients_ciber <- intersect(colnames(expr_tumor_ciber), survival_patient_ciber$patient_id)
length(common_patients_ciber)   # expect 541

expr_tumor_ciber <- expr_tumor_ciber[, common_patients_ciber]

survival_patient_ciber <- survival_patient_ciber[
  survival_patient_ciber$patient_id %in% common_patients_ciber,
]
survival_patient_ciber <- survival_patient_ciber[
  match(colnames(expr_tumor_ciber), survival_patient_ciber$patient_id),
]

stopifnot(all(colnames(expr_tumor_ciber) == survival_patient_ciber$patient_id))

# FILTER TO LM22 SIGNATURE GENES, WRITE CIBERSORT INPUT
lm22_ciber <- read.table("LM22.txt", header = TRUE, sep = "\t", row.names = 1)
dim(lm22_ciber)   # expect 547 x 22 (LM22 standard reference)

common_genes_ciber <- intersect(rownames(expr_tumor_ciber), rownames(lm22_ciber))
length(common_genes_ciber)

tumor_filtered_ciber <- expr_tumor_ciber[common_genes_ciber, ]
dim(tumor_filtered_ciber)

write.table(
  tumor_filtered_ciber,
  file = "CIBERSORT_tumor_filtered_updated1.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

# RUN CIBERSORT (standalone CIBERSORT.R implementation, LM22 reference)
source("CIBERSORT.R")

cibersort_res_risk <- CIBERSORT(
  "LM22.txt",
  "CIBERSORT_tumor_filtered_updated1.txt",
  perm = 100,
  QN = FALSE   # No quantile normalization -- correct choice for RNA-seq data
)
dim(cibersort_res_risk)   # expect 541 x 25 (22 cell types + P-value, Correlation, RMSE)

saveRDS(cibersort_res_risk, "results1/cibersort_res_risk.rds")

cibersort_res_risk <- readRDS("results1/cibersort_res_risk.rds")

cibersort_df_risk <- as.data.frame(cibersort_res_risk)

# DECONVOLUTION QUALITY CHECK
table(cibersort_df_risk$`P-value` < 0.05)
mean(cibersort_df_risk$`P-value` < 0.05)   # fraction of patients with significant deconvolution
summary(cibersort_df_risk$RMSE)
summary(cibersort_df_risk$Correlation)

cibersort_df_risk$patient_id <- rownames(cibersort_df_risk)

# MERGE WITH RISK METADATA
ciber_merged_risk <- merge(survival_patient_ciber, cibersort_df_risk, by = "patient_id")
dim(ciber_merged_risk)   # expect 541 x (9 + 25 - 1 dup patient_id col)
table(ciber_merged_risk$risk_group)
sum(is.na(ciber_merged_risk))   # expect 0

saveRDS(ciber_merged_risk, "results1/ciber_merged_risk.rds")

immune_cols_risk <- colnames(cibersort_df_risk)[1:22]
length(immune_cols_risk)   # expect 22

# RISK-GROUP COMPARISON (Wilcoxon + effect size + FDR) -- ALL 22 CELLS
pvals_risk <- sapply(immune_cols_risk, function(cell) {
  wilcox.test(ciber_merged_risk[[cell]] ~ ciber_merged_risk$risk_group)$p.value
})
p_adj_risk <- p.adjust(pvals_risk, method = "BH")

effect_sizes <- sapply(immune_cols_risk, function(cell) {
  rstatix::wilcox_effsize(
    data = ciber_merged_risk,
    formula = as.formula(paste0("`", cell, "` ~ risk_group"))
  )$effsize
})
effect_sizes <- unname(effect_sizes)

results_risk <- data.frame(
  Cell = immune_cols_risk,
  Effect_size = effect_sizes,
  P_value = pvals_risk,
  FDR = p_adj_risk,
  row.names = NULL
) %>% arrange(FDR)

results_risk
cat("Significant immune cells (FDR < 0.05):", sum(results_risk$FDR < 0.05), "\n")

write.csv(results_risk, "results1/CIBERSORT_group_comparison.csv", row.names = FALSE)

# CORRELATION WITH CONTINUOUS RISK SCORE (Spearman) -- ALL 22 CELLS
cor_results_risk <- data.frame()
for (cell in immune_cols_risk) {
  test <- cor.test(ciber_merged_risk[[cell]], ciber_merged_risk$risk_score, method = "spearman")
  cor_results_risk <- rbind(cor_results_risk,
                            data.frame(Cell = cell, Correlation = test$estimate, P_value = test$p.value))
}
cor_results_risk$FDR <- p.adjust(cor_results_risk$P_value, method = "BH")
cor_results_risk <- cor_results_risk[order(cor_results_risk$FDR), ]

cor_results_risk
cat("Risk-score correlated immune cells (FDR < 0.05):", sum(cor_results_risk$FDR < 0.05), "\n")

write.csv(cor_results_risk, "results1/CIBERSORT_risk_correlations.csv", row.names = FALSE)

# OVERLAP -- significant in BOTH group comparison AND score correlation
sig_group_cells <- results_risk$Cell[results_risk$FDR < 0.05]
sig_cor_cells <- cor_results_risk$Cell[cor_results_risk$FDR < 0.05]
overlap_cells <- intersect(sig_group_cells, sig_cor_cells)

length(overlap_cells)
overlap_cells

# IMMUNE HEATMAP (overlap cells, patients ordered by risk score)
col_low            <- "#3B6FB6"   # blue   -- low risk
col_high           <- "#D73027"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

col_risk_low  <- "#6BBFB0"   # seafoam -- RiskGroup annotation, low risk
col_risk_high <- "salmon"    # RiskGroup annotation, high risk

mat_risk <- cibersort_df_risk[, overlap_cells]
rownames(mat_risk) <- cibersort_df_risk$patient_id

order_idx_risk <- order(ciber_merged_risk$risk_score)
mat_risk <- mat_risk[order_idx_risk, ]
mat_risk <- t(mat_risk)

annotation_col_risk <- data.frame(
  RiskGroup = ciber_merged_risk$risk_group[order_idx_risk],
  RiskScore = ciber_merged_risk$risk_score[order_idx_risk]
)
rownames(annotation_col_risk) <- colnames(mat_risk)

ann_colors <- list(
  RiskGroup = c(Low = col_risk_low, High = col_risk_high),
  RiskScore = colorRampPalette(c("white", "maroon"))(100)
)
heatmap_colors_ciber <- colorRampPalette(c(col_low, col_diverging_mid, col_high))(100)

pheatmap(
  mat_risk, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
  annotation_col = annotation_col_risk, annotation_colors = ann_colors,
  show_colnames = FALSE, fontsize_row = 9, border_color = NA,
  color = heatmap_colors_ciber,
  breaks = seq(-2.5, 2.5, length.out = 101),
  main = "Immune Landscape Ordered by Risk Score"
)

dev.copy(pdf, "plots1/IMMUNE HEATMAP.pdf", width = 8, height = 5)
dev.off()
dev.copy(png, filename = "plots1/IMMUNE HEATMAP.png",
         width = 8, height = 5, units = "in", res = 300)
dev.off()

# IMMUNE CELL DIFFERENCES ACROSS RISK GROUPS (violin + box)

# NOTE: this uses a manually-ordered, hardcoded list of 8 cells for
cells_use_risk <- c(
  "Dendritic cells activated",
  "T cells regulatory (Tregs)",
  "T cells CD8",
  "Dendritic cells resting",
  "T cells CD4 memory activated",
  "NK cells resting",
  "Plasma cells"
)

col_low            <- "#6BBFB0"   # blue   -- low risk
col_high           <- "salmon"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

immune_plot_data_risk <- ciber_merged_risk %>%
  dplyr::select(risk_group, all_of(cells_use_risk)) %>%
  pivot_longer(-risk_group, names_to = "Cell", values_to = "Fraction")

immune_plot_data_risk$Cell <- factor(immune_plot_data_risk$Cell, levels = cells_use_risk)

ggplot(immune_plot_data_risk, aes(x = risk_group, y = Fraction, fill = risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.65, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.5, alpha = 0.9) +
  geom_jitter(width = 0.12, alpha = 0.1, size = 0.9, color = "black") +
  facet_wrap(~Cell, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.y.npc = 0.92, label.x = 1.25, size = 3) +
  scale_fill_manual(values = c(Low = col_low, High = col_high)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "plain", size = 7),
        axis.text = element_text(color = "black")) +
  labs(x = "Risk Group", y = "Estimated Cell Fraction")

dev.copy(pdf, "plots1/IMMUNE CELL DIFFERENCES ACROSS RISK GROUPS.pdf", width = 7, height = 5)
dev.off()
dev.copy(png, filename = "plots1/IMMUNE CELL DIFFERENCES ACROSS RISK GROUPS.png",
         width = 7, height = 5, units = "in", res = 300)
dev.off()

# CORRELATION SCATTERPLOTS (overlap cells vs continuous risk_score)
top_cor_cells_risk <- overlap_cells
colors_use <- col_categorical_8[((seq_along(top_cor_cells_risk) - 1) %% length(col_categorical_8)) + 1]

plot_list <- list()
for (i in seq_along(top_cor_cells_risk)) {
  cell <- top_cor_cells_risk[i]
  p <- ggplot(ciber_merged_risk, aes(x = risk_score, y = .data[[cell]])) +
    geom_point(color = "grey30", alpha = 0.25, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = colors_use[i], fill = colors_use[i],
                alpha = 0.18, linewidth = 1.3) +
    stat_cor(method = "spearman", size = 3.5, label.x.npc = "left", label.y.npc = "top") +
    theme_classic(base_size = 12) +
    labs(title = cell, x = "Risk Score", y = NULL)
  plot_list[[i]] <- p
}

ggarrange(plotlist = plot_list, ncol = 2, nrow = ceiling(length(plot_list) / 2))

dev.copy(pdf, "plots1/CORRELATION SCATTERPLOTS.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/CORRELATION SCATTERPLOTS.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# IMMUNE CELL CORRELATION MATRIX (all 22 cells, exploratory)
library(corrplot)
corr_mat <- cor(cibersort_df_risk[, immune_cols_risk], method = "spearman")

corrplot(
  corr_mat, method = "color", type = "upper", order = "hclust",
  tl.col = "black", tl.cex = 0.8, cl.cex = 0.8, cl.ratio = 0.2,
  col = colorRampPalette(c(col_high, col_diverging_mid, col_low))(100),
  mar = c(0, 0, 1, 0)
)

dev.copy(pdf, "plots1/Immune Cell Correlation Matrix.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/Immune Cell Correlation Matrix.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

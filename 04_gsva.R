# 04_gsva_WL.R
#
# GSVA of the 50 MSigDB Hallmark pathways and limma comparison of High vs Low risk groups.
# Input: results1/expr_tumor_vst_genesymbol.rds, results1/full_survival_with_risk.rds  |  Output: results1/, plots1/

library(GSVA)
library(msigdbr)
library(limma)
library(dplyr)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)
dir.create("plots1", showWarnings = FALSE, recursive = TRUE)

# COLOR PALETTE -- same muted blue/orange pair as 03a, applied consistently
# so "risk group" means the same color across the entire figure set.
col_low            <- "#3B6FB6"   # blue   -- low risk
col_high           <- "#D73027"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

# LOAD FULL VST EXPRESSION MATRIX + RISK-SCORED SURVIVAL DATA
expr_tumor_gsva <- readRDS("results1/expr_tumor_vst_genesymbol.rds")
survival_patient_gsva <- readRDS("results1/full_survival_with_risk.rds")
dim(expr_tumor_gsva)          # expect 19,616 x 545
dim(survival_patient_gsva)    # expect 541 x n (n includes risk_score, risk_group)

# Match patient IDs (VST matrix columns are full TCGA barcodes; trim to 12-char ID)
colnames(expr_tumor_gsva) <- substr(colnames(expr_tumor_gsva), 1, 12)

# SUBSET TO COMMON PATIENTS (VST matrix has 545; risk data has 541)
common_patients <- intersect(colnames(expr_tumor_gsva), survival_patient_gsva$patient_id)
length(common_patients)   # expect 541

expr_tumor_gsva <- expr_tumor_gsva[, common_patients]

survival_patient_gsva <- survival_patient_gsva %>%
  filter(patient_id %in% common_patients)

survival_patient_gsva <- survival_patient_gsva[
  match(colnames(expr_tumor_gsva), survival_patient_gsva$patient_id),
]

stopifnot(all(colnames(expr_tumor_gsva) == survival_patient_gsva$patient_id))

# LOAD HALLMARK GENE SETS (MSigDB, via msigdbr)
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
length(hallmark_list)   # expect 50

# RUN GSVA (Gaussian kernel, current gsvaParam()/gsva() API)
gsva_param <- gsvaParam(
  exprData = expr_tumor_gsva,
  geneSets = hallmark_list,
  kcdf = "Gaussian"
)

gsva_res <- gsva(gsva_param)
dim(gsva_res)   # expect 50 x 541

# MERGE GSVA SCORES WITH RISK METADATA
gsva_df <- as.data.frame(t(gsva_res))
gsva_df$patient_id <- rownames(gsva_df)

meta_df <- survival_patient_gsva %>%
  dplyr::select(patient_id, risk_score, risk_group)

final_gsva_df <- merge(meta_df, gsva_df, by = "patient_id")
dim(final_gsva_df)   # expect 541 x 53 (patient_id + risk_score + risk_group + 50 pathways)

saveRDS(final_gsva_df, "results1/final_gsva_df.rds")

# DIFFERENTIAL PATHWAY ACTIVITY (High vs Low risk, limma)
final_gsva_df <- readRDS("results1/final_gsva_df.rds")
group <- factor(final_gsva_df$risk_group, levels = c("Low", "High"))
gsva_matrix <- t(as.matrix(final_gsva_df[, -(1:3)]))   # pathways x patients

design <- model.matrix(~group)   # standard treatment-contrast design (intercept + groupHigh)

fit <- lmFit(gsva_matrix, design)
fit <- eBayes(fit)

pathway_results <- topTable(fit, coef = "groupHigh", number = Inf)
pathway_results$Pathway <- rownames(pathway_results)
pathway_results <- pathway_results %>% arrange(adj.P.Val)

sig_pathways <- pathway_results %>% filter(adj.P.Val < 0.05)
nrow(sig_pathways)

write.csv(pathway_results, "results1/pathway_results.csv", row.names = FALSE)
saveRDS(pathway_results, "results1/pathway_results.rds")

# VISUALIZATION -- full significant-pathway barplot
col_low            <- "#3B6FB6"   # blue   -- low risk
col_high           <- "#D73027"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

library(ggplot2)

sig_pathways$Pathway_clean <- gsub("HALLMARK_", "", sig_pathways$Pathway)
sig_pathways$Pathway_clean <- gsub("_", " ", sig_pathways$Pathway_clean)
sig_pathways <- sig_pathways %>% arrange(logFC)

ggplot(sig_pathways, aes(x = logFC, y = reorder(Pathway_clean, logFC), fill = logFC)) +
  geom_col(width = .8) +
  geom_vline(xintercept = 0) +
  scale_fill_gradient2(
    low = col_low, mid = col_diverging_mid, high = col_high, midpoint = 0,
    limits = c(-0.3, 0.3), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
    labels = c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3")
  ) +
  theme_bw() +
  labs(x = "logFC (High vs Low)", y = "")

dev.copy(pdf, "plots1/GSVA BARPLOT (ALL SIGNIFICANT PATHWAYS).pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/GSVA BARPLOT (ALL SIGNIFICANT PATHWAYS).png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

# REPRESENTATIVE VIOLIN PLOTS (4 headline pathways)
col_low            <- "#6BBFB0"   # blue   -- low risk
col_high           <- "salmon"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

library(tidyr)
library(ggpubr)

panel_df <- final_gsva_df %>%
  dplyr::select(
    risk_group,
    HALLMARK_ALLOGRAFT_REJECTION, HALLMARK_INTERFERON_GAMMA_RESPONSE,
    HALLMARK_G2M_CHECKPOINT, HALLMARK_E2F_TARGETS
  )

panel_long <- panel_df %>%
  pivot_longer(cols = -risk_group, names_to = "Pathway", values_to = "GSVA_score")

panel_long$Pathway <- factor(
  panel_long$Pathway,
  levels = c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
             "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS"),
  labels = c("Allograft Rejection", "IFN-γ Response", "G2M Checkpoint", "E2F Targets")
)

ggplot(panel_long, aes(risk_group, GSVA_score, fill = risk_group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.45, alpha = 0.9) +
  facet_wrap(~Pathway, ncol = 2) +
  stat_compare_means(method = "wilcox.test", label = "p", size = 3.5,
                     label.y.npc = 1.0, label.x = 1.35) +
  scale_fill_manual(values = c(Low = col_low, High = col_high)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(face = "plain", size = 10),
        axis.title.y = element_text(face = "plain", margin = margin(r = 12)),
        axis.text = element_text(color = "black"))

dev.copy(pdf, "plots1/GSVA Representative Pathways.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/GSVA Representative Pathways.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

# TOP-PATHWAY HEATMAP (all significant pathways, patients ordered by risk)
col_low            <- "#3B6FB6"   # blue   -- low risk
col_high           <- "#D73027"   # orange -- high risk
col_diverging_mid  <- "white"   # neutral midpoint for gradient/heatmap ramps

col_risk_low  <- "#6BBFB0"   # seafoam -- RiskGroup annotation, low risk
col_risk_high <- "salmon"    # RiskGroup annotation, high risk

library(pheatmap)

ordered_patients <- final_gsva_df %>% arrange(risk_group, risk_score)
ordered_sample_ids <- ordered_patients$patient_id

top_pathway_names <- sig_pathways$Pathway

heatmap_data <- as.matrix(final_gsva_df[, top_pathway_names])
rownames(heatmap_data) <- final_gsva_df$patient_id
heatmap_data <- t(heatmap_data)
heatmap_data <- heatmap_data[, ordered_sample_ids]

rownames(heatmap_data) <- gsub("_", " ", gsub("HALLMARK_", "", rownames(heatmap_data)))

annotation_risk <- data.frame(
  RiskGroup = ordered_patients$risk_group,
  RiskScore = ordered_patients$risk_score
)
rownames(annotation_risk) <- ordered_sample_ids

heatmap_colors <- colorRampPalette(c(col_low, col_diverging_mid, col_high))(100)

pheatmap(
  heatmap_data, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
  show_colnames = FALSE, annotation_col = annotation_risk,
  annotation_colors = list(
    RiskGroup = c(Low = col_risk_low, High = col_risk_high),
    RiskScore = colorRampPalette(c("white", "maroon"))(100)
  ),
  fontsize_row = 8, color = heatmap_colors,
  main = "GSVA Pathway Landscape Ordered by Risk"
)

dev.copy(pdf, "plots1/GSVA HEATMAP (ALL SIGNIFICANT PATHWAYS).pdf", width = 10, height = 6)
dev.off()
dev.copy(png, filename = "plots1/GSVA HEATMAP (ALL SIGNIFICANT PATHWAYS).png",
         width = 10, height = 6, units = "in", res = 300)
dev.off()

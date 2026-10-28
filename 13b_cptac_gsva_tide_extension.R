# 13b_cptac_gsva_tide_extension_WL.R
#
# CPTAC-UCEC extension: GSVA Hallmark pathways, TIDE and IPS by risk group, plus the corresponding figures.
# Input: cptac_validation/ (from 13), results1/pathway_results.rds  |  Output: cptac_validation/, plots1/cptac_figures/

library(GSVA)
library(msigdbr)
library(limma)
library(dplyr)
library(reticulate)

set.seed(123)

# GSVA HALLMARK PATHWAYS ON CPTAC
cptac_expr_raw <- readRDS("cptac_validation/cptac_expr_raw.rds")
cptac_merged_cp <- readRDS("cptac_validation/cptac_final_validation_data.rds")

# Build gene x patient matrix from CPTAC's RSEM data
cptac_expr_mat_full <- cptac_expr_raw %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  filter(!is.na(Hugo_Symbol), Hugo_Symbol != "")
rownames(cptac_expr_mat_full) <- cptac_expr_mat_full$Hugo_Symbol
cptac_expr_mat_full <- as.matrix(
  cptac_expr_mat_full %>% dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id)
)

# RSEM values are typically already non-negative counts/TPM-like; log-transform
# for GSVA, consistent with standard practice (TCGA side used VST, an
# analogous variance-stabilized log-scale transform)
cptac_expr_log <- log2(cptac_expr_mat_full + 1)

cat("[CHECK] CPTAC expression matrix for GSVA:", dim(cptac_expr_log), "\n")

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
stopifnot(length(hallmark_list) == 50)

gsva_param_cptac <- gsvaParam(
  exprData = cptac_expr_log,
  geneSets = hallmark_list,
  kcdf = "Gaussian"
)
gsva_res_cptac <- gsva(gsva_param_cptac)
cat("[CHECK] CPTAC GSVA result dim:", dim(gsva_res_cptac), "\n")   # expect 50 x ~95

gsva_df_cptac <- as.data.frame(t(gsva_res_cptac))
gsva_df_cptac$patient_id <- rownames(gsva_df_cptac)

cptac_gsva_merged <- merge(cptac_merged_cp, gsva_df_cptac, by = "patient_id")
stopifnot(nrow(cptac_gsva_merged) == 95)

pathway_results_tcga <- readRDS("results1/pathway_results.rds")
tcga_sig_pathways <- pathway_results_tcga$Pathway[pathway_results_tcga$adj.P.Val < 0.05]
cat("[CHECK]", length(tcga_sig_pathways), "TCGA Hallmark pathways significant (FDR<0.05) -- testing these for replication in CPTAC\n")
missing_pathways_cptac <- setdiff(tcga_sig_pathways, colnames(cptac_gsva_merged))
if (length(missing_pathways_cptac) > 0) {
  cat("[CHECK] TCGA-significant pathways NOT found in CPTAC's GSVA output:",
      paste(missing_pathways_cptac, collapse = ", "), "\n")
}
stopifnot(all(tcga_sig_pathways %in% colnames(cptac_gsva_merged)))

run_pathway_test_cptac <- function(data, var) {
  formula_unadj <- as.formula(paste0("`", var, "` ~ risk_group"))
  wt <- wilcox.test(formula_unadj, data = data)
  formula_adj <- as.formula(paste0("`", var, "` ~ risk_group + GENOMICS_SUBTYPE"))
  lm_adj <- lm(formula_adj, data = data)
  adj_coef <- summary(lm_adj)$coefficients["risk_groupHigh", ]
  data.frame(
    pathway      = var,
    unadj_p      = wt$p.value,
    logFC        = median(data[[var]][data$risk_group == "High"], na.rm = TRUE) -
      median(data[[var]][data$risk_group == "Low"],  na.rm = TRUE),
    adj_estimate = adj_coef["Estimate"],
    adj_p        = adj_coef["Pr(>|t|)"]
  )
}

cptac_pathway_results <- do.call(rbind, lapply(tcga_sig_pathways, run_pathway_test_cptac, data = cptac_gsva_merged))
cptac_pathway_results$unadj_p_fdr <- p.adjust(cptac_pathway_results$unadj_p, method = "fdr")
cptac_pathway_results$adj_p_fdr   <- p.adjust(cptac_pathway_results$adj_p, method = "fdr")
cptac_pathway_results <- cptac_pathway_results %>% arrange(adj_p_fdr)

cat("\n=== CPTAC: GSVA pathway replication (of TCGA's FDR-significant pathways) ===\n")
print(cptac_pathway_results, digits = 4)

n_same_direction <- sum(sign(cptac_pathway_results$logFC) ==
                          sign(pathway_results_tcga$logFC[match(cptac_pathway_results$pathway, pathway_results_tcga$Pathway)]))
n_replicated <- sum(cptac_pathway_results$adj_p_fdr < 0.05)

cat("\n[CHECK] Pathways with same direction as TCGA:", n_same_direction, "/", length(tcga_sig_pathways), "\n")
cat("[CHECK] Pathways surviving FDR+subtype adjustment in CPTAC:", n_replicated, "/", length(tcga_sig_pathways), "\n")

write.csv(cptac_pathway_results, "cptac_validation/cptac_gsva_pathway_results.csv", row.names = FALSE)
saveRDS(cptac_gsva_merged, "cptac_validation/cptac_gsva_merged.rds")

# Median-center CPTAC expression (TIDE's documented requirement)
cptac_gene_medians <- apply(cptac_expr_log, 1, median, na.rm = TRUE)
cptac_expr_tide_centered <- sweep(cptac_expr_log, 1, cptac_gene_medians, FUN = "-")

med_check_cptac <- summary(apply(cptac_expr_tide_centered, 1, median))
stopifnot(abs(med_check_cptac["Median"]) < 1e-8)
cat("\n[CHECK] CPTAC gene-wise median after centering (expect ~0):", med_check_cptac["Median"], "\n")

cptac_tide_export <- data.frame(GeneSymbol = rownames(cptac_expr_tide_centered),
                                cptac_expr_tide_centered, check.names = FALSE)
write.table(cptac_tide_export, file = "cptac_validation/tide_input_cptac_n95.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

Sys.setenv(CONDA_SOLVER = "classic")
reticulate::use_condaenv("tide_env", required = TRUE)

py_run_string("
import pandas as pd
from tidepy.pred import TIDE
 
df = pd.read_csv('cptac_validation/tide_input_cptac_n95.txt', sep='\t', index_col=0)
print('CPTAC TIDE input shape:', df.shape)
 
result = TIDE(df, cancer='Other', pretreat=False, ignore_norm=True, vthres=0.)
result.to_csv('cptac_validation/tide_scores_cptac_n95.csv')
")

cptac_tide_scores <- read.csv("cptac_validation/tide_scores_cptac_n95.csv", row.names = 1)
cptac_tide_scores$patient_id <- rownames(cptac_tide_scores)
stopifnot(nrow(cptac_tide_scores) == 95)
cat("[CHECK] CPTAC TIDE scores:", nrow(cptac_tide_scores), "patients\n")

cptac_tide_merged <- merge(cptac_merged_cp, cptac_tide_scores, by = "patient_id")
stopifnot(nrow(cptac_tide_merged) == 95)

cat("\n=== CPTAC: TIDE unadjusted ===\n")
print(wilcox.test(TIDE ~ risk_group, data = cptac_tide_merged))
cat("\n=== CPTAC: TIDE subtype-adjusted ===\n")
tide_adj_coef_cptac <- summary(lm(TIDE ~ risk_group + GENOMICS_SUBTYPE, data = cptac_tide_merged))$coefficients["risk_groupHigh", ]
print(tide_adj_coef_cptac)

# Responder rate comparison
cptac_responder_table <- table(cptac_tide_merged$risk_group, cptac_tide_merged$Responder)
cat("\n[CHECK] CPTAC responder rate by risk group:\n")
print(prop.table(cptac_responder_table, margin = 1))
print(chisq.test(cptac_responder_table))

dir.create("cptac_validation/ips_calc", showWarnings = FALSE)
cptac_ips_export <- data.frame(GeneSymbol = rownames(cptac_expr_mat_full),
                               cptac_expr_mat_full, check.names = FALSE)
write.table(cptac_ips_export, file = "cptac_validation/ips_calc/EXPR.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

stopifnot(file.exists("ips_calc/IPS_genes.txt"))
file.copy("ips_calc/IPS_genes.txt", "cptac_validation/ips_calc/IPS_genes.txt", overwrite = TRUE)

ipsmap <- function(x) {
  if (x <= 0) return(0)
  if (x >= 3) return(10)
  return(round(x * 10 / 3, digits = 0))
}

gene_expression_cptac <- read.table("cptac_validation/ips_calc/EXPR.txt", row.names = 1,
                                    header = TRUE, sep = "\t", check.names = FALSE)
IPSG <- read.table("cptac_validation/ips_calc/IPS_genes.txt", header = TRUE, sep = "\t", check.names = FALSE)
unique_ips_genes <- as.vector(unique(IPSG$NAME))
sample_names_cptac <- names(gene_expression_cptac)

missing_genes_cptac <- setdiff(as.vector(IPSG$GENE), rownames(gene_expression_cptac))
cat("\n[CHECK] IPS genes missing from CPTAC expression matrix:",
    if(length(missing_genes_cptac)==0) "NONE" else paste(missing_genes_cptac, collapse=", "), "\n")

IPS <- MHC <- CP <- EC <- SC <- AZ <- numeric(length(sample_names_cptac))
for (i in seq_along(sample_names_cptac)) {
  GE  <- gene_expression_cptac[[i]]
  Z1  <- (gene_expression_cptac[as.vector(IPSG$GENE), i] - mean(GE)) / sd(GE)
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

cptac_ips_df <- data.frame(patient_id = sample_names_cptac, MHC = MHC, EC = EC, SC = SC, CP = CP, AZ = AZ, IPS = IPS)
stopifnot(nrow(cptac_ips_df) == 95, all(cptac_ips_df$IPS >= 0 & cptac_ips_df$IPS <= 10))
cat("[CHECK] CPTAC IPS scores: n=", nrow(cptac_ips_df), " range=[", min(cptac_ips_df$IPS), ",", max(cptac_ips_df$IPS), "]\n")

cptac_ips_merged <- merge(cptac_merged_cp, cptac_ips_df, by = "patient_id")
cat("\n=== CPTAC: IPS unadjusted ===\n")
print(wilcox.test(IPS ~ risk_group, data = cptac_ips_merged))
cat("\n=== CPTAC: IPS subtype-adjusted ===\n")
print(summary(lm(IPS ~ risk_group + GENOMICS_SUBTYPE, data = cptac_ips_merged))$coefficients["risk_groupHigh", ])

# SAVE ALL RESULTS
cptac_full_extension <- cptac_tide_merged %>%
  merge(cptac_ips_df, by = "patient_id") %>%
  merge(gsva_df_cptac, by = "patient_id")

saveRDS(cptac_full_extension, "cptac_validation/cptac_full_extension.rds")
cptac_full_extension <- readRDS("cptac_validation/cptac_full_extension.rds")   # checkpoint reload
write.csv(cptac_pathway_results, "cptac_validation/cptac_gsva_pathway_results.csv", row.names = FALSE)

cat("\n[SAVED] cptac_validation/cptac_full_extension.rds\n")
cat("[SAVED] cptac_validation/cptac_gsva_pathway_results.csv\n")

library(ggplot2)
library(dplyr)

dir.create("plots1/cptac_figures", showWarnings = FALSE, recursive = TRUE)

palette_risk <- c("Low" = "#6BBFB0", "High" = "salmon")

# PLOT 1: CPTAC TIDE score by risk group
cptac_tide_merged <- readRDS("cptac_validation/cptac_full_extension.rds")

tide_adj_p_label <- sprintf("Subtype-adjusted p = %.4g", tide_adj_coef_cptac["Pr(>|t|)"])

p_cptac_tide <- ggplot(cptac_tide_merged, aes(x = risk_group, y = TIDE, fill = risk_group)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.size = 0.6) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "TIDE Score",
       title = "CPTAC-UCEC: TIDE Score by Risk Group") +
  annotate("text", x = 1.5, y = max(cptac_tide_merged$TIDE, na.rm = TRUE) * 1.05,
           label = tide_adj_p_label, size = 3, fontface = "italic") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_cptac_tide

ggsave("plots1/cptac_figures/CPTAC_TIDE.png", p_cptac_tide,
       width = 5.5, height = 4.5, dpi = 300, bg = "white")

# PLOT 2: CPTAC immunotherapy responder rate by risk group

# Responder column, same "True"/"False" string coding as TCGA
cptac_responder_summary <- cptac_tide_merged %>%
  group_by(risk_group, Responder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(risk_group) %>%
  mutate(pct = 100 * n / sum(n))

p_cptac_responder <- ggplot(cptac_responder_summary %>% filter(Responder == "True"),
                            aes(x = risk_group, y = pct, fill = risk_group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = palette_risk) +
  labs(title = "CPTAC-UCEC: Predicted Immunotherapy Responder Rate",
       x = "", y = "Responder Rate (%)") +
  ylim(0, max(cptac_responder_summary$pct[cptac_responder_summary$Responder == "True"]) * 1.2) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        legend.position = "none",
        axis.text = element_text(color = "black"))

ggsave("plots1/cptac_figures/CPTAC_Responder_Rate.png", p_cptac_responder,
       width = 5, height = 5, dpi = 300, bg = "white")

# PLOT 3: GSVA pathway replication (TCGA vs. CPTAC, effect sizes side by side)
pathway_results_tcga <- readRDS("results1/pathway_results.rds") %>%
  filter(adj.P.Val < 0.05) %>%
  rename(pathway = Pathway) %>%
  mutate(cohort = "TCGA-UCEC (n=541)", adj_p_fdr = adj.P.Val)

cptac_pathway_results <- read.csv("cptac_validation/cptac_gsva_pathway_results.csv") %>%
  mutate(cohort = "CPTAC-UCEC (n=95)")

# Match columns needed for plotting
tcga_plot_data <- pathway_results_tcga %>% dplyr::select(pathway, logFC, adj_p_fdr, cohort)
cptac_plot_data <- cptac_pathway_results %>% dplyr::select(pathway, logFC, adj_p_fdr, cohort)

combined_pathways <- bind_rows(tcga_plot_data, cptac_plot_data) %>%
  mutate(
    sig = ifelse(adj_p_fdr < 0.05, "FDR < 0.05", "n.s."),
    cohort = factor(cohort, levels = c("TCGA-UCEC (n=541)", "CPTAC-UCEC (n=95)")),
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean)
  )

# Order by TCGA logFC for consistent, readable ordering
gene_order <- tcga_plot_data %>% arrange(logFC) %>% pull(pathway)
combined_pathways$pathway <- factor(combined_pathways$pathway, levels = gene_order)
combined_pathways$pathway_clean <- factor(combined_pathways$pathway_clean,
                                          levels = gsub("_", " ", gsub("HALLMARK_", "", gene_order)))

p_gsva_replication <- ggplot(combined_pathways, aes(x = logFC, y = pathway_clean, color = cohort, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("TCGA-UCEC (n=541)" = "#3B6FB6", "CPTAC-UCEC (n=95)" = "#D73027"), name = "Cohort") +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "n.s." = 1), name = NULL) +
  labs(x = "Log Fold Change (High- vs. Low-Risk)", y = NULL,
       title = "GSVA Hallmark Pathway Replication: TCGA vs. CPTAC-UCEC") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top",
        axis.text.y = element_text(size = 7))

ggsave("plots1/cptac_figures/CPTAC_GSVA_Replication.png", p_gsva_replication,
       width = 9, height = 8, dpi = 300, bg = "white")

cat("[SAVED] 3 CPTAC plots to plots1/cptac_figures/\n")

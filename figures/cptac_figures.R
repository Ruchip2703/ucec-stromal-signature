# cptac_figures_WL.R
#
# CPTAC-UCEC replication figure: checkpoint gene effects (TCGA vs CPTAC) and TME scores by risk group.
# Input: results1/checkpoint_gene_results.csv, cptac_validation/  |  Output: plots1/cptac_figures/

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

dir.create("plots1/cptac_figures", showWarnings = FALSE, recursive = TRUE)

palette_risk <- c("Low" = "#6BBFB0", "High" = "salmon")

# PANEL A: TCGA VS. CPTAC CHECKPOINT GENE EFFECT SIZES, SIDE BY SIDE
tcga_checkpoint <- read.csv("results1/checkpoint_gene_results.csv") %>%
  rename(gene = variable) %>%
  mutate(cohort = "TCGA-UCEC (n=541)")

cptac_checkpoint <- read.csv("cptac_validation/cptac_checkpoint_validation_results.csv") %>%
  rename(gene = variable) %>%
  mutate(cohort = "CPTAC-UCEC (n=95)")

stopifnot(nrow(tcga_checkpoint) == 10, nrow(cptac_checkpoint) == 10)

combined_checkpoint <- bind_rows(tcga_checkpoint, cptac_checkpoint) %>%
  mutate(
    sig = ifelse(adj_p_fdr < 0.05, "FDR < 0.05", "n.s."),
    cohort = factor(cohort, levels = c("TCGA-UCEC (n=541)", "CPTAC-UCEC (n=95)"))
  )

# Order genes by TCGA effect size for consistent, readable ordering
gene_order <- tcga_checkpoint %>% arrange(adj_estimate) %>% pull(gene)
combined_checkpoint$gene <- factor(combined_checkpoint$gene, levels = gene_order)

common_names <- c(
  PDCD1 = "PDCD1 (PD-1)", CTLA4 = "CTLA4", LAG3 = "LAG3", TIGIT = "TIGIT",
  HAVCR2 = "HAVCR2 (TIM-3)", CD274 = "CD274 (PD-L1)", PDCD1LG2 = "PDCD1LG2 (PD-L2)",
  CD27 = "CD27", ICOS = "ICOS", IDO1 = "IDO1"
)
combined_checkpoint$gene_label <- common_names[as.character(combined_checkpoint$gene)]
combined_checkpoint$gene_label <- factor(combined_checkpoint$gene_label, levels = common_names[gene_order])

panel_A <- ggplot(combined_checkpoint, aes(x = adj_estimate, y = gene_label, color = cohort, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("TCGA-UCEC (n=541)" = "#2166AC", "CPTAC-UCEC (n=95)" = "#B2182B"), name = "Cohort") +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "n.s." = 1), name = NULL) +
  labs(x = "Subtype-Adjusted Effect Size (High- vs. Low-Risk)", y = NULL,
       title = "Immune Checkpoint Gene Replication: TCGA vs. CPTAC-UCEC")+
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

panel_A

# PANEL B: CPTAC TME SCORES BY RISK GROUP
cptac_merged <- readRDS("cptac_validation/cptac_final_validation_data.rds")
stopifnot(all(c("EPIC_CAFS", "CIBERSORT_MACROPHAGES _M2", "ESTIMATE_STROMALSCORE",
                "ESTIMATE_IMMUNESCORE", "risk_group") %in% colnames(cptac_merged)))

tme_long <- cptac_merged %>%
  dplyr::select(risk_group, EPIC_CAFS, `CIBERSORT_MACROPHAGES _M2`,
                ESTIMATE_STROMALSCORE, ESTIMATE_IMMUNESCORE) %>%
  rename(`CAF (EPIC)` = EPIC_CAFS,
         `M2 Macrophages (CIBERSORT)` = `CIBERSORT_MACROPHAGES _M2`,
         `Stromal Score (ESTIMATE)` = ESTIMATE_STROMALSCORE,
         `Immune Score (ESTIMATE)` = ESTIMATE_IMMUNESCORE) %>%
  pivot_longer(cols = -risk_group, names_to = "score_type", values_to = "value")

panel_B <- ggplot(tme_long, aes(x = risk_group, y = value, fill = risk_group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.size = 0.5) +
  facet_wrap(~ score_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "Score",
       title = "Independently-Derived TME Scores by Risk Group (CPTAC-UCEC)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), strip.text = element_text(size = 8.5))

panel_B

# SAVE INDIVIDUAL PANELS (for slide decks / standalone use)
ggsave("plots1/cptac_figures/cptac_panelA_checkpoints.png", panel_A, width = 10, height = 5.5, dpi = 300, units = "in")
ggsave("plots1/cptac_figures/cptac_panelB_TME.png",         panel_B, width = 10, height = 3.5, dpi = 300, units = "in")
cat("[SAVED] 2 individual CPTAC panels to plots1/cptac_figures/\n")

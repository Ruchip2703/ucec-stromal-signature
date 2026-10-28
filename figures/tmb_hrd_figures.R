# tmb_hrd_figures_WL.R
#
# Molecular subtype, TMB and HRD figure panels.
# Input: results1/risk_subtype_tmb_hrd.rds  |  Output: plots1/tmb_hrd_figures_clean/

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

dir.create("plots1/tmb_hrd_figures_clean", showWarnings = FALSE, recursive = TRUE)

risk_full <- readRDS("results1/risk_subtype_tmb_hrd.rds")
stopifnot(all(c("risk_group", "Subtype_Selected", "TMB", "HRD_score") %in% colnames(risk_full)))

palette_risk <- c("Low" = "#6BBFB0", "High" = "salmon")

subtype_labels <- c(
  "UCEC.CN_HIGH" = "CN-high", "UCEC.CN_LOW" = "CN-low", "UCEC.MSI" = "MSI",
  "UCEC.POLE" = "POLE", "UCEC.NA" = "No subtype\ncall"
)
subtype_order <- c("UCEC.POLE", "UCEC.MSI", "UCEC.CN_LOW", "UCEC.CN_HIGH", "UCEC.NA")
risk_full$Subtype_Selected <- factor(risk_full$Subtype_Selected, levels = subtype_order)

# PANEL A: SUBTYPE COMPOSITION BY RISK GROUP
subtype_data <- risk_full %>%
  filter(!is.na(Subtype_Selected)) %>%
  count(risk_group, Subtype_Selected) %>%
  group_by(risk_group) %>%
  mutate(pct = n / sum(n) * 100)

panel_A <- ggplot(subtype_data, aes(x = risk_group, y = pct, fill = Subtype_Selected)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2", name = "Molecular\nSubtype", labels = subtype_labels) +
  labs(x = NULL, y = "Percent of Patients", title = "Molecular Subtype by Risk Group") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_blank())

panel_A

# PANEL B: TMB BY RISK GROUP, FACETED BY SUBTYPE
panel_B <- risk_full %>%
  filter(!is.na(TMB), !is.na(Subtype_Selected)) %>%
  ggplot(aes(x = risk_group, y = TMB, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  facet_wrap(~ Subtype_Selected, nrow = 1, labeller = as_labeller(subtype_labels)) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  scale_y_log10(labels = label_number()) +
  labs(x = NULL, y = "TMB (mut/Mb, log scale)",
       title = "TMB by Risk Group, Within Molecular Subtype") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

panel_B

# PANEL C: HRD BY RISK GROUP, FACETED BY SUBTYPE
panel_C <- risk_full %>%
  filter(!is.na(HRD_score), !is.na(Subtype_Selected)) %>%
  ggplot(aes(x = risk_group, y = HRD_score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.85) +
  facet_wrap(~ Subtype_Selected, nrow = 1, labeller = as_labeller(subtype_labels)) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "HRD Score",
       title = "HRD Score by Risk Group, Within Molecular Subtype") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

panel_C

# PANEL D: UNADJUSTED VS. SUBTYPE-ADJUSTED SIGNIFICANCE SUMMARY
# (TMB before HRD, matching Results narrative order)
p_vals <- c(0.0070, 0.972, 3.83e-16, 0.105)

summary_stats <- data.frame(
  marker = factor(rep(c("TMB", "HRD"), each = 2), levels = c("TMB", "HRD")),
  model  = rep(c("Unadjusted", "Subtype-adjusted"), 2),
  neg_log10_p = -log10(p_vals)
)
summary_stats$model <- factor(summary_stats$model, levels = c("Unadjusted", "Subtype-adjusted"))
summary_stats$p_label <- sprintf("p = %.3g", p_vals)

panel_D <- ggplot(summary_stats, aes(x = marker, y = neg_log10_p, fill = model)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = p_label), position = position_dodge(0.7),
            vjust = -0.4, size = 2.9, color = "gray20") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  annotate("text", x = 2.53, y = -log10(0.05) + 0.3, label = "p = 0.05",
           size = 3, color = "gray40", hjust = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_fill_manual(values = c("Unadjusted" = "#4C72B0", "Subtype-adjusted" = "#C44E52"), name = NULL) +
  labs(x = NULL, y = expression(-log[10](p)),
       title = "Significance Before and After Subtype Adjustment") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")

panel_D

# SAVE INDIVIDUAL PANELS (for slide decks / standalone use)
ggsave("plots1/tmb_hrd_figures_clean/A3_panelA_subtype.png",      panel_A, width = 6,  height = 5,   dpi = 300, units = "in")
ggsave("plots1/tmb_hrd_figures_clean/A3_panelB_TMB.png",           panel_B, width = 10, height = 3.5, dpi = 300, units = "in")
ggsave("plots1/tmb_hrd_figures_clean/A3_panelC_HRD.png",           panel_C, width = 10, height = 3.5, dpi = 300, units = "in")
ggsave("plots1/tmb_hrd_figures_clean/A3_panelD_significance.png",  panel_D, width = 7,  height = 7,   dpi = 300, units = "in")
cat("[SAVED] 4 individual A3 panels to plots1/tmb_hrd_figures_clean/\n")

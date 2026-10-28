# tide_figures_WL.R
#
# TIDE, dysfunction/exclusion, cell-type, mechanism, checkpoint, responder-rate and IPS panels.
# Input: results1/survival_with_checkpoints.rds, results1/checkpoint_gene_results.csv, results1/survival_with_tide_ips.rds  |  Output: plots1/tide_figures_clean/

library(ggplot2)
library(dplyr)
library(tidyr)

survival_cp <- readRDS("results1/survival_with_checkpoints.rds")
checkpoint_results <- read.csv("results1/checkpoint_gene_results.csv") %>% rename(gene = variable)

col_low  <- "#6BBFB0"   # blue  -- Low risk
col_high <- "salmon"   # red   -- High risk
palette_risk <- c("Low" = col_low, "High" = col_high)

dir.create("plots1/tide_figures_clean", showWarnings = FALSE, recursive = TRUE)

# TIDE Score
cat("[CHECK] Subtype_Selected present in survival_cp:", "Subtype_Selected" %in% colnames(survival_cp), "\n")
tide_adj_coef <- summary(lm(TIDE ~ risk_group + Subtype_Selected, data = survival_cp))$coefficients["risk_groupHigh", ]
cat("[CHECK] Subtype-adjusted TIDE p-value (used in panel annotation):", tide_adj_coef["Pr(>|t|)"], "\n")
tide_p_label <- sprintf("Subtype-adjusted p = %.3g", tide_adj_coef["Pr(>|t|)"])

panel_TIDE <- ggplot(survival_cp, aes(x = risk_group, y = TIDE, fill = risk_group)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.size = 0.6) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "TIDE Score", title = "TIDE Score by Risk Group") +
  annotate("text", x = 1.5, y = max(survival_cp$TIDE, na.rm = TRUE) * 1.05,
           label = tide_p_label, size = 3, fontface = "italic") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

panel_TIDE
ggsave("plots1/tide_figures_clean/TIDE.png", panel_TIDE, width = 5.5, height = 4.5, dpi = 300, bg = "white")

# Dysfunction vs. Exclusion
db_long <- survival_cp %>%
  dplyr::select(risk_group, Dysfunction, Exclusion) %>%
  pivot_longer(cols = c(Dysfunction, Exclusion), names_to = "component", values_to = "score")

panel_DysExc <- ggplot(db_long, aes(x = component, y = score, fill = risk_group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.size = 0.5, position = position_dodge(0.7)) +
  scale_fill_manual(values = palette_risk, name = "Risk Group") +
  labs(x = NULL, y = "Score", title = "Dysfunction vs. Exclusion by Risk Group") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

panel_DysExc

ggsave("plots1/tide_figures_clean/DysfunctionExclusion.png", panel_DysExc, width = 5.5, height = 4.5, dpi = 300, bg = "white")

# Exclusion-Driving Cell Types
cell_long <- survival_cp %>%
  dplyr::select(risk_group, CAF, TAM.M2, MDSC) %>%
  pivot_longer(cols = c(CAF, TAM.M2, MDSC), names_to = "cell_type", values_to = "score") %>%
  mutate(cell_type = factor(cell_type, levels = c("CAF", "TAM.M2", "MDSC"),
                            labels = c("CAF", "M2 Macrophages", "MDSC")))

panel_CellTypes <- ggplot(cell_long, aes(x = cell_type, y = score, fill = risk_group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.size = 0.5, position = position_dodge(0.7)) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "TIDE Cell-Type Score",
       title = "Exclusion-Driving Cell Types by Risk Group") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

panel_CellTypes

ggsave("plots1/tide_figures_clean/CellTypes.png", panel_CellTypes, width = 5.5, height = 4.5, dpi = 300, bg = "white")

# Primary Mechanism
mech_data <- survival_cp %>%
  filter(!is.na(CTL.flag)) %>%
  mutate(Mechanism = ifelse(CTL.flag, "Dysfunction-primary", "Exclusion-primary")) %>%
  count(risk_group, Mechanism) %>%
  group_by(risk_group) %>%
  mutate(pct = n / sum(n) * 100)

panel_Mechanism <- ggplot(mech_data, aes(x = risk_group, y = pct, fill = Mechanism)) +
  geom_col(width = 0.6, color = "white") +
  scale_fill_manual(values = c("Exclusion-primary" = col_high, "Dysfunction-primary" = col_low), name = NULL) +
  labs(x = NULL, y = "Percent of Patients",
       title = "Primary Immune-Evasion Mechanism by Risk Group") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

panel_Mechanism

ggsave("plots1/tide_figures_clean/Mechanism.png", panel_Mechanism, width = 5.5, height = 4.5, dpi = 300, bg = "white")

# Checkpoint Gene Forest Plot
common_names <- c(
  PDCD1 = "PDCD1 (PD-1)", CTLA4 = "CTLA4", LAG3 = "LAG3", TIGIT = "TIGIT",
  HAVCR2 = "HAVCR2 (TIM-3)", CD274 = "CD274 (PD-L1)", PDCD1LG2 = "PDCD1LG2 (PD-L2)",
  CD27 = "CD27", ICOS = "ICOS", IDO1 = "IDO1"
)

checkpoint_results <- checkpoint_results %>%
  mutate(sig_label = ifelse(adj_p_fdr < 0.001, "FDR < 0.001",
                            ifelse(adj_p_fdr < 0.05, "FDR < 0.05", "n.s.")))
checkpoint_results$gene_label <- common_names[as.character(checkpoint_results$gene)]
checkpoint_results$gene_label <- reorder(checkpoint_results$gene_label, checkpoint_results$adj_estimate)

# NOTE: significance-tier colors here (FDR<0.001 / FDR<0.05 / n.s.) are a
panel_Checkpoints <- ggplot(checkpoint_results, aes(x = adj_estimate, y = gene_label, color = sig_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(-1.5, 0, by = 0.25)) +
  scale_color_manual(values = c("FDR < 0.001" = "orange", "FDR < 0.05" = "#C44E52", "n.s." = "gray60"), name = NULL) +
  labs(x = "Subtype-Adjusted Effect (risk_groupHigh)", y = NULL,
       title = "Immune Checkpoint Gene Expression: High- vs. Low-Risk") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

panel_Checkpoints

ggsave("plots1/tide_figures_clean/Checkpoints.png", panel_Checkpoints, width = 9, height = 5, dpi = 300, bg = "white")

cat("[SAVED] 5 clean individual panels to plots1/A2_figures_clean/\n")

library(ggplot2)
library(dplyr)

survival_full <- readRDS("results1/survival_with_tide_ips.rds")

# Check what values these columns actually contain
str(survival_full$Responder)
table(survival_full$Responder, useNA = "always")

str(survival_full$risk_group)
table(survival_full$risk_group, useNA = "always")

# PLOT 1: Immunotherapy responder rate by risk group (bar chart)
responder_summary <- survival_full %>%
  group_by(risk_group, Responder) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(risk_group) %>%
  mutate(pct = 100 * n / sum(n))

p_responder <- ggplot(responder_summary %>% filter(Responder == "True"),
                      aes(x = risk_group, y = pct, fill = risk_group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = palette_risk) +
  labs(title = "Predicted Immunotherapy Responder Rate",
       x = "", y = "Responder Rate (%)") +
  ylim(0, max(responder_summary$pct[responder_summary$Responder == "True"]) * 1.2) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text = element_text(color = "black"))

p_responder

ggsave("plots1/tide_figures_clean/Responder_Rate_Barplot.png", p_responder, width = 5, height = 5, dpi = 300, bg = "white")

# PLOT 2: IPS boxplot by risk group
cat("[CHECK] Subtype_Selected present in survival_full:", "Subtype_Selected" %in% colnames(survival_full), "\n")
ips_adj_coef <- summary(lm(IPS ~ risk_group + Subtype_Selected, data = survival_full))$coefficients["risk_groupHigh", ]
cat("[CHECK] Subtype-adjusted IPS p-value (used in panel annotation):", ips_adj_coef["Pr(>|t|)"], "\n")
ips_p_label <- sprintf("Subtype-adjusted p = %.3g", ips_adj_coef["Pr(>|t|)"])

p_ips <- ggplot(survival_full, aes(x = risk_group, y = IPS, fill = risk_group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = palette_risk) +
  labs(title = "Immunophenoscore (IPS) by Risk Group",
       x = "", y = "IPS") +
  annotate("text", x = 1.5, y = max(survival_full$IPS, na.rm = TRUE) * 1.05,
           label = ips_p_label, size = 3, fontface = "italic") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text = element_text(color = "black"))

p_ips

ggsave("plots1/tide_figures_clean/IPS_Boxplot.png", p_ips, width = 5, height = 5, dpi = 300, bg = "white")

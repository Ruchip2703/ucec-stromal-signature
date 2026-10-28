# drug_sensitivity_figures_WL.R
#
# Drug sensitivity figure: percent change in predicted IC50 for 13 drugs and representative distributions.
# Input: results1/survival_with_drug_sensitivity.rds, results1/drug_sensitivity_results.csv  |  Output: plots1/drug_figures_clean/

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

dir.create("plots1/drug_figures_clean", showWarnings = FALSE, recursive = TRUE)

survival_drug <- readRDS("results1/survival_with_drug_sensitivity.rds")
drug_results  <- read.csv("results1/drug_sensitivity_results.csv")

stopifnot(nrow(drug_results) == 13)
stopifnot(all(drug_results$adj_estimate > 0))
cat("[CHECK] All 13 drugs show adj_estimate > 0 (High-risk = more resistant): TRUE\n")

palette_risk <- c("Low" = "#6BBFB0", "High" = "salmon")

# PANEL A: FOREST PLOT, ALL 13 DRUGS, PERCENT CHANGE (SCALE-INVARIANT)
drug_class <- c(
  Paclitaxel_1080 = "Taxane", Docetaxel_1007 = "Taxane", Docetaxel_1819 = "Taxane",
  Cisplatin_1005 = "Platinum", Oxaliplatin_1089 = "Platinum", Oxaliplatin_1806 = "Platinum",
  Epirubicin_1511 = "Anthracycline",
  Olaparib_1017 = "PARP inhibitor", Talazoparib_1259 = "PARP inhibitor", Niraparib_1177 = "PARP inhibitor",
  Tamoxifen_1199 = "Hormonal", Fulvestrant_1200 = "Hormonal", Fulvestrant_1816 = "Hormonal"
)
drug_label <- c(
  Paclitaxel_1080 = "Paclitaxel", Docetaxel_1007 = "Docetaxel (screen 1)",
  Docetaxel_1819 = "Docetaxel (screen 2)", Cisplatin_1005 = "Cisplatin",
  Oxaliplatin_1089 = "Oxaliplatin (screen 1)", Oxaliplatin_1806 = "Oxaliplatin (screen 2)",
  Epirubicin_1511 = "Epirubicin", Olaparib_1017 = "Olaparib", Talazoparib_1259 = "Talazoparib",
  Niraparib_1177 = "Niraparib", Tamoxifen_1199 = "Tamoxifen",
  Fulvestrant_1200 = "Fulvestrant (screen 1)", Fulvestrant_1816 = "Fulvestrant (screen 2)"
)

stopifnot(all(drug_results$drug %in% names(drug_class)))

drug_results <- drug_results %>%
  mutate(
    pct_change = (median_High - median_Low) / median_Low * 100,
    class = drug_class[drug], label = drug_label[drug],
    label = reorder(label, pct_change)
  )

panel_A <- ggplot(drug_results, aes(x = pct_change, y = label, color = class)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.8) +
  scale_color_brewer(palette = "Set1", name = "Drug Class") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Percent Change in Median Predicted IC50 (High- vs. Low-Risk)", y = NULL,
       title = "Predicted Drug Sensitivity: High- vs. Low-Risk (all higher in High-risk)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

panel_A

# PANEL B: RAW DISTRIBUTIONS, ONE REPRESENTATIVE DRUG PER CLASS (5 CLASSES)
representative_drugs <- c("Paclitaxel_1080", "Cisplatin_1005", "Epirubicin_1511",
                          "Olaparib_1017", "Tamoxifen_1199")
rep_labels <- c(Paclitaxel_1080 = "Paclitaxel\n(Taxane)", Cisplatin_1005 = "Cisplatin\n(Platinum)",
                Epirubicin_1511 = "Epirubicin\n(Anthracycline)", Olaparib_1017 = "Olaparib\n(PARP inhibitor)",
                Tamoxifen_1199 = "Tamoxifen\n(Hormonal)")

stopifnot(all(representative_drugs %in% colnames(survival_drug)))

rep_data <- survival_drug %>%
  dplyr::select(risk_group, all_of(representative_drugs)) %>%
  pivot_longer(cols = -risk_group, names_to = "drug", values_to = "IC50") %>%
  mutate(drug_label = rep_labels[drug], drug_label = factor(drug_label, levels = rep_labels[representative_drugs]))

panel_B <- ggplot(rep_data, aes(x = risk_group, y = IC50, fill = risk_group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.size = 0.5) +
  facet_wrap(~ drug_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = palette_risk, guide = "none") +
  labs(x = NULL, y = "Predicted IC50", title = "Representative Drugs, One Per Class") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), strip.text = element_text(size = 9))

panel_B

# SAVE INDIVIDUAL PANELS (for slide decks / standalone use)
ggsave("plots1/drug_figures_clean/drug_panelA_forest.png",   panel_A, width = 10, height = 6, dpi = 300, units = "in")
ggsave("plots1/drug_figures_clean/drug_panelB_boxplots.png", panel_B, width = 10, height = 4, dpi = 300, units = "in")
cat("[SAVED] 2 individual drug panels to plots1/drug_figures_clean/\n")

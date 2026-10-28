# 12_crispr_check_WL.R
#
# DepMap CRISPR (Chronos) essentiality of the 9 signature genes, pan-cancer and endometrial/uterine lines.
# Output: results1/depmap_*.csv, plots1/depmap_figures/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("depmap", quietly = TRUE)) BiocManager::install("depmap")
if (!requireNamespace("ExperimentHub", quietly = TRUE)) BiocManager::install("ExperimentHub")

library(depmap)
library(ExperimentHub)
library(dplyr)
library(ggplot2)

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)
dir.create("plots1/depmap_figures", showWarnings = FALSE, recursive = TRUE)

col_high <- "#eb6834"

# LOAD CRISPR GENE EFFECT (CHRONOS) DATA + CELL LINE METADATA
# First run downloads via ExperimentHub's automatic cache -- may take a few
# minutes depending on connection, but is a standard, reliable Bioconductor
# download (not the ad-hoc OSF situation from GDSC2).
crispr_data <- depmap_crispr()      # gene, cell_line, dependency (Chronos score), etc.
metadata    <- depmap_metadata()    # cell line lineage/disease annotations

cat("[CHECK] CRISPR dataset dimensions:", paste(dim(crispr_data), collapse = " x "), "\n")
cat("[CHECK] Metadata dimensions:", paste(dim(metadata), collapse = " x "), "\n")

# FILTER TO THE 9 SIGNATURE GENES
sig_genes <- c("TIMP3", "PTGIS", "PPP1R16B", "DIO2", "CTSW", "ZNF683", "GZMK", "IHH", "SPN")

gene_crispr <- crispr_data %>%
  filter(gene_name %in% sig_genes) %>%
  left_join(metadata %>% dplyr::select(depmap_id, lineage, lineage_subtype, primary_disease),
            by = "depmap_id")

cat("[CHECK] Signature genes found in CRISPR data:", paste(unique(gene_crispr$gene_name), collapse = ", "), "\n")
missing_genes <- setdiff(sig_genes, unique(gene_crispr$gene_name))
if (length(missing_genes) > 0) cat("[CHECK] NOT found in DepMap:", paste(missing_genes, collapse = ", "), "\n")

# PAN-CANCER SUMMARY (context -- how essential is each gene, overall?)
pan_cancer_summary <- gene_crispr %>%
  filter(!is.na(dependency)) %>%
  group_by(gene_name) %>%
  summarise(
    n_cell_lines = n(),
    median_chronos = median(dependency, na.rm = TRUE),
    pct_essential_minus0.5 = round(100 * mean(dependency <= -0.5, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(median_chronos)

cat("\n=== Pan-cancer CRISPR essentiality summary (all cell lines) ===\n")
print(pan_cancer_summary, digits = 3)

# ENDOMETRIAL/UTERINE-SPECIFIC SUBSET
cat("\n[CHECK] Lineage values containing 'uter' or 'endometr':\n")
print(unique(metadata$lineage[grepl("uter|endometr", metadata$lineage, ignore.case = TRUE)]))
print(unique(metadata$lineage_subtype[grepl("uter|endometr", metadata$lineage_subtype, ignore.case = TRUE)]))

ec_crispr <- gene_crispr %>%
  filter(grepl("uter|endometr", lineage, ignore.case = TRUE) |
           grepl("uter|endometr", lineage_subtype, ignore.case = TRUE))

cat("\n[CHECK] Endometrial/uterine cell lines found:", length(unique(ec_crispr$cell_line)), "\n")
print(unique(ec_crispr$cell_line))

ec_summary <- ec_crispr %>%
  filter(!is.na(dependency)) %>%
  group_by(gene_name) %>%
  summarise(
    n_cell_lines = n(),
    mean_chronos = round(mean(dependency, na.rm = TRUE), 3),
    min_chronos  = round(min(dependency, na.rm = TRUE), 3),
    max_chronos  = round(max(dependency, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(mean_chronos)

cat("\n=== Endometrial/uterine-specific CRISPR essentiality ===\n")
print(ec_summary, digits = 3)

# VISUALIZE
plot_data <- gene_crispr %>% filter(!is.na(dependency))
plot_data$is_EC <- grepl("uter|endometr", plot_data$lineage, ignore.case = TRUE) |
  grepl("uter|endometr", plot_data$lineage_subtype, ignore.case = TRUE)

p <- ggplot(plot_data, aes(x = reorder(gene_name, dependency, median), y = dependency)) +
  geom_boxplot(outlier.size = 0.5, fill = "gray90") +
  geom_jitter(data = filter(plot_data, is_EC), aes(color = "Endometrial/Uterine line"),
              width = 0.15, size = 2) +
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Endometrial/Uterine line" = col_high), name = NULL) +
  labs(x = NULL, y = "CRISPR Gene Effect (Chronos Score)",
       title = "Essentiality of the 9 Signature Genes Across Cancer Cell Lines",
       subtitle = "Dashed red line = conventional essentiality threshold (-0.5)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")+
  theme(plot.title = element_text(face = "bold", size = 11))

p

ggsave("plots1/depmap_figures/DepMap_signature_gene_essentiality.png", p,
       width = 9, height = 6, dpi = 300, units = "in")
ggsave("plots1/depmap_figures/DepMap_signature_gene_essentiality.pdf", p,
       width = 9, height = 6, units = "in")

cat("\n[SAVED] plots1/depmap_figures/DepMap_signature_gene_essentiality.png (+ .pdf)\n")

saveRDS(gene_crispr, "results1/depmap_signature_gene_crispr.rds")
gene_crispr <- readRDS("results1/depmap_signature_gene_crispr.rds")   # checkpoint reload
write.csv(pan_cancer_summary, "results1/depmap_pancancer_summary.csv", row.names = FALSE)
write.csv(ec_summary, "results1/depmap_EC_summary.csv", row.names = FALSE)

# VERIFICATION: CRISPR essentiality numbers
pan_cancer_summary <- read.csv("results1/depmap_pancancer_summary.csv")
ec_summary <- read.csv("results1/depmap_EC_summary.csv")

cat("=== Pan-cancer summary ===\n")
print(pan_cancer_summary, digits = 4)

cat("\n=== Endometrial/uterine-specific summary ===\n")
print(ec_summary, digits = 4)

# Check 1: n_cell_lines consistency across genes (pan-cancer)
cat("\n[VERIFY] Pan-cancer n_cell_lines per gene:\n")
print(pan_cancer_summary[, c("gene_name", "n_cell_lines")])
cat("[VERIFY] Range of n_cell_lines across genes:", range(pan_cancer_summary$n_cell_lines), "\n")

# Check 2: n_cell_lines consistency (EC-specific)
cat("\n[VERIFY] EC-specific n_cell_lines per gene:\n")
print(ec_summary[, c("gene_name", "n_cell_lines")])

# Check 3: unique EC cell line count directly
ec_crispr_check <- gene_crispr %>%
  dplyr::filter(grepl("uter|endometr", lineage, ignore.case = TRUE) |
                  grepl("uter|endometr", lineage_subtype, ignore.case = TRUE))
cat("\n[VERIFY] Total unique EC/uterine cell lines (direct count):",
    length(unique(ec_crispr_check$cell_line)), "\n")

# Check 4: pct essential below -0.5 threshold, all 9 genes
cat("\n[VERIFY] % essential (Chronos <= -0.5) per gene, pan-cancer:\n")
print(pan_cancer_summary[, c("gene_name", "pct_essential_minus0.5")])

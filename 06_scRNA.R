# 06_scRNA_WL.R
#
# Single-cell validation (GSE173682): QC, SCTransform, clustering, annotation, module scores of the signature and patient-level pseudobulk tests.
# Input: raw GEO files in GSE173682_raw/ (GSE173682_RAW.tar)  |  Output: results/, results1/, plots1/scRNA/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  BiocManager::install("glmGamPoi", update = FALSE, ask = FALSE)
}
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(Seurat)
library(R.utils)
library(glmGamPoi)
library(future)
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 3000 * 1024^2)
plan("sequential")

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)
dir.create("plots1/scRNA", showWarnings = FALSE, recursive = TRUE)

# Project root: run this script from the project root directory
BASE_DIR <- "."

# EXTRACT AND PREPARE RAW GEO FILES
# DATA ACQUISITION: raw 10X files (barcodes/features/matrix per sample) for
# GSE173682 downloaded from GEO
# as GSE173682_RAW.tar before running this script.
getwd()
setwd(BASE_DIR)
untar(file.path(BASE_DIR, "GSE173682_raw/GSE173682_RAW.tar"),
      exdir = file.path(BASE_DIR, "GSE173682_raw"))
list.files(file.path(BASE_DIR, "GSE173682_raw"))

files <- list.files(file.path(BASE_DIR, "GSE173682_raw"), full.names = TRUE)
rna_files <- files[grepl("matrix|features|barcodes", files)]

sapply(rna_files, gunzip, overwrite = TRUE)
list.files(file.path(BASE_DIR, "GSE173682_raw"))

# ORGANIZE FILES INTO SAMPLE-WISE DIRECTORIES
base_dir <- file.path(BASE_DIR, "GSE173682_raw")
files <- list.files(base_dir, full.names = TRUE)

# Extract sample IDs (GSM numbers)
sample_ids <- unique(sub("_(barcodes|features|matrix).*", "", basename(files)))

for (sid in sample_ids) {
  
  sample_dir <- file.path(base_dir, sid)
  dir.create(sample_dir, showWarnings = FALSE)
  
  sample_files <- files[grepl(sid, files) & grepl("matrix|features|barcodes", files)]
  
  file.rename(sample_files,
              file.path(sample_dir, basename(sample_files)))
}

list.files(base_dir)
list.files(file.path(base_dir, "GSM5276933"))

# STANDARDIZE FILE NAMES FOR 10X FORMAT
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

for (dir in sample_dirs) {
  
  files <- list.files(dir, full.names = TRUE)
  
  file.rename(files[grepl("barcodes", files)],
              file.path(dir, "barcodes.tsv"))
  
  file.rename(files[grepl("features", files)],
              file.path(dir, "features.tsv"))
  
  file.rename(files[grepl("matrix", files)],
              file.path(dir, "matrix.mtx"))
}
list.files(file.path(base_dir, "GSM5276933"))

# LOAD DATA INTO SEURAT OBJECTS (PER SAMPLE)
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

seurat_list <- list()

for (dir in sample_dirs) {
  
  sample_name <- basename(dir)
  
  matrix <- ReadMtx(
    mtx = file.path(dir, "matrix.mtx"),
    features = file.path(dir, "features.tsv"),
    cells = file.path(dir, "barcodes.tsv")
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = matrix,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  seurat_obj$sample <- sample_name   # patient/GSM identifier -- used later
  
  seurat_list[[sample_name]] <- seurat_obj
}
length(seurat_list)
seurat_list[[1]]

# MERGE ALL SAMPLES INTO A SINGLE OBJECT
GSE173682_combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE173682_combined"
)
# SANITY CHECKS
GSE173682_combined
table(GSE173682_combined$sample)
head(colnames(GSE173682_combined))

# QUALITY CONTROL METRICS
GSE173682_combined[["percent.mt"]] <- PercentageFeatureSet(
  GSE173682_combined,  pattern = "^MT-")
head(GSE173682_combined@meta.data)

VlnPlot(
  GSE173682_combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
) &
  theme_classic(base_size = 13)

dev.copy(pdf, "plots1/scRNA/07_QC_metrics_violin.pdf", width = 10, height = 4)
dev.off()
dev.copy(png, filename = "plots1/scRNA/07_QC_metrics_violin.png",
         width = 10, height = 4, units = "in", res = 300)
dev.off()

# QUALITY CONTROL FILTERING
GSE173682_filtered <- subset(
  GSE173682_combined,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 7000 &
    percent.mt < 12
)

GSE173682_filtered
table(GSE173682_filtered$sample)

# NORMALIZATION (SCTransform)
GSE173682_filtered <- SCTransform(
  GSE173682_filtered,
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

# PCA
GSE173682_filtered <- RunPCA(GSE173682_filtered, verbose = FALSE)
GSE173682_filtered[["pca"]]
saveRDS(GSE173682_filtered, file = "results/GSE173682_postSCT_PCA.rds")
GSE173682_filtered <- readRDS("results/GSE173682_postSCT_PCA.rds")

pdf("plots1/scRNA/ElbowPlot.pdf", width = 6, height = 4)

ElbowPlot(GSE173682_filtered)

dev.off()

ElbowPlot(GSE173682_filtered)
dev.copy(png, filename = "plots1/scRNA/ElbowPlot.png",
         width = 6, height = 4, units = "in", res = 300)
dev.off()

# Elbow plot indicated that the first 15 principal components
# captured the majority of biologically meaningful variation
# and were retained for downstream analyses.
DimPlot(
  GSE173682_filtered,
  reduction = "pca",
  group.by = "sample"
)
dev.copy(pdf, "plots1/scRNA/10_PCA_by_sample.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/10_PCA_by_sample.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# FIND NEIGHBORS (BUILD CELL SIMILARITY GRAPH)
GSE173682_neighbors <- FindNeighbors(
  GSE173682_filtered,
  dims = 1:15
)
saveRDS(GSE173682_neighbors, file = "results/GSE173682_neighbors.rds",
        compress = FALSE)

# CLUSTERING (IDENTIFY CELL POPULATIONS)
GSE173682_neighbors <- readRDS("results/GSE173682_neighbors.rds")
GSE173682_clustered <- FindClusters(
  GSE173682_neighbors,
  resolution = 0.5
)
# Shows how many cells per cluster
table(GSE173682_clustered$seurat_clusters)
cat(
  "Number of clusters:",
  length(unique(Idents(GSE173682_clustered))),
  "\n"
)
saveRDS(GSE173682_clustered, file = "results/GSE173682_clustered.rds", compress = FALSE)

DimPlot(
  GSE173682_clustered,
  label = TRUE,
  repel = TRUE
)
dev.copy(pdf, "plots1/scRNA/12_Clusters_labeled.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/12_Clusters_labeled.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

# UMAP (2D VISUALIZATION OF CELL CLUSTERS)
GSE173682_clustered <- readRDS("results/GSE173682_clustered.rds")

GSE173682_umap <- RunUMAP(
  GSE173682_clustered,
  dims = 1:15
)

# Plot clusters with labels
DimPlot(
  GSE173682_umap,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  NoLegend() +
  theme_classic(base_size = 13)

dev.copy(pdf, "plots1/scRNA/13_UMAP_labeled.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/13_UMAP_labeled.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

DimPlot(
  GSE173682_umap,
  group.by="sample"
)
dev.copy(pdf, "plots1/scRNA/13_UMAP_by_sample.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/13_UMAP_by_sample.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

table(Idents(GSE173682_umap))

saveRDS(GSE173682_umap, file = "results/GSE173682_umap.rds", compress = FALSE)

GSE173682_umap <- readRDS("results/GSE173682_umap.rds")

# FIND MARKER GENES (IDENTIFY CLUSTER-SPECIFIC GENES)

# PREP SCT FOR MARKER ANALYSIS
GSE173682_prepped <- PrepSCTFindMarkers(GSE173682_umap)
saveRDS(GSE173682_prepped, file = "results/GSE173682_prepped.rds",
        compress = FALSE)

# FIND CLUSTER MARKERS
GSE173682_prepped <- readRDS("results/GSE173682_prepped.rds")

GSE173682_markers <- FindAllMarkers(
  GSE173682_prepped,
  only.pos = TRUE,        # only positive markers
  min.pct = 0.25,         # gene expressed in >=25% cells
  logfc.threshold = 0.25, # minimum log fold change
  test.use = "wilcox"
)

cat(
  "Number of marker genes:",
  nrow(GSE173682_markers),
  "\n"
)

cat(
  "Clusters:",
  length(unique(GSE173682_markers$cluster)),
  "\n"
)

saveRDS(GSE173682_markers, file = "results/GSE173682_markers.rds", compress = FALSE)

# TOP MARKERS
GSE173682_markers <- readRDS("results/GSE173682_markers.rds")
top10_markers <- GSE173682_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

head(top10_markers)

top5_all_clusters <- GSE173682_markers %>%
  group_by(cluster) %>%
  filter(
    pct.1 > 0.3 &
      pct.2 < 0.1
  ) %>%
  slice_max(
    n=5,
    order_by=avg_log2FC
  )

print(
  top5_all_clusters,
  n=150
)

write.csv(
  top10_markers,
  "results/top10_cluster_markers.csv",
  row.names = FALSE
)

write.csv(
  top5_all_clusters,
  "results/top5_specific_markers.csv",
  row.names = FALSE
)

top_marker <- GSE173682_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 1)

for(i in seq(1, length(top_marker$gene), by = 4)){
  
  p <- FeaturePlot(
    GSE173682_umap,
    features = top_marker$gene[i:min(i+3, length(top_marker$gene))]
  )
  
  ggsave(
    filename = paste0("plots1/scRNA/FeaturePlot_", i, ".pdf"),
    plot = p,
    width = 10,
    height = 8
  )
  ggsave(
    filename = paste0("plots1/scRNA/FeaturePlot_", i, ".png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# DOTPLOT
DotPlot(
  GSE173682_umap,
  features = unique(top5_all_clusters$gene)
) +
  
  RotatedAxis() +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  ) +
  
  ggtitle("Top Marker Genes across Annotated Cell Clusters")

dev.copy(pdf, "plots1/scRNA/17_DotPlot_top_markers.pdf", width = 10, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/17_DotPlot_top_markers.png",
         width = 10, height = 6, units = "in", res = 300)
dev.off()

FeaturePlot(
  GSE173682_umap,
  features = c(
    "EPCAM",   # tumor epithelial
    "CD3D",    # T cells
    "CD8A",    # CD8 T cells
    "CD68",    # macrophages
    "COL1A1",  # fibroblasts
    "VWF",     # endothelial
    "MS4A1",   # B cells
    "JCHAIN"   # plasma cells
  )
)
dev.copy(pdf, "plots1/scRNA/18_FeaturePlot_canonical_markers.pdf", width = 10, height = 8)
dev.off()
dev.copy(png, filename = "plots1/scRNA/18_FeaturePlot_canonical_markers.png",
         width = 10, height = 8, units = "in", res = 300)
dev.off()

# CLUSTER ANNOTATION
cluster_annotations <- c(
  "0"="CAF ECM",
  "1"="Macrophages M2",
  "2"="CAF inflammatory",
  "3"="Endothelial",
  "4"="Myofibroblasts",
  "5"="T NK cells",
  "6"="Tumor cells",
  "7"="Pericytes",
  "8"="CD8 T cells",
  "9"="Naive T cells",
  "10"="Tumor secretory",
  "11"="Tumor cells",
  "12"="T cells",
  "13"="Cycling cells",
  "14"="CAF ECM",
  "15"="Tumor epithelial",
  "16"="Tumor cells",
  "17"="Tumor stem like",
  "18"="Secretory epithelial",
  "19"="Activated CAF",
  "20"="WNT CAF",
  "21"="Ciliated epithelial",
  "22"="Tumor cells",
  "23"="Plasma cells",
  "24"="Dendritic cells",
  "25"="Activated T cells",
  "26"="Cycling tumor cells",
  "27"="B cells"
)

stopifnot(
  length(cluster_annotations) ==
    length(levels(GSE173682_umap))
)

# ADD CELL TYPE ANNOTATION TO OBJECT
GSE173682_annotated_v1 <- GSE173682_umap
head(colnames(GSE173682_annotated_v1))

clusters <- as.character(Idents(GSE173682_annotated_v1))

cell_types <- cluster_annotations[clusters]

# assign cell names
names(cell_types) <- colnames(GSE173682_annotated_v1)

GSE173682_annotated_v1 <- AddMetaData(
  GSE173682_annotated_v1,
  metadata = cell_types,
  col.name = "cell_type"
)

table(GSE173682_annotated_v1$cell_type)
table(
  Idents(GSE173682_annotated_v1),
  GSE173682_annotated_v1$cell_type
)

DimPlot(
  GSE173682_annotated_v1,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5
) + NoLegend()

dev.copy(pdf, "plots1/scRNA/20_DimPlot_annotated_celltype.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/20_DimPlot_annotated_celltype.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

DotPlot(
  GSE173682_annotated_v1,
  features = c(
    "EPCAM","KRT8","KRT18",          # epithelial/tumor lineage markers
    "COL1A1","POSTN","FAP",          # CAF
    "CD3D","CD8A","CCR7",            # T cells
    "MS4A1","JCHAIN",                # B/plasma
    "CD68","CD163",                  # macrophages
    "VWF","PECAM1",                  # endothelial
    "HLA-DRA"                        # Dendritic
    
  ),
  group.by = "cell_type"
) + RotatedAxis()

dev.copy(pdf, "plots1/scRNA/20_DotPlot_canonical_markers_by_celltype.pdf", width = 9, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/20_DotPlot_canonical_markers_by_celltype.png",
         width = 9, height = 5, units = "in", res = 300)
dev.off()

saveRDS(GSE173682_annotated_v1, file = "results/GSE173682_annotated_v1.rds", compress = FALSE)
GSE173682_annotated_v1 <- readRDS("results/GSE173682_annotated_v1.rds")

# CHECK MODULE SCORE
lasso_genes <- readRDS("results1/lasso_selected_genes_train.rds")
lasso_genes

present_genes <- lasso_genes[
  lasso_genes %in% rownames(GSE173682_annotated_v1)
]

present_genes
length(present_genes)
cat(
  "Genes detected:",
  length(present_genes),
  "of",
  length(lasso_genes),
  "\n"
)

GSE173682_annotated_v1_sc <- AddModuleScore(
  GSE173682_annotated_v1,
  features = list(present_genes),
  name = "RiskScore"
)

colnames(
  GSE173682_annotated_v1_sc@meta.data
)

summary(
  GSE173682_annotated_v1_sc$RiskScore1
)

hist(
  GSE173682_annotated_v1_sc$RiskScore1,
  breaks = 50,
  main = "Distribution of Signature Scores"
)
dev.copy(pdf, "plots1/scRNA/21_RiskScore_histogram.pdf", width = 6, height = 4)
dev.off()
dev.copy(png, filename = "plots1/scRNA/21_RiskScore_histogram.png",
         width = 6, height = 4, units = "in", res = 300)
dev.off()

# VISUALIZE FULL SIGNATURE SCORE
FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "RiskScore1",
  reduction = "umap",
  cols = c("lightgrey", "red"),
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster = FALSE
)
dev.copy(pdf, "plots1/scRNA/22_FeaturePlot_RiskScore1.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/22_FeaturePlot_RiskScore1.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# VIOLIN PLOT ACROSS CELL TYPES
VlnPlot(
  GSE173682_annotated_v1_sc,
  features = "RiskScore1",
  group.by = "cell_type",
  pt.size = 0
) +
  
  theme_classic(base_size = 15) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5, size = 15,
      face = "bold"
    )
  )

dev.copy(pdf, "plots1/scRNA/23_VlnPlot_RiskScore1_by_celltype.pdf", width = 9, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/23_VlnPlot_RiskScore1_by_celltype.png",
         width = 9, height = 5, units = "in", res = 300)
dev.off()

# AVERAGE SCORE PER CELL TYPE
avg_scores <-
  GSE173682_annotated_v1_sc@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    mean_score = mean(RiskScore1),
    median_score = median(RiskScore1),
    n_cells = n()
  ) %>%
  arrange(desc(mean_score))

kruskal_result <- kruskal.test(
  RiskScore1 ~ cell_type,
  data = GSE173682_annotated_v1_sc@meta.data
)

kruskal_result

# BARPLOT OF CELL TYPE SCORES
ggplot(
  avg_scores,
  aes(
    x = reorder(cell_type, mean_score),
    y = mean_score,
    fill = mean_score
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  
  coord_flip() +
  
  theme_classic(base_size = 16) +
  
  labs(
    title = "Average Prognostic Signature Score Across Cell Types",
    x = "",
    y = "Mean Signature Score"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.position = "none"
  )

dev.copy(pdf, "plots1/scRNA/25_Barplot_avg_scores.pdf", width = 8, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/25_Barplot_avg_scores.png",
         width = 8, height = 6, units = "in", res = 300)
dev.off()

# DEFINE RISK AND PROTECTIVE GENES
risk_genes <- c(
  "PTGIS",
  "TIMP3"
)

protective_genes <- c(
  "CTSW",
  "DIO2",
  "GZMK",
  "IHH",
  "PPP1R16B",
  "SPN",
  "ZNF683"
)

# KEEP ONLY PRESENT GENES
risk_genes <- risk_genes[
  risk_genes %in%
    rownames(GSE173682_annotated_v1_sc)
]

protective_genes <- protective_genes[
  protective_genes %in%
    rownames(GSE173682_annotated_v1_sc)
]

risk_genes
protective_genes

cat(
  "Risk genes detected:",
  length(risk_genes),
  "\n"
)

cat(
  "Protective genes detected:",
  length(protective_genes),
  "\n"
)

# ADD RISK PROGRAM SCORE
GSE173682_annotated_v1_sc <- AddModuleScore(
  GSE173682_annotated_v1_sc,
  features = list(risk_genes),
  name = "RiskProgram"
)

colnames(GSE173682_annotated_v1_sc@meta.data)

summary(
  GSE173682_annotated_v1_sc$RiskProgram1
)

# ADD PROTECTIVE PROGRAM SCORE
GSE173682_annotated_v1_sc <- AddModuleScore(
  GSE173682_annotated_v1_sc,
  features = list(protective_genes),
  name = "ProtectiveProgram"
)

colnames(GSE173682_annotated_v1_sc@meta.data)

summary(
  GSE173682_annotated_v1_sc$ProtectiveProgram1
)

# STATISTICAL TESTS FOR RISK AND PROTECTIVE PROGRAMS
risk_program_test <- kruskal.test(
  RiskProgram1 ~ cell_type,
  data = GSE173682_annotated_v1_sc@meta.data
)
risk_program_test

protective_program_test <- kruskal.test(
  ProtectiveProgram1 ~ cell_type,
  data = GSE173682_annotated_v1_sc@meta.data
)
protective_program_test

saveRDS(
  list(risk = risk_program_test, protective = protective_program_test),
  file = "results1/scrna_subprogram_kruskal_tests_PERCELL_NAIVE.rds"
)

naive_kruskal_tests <- readRDS("results1/scrna_subprogram_kruskal_tests_PERCELL_NAIVE.rds")

pseudobulk_scores <- GSE173682_annotated_v1_sc@meta.data %>%
  group_by(sample, cell_type) %>%
  summarise(
    n_cells             = n(),
    mean_RiskScore      = mean(RiskScore1),
    mean_RiskProgram     = mean(RiskProgram1),
    mean_ProtectiveProgram = mean(ProtectiveProgram1),
    .groups = "drop"
  )

saveRDS(pseudobulk_scores, "results1/scrna_pseudobulk_scores_by_patient_celltype.rds")
write.csv(pseudobulk_scores, "results1/scrna_pseudobulk_scores_by_patient_celltype.csv",
          row.names = FALSE)

pseudobulk_scores <- readRDS("results1/scrna_pseudobulk_scores_by_patient_celltype.rds")

# How many patients contribute to each cell type -- check before trusting
# a given cell type's pseudobulk comparison.
patients_per_celltype <- pseudobulk_scores %>%
  group_by(cell_type) %>%
  summarise(n_patients = n_distinct(sample)) %>%
  arrange(n_patients)
cat("\n=== Patients contributing per cell type (pseudobulk N) ===\n")
print(patients_per_celltype, n = Inf)

risk_program_test_pseudobulk <- kruskal.test(
  mean_RiskProgram ~ cell_type,
  data = pseudobulk_scores
)
cat("\n=== Risk program, PATIENT-LEVEL pseudobulk Kruskal-Wallis ===\n")
print(risk_program_test_pseudobulk)

protective_program_test_pseudobulk <- kruskal.test(
  mean_ProtectiveProgram ~ cell_type,
  data = pseudobulk_scores
)
cat("\n=== Protective program, PATIENT-LEVEL pseudobulk Kruskal-Wallis ===\n")
print(protective_program_test_pseudobulk)

kruskal_result_pseudobulk <- kruskal.test(
  mean_RiskScore ~ cell_type,
  data = pseudobulk_scores
)
cat("\n=== Combined signature score, PATIENT-LEVEL pseudobulk Kruskal-Wallis ===\n")
print(kruskal_result_pseudobulk)

saveRDS(
  list(
    risk_pseudobulk       = risk_program_test_pseudobulk,
    protective_pseudobulk = protective_program_test_pseudobulk,
    combined_pseudobulk   = kruskal_result_pseudobulk,
    patients_per_celltype = patients_per_celltype
  ),
  file = "results1/scrna_subprogram_kruskal_tests_PSEUDOBULK.rds"
)

pseudobulk_kruskal_tests <- readRDS("results1/scrna_subprogram_kruskal_tests_PSEUDOBULK.rds")

# VISUALIZE RISK PROGRAM
GSE173682_annotated_v1_sc$`Risk-Associated Program` <-
  GSE173682_annotated_v1_sc$RiskProgram1

FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "Risk-Associated Program",
  reduction = "umap",
  cols = c("lightgrey", "#D73027"),
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  theme_classic(base_size = 15) &
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )
dev.copy(pdf, "plots1/scRNA/31_FeaturePlot_RiskProgram.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/31_FeaturePlot_RiskProgram.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# VISUALIZE PROTECTIVE PROGRAM
GSE173682_annotated_v1_sc$`Protective Immune Program` <-
  GSE173682_annotated_v1_sc$ProtectiveProgram1

FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "Protective Immune Program",
  reduction = "umap",
  cols = c("lightgrey", "#4575B4"),
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  theme_classic(base_size = 15) &
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )
dev.copy(pdf, "plots1/scRNA/32_FeaturePlot_ProtectiveProgram.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/32_FeaturePlot_ProtectiveProgram.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# COMPARE PROGRAMS ACROSS CELL TYPES
options(repr.plot.width = 18, repr.plot.height = 7)
VlnPlot(
  GSE173682_annotated_v1_sc,
  features = c(
    "Risk-Associated Program",
    "Protective Immune Program"
  ),
  group.by = "cell_type",
  pt.size = 0
) &
  
  theme_classic(base_size = 15) &
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5, size = 15,
      face = "bold"
    )
  )

dev.copy(pdf, "plots1/scRNA/33_VlnPlot_RiskProtective_by_celltype.pdf", width = 12, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/33_VlnPlot_RiskProtective_by_celltype.png",
         width = 12, height = 5, units = "in", res = 300)
dev.off()

# DOTPLOT OF SIGNATURE GENES
DotPlot(
  GSE173682_annotated_v1_sc,
  features = present_genes,
  group.by = "cell_type"
) +
  RotatedAxis()

dev.copy(pdf, "plots1/scRNA/34_DotPlot_signature_genes.pdf", width = 8, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/34_DotPlot_signature_genes.png",
         width = 8, height = 5, units = "in", res = 300)
dev.off()

# SAVE FINAL OBJECT
saveRDS(
  GSE173682_annotated_v1_sc,
  file = "results1/GSE173682_signature_validation.rds",
  compress = FALSE
)

GSE173682_annotated_v1_sc <- readRDS("results1/GSE173682_signature_validation.rds")

stopifnot(nrow(GSE173682_annotated_v1_sc@meta.data) == 65221)
stopifnot(length(unique(GSE173682_annotated_v1_sc$seurat_clusters)) == 28)
stopifnot(length(unique(GSE173682_annotated_v1_sc$cell_type)) == 24)
cat("[CHECK] Final cells:", nrow(GSE173682_annotated_v1_sc@meta.data),
    "| clusters:", length(unique(GSE173682_annotated_v1_sc$seurat_clusters)),
    "| cell types:", length(unique(GSE173682_annotated_v1_sc$cell_type)), "\n")

# FINAL PLOTTING

# Annotated UMAP
DimPlot(
  GSE173682_annotated_v1_sc,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    
    legend.position = "none"
  ) +
  
  ggtitle("Cell-Type Annotation of UCEC Single-Cell Landscape")

dev.copy(pdf, "plots1/scRNA/FINAL_Celltype_Annotation_UMAP.pdf", width = 7, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Celltype_Annotation_UMAP.png",
         width = 7, height = 6, units = "in", res = 300)
dev.off()

DotPlot(
  GSE173682_annotated_v1_sc,
  features = c(
    "EPCAM","KRT8","KRT18",
    "COL1A1","POSTN","ACTA2",
    "CD3D","CD8A","CCR7",
    "MS4A1","JCHAIN",
    "CD68","CD163",
    "VWF","PECAM1"
  ),
  group.by = "cell_type"
) +
  
  RotatedAxis() +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  ) +
  
  ggtitle("Canonical Marker Expression Across Cell Types")

dev.copy(pdf, "plots1/scRNA/FINAL_Canonical_Marker_DotPlot.pdf", width = 9, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Canonical_Marker_DotPlot.png",
         width = 9, height = 5, units = "in", res = 300)
dev.off()

# Combined 9-Gene Signature FeaturePlot
FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "RiskScore1",
  reduction = "umap",
  cols = c("lightgrey", "red"),
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  
  theme_classic(base_size = 13) &
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  )&
  
  ggtitle("Cellular Distribution of the 9-Gene Prognostic Signature")

dev.copy(pdf, "plots1/scRNA/FINAL_Cellular_Distribution_9Gene_Signature.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Cellular_Distribution_9Gene_Signature.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# Mean Signature Score Barplot
ggplot(
  avg_scores,
  aes(
    x = reorder(cell_type, mean_score),
    y = mean_score,
    fill = mean_score
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  
  coord_flip() +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Average 9-Gene Prognostic Signature Score Across Cell Types",
    subtitle = paste0(
      "Per-cell Kruskal-Wallis p < 2.2e-16 (naive, exploratory) | ",
      "Patient-level pseudobulk p = ",
      signif(kruskal_result_pseudobulk$p.value, 3)
    ),
    x = "",
    y = "Mean Signature Score"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0.5
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    legend.position = "none"
  )

dev.copy(pdf, "plots1/scRNA/FINAL_Mean_Signature_Score_Barplot.pdf", width = 8, height = 6)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Mean_Signature_Score_Barplot.png",
         width = 8, height = 6, units = "in", res = 300)
dev.off()

kruskal_result <- kruskal.test(
  RiskScore1 ~ cell_type,
  data = GSE173682_annotated_v1_sc@meta.data
)

kruskal_result

# Risk vs Protective FeaturePlots
# Risk-associated
FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "Risk-Associated Program",
  reduction = "umap",
  cols = c("lightgrey", "#D73027"),
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  
  theme_classic(base_size = 13) &
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  )

dev.copy(pdf, "plots1/scRNA/FINAL_RiskAssociated_Program_FeaturePlot.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_RiskAssociated_Program_FeaturePlot.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# Protective immune
FeaturePlot(
  GSE173682_annotated_v1_sc,
  features = "Protective Immune Program",
  reduction = "umap",
  cols = c("lightgrey", "#4575B4"),
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  
  theme_classic(base_size = 13) &
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  )

dev.copy(pdf, "plots1/scRNA/FINAL_ProtectiveImmune_Program_FeaturePlot.pdf", width = 6, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_ProtectiveImmune_Program_FeaturePlot.png",
         width = 6, height = 5, units = "in", res = 300)
dev.off()

# Risk vs Protective Violin Plots
VlnPlot(
  GSE173682_annotated_v1_sc,
  features = c(
    "Risk-Associated Program",
    "Protective Immune Program"
  ),
  group.by = "cell_type",
  pt.size = 0
) &
  
  theme_classic(base_size = 13) &
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    )
  )

dev.copy(pdf, "plots1/scRNA/FINAL_Risk_vs_Protective_Violin.pdf", width = 12, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Risk_vs_Protective_Violin.png",
         width = 12, height = 5, units = "in", res = 300)
dev.off()

# Grouped Biology Comparison
GSE173682_annotated_v1_sc$major_group <- case_when(
  
  GSE173682_annotated_v1_sc$cell_type %in% c(
    "CD8 T cells",
    "T cells",
    "T NK cells",
    "Activated T cells",
    "Naive T cells",
    "B cells",
    "Plasma cells",
    "Dendritic cells",
    "Macrophages M2"
  ) ~ "Immune",
  
  GSE173682_annotated_v1_sc$cell_type %in% c(
    "CAF ECM",
    "CAF inflammatory",
    "Activated CAF",
    "WNT CAF",
    "Myofibroblasts",
    "Pericytes",
    "Endothelial"
  ) ~ "Stromal",
  
  TRUE ~ "Tumor"
)

kruskal.test(
  RiskScore1 ~ major_group,
  data = GSE173682_annotated_v1_sc@meta.data
)

major_group_lookup <- GSE173682_annotated_v1_sc@meta.data %>%
  dplyr::select(cell_type, major_group) %>%
  distinct()

pseudobulk_major_group <- pseudobulk_scores %>%
  left_join(major_group_lookup, by = "cell_type")

kruskal_major_group_pseudobulk <- kruskal.test(
  mean_RiskScore ~ major_group,
  data = pseudobulk_major_group
)
cat("\n=== major_group, PATIENT-LEVEL pseudobulk Kruskal-Wallis ===\n")
print(kruskal_major_group_pseudobulk)

ggplot(
  GSE173682_annotated_v1_sc@meta.data,
  aes(
    x = major_group,
    y = RiskScore1,
    fill = major_group
  )
) +
  
  geom_violin(trim = FALSE) +
  
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white"
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Distribution of Prognostic Signature Scores Across Major Cellular Compartments",
    subtitle = paste0(
      "Patient-level pseudobulk Kruskal-Wallis p = ",
      signif(kruskal_major_group_pseudobulk$p.value, 3)
    ),
    x = "",
    y = "9-Gene Prognostic Signature Score"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.1,
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0.1
    ),
    legend.position = "none"
  )

dev.copy(pdf, "plots1/scRNA/FINAL_Major_Cellular_Compartments_Violin.pdf", width = 7, height = 5)
dev.off()
dev.copy(png, filename = "plots1/scRNA/FINAL_Major_Cellular_Compartments_Violin.png",
         width = 7, height = 5, units = "in", res = 300)
dev.off()

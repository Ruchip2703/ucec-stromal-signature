# 01_mrna_raw_counts.R
#
# TCGA-UCEC: download STAR raw counts from GDC, keep primary tumours, keep one vial per patient (A > B > C) and strip Ensembl version suffixes.
# Output: results/mrna_counts.rds, results/mrna_counts_tumor.rds, results/mrna_counts_unique.rds

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("results", showWarnings = FALSE)

# DOWNLOAD RAW STAR-COUNTS DATA
query <- GDCquery(
  project       = "TCGA-UCEC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 20)

data <- GDCprepare(query)
saveRDS(data, "results/data from SummarizedExp.rds")

data <- readRDS("results/data from SummarizedExp.rds")

# SANITY CHECK: data object is a valid SummarizedExperiment
stopifnot(inherits(data, "SummarizedExperiment"))
cat("[CHECK] colData columns:", paste(head(colnames(colData(data)), 10), collapse = ", "), "...\n")

# EXTRACT RAW COUNT MATRIX
mrna_counts <- assay(data)
cat("[CHECK] mrna_counts dim:", paste(dim(mrna_counts), collapse = " x "))
cat("[CHECK] Value range sanity (raw counts should be non-negative integers):",
    min(mrna_counts, na.rm = TRUE), "to", max(mrna_counts, na.rm = TRUE), "\n")
stopifnot(min(mrna_counts, na.rm = TRUE) >= 0)

saveRDS(mrna_counts, "results/mrna_counts.rds")

# SAMPLE-TYPE BREAKDOWN, KEEP TUMOR SAMPLES ONLY
mrna_counts <- readRDS("results/mrna_counts.rds")
sample_types <- substr(colnames(mrna_counts), 14, 15)
cat("[CHECK] Sample type breakdown:\n"); print(table(sample_types))

tumor_id <- sample_types == "01"
mrna_counts_tumor <- mrna_counts[, tumor_id]

cat("[CHECK] Tumor-only samples:", ncol(mrna_counts_tumor), "\n")
stopifnot(ncol(mrna_counts_tumor) > 500)   # sanity floor, catches gross filtering errors

saveRDS(mrna_counts_tumor, "results/mrna_counts_tumor.rds")

# VIAL-PRIORITY DEDUPLICATION (A > B > C), ONE SAMPLE PER PATIENT
barcodes <- colnames(mrna_counts_tumor)
patient_ids <- substr(barcodes, 1, 12)

table_counts <- table(patient_ids)
duplicate_patients <- names(table_counts[table_counts > 1])
cat("[CHECK] Patients with multiple tumor aliquots:", length(duplicate_patients), "\n")

sample_info <- data.frame(
  barcode = barcodes,
  patient = patient_ids,
  vial    = substr(barcodes, 16, 16),
  stringsAsFactors = FALSE
)
sample_info$vial_priority <- match(sample_info$vial, c("A", "B", "C"))

sample_info <- sample_info[order(sample_info$patient, sample_info$vial_priority), ]
sample_info_unique <- sample_info[!duplicated(sample_info$patient), ]

# SANITY CHECK: dedup produced exactly one row per unique patient
stopifnot(nrow(sample_info_unique) == length(unique(sample_info$patient)))
stopifnot(sum(duplicated(sample_info_unique$patient)) == 0)
cat("[CHECK] Unique patients after vial-priority dedup:", nrow(sample_info_unique))

mrna_counts_unique <- mrna_counts_tumor[, sample_info_unique$barcode]
stopifnot(ncol(mrna_counts_unique) == nrow(sample_info_unique))

# ENSEMBL ID CLEANUP -- STRIP VERSION SUFFIX, COLLAPSE DUPLICATES
rownames(mrna_counts_unique) <- sub("\\..*", "", rownames(mrna_counts_unique))
n_dup_genes_before <- sum(duplicated(rownames(mrna_counts_unique)))
cat("[CHECK] Duplicate gene IDs after version-stripping:", n_dup_genes_before, "\n")

mrna_counts_unique <- rowsum(mrna_counts_unique, group = rownames(mrna_counts_unique))

# SANITY CHECK: collapse worked, no duplicate rownames remain
stopifnot(sum(duplicated(rownames(mrna_counts_unique))) == 0)
cat("[CHECK] Final matrix:", paste(dim(mrna_counts_unique), collapse = " x "))

saveRDS(mrna_counts_unique, "results/mrna_counts_unique.rds")
cat("\n[SAVED] results/mrna_counts_unique.rds\n")

# 02_stromal_deg_WL.R
#
# TCGA-UCEC: filtering and VST normalisation, ESTIMATE scores, 75/25 train/test split (seed 123) and High vs Low stromal limma DEGs from the training set only.
# Input: results/mrna_counts.rds, mrna_counts_unique.rds, clinical_patient.rds, survival_patient.rds  |  Output: results1/

library(dplyr)
library(org.Hs.eg.db)
library(DESeq2)
library(edgeR)
library(limma)
library(estimate)   # filterCommonGenes(), estimateScore()

set.seed(123)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

# BUILD COMBINED TUMOR + NORMAL COUNT MATRIX
mrna_counts         <- readRDS("results/mrna_counts.rds")
mrna_counts_unique  <- readRDS("results/mrna_counts_unique.rds")

normal_ids    <- substr(colnames(mrna_counts), 14, 15) == "11"
normal_counts <- mrna_counts[, normal_ids]
rownames(normal_counts) <- sub("\\..*", "", rownames(normal_counts))

normal_counts_aligned <- normal_counts[rownames(mrna_counts_unique), ]

deg_counts <- cbind(mrna_counts_unique, normal_counts_aligned)
dim(deg_counts)   # expect 60,616 genes x 580 samples (545 tumor + 35 normal)

condition <- c(
  rep("Tumor",  ncol(mrna_counts_unique)),
  rep("Normal", ncol(normal_counts_aligned))
)
table(condition)   # expect Normal=35, Tumor=545

# DESeq2 SETUP, CPM FILTERING, SIZE-FACTOR/DISPERSION FIT
# The DESeq() fit is required infrastructure for vst(..., blind=FALSE) later
# (dispersions must already be estimated). Its Tumor-vs-Normal test results
# are not used or reported anywhere downstream. This normalization step uses
# no survival/outcome/grouping labels, so it is NOT a feature-selection step
# and does not need to be restricted to training data.
col_data <- data.frame(
  row.names = colnames(deg_counts),
  condition = factor(condition)
)

storage.mode(deg_counts) <- "numeric"

dds <- DESeqDataSetFromMatrix(
  countData = deg_counts,
  colData   = col_data,
  design    = ~condition
)
dim(dds)   # expect 60,616 x 580

# CPM filter: keep genes with CPM >= 1 in at least 10 of the 580 samples
keep <- rowSums(cpm(counts(dds)) >= 1) >= 10
dds  <- dds[keep, ]
dim(dds)   # expect 23,389 x 580

dds <- DESeq(dds)
saveRDS(dds, "results1/dds_after_DEG.rds")

# SUBSET TO TUMOR-ONLY, VARIANCE-STABILIZING TRANSFORM
dds <- readRDS("results1/dds_after_DEG.rds")
dds_tumor <- dds[, dds$condition == "Tumor"]
design(dds_tumor) <- ~1

vsd      <- vst(dds_tumor, blind = FALSE)
expr_tum <- assay(vsd)
dim(expr_tum)   # expect 23,389 x 545

# ENSEMBL -> GENE SYMBOL MAPPING + DUPLICATE COLLAPSE
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(expr_tum),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)
table(is.na(gene_symbols))   # expect FALSE=19,642, TRUE=3,747

expr_tum_df <- as.data.frame(expr_tum)
expr_tum_df$GeneSymbol <- gene_symbols

# Drop unmapped genes
expr_tum_df <- expr_tum_df[!is.na(expr_tum_df$GeneSymbol) &
                             expr_tum_df$GeneSymbol != "", ]
dim(expr_tum_df)   # expect 19,642 x 546

# Collapse duplicate gene symbols: keep the highest-variance row per symbol
gene_variance <- rowVars(as.matrix(expr_tum_df[, -ncol(expr_tum_df)]))
expr_tum_df$variance <- gene_variance
expr_tum_df <- expr_tum_df[order(expr_tum_df$GeneSymbol, -expr_tum_df$variance), ]
expr_tum_df <- expr_tum_df[!duplicated(expr_tum_df$GeneSymbol), ]
expr_tum_df$variance <- NULL

rownames(expr_tum_df) <- expr_tum_df$GeneSymbol
expr_tum_df$GeneSymbol <- NULL
expr_tum_df <- as.matrix(expr_tum_df)
dim(expr_tum_df)   # expect 19,616 x 545

saveRDS(expr_tum_df, "results1/expr_tumor_vst_genesymbol.rds")

# DETERMINE ANALYTIC COHORT (n=541) + TRAIN/TEST SPLIT
clinical_patient <- readRDS("results/clinical_patient.rds")   # raw input
survival_patient  <- readRDS("results/survival_patient.rds")  # raw input
survival_patient$patient_id <- substr(survival_patient$sample, 1, 12)

expr_patient_ids <- substr(colnames(expr_tum_df), 1, 12)

common_patients <- Reduce(
  intersect,
  list(
    expr_patient_ids,
    survival_patient$patient_id,
    clinical_patient$patient_id
  )
)
cat("Common patients across all datasets:", length(common_patients), "\n")
stopifnot(length(common_patients) == 541)

set.seed(123)
train_index    <- sample(1:length(common_patients), size = 0.75 * length(common_patients))
train_patients <- common_patients[train_index]
test_patients  <- common_patients[-train_index]
stopifnot(length(train_patients) == 405, length(test_patients) == 136)

saveRDS(
  list(
    common_patients = common_patients,
    train_patients  = train_patients,
    test_patients   = test_patients
  ),
  "results1/analytic_cohort_train_test_split.rds"
)
cat("[SAVED] results1/analytic_cohort_train_test_split.rds\n\n")

# RUN ESTIMATE (STROMAL / IMMUNE / ESTIMATE SCORES + TUMOR PURITY)
estimate_input <- data.frame(
  NAME        = rownames(expr_tum_df),
  Description = rownames(expr_tum_df),
  expr_tum_df
)
write.table(
  estimate_input,
  file      = "estimate_input_tumor1.gct",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

filterCommonGenes(
  input.f = "estimate_input_tumor1.gct",
  output.f = "estimate_tumor_filtered1.gct",
  id = "GeneSymbol"
)

estimateScore(
  input.ds  = "estimate_tumor_filtered1.gct",
  output.ds = "estimate_scores_tum1.gct"
)

tum_est_scores <- read.table(
  "estimate_scores_tum1.gct",
  skip = 2, header = TRUE, sep = "\t", check.names = FALSE
)

tum_est_df <- as.data.frame(t(tum_est_scores[, -c(1, 2)]))
colnames(tum_est_df) <- tum_est_scores$NAME
tum_est_df$sample <- rownames(tum_est_df)
tum_est_df <- tum_est_df[-1, ]
dim(tum_est_df)   # expect 545 x 5 (StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity, sample)

saveRDS(tum_est_df, "results1/tumor_estimate_scores.rds")

# STROMAL GROUPING -- MEDIAN CUTOFF FROM TRAINING SET ONLY
tum_est_df <- readRDS("results1/tumor_estimate_scores.rds")
tum_est_df$sample <- gsub("\\.", "-", tum_est_df$sample)
tum_est_df$patient_id <- substr(tum_est_df$sample, 1, 12)
tum_est_df$StromalScore <- as.numeric(tum_est_df$StromalScore)

split <- readRDS("results1/analytic_cohort_train_test_split.rds")

train_stromal_scores <- tum_est_df$StromalScore[tum_est_df$patient_id %in% split$train_patients]
median_stromal_train <- median(train_stromal_scores)
cat("[CHECK] Training-set median StromalScore (n=", length(train_stromal_scores),
    "):", median_stromal_train, "\n")

tum_est_df$stromal_group <- ifelse(
  tum_est_df$StromalScore >= median_stromal_train, "High", "Low"
)
table(tum_est_df$stromal_group)

# Align ESTIMATE score table to expression matrix column order
tum_est_df <- tum_est_df[match(colnames(expr_tum_df), tum_est_df$sample), ]
stopifnot(all(colnames(expr_tum_df) == tum_est_df$sample))

saveRDS(tum_est_df, "results1/tumor_estimate_scores_grouped.rds")

# STROMAL HIGH vs LOW DIFFERENTIAL EXPRESSION (limma)
train_sample_ids <- tum_est_df$sample[tum_est_df$patient_id %in% split$train_patients]
train_cols        <- colnames(expr_tum_df) %in% train_sample_ids

expr_tum_train <- expr_tum_df[, train_cols]
tum_est_train  <- tum_est_df[match(colnames(expr_tum_train), tum_est_df$sample), ]
stopifnot(ncol(expr_tum_train) == 405)
stopifnot(all(colnames(expr_tum_train) == tum_est_train$sample))

group_stromal <- factor(tum_est_train$stromal_group, levels = c("Low", "High"))

design_stromal <- model.matrix(~0 + group_stromal)
colnames(design_stromal) <- c("Low", "High")

fit  <- lmFit(expr_tum_train, design_stromal)
contrast_matrix <- makeContrasts(High - Low, levels = design_stromal)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

stromal_deg <- topTable(fit2, adjust.method = "BH", number = Inf)
dim(stromal_deg)   # expect 19,616 x 6 (same gene universe; test statistics
saveRDS(stromal_deg, "results1/stromal_deg_full.rds")

stromal_sig_deg <- stromal_deg %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)
dim(stromal_sig_deg)
saveRDS(stromal_sig_deg, "results1/stromal_sig_deg.rds")

# EXPRESSION MATRIX FOR CANDIDATE STROMAL DEGs
cohort_sample_ids <- tum_est_df$sample[tum_est_df$patient_id %in% split$common_patients]

stromal_deg_expr <- expr_tum_df[rownames(expr_tum_df) %in% rownames(stromal_sig_deg), cohort_sample_ids]
stromal_deg_expr <- as.data.frame(t(stromal_deg_expr))
rownames(stromal_deg_expr) <- substr(rownames(stromal_deg_expr), 1, 12)
dim(stromal_deg_expr)

saveRDS(stromal_deg_expr, "results1/stromal_deg_expr_matrix.rds")

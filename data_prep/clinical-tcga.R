# clinical-tcga.R
#
# TCGA-UCEC clinical data: keep primary tumour samples, one aliquot per patient (01A > 01B > 01C).
# Input: data/TCGA-UCEC clinical.xlsx (UCSC Xena)  |  Output: results/clinical_patient.rds

options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)

clinical_tcga <- read_excel(
  "data/TCGA-UCEC clinical.xlsx"
)

# Basic inspection
dim(clinical_tcga)
colnames(clinical_tcga)
head(clinical_tcga)

# Keep only primary tumor samples
library(dplyr)

clinical_tcga_tumor <- clinical_tcga %>%
  filter(grepl("-01", sample))

# Check dimensions
dim(clinical_tcga_tumor)
colnames(clinical_tcga_tumor)
head(clinical_tcga_tumor)

# Confirm only primary tumors
table(substr(clinical_tcga_tumor$sample, 14, 15))

# Ensure no normals or recurrent tumors
any(grepl("-11", clinical_tcga_tumor$sample))
any(grepl("-02", clinical_tcga_tumor$sample))

clinical_tcga_tumor <- clinical_tcga_tumor %>%
  mutate(
    patient_id = substr(sample, 1, 12),
    aliquot = substr(sample, 14, 16),
    aliquot_priority = case_when(
      aliquot == "01A" ~ 1,
      aliquot == "01B" ~ 2,
      aliquot == "01C" ~ 3,
      TRUE ~ 99
    )
  )

clinical_tcga_tumor %>%
  dplyr::count(patient_id) %>%
  dplyr::filter(n > 1)

clinical_patient <- clinical_tcga_tumor %>%
  dplyr::arrange(patient_id, aliquot_priority) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::slice(1) %>%
  ungroup()

any(duplicated(clinical_patient$patient_id))  

saveRDS(clinical_patient, "results/clinical_patient.rds")

dim(clinical_patient)

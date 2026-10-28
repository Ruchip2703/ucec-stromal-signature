# survival_tcga.R
#
# TCGA-UCEC survival data: keep primary tumour samples, one aliquot per patient (01A > 01B > 01C).
# Input: data/TCGA-UCEC survival.xlsx (UCSC Xena)  |  Output: results/survival_patient.rds

library(readxl)
library(dplyr)

# Load TCGA-UCEC survival data
survival_tcga <- read_excel(
  "data/TCGA-UCEC survival.xlsx")

# Basic inspection
dim(survival_tcga)
colnames(survival_tcga)
head(survival_tcga)

# Number of unique patients
length(unique(survival_tcga$`_PATIENT`))

# Number of duplicated patients
sum(duplicated(survival_tcga$`_PATIENT`))

# List duplicated patients
unique(survival_tcga$`_PATIENT`[
  duplicated(survival_tcga$`_PATIENT`)
])

# Keep only primary tumour samples
survival_tcga_primary <- survival_tcga %>%
  filter(substr(sample, 14, 15) == "01")

# Check dimensions
dim(survival_tcga_primary)

# Confirm only primary tumour samples remain
table(substr(survival_tcga_primary$sample, 14, 16))

# Ensure no recurrent or normal samples remain
any(substr(survival_tcga_primary$sample, 14, 15) == "02")
any(substr(survival_tcga_primary$sample, 14, 15) == "11")

# Create patient IDs and aliquot priority
survival_tcga_primary <- survival_tcga_primary %>%
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

# Check duplicated patients after filtering
survival_tcga_primary %>%
  dplyr::count(patient_id) %>%
  dplyr::filter(n > 1)

# Retain one aliquot per patient (01A > 01B > 01C)
survival_patient <- survival_tcga_primary %>%
  arrange(patient_id, aliquot_priority) %>%
  distinct(patient_id, .keep_all = TRUE)

# Final quality control
dim(survival_patient)

length(unique(survival_patient$patient_id))

sum(duplicated(survival_patient$patient_id))

head(survival_patient)

# Save processed survival dataset
saveRDS(
  survival_patient,
  "results/survival_patient.rds"
)

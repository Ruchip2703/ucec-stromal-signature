# Table1_WL.R
#
# Table 1: clinical characteristics by risk group.
# Input: results/clinical_patient.rds, results1/full_survival_with_risk.rds, results1/risk_subtype_tmb_hrd.rds, data/TCGA Firehose.xlsx (cBioPortal clinical export)  |  Output: results1/Table1_Clinical_Characteristics.docx

if (!requireNamespace("gtsummary", quietly = TRUE)) install.packages("gtsummary")
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")

library(dplyr)
library(gtsummary)
library(readxl)
library(flextable)
library(officer)

dir.create("results1", showWarnings = FALSE, recursive = TRUE)

clinical_patient <- readRDS("results/clinical_patient.rds")
survival_full <- readRDS("results1/full_survival_with_risk.rds")

# PULL BMI DIRECTLY FROM THE FIREHOSE FILE
firehose_raw <- read_excel("data/TCGA Firehose.xlsx", sheet = "ucec_tcga_clinical_data")

dim(firehose_raw)
colnames(firehose_raw)

stopifnot(all(c("Patient ID", "BMI (kg/m2)") %in% colnames(firehose_raw)))

bmi_data <- firehose_raw %>%
  dplyr::select(patient_id = `Patient ID`, BMI = `BMI (kg/m2)`) %>%
  mutate(BMI = as.numeric(BMI)) %>%
  distinct(patient_id, .keep_all = TRUE)   # one row per patient (file has per-sample duplicates)

cat("[CHECK] Unique patients with BMI pulled:", nrow(bmi_data), "\n")
cat("[CHECK] Non-missing BMI:", sum(!is.na(bmi_data$BMI)), "\n")

risk_full <- readRDS("results1/risk_subtype_tmb_hrd.rds")
id_col <- intersect(c("patient_id", "tcga_patient"), colnames(risk_full))[1]
stopifnot(!is.na(id_col), "Subtype_Selected" %in% colnames(risk_full))

subtype_map <- risk_full %>%
  dplyr::transmute(patient_id = .data[[id_col]], Subtype_Selected) %>%
  distinct(patient_id, .keep_all = TRUE)

# PanCancer Atlas labels -> clinical names. "UCEC.NA" is a literal string in the
# Atlas (tumor with no assigned subtype), NOT an R NA, so it is kept as its own level.
subtype_labels <- c(
  "UCEC.POLE"    = "POLE",
  "UCEC.MSI"     = "MSI (MMRd)",
  "UCEC.CN_LOW"  = "CN-low (NSMP)",
  "UCEC.CN_HIGH" = "CN-high (p53abn)",
  "UCEC.NA"      = "No subtype call"
)
stopifnot(all(na.omit(unique(as.character(subtype_map$Subtype_Selected))) %in% names(subtype_labels)))

# BUILD TABLE 1 DATA
table1_data <- clinical_patient %>%
  dplyr::select(
    patient_id,
    Age               = age_at_index.demographic,
    Race              = race.demographic,
    Ethnicity         = ethnicity.demographic,
    `FIGO Stage`      = figo_stage.diagnoses,
    `Tumor Grade`     = tumor_grade.diagnoses,
    `Histologic Type` = primary_diagnosis.diagnoses,
    `Vital Status`    = vital_status.demographic
  ) %>%
  left_join(bmi_data, by = "patient_id") %>%
  inner_join(survival_full %>% dplyr::select(patient_id, risk_group), by = "patient_id") %>%
  left_join(subtype_map, by = "patient_id") %>%
  mutate(`Molecular Subtype` = factor(unname(subtype_labels[as.character(Subtype_Selected)]),
                                      levels = unname(subtype_labels))) %>%
  dplyr::select(-Subtype_Selected)

stopifnot(nrow(table1_data) == 541)
cat("[CHECK] Table 1 cohort: n =", nrow(table1_data), "\n")

cat("[CHECK] Molecular subtype by risk group (NA = not in PanCancer Atlas merge):\n")
print(table(table1_data$`Molecular Subtype`, table1_data$risk_group, useNA = "ifany"))

# Filter implausible BMI values (standard adult range 10-70)
cat("[CHECK] BMI values outside plausible range (10-70):\n")
print(table1_data$patient_id[which(table1_data$BMI < 10 | table1_data$BMI > 70)])
print(table1_data$BMI[which(table1_data$BMI < 10 | table1_data$BMI > 70)])
table1_data$BMI[table1_data$BMI < 10 | table1_data$BMI > 70] <- NA
cat("[CHECK] BMI summary after filtering:\n"); print(summary(table1_data$BMI))

cat("[CHECK] Missingness per variable:\n")
print(sapply(table1_data, function(x) sum(is.na(x) | x %in% c("not reported", "Not Reported", "unknown"))))

table1_data$`FIGO Stage Group` <- ifelse(
  table1_data$`FIGO Stage` %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                  "Stage IIIC1", "Stage IIIC2", "Stage IV", "Stage IVA", "Stage IVB"),
  "Advanced (III-IV)", "Early (I-II)"
)

table1_data <- table1_data %>% dplyr::relocate(`Molecular Subtype`, .after = `FIGO Stage Group`)

# BUILD AND EXPORT TABLE 1
# Tumor Grade dropped: 100% "not reported", uninformative, breaks add_p()
table1 <- table1_data %>%
  dplyr::select(-patient_id, -`FIGO Stage`, -`Tumor Grade`) %>%
  tbl_summary(
    by = risk_group,
    statistic = list(Age ~ "{mean} ({sd})", BMI ~ "{mean} ({sd})"),
    digits = list(BMI ~ 1),
    missing = "ifany",
    missing_text = "Not reported"
  ) %>%
  add_p() %>%
  add_overall() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_caption("**Table 1. Clinical Characteristics by Risk Group**") %>%
  modify_footnote(label ~ paste("Molecular subtype per TCGA PanCancer Atlas: MSI = MMRd, CN-low = NSMP,",
                                "CN-high = p53abn; 'No subtype call' = tumors without an assigned subtype;",
                                "'Not reported' = patients absent from the PanCancer Atlas subtype table.")) %>%
  bold_labels()

table1

table1_flex <- as_flex_table(table1)
save_as_docx(table1_flex, path = "results1/Table1_Clinical_Characteristics.docx")
cat("[SAVED] results1/Table1_Clinical_Characteristics.docx\n")

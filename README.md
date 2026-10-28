# Stromal-associated 9-gene prognostic signature in endometrial carcinoma (TCGA-UCEC)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23722384.svg)]
(https://doi.org/10.5281/zenodo.23722384)

Analysis code for a stromal-associated 9-gene prognostic signature in uterine corpus endometrial carcinoma. The signature was derived in TCGA-UCEC (ESTIMATE-based stromal stratification, limma, univariate Cox and LASSO-Cox) and replicated in CPTAC-UCEC.

- Risk genes: TIMP3, PTGIS
- Protective genes: IHH, PPP1R16B, CTSW, ZNF683, DIO2, GZMK, SPN

## How to run

1. Open R in the project root (all paths are relative to it).
2. Create the folders `data/`, `results/`, `results1/` and `plots1/` (most are created automatically by the scripts).
3. Put the required input files in `data/` (see Data sources).
4. Run the scripts in the order below. Each script reads the outputs of the earlier ones.

Random seed: `set.seed(123)` throughout. Analytic cohort: n = 541 patients, split into 405 training and 136 test patients.

## Scripts and order

| Step | Script | Purpose |
|---|---|---|
| 0 | `clinical-tcga.R`, `survival_tcga.R` | Clinical and survival tables, one primary-tumour aliquot per patient |
| 1 | `01_mrna_raw_counts.R` | GDC STAR raw counts, tumour filtering, de-duplication |
| 2 | `02_stromal_deg.R` | VST, ESTIMATE, train/test split, stromal DEGs (training set only) |
| 3 | `03_cox_lasso.R` | Univariate Cox, LASSO-Cox, risk score, Kaplan-Meier, ROC, multivariable Cox |
| 4 | `04_gsva.R` | GSVA Hallmark pathways by risk group |
| 5 | `05_cibersort.R` | CIBERSORT (LM22) immune deconvolution |
| 6 | `06_scRNA.R` | Single-cell validation (GSE173682) |
| 7 | `07_tmb_hrd.R` | Molecular subtype, TMB and HRD |
| 8 | `08_tide_ips.R` | TIDE (TIDEpy) and Immunophenoscore |
| 9 | `09_checkpoint_exclusion.R` | Checkpoint genes, dysfunction vs exclusion |
| 10 | `10_oncopredict.R` | Predicted drug sensitivity (GDSC2) |
| 11 | `11_mediation.R` | Mediation analysis |
| 12 | `12_crispr_check.R` | DepMap CRISPR essentiality |
| 13 | `13_cptac_external_validation.R`, `13b_cptac_gsva_tide_extension_WL.R` | CPTAC-UCEC replication |
| T | `Table1.R` | Table 1 |
| F | `03a_stromal_lasso_figures.R`, `tmb_hrd_figures.R`, `tide_figures.R`, `drug_sensitivity_figures.R`, `cptac_figures.R` | Individual figure panels |
| F | `paper_fig.R` | Assembles the multi-panel figures (run last) |

## Data sources (not included in this repository)

| Data | Source | Used in |
|---|---|---|
| TCGA-UCEC RNA-seq (STAR counts) and MAF | GDC Data Release 45.0 (11 Feb 2026), via TCGAbiolinks | 01, 07 |
| TCGA-UCEC clinical and survival tables | UCSC Xena -> `data/TCGA-UCEC clinical.xlsx`, `data/TCGA-UCEC survival.xlsx` | 00 |
| TCGA clinical export including BMI | cBioPortal -> `data/TCGA Firehose.xlsx` | Table 1 |
| Molecular subtypes | TCGA PanCancer Atlas (`TCGAbiolinks::PanCancerAtlas_subtypes()`) | 03, 07 |
| HRD scores | Marquard et al. 2015, Biomark Res (Additional file 8), downloaded automatically | 07 |
| CPTAC-UCEC | cBioPortal study `ucec_cptac_2020` -> `ucec_cptac_2020/ucec_cptac_2020/` | 13, 13b |
| Single-cell RNA-seq | GEO GSE173682 (`GSE173682_RAW.tar` in `GSE173682_raw/`) | 06 |
| GDSC2 training data | https://osf.io/c6tfx/ (`DataFiles.zip`, unzip into `oncopredict_data/`) | 10 |
| IPS reference files | github.com/icbi-lab/Immunophenogram (downloaded automatically) | 08 |
| DepMap CRISPR (Chronos) | `depmap` Bioconductor package via ExperimentHub | 12 |

## CIBERSORT

`CIBERSORT.R` and the LM22 signature matrix (`LM22.txt`) are distributed under the CIBERSORT licence and are not redistributed here. Register at https://cibersortx.stanford.edu, download both files and place them in the project root before running `05_cibersort.R`.

## Software

- R 4.5.1. Key packages: TCGAbiolinks 2.38.0, DESeq2, limma, glmnet, survival, GSVA 2.4.9, Seurat, oncoPredict, mediation, rms, dcurves.
- TIDE is run with `tidepy` in a Python 3.11 conda environment (`tide_env`), created automatically by `08_tide_ips.R` through `reticulate`.
- Exact package versions are recorded in `sessionInfo.txt`.

## Outputs

Data objects and tables are written to `results/` and `results1/`, and figures to `plots1/`. These folders are not tracked in this repository.

## Licence and citation

Code is released under the MIT licence. If you use this code, please cite the associated manuscript and this repository.
Archived at Zenodo: https://doi.org/10.5281/zenodo.23722384

## Contact

Ruchi Punjabi, ruchipunjabi27@gmail.com

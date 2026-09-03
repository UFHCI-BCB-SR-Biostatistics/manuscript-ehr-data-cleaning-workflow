# EHR Data Cleaning and Processing Workflow

## Overview

This repository contains R code developed to support the data processing and cleaning workflow described in the manuscript **“Managing and Cleaning Large-Scale Electronic Health Record Data.”**

The code documents a reproducible workflow for preparing large-scale electronic health record (EHR) data for analysis. The primary data-processing workflow consists of four sequential steps: cohort construction, laboratory biomarker identification and standardization, biomarker condensing and monthly aggregation, and biomarker merging.

Additional R scripts used to conduct the statistical analyses presented in the manuscript are provided separately in the `Analysis/` directory.

All data processing and statistical analysis code included in this repository was implemented in **R**.

## Data Availability

The patient-level EHR data used in this study are not included in this repository because the underlying data are restricted.

Access to the source data is subject to the requirements and approvals of the original data provider. Accordingly, this repository contains the code used to implement the data-processing and analysis workflow but does not contain patient-level data.

Local directory structures and source filenames have been replaced with generic placeholders in the public code. Users with appropriate access to comparable data should specify the corresponding local paths and filenames before running the scripts.

## Repository Structure

The repository is organized as follows:

```text
.
├── Code/
│   ├── Demographics and Diagnosis/
│   │   └── demographics_diagnosis_cleaning.R
│   │
│   ├── Laboratory Biomarkers/
│   │   ├── 01_obtaining_biomarkers.R
│   │   ├── 02_condensing_biomarkers.R
│   │   └── 03_merging_biomarkers.R
│   │
│   └── Analysis/
│       └── merge_biomarkers_demographics.R
│
├── Analysis/
│   ├── [FIRST_ANALYSIS_FILE].R
│   └── [SECOND_ANALYSIS_FILE].R
│
└── README.md
```

The `Code/` directory contains scripts used for cohort construction, EHR data cleaning and processing, and preparation of the final analysis-ready dataset. The separate `Analysis/` directory contains the R scripts used to conduct the statistical analyses presented in the manuscript.

## Primary Data-Processing Workflow

The primary data-processing workflow consists of the following four sequential steps:

### 1. Cohort Construction

`Code/Demographics and Diagnosis/demographics_diagnosis_cleaning.R`

Constructs the study cohort using demographic and diagnosis data. The script identifies liver disease and hepatocellular carcinoma (HCC) diagnoses, determines the first qualifying diagnosis dates, applies cohort eligibility criteria, and performs quality-control checks.

### 2. Laboratory Biomarker Identification and Standardization

`Code/Laboratory Biomarkers/01_obtaining_biomarkers.R`

Identifies laboratory records corresponding to the biomarkers evaluated in the study and standardizes laboratory names and measurement units where appropriate. The script also performs quality-control checks on laboratory names, measurement units, and biomarker value distributions.

The biomarkers processed include:

- AST
- ALT
- Albumin
- Total bilirubin
- Total cholesterol
- LDL cholesterol
- HDL cholesterol
- Triglycerides
- Hemoglobin
- INR
- Platelets
- Alpha-fetoprotein
- Neutrophils
- Lymphocytes

### 3. Biomarker Condensing and Monthly Aggregation

`Code/Laboratory Biomarkers/02_condensing_biomarkers.R`

Restricts biomarker measurements to the study-specific observation periods and applies biomarker-specific value thresholds.

For patients with HCC, measurements are evaluated between the first liver disease diagnosis and first HCC diagnosis. For patients without HCC, measurements following the first liver disease diagnosis are evaluated.

When multiple measurements of the same biomarker are available within a calendar month, measurements are summarized to create monthly biomarker values.

This script also implements the lipid-component calculation used to derive a missing lipid value when the other required lipid components are available.

### 4. Biomarker Merging

`Code/Laboratory Biomarkers/03_merging_biomarkers.R`

Combines the processed monthly biomarker datasets and reshapes the data into a subject-month format, with biomarker measurements represented as separate variables for subsequent analysis.

## Final Data Preparation

`Code/Analysis/merge_biomarkers_demographics.R`

Following completion of the primary data-processing workflow, the processed biomarker data are combined with the demographic and diagnosis data to create the final analysis-ready dataset.

## Statistical Analysis Files

The separate `Analysis/` directory contains the R scripts used to conduct the statistical analyses presented in the manuscript. These scripts are provided separately from the primary data-processing workflow to distinguish the EHR data-cleaning and preparation procedures from the subsequent statistical analyses.

## Reproducibility and Use with Other EHR Data

This repository is intended to document the major data-processing procedures used to prepare the EHR data for analysis and to improve transparency and reproducibility of the study methodology.

The workflow reflects the structure and characteristics of the EHR data used in this study. EHR datasets may differ across data sources and institutions with respect to data structure, variable definitions, coding practices, laboratory naming conventions, measurement units, completeness, and data quality. Therefore, although the general workflow may be applicable to other EHR data sources, individual processing steps and code may require adaptation before being applied to a different EHR dataset.

Because the underlying patient-level data cannot be publicly distributed, the repository does not provide a fully executable replication dataset. Instead, the code documents the processing logic used for cohort construction, data cleaning and standardization, biomarker processing, quality control, final data preparation, and statistical analysis.

## Citation

If using or adapting code from this repository, please cite the associated manuscript:

**Managing and Cleaning Large-Scale Electronic Health Record Data**

Full citation information will be added following publication.

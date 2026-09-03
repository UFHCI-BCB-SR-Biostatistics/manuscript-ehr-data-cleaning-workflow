# EHR Data Cleaning and Processing Workflow

## Overview

This repository contains R code developed to support the data processing and cleaning workflow described in the manuscript **“Managing and Cleaning Large-Scale Electronic Health Record Data.”**

The code documents a reproducible workflow for preparing large-scale electronic health record (EHR) data for analysis. The workflow includes cohort construction, BMI data processing, laboratory biomarker identification and standardization, longitudinal biomarker processing, monthly aggregation of biomarker measurements, and creation of the final analytic dataset.

All data processing and analysis code included in this repository was implemented in **R**.

## Data Availability

The patient-level EHR data used in this study are not included in this repository because the underlying data are restricted.

Access to the source data is subject to the requirements and approvals of the original data provider. Accordingly, this repository contains the code used to implement the data-processing workflow but does not contain patient-level data.

Local directory structures and source filenames have been replaced with generic placeholders in the public code. Users with appropriate access to comparable data should specify the corresponding local paths and filenames before running the scripts.

## Repository Structure

The code is organized into three main directories:

- `Code/Demographics and Diagnosis/`
- `Code/Laboratory Biomarkers/`
- `Code/Analysis/`

### Demographics and Diagnosis

#### `demographics_diagnosis_cleaning.R`

Constructs the study cohort using demographic and diagnosis data. The script identifies liver disease and hepatocellular carcinoma (HCC) diagnoses, determines the first qualifying diagnosis dates, applies cohort eligibility criteria, and performs quality-control checks.

#### `bmi_data_processing.R`

Processes BMI-related data after construction of the study cohort. Available height, weight, and BMI measurements are evaluated in relation to the first liver disease diagnosis date and used to derive BMI measurements for the analytic cohort.

BMI processing is performed after cohort construction so that BMI measurements can be evaluated for patients meeting the study cohort criteria and in relation to the relevant diagnosis dates.

### Laboratory Biomarkers

Laboratory biomarker processing is performed sequentially using three scripts.

#### `01_obtaining_biomarkers.R`

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

#### `02_condensing_biomarkers.R`

Restricts biomarker measurements to the study-specific observation periods and applies biomarker-specific value thresholds.

For patients with HCC, measurements are evaluated between the first liver disease diagnosis and first HCC diagnosis. For patients without HCC, measurements following the first liver disease diagnosis are evaluated.

When multiple measurements of the same biomarker are available within a calendar month, measurements are summarized to create monthly biomarker values.

This script also implements the lipid-component calculation used to derive a missing lipid value when the other required lipid components are available.

#### `03_merging_biomarkers.R`

Combines the processed monthly biomarker datasets and reshapes the data into a subject-month format, with biomarker measurements represented as separate variables for subsequent analysis.

### Analysis

#### `merge_biomarkers_demographics.R`

Combines the processed biomarker data with the demographic and diagnosis data to create the dataset used for subsequent analyses.

## Workflow

The general data-processing sequence is:

1. `Code/Demographics and Diagnosis/demographics_diagnosis_cleaning.R`
   - Construct the study cohort from demographic and diagnosis data.

2. `Code/Demographics and Diagnosis/bmi_data_processing.R`
   - Process BMI data for the constructed cohort.

3. `Code/Laboratory Biomarkers/01_obtaining_biomarkers.R`
   - Identify and standardize laboratory biomarker records.

4. `Code/Laboratory Biomarkers/02_condensing_biomarkers.R`
   - Apply study-period and value criteria and aggregate biomarker measurements by month.

5. `Code/Laboratory Biomarkers/03_merging_biomarkers.R`
   - Merge monthly biomarker measurements and reshape the biomarker data.

6. `Code/Analysis/merge_biomarkers_demographics.R`
   - Merge processed biomarker data with demographic and diagnosis data for analysis.

## Reproducibility and Use with Other EHR Data

This repository is intended to document the major data-processing procedures used to prepare the EHR data for analysis and to improve transparency and reproducibility of the study methodology.

The workflow reflects the structure and characteristics of the EHR data used in this study. EHR datasets may differ across data sources and institutions with respect to data structure, variable definitions, coding practices, laboratory naming conventions, measurement units, completeness, and data quality. Therefore, although the general workflow may be applicable to other EHR data sources, individual processing steps and code may require adaptation before being applied to a different EHR dataset.

Because the underlying patient-level data cannot be publicly distributed, the repository does not provide a fully executable replication dataset. Instead, the code documents the processing logic used for cohort construction, data cleaning and standardization, biomarker processing, quality control, and preparation of analysis-ready data.

## Citation

If using or adapting code from this repository, please cite the associated manuscript:

**Managing and Cleaning Large-Scale Electronic Health Record Data**

Full citation information will be added following publication.

# ==============================================================================
# OneFlorida+ Laboratory Biomarker Processing
# ==============================================================================

# ---- Load libraries -----------------------------------------------------------

library(tidyverse)
library(here)


# ---- Project paths ------------------------------------------------------------

# Patient-level EHR data are not included in the public repository.
# Replace the placeholder file name with the appropriate source file name.

raw_lab_file <- here(
  "Data",
  "Raw OneFL Data",
  "laboratory_data.csv"
)

# Output directory used for biomarker-specific datasets
raw_biomarker_dir <- here(
  "Data",
  "Biomarker Data",
  "UPDATED Raw Biomarker Data"
)

# Set to TRUE only when files should be written to disk
SAVE_OUTPUTS <- FALSE


# ---- Load laboratory data -----------------------------------------------------

lab <- read_csv(
  raw_lab_file,
  show_col_types = FALSE
)

# Remove automatically generated index column if present
if (names(lab)[1] %in% c("X", "...1")) {
  lab <- lab %>%
    select(-1)
}


# ==============================================================================
# Helper functions
# ==============================================================================

# Variables retained for each biomarker
lab_variables <- c(
  "Subject_ID",
  "RAW_LAB_NAME",
  "LAB_ORDER_DATE",
  "RESULT_DATE",
  "RESULT_NUM",
  "RESULT_UNIT"
)


# Extract a biomarker based on predefined laboratory names
get_biomarker <- function(data, biomarker_names) {
  data %>%
    filter(RAW_LAB_NAME %in% biomarker_names) %>%
    select(all_of(lab_variables))
}


# Summarize raw laboratory names
check_lab_names <- function(data) {
  data %>%
    count(RAW_LAB_NAME, sort = TRUE)
}


# Summarize measurement units
check_units <- function(data) {
  data %>%
    count(RESULT_UNIT, sort = TRUE)
}


# Calculate 99th percentile
get_99th_percentile <- function(data) {
  quantile(
    data$RESULT_NUM,
    probs = 0.99,
    na.rm = TRUE
  )
}


# Save biomarker-specific dataset
save_biomarker <- function(data, file_name) {
  
  if (SAVE_OUTPUTS) {
    
    dir.create(
      raw_biomarker_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write_csv(
      data,
      file.path(raw_biomarker_dir, file_name)
    )
  }
}


# ==============================================================================
# AST
# ==============================================================================

ast_names <- c(
  "AST",
  "AST (POC)",
  "AST (SGOT)",
  "AST (SGOT) P5P",
  "AST:",
  "Aspartate Aminotransferase (AST)",
  "ASPARTATE AMINOTRANSFERASE-PLASM"
)

ast1 <- get_biomarker(
  lab,
  ast_names
)

# Check original unit distribution
ast_original_unit_counts <- check_units(ast1)

# Standardize measurement unit
ast1 <- ast1 %>%
  mutate(
    RESULT_UNIT = "U/L"
  )

# Quality-control checks
ast_name_counts <- check_lab_names(ast1)
ast_unit_counts <- check_units(ast1)

summary(ast1$RESULT_NUM)
ast_99 <- get_99th_percentile(ast1)

save_biomarker(
  ast1,
  "All Raw AST Data.csv"
)


# ==============================================================================
# Albumin
# ==============================================================================

albumin_names <- c(
  "ALBUMIN",
  "ALBUMIN, SERUM",
  "ALBUMIN:",
  "Albumin",
  "Albumin Level",
  "Albumin Lvl",
  "Albumin, Serum",
  "Serum Albumin",
  "Albumin, HFPA",
  "\"Albumin, HFPA\"\"\"\"\""
)

alb1 <- get_biomarker(
  lab,
  albumin_names
)

# Examine additional names containing "Alb"
alblook <- lab %>%
  filter(
    str_detect(
      RAW_LAB_NAME,
      regex("Alb", ignore_case = TRUE)
    )
  ) %>%
  count(
    RAW_LAB_NAME,
    sort = TRUE
  )

# Check original units
alb_original_unit_counts <- check_units(alb1)

# Remove blank units and convert mg/dL to g/dL
alb1 <- alb1 %>%
  filter(
    !is.na(RESULT_UNIT),
    RESULT_UNIT != ""
  ) %>%
  mutate(
    RESULT_NUM = case_when(
      RESULT_UNIT == "mg/dL" ~ RESULT_NUM / 1000,
      TRUE ~ RESULT_NUM
    ),
    RESULT_UNIT = "g/dL"
  )

# Quality-control checks
alb_name_counts <- check_lab_names(alb1)
alb_unit_counts <- check_units(alb1)

summary(alb1$RESULT_NUM)
alb_99 <- get_99th_percentile(alb1)

save_biomarker(
  alb1,
  "All Raw Albumin Data.csv"
)


# ==============================================================================
# ALT
# ==============================================================================

alt_names <- c(
  "ALT",
  "ALT (POC)",
  "ALT (SGPT)",
  "ALT (SGPT) P5P",
  "ALT/SGPT",
  "ALT:"
)

alt1 <- get_biomarker(
  lab,
  alt_names
)

# Check original names and units
alt_original_name_counts <- check_lab_names(alt1)
alt_original_unit_counts <- check_units(alt1)

# Standardize unit
alt1 <- alt1 %>%
  mutate(
    RESULT_UNIT = "U/L"
  )

# Quality-control checks
alt_name_counts <- check_lab_names(alt1)
alt_unit_counts <- check_units(alt1)

summary(alt1$RESULT_NUM)
alt_99 <- get_99th_percentile(alt1)

save_biomarker(
  alt1,
  "All Raw ALT Data.csv"
)


# ==============================================================================
# Total Bilirubin
# ==============================================================================

bilirubin_names <- c(
  "BILIRUBIN",
  "BILIRUBIN  TOTAL",
  "BILIRUBIN TOTAL",
  "BILIRUBIN, TOTAL",
  "TOTAL BILI:",
  "TOTAL BILIRUBIN",
  "Total Bilirubin, HFPA",
  "Bilirubin",
  "-- Bilirubin",
  "Bili Total",
  "Bilirubin Total",
  "Bilirubin, Total",
  "Total Bilirubin",
  "\"Total Bilirubin, HFPA\"\"\"\"\""
)

bill1 <- get_biomarker(
  lab,
  bilirubin_names
)

# Examine additional names containing "Bili"
blook <- lab %>%
  filter(
    str_detect(
      RAW_LAB_NAME,
      regex("Bili", ignore_case = TRUE)
    )
  )

bill_name_counts <- check_lab_names(bill1)
bill_original_unit_counts <- check_units(bill1)

# Remove unusable unit labels and standardize retained measurements
bill1 <- bill1 %>%
  filter(
    !is.na(RESULT_UNIT),
    !RESULT_UNIT %in% c("UN", "")
  ) %>%
  mutate(
    RESULT_UNIT = "mg/dL"
  )

bill_unit_counts <- check_units(bill1)

summary(bill1$RESULT_NUM)
bill_99 <- get_99th_percentile(bill1)

save_biomarker(
  bill1,
  "All Raw Bilirubin Data.csv"
)


# ==============================================================================
# Total Cholesterol
# ==============================================================================

total_cholesterol_names <- c(
  "Cholesterol",
  "Cholesterol Total",
  "Cholesterol, Total",
  "CHOLESTEROL",
  "CHOLESTEROL  TOTAL",
  "CHOLESTEROL TOTAL",
  "NMR Total Cholesterol"
)

tchol1 <- get_biomarker(
  lab,
  total_cholesterol_names
)

# Check original names and units
tchol_name_counts <- check_lab_names(tchol1)
tchol_original_unit_counts <- check_units(tchol1)

# Standardize unit
tchol1 <- tchol1 %>%
  mutate(
    RESULT_UNIT = "mg/dL"
  )

tchol_unit_counts <- check_units(tchol1)

summary(tchol1$RESULT_NUM)
tchol_99 <- get_99th_percentile(tchol1)

save_biomarker(
  tchol1,
  "All Raw Total Cholesterol Data.csv"
)


# ==============================================================================
# LDL Cholesterol
# ==============================================================================

ldl_names <- c(
  "CALCULATED LDL",
  "Cholesterol LDL",
  "CHOLESTEROL, LDL",
  "DIRECT LDL",
  "LDL",
  "LDL-CHOLESTEROL",
  "LDL CALC",
  "LDL CHOL CALC (NIH)",
  "LDL Cholesterol",
  "LDL CHOLESTEROL",
  "LDL Cholesterol Calc",
  "LDL CHOLESTEROL CALC",
  "LDL CHOLESTEROL CALCULATED",
  "LDL Cholesterol Direct",
  "LDL Cholesterol Total",
  "LDL CHOLSTEROL",
  "LDL Direct",
  "LDL DIRECT",
  "Low Density Lipoprotein Cholesterol"
)

ldl1 <- get_biomarker(
  lab,
  ldl_names
)

# Check original names and units
ldl_name_counts <- check_lab_names(ldl1)
ldl_original_unit_counts <- check_units(ldl1)

# Standardize unit
ldl1 <- ldl1 %>%
  mutate(
    RESULT_UNIT = "mg/dL"
  )

ldl_unit_counts <- check_units(ldl1)

summary(ldl1$RESULT_NUM)
ldl_99 <- get_99th_percentile(ldl1)

save_biomarker(
  ldl1,
  "All Raw LDL Cholesterol Data.csv"
)


# ==============================================================================
# HDL Cholesterol
# ==============================================================================

hdl_names <- c(
  "CHOLESTEROL, HDL",
  "Cholesterol HDL",
  "HDL",
  "HDL CHOLESTEROL",
  "HDL Cholesterol",
  "HIGH DENSITY CHOLESTEROL"
)

hdl1 <- get_biomarker(
  lab,
  hdl_names
)

# Check original names and units
hdl_name_counts <- check_lab_names(hdl1)
hdl_original_unit_counts <- check_units(hdl1)

# Standardize unit
hdl1 <- hdl1 %>%
  mutate(
    RESULT_UNIT = "mg/dL"
  )

hdl_unit_counts <- check_units(hdl1)

summary(hdl1$RESULT_NUM)
hdl_99 <- get_99th_percentile(hdl1)

save_biomarker(
  hdl1,
  "All Raw HDL Cholesterol Data.csv"
)


# ==============================================================================
# Hemoglobin
# ==============================================================================

hemoglobin_names <- c(
  "HGB",
  "TOTAL HGB",
  "HEMOGLOBIN",
  "HEMOGLOBIN  ARTERIAL",
  "HEMOGLOBIN  MIXED VENOUS",
  "HEMOGLOBIN  VENOUS",
  "HEMOGLOBIN (I-STAT)",
  "HEMOGLOBIN (POINT OF CARE)",
  "HEMOGLOBIN BG",
  "HEMOGLOBIN-BG",
  "Hemoglobin (Calc), ISTAT POC",
  "POCT Hemoglobin calc",
  "ABG Hgb (Calc)",
  "Hgb (Calc)  POC",
  "Hgb POC IStat",
  "\"Hemoglobin (Calc), ISTAT POC\"\"\"\"\""
)

hemo1 <- get_biomarker(
  lab,
  hemoglobin_names
)

# Examine additional names containing "Hemo"
hlook <- lab %>%
  filter(
    str_detect(
      RAW_LAB_NAME,
      regex("Hemo", ignore_case = TRUE)
    )
  )

hemo_candidate_names <- hlook %>%
  count(
    RAW_LAB_NAME,
    sort = TRUE
  )

# Check original names and units
hemo_name_counts <- check_lab_names(hemo1)
hemo_original_unit_counts <- check_units(hemo1)

# Standardize unit
hemo1 <- hemo1 %>%
  mutate(
    RESULT_UNIT = "g/dL"
  )

hemo_unit_counts <- check_units(hemo1)

summary(hemo1$RESULT_NUM)
hemo_99 <- get_99th_percentile(hemo1)

save_biomarker(
  hemo1,
  "All Raw Hemoglobin Data.csv"
)


# ==============================================================================
# INR
# ==============================================================================

inr_names <- c(
  "INR",
  "INR POC",
  "INR POC IStat",
  "ISTAT INR",
  "LWB INR",
  "PROTHROMBIN TIME INR (POC)",
  "PT INR",
  "PT-INR",
  "Whole Blood INR"
)

inr1 <- get_biomarker(
  lab,
  inr_names
)

# Check original names and units
inr_name_counts <- check_lab_names(inr1)
inr_original_unit_counts <- check_units(inr1)

# Standardize unit
inr1 <- inr1 %>%
  mutate(
    RESULT_UNIT = "{ratio}"
  )

inr_unit_counts <- check_units(inr1)

summary(inr1$RESULT_NUM)
inr_99 <- get_99th_percentile(inr1)

save_biomarker(
  inr1,
  "All Raw INR Data.csv"
)


# ==============================================================================
# Platelets
# ==============================================================================

platelet_names <- c(
  "PLATELET  CITRATED BLOOD",
  "PLATELET COUNT",
  "PLATELET COUNT(SODIUM CITRATE)",
  "PLATELET ESTIMATE",
  "PLATELET ESTIMATION",
  "PLATELETS",
  "PLATELETS.",
  "Platelet",
  "Platelet Cnt",
  "Platelet Count",
  "Platelet Count (k/mm3)",
  "Platelets"
)

pl1 <- get_biomarker(
  lab,
  platelet_names
)

# Check original names and units
pl_name_counts <- check_lab_names(pl1)
pl_original_unit_counts <- check_units(pl1)

# Standardize unit
pl1 <- pl1 %>%
  mutate(
    RESULT_UNIT = "10*3/uL"
  )

pl_unit_counts <- check_units(pl1)

summary(pl1$RESULT_NUM)
pl_99 <- get_99th_percentile(pl1)

save_biomarker(
  pl1,
  "All Raw Platelets Data.csv"
)


# ==============================================================================
# Triglycerides
# ==============================================================================

triglyceride_names <- c(
  "Triglycerides"
)

trig1 <- get_biomarker(
  lab,
  triglyceride_names
)

# Check original names and units
trig_name_counts <- check_lab_names(trig1)
trig_original_unit_counts <- check_units(trig1)

# Standardize unit
trig1 <- trig1 %>%
  mutate(
    RESULT_UNIT = "mg/dL"
  )

trig_unit_counts <- check_units(trig1)

summary(trig1$RESULT_NUM)
trig_99 <- get_99th_percentile(trig1)

save_biomarker(
  trig1,
  "All Raw Triglycerides Data.csv"
)


# ==============================================================================
# Alpha-Fetoprotein
# ==============================================================================

afp_names <- c(
  "ALPHA FETO PROT",
  "ALPHA FETOPROTEIN,$TUMOR MARKER",
  "ALPHA-FETOPROTEIN",
  "Alpha Fetoprotein"
)

alp1 <- get_biomarker(
  lab,
  afp_names
)

# Check original names and units
afp_name_counts <- check_lab_names(alp1)
afp_original_unit_counts <- check_units(alp1)

# Standardize unit
alp1 <- alp1 %>%
  mutate(
    RESULT_UNIT = "ng/mL"
  )

afp_unit_counts <- check_units(alp1)

summary(alp1$RESULT_NUM)
afp_99 <- get_99th_percentile(alp1)

save_biomarker(
  alp1,
  "All Raw Alpha-Feto Protein Data.csv"
)


# ==============================================================================
# Neutrophils
# ==============================================================================

neutrophil_names <- c(
  "ABSOLUTE NEUT COUNT",
  "ABSOLUTE NEUTROPHIL COUNT",
  "ABSOLUTE NEUTROPHILS",
  "Abs Neutrophil Cnt",
  "Abs Neuts Manual",
  "Absolute Count Neutrophils",
  "NEUT ABSOLUTE",
  "NEUTROPHIL ABS MAN",
  "NEUTROPHILS",
  "NEUTROPHILS (ABSOLUTE)",
  "NEUTROPHILS ABSOLUTE COUNT",
  "NEUTROS ABS",
  "Neut Abs",
  "Neut Abs Man",
  "Neutrophil Ab",
  "Neutrophils",
  "Neutrophils (Absolute)",
  "Neuts Absolute",
  "Neuts Auto"
)

neut1 <- get_biomarker(
  lab,
  neutrophil_names
)

# Check original names and units
neut_name_counts <- check_lab_names(neut1)
neut_original_unit_counts <- check_units(neut1)

# Retain absolute-count measurements only.
# Percentage and unknown-unit records are excluded.
# Counts reported per uL are converted to 10^3/uL.
neut1 <- neut1 %>%
  filter(
    !is.na(RESULT_UNIT),
    !RESULT_UNIT %in% c("%", "UN")
  ) %>%
  mutate(
    RESULT_NUM = case_when(
      RESULT_UNIT %in% c("/uL", "{cells}/uL") ~ RESULT_NUM / 1000,
      TRUE ~ RESULT_NUM
    ),
    RESULT_UNIT = "10*3/uL"
  )

neut_unit_counts <- check_units(neut1)

summary(neut1$RESULT_NUM)
neut_99 <- get_99th_percentile(neut1)

save_biomarker(
  neut1,
  "All Raw Neutrophils Data.csv"
)


# ==============================================================================
# Lymphocytes
# ==============================================================================

lymphocyte_names <- c(
  "Abs Lymph Cnt",
  "Abs Lymphocyte Cnt",
  "Abs Lymphs Manual",
  "Absolute Count Lymphocytes",
  "LYMPH ABSOLUTE",
  "LYMPHOCYTES",
  "LYMPHOCYTES ABSOLUTE COUNT",
  "LYMPHS",
  "LYMPHS (ABSOLUTE)",
  "LYMPHS ABS",
  "LYMPHS ABS MAN",
  "Lymph Absolute",
  "Lymphocytes",
  "Lymphocytes (Absolute)",
  "Lymphs",
  "Lymphs (Absolute)",
  "Lymphs Absolute",
  "Lymphs Auto",
  "Lymphs Man",
  "Manual Count Lymphocyte"
)

lymph1 <- get_biomarker(
  lab,
  lymphocyte_names
)

# Check original names and units
lymph_name_counts <- check_lab_names(lymph1)
lymph_original_unit_counts <- check_units(lymph1)

# Retain absolute-count measurements only.
# Percentage and unknown-unit records are excluded.
# Counts reported as cells/uL are converted to 10^3/uL.
lymph1 <- lymph1 %>%
  filter(
    !is.na(RESULT_UNIT),
    !RESULT_UNIT %in% c("%", "UN")
  ) %>%
  mutate(
    RESULT_NUM = case_when(
      RESULT_UNIT == "{cells}/uL" ~ RESULT_NUM / 1000,
      TRUE ~ RESULT_NUM
    ),
    RESULT_UNIT = "10*3/uL"
  )

lymph_unit_counts <- check_units(lymph1)

summary(lymph1$RESULT_NUM)
lymph_99 <- get_99th_percentile(lymph1)

save_biomarker(
  lymph1,
  "All Raw Lymphocytes Data.csv"
)


# ==============================================================================
# Summary of 99th Percentiles
# ==============================================================================

biomarker_99th_percentiles <- tibble(
  biomarker = c(
    "AST",
    "Albumin",
    "ALT",
    "Total Bilirubin",
    "Total Cholesterol",
    "LDL Cholesterol",
    "HDL Cholesterol",
    "Hemoglobin",
    "INR",
    "Platelets",
    "Triglycerides",
    "Alpha-Fetoprotein",
    "Neutrophils",
    "Lymphocytes"
  ),
  percentile_99 = c(
    ast_99,
    alb_99,
    alt_99,
    bill_99,
    tchol_99,
    ldl_99,
    hdl_99,
    hemo_99,
    inr_99,
    pl_99,
    trig_99,
    afp_99,
    neut_99,
    lymph_99
  )
)

biomarker_99th_percentiles
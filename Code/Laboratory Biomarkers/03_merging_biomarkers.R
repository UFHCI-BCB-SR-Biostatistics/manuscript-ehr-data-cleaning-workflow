# ==============================================================================
# Merge Monthly Biomarker Files
# ==============================================================================

# ---- Load libraries -----------------------------------------------------------

library(tidyverse)
library(here)


# ---- Project paths ------------------------------------------------------------

# Restricted patient-level data are not included in this repository.
# Specify the appropriate local input and output paths before running this script.

mean_biomarker_dir <- here(
  "SPECIFY_PATH_TO_MONTHLY_BIOMARKER_DATA"
)

merged_biomarker_dir <- here(
  "SPECIFY_OUTPUT_PATH_TO_MERGED_BIOMARKER_DATA"
)

dir.create(
  merged_biomarker_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Identify biomarker files -------------------------------------------------

# Locate monthly biomarker CSV files.

files <- list.files(
  path = mean_biomarker_dir,
  pattern = "_monthly_data\\.csv$",
  full.names = TRUE
)


# ---- Merge biomarker files ----------------------------------------------------

combined_biomarkers <- files %>%
  lapply(
    read_csv,
    show_col_types = FALSE
  ) %>%
  bind_rows() %>%
  dplyr::select(-1) %>%
  arrange(
    Subject_ID,
    RESULT_MONTH
  )


# ---- Create subject-month identifier -----------------------------------------

combined_biomarkers <- combined_biomarkers %>%
  group_by(
    Subject_ID,
    RESULT_MONTH
  ) %>%
  mutate(
    row_id = cur_group_id()
  ) %>%
  ungroup()


# ---- Reshape biomarker data to wide format -----------------------------------

biomarkers_wide <- combined_biomarkers %>%
  dplyr::select(
    -RESULT_DATE
  ) %>%
  pivot_wider(
    id_cols = c(
      Subject_ID,
      row_id,
      RESULT_MONTH
    ),
    names_from = LAB_TYPE,
    values_from = c(
      RESULT_UNIT,
      mean_RESULT_NUM
    ),
    names_glue = "{ifelse(
      .value == 'mean_RESULT_NUM',
      paste0('mean_', LAB_TYPE, '_RESULT_NUM'),
      paste0(LAB_TYPE, '_', .value)
    )}"
  ) %>%
  arrange(
    Subject_ID,
    RESULT_MONTH
  ) %>%
  dplyr::select(
    -row_id
  )


# ---- Check for duplicate subject-month records -------------------------------

duplicate_subject_month <- biomarkers_wide %>%
  count(
    Subject_ID,
    RESULT_MONTH
  ) %>%
  filter(
    n > 1
  )


# ---- Save merged biomarker dataset -------------------------------------------

write_csv(
  biomarkers_wide,
  file.path(
    merged_biomarker_dir,
    "merged_biomarker_data.csv"
  )
)
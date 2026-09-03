# ==============================================================================
# Merge Biomarker Data with Demographic Variables
# ==============================================================================

# ---- Load libraries -----------------------------------------------------------

library(tidyverse)
library(here)


# ---- Project paths ------------------------------------------------------------

# Restricted patient-level data are not included in this repository.
# Specify the appropriate local paths and file names before running this script.

biomarker_file <- here(
  "SPECIFY_PATH",
  "biomarker_data.csv"
)

demographic_file <- here(
  "SPECIFY_PATH",
  "demographic_data.csv"
)

output_dir <- here(
  "SPECIFY_OUTPUT_PATH"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---- Load data ----------------------------------------------------------------

biomarkers <- read_csv(
  biomarker_file,
  show_col_types = FALSE
)

demographics <- read_csv(
  demographic_file,
  show_col_types = FALSE
)


# ---- Select demographic variables --------------------------------------------

demo_vars <- demographics %>%
  dplyr::select(
    Subject_ID,
    within_period,
    First_HCC_MONTH,
    First_HCC_YEAR
  )


# ---- Merge demographic variables with biomarker data -------------------------

final_data <- biomarkers %>%
  left_join(
    demo_vars,
    by = "Subject_ID"
  )


# ---- Prepare date variables ---------------------------------------------------

final_data <- final_data %>%
  mutate(
    RESULT_MONTH = as.Date(RESULT_MONTH)
  ) %>%
  arrange(
    Subject_ID,
    RESULT_MONTH
  ) %>%
  relocate(
    RESULT_MONTH,
    .after = First_HCC_YEAR
  )


# ---- Check for duplicate subject-month records -------------------------------

duplicate_subject_month <- final_data %>%
  count(
    Subject_ID,
    RESULT_MONTH
  ) %>%
  filter(
    n > 1
  )


# ---- Check number of unique subjects -----------------------------------------

n_distinct(
  final_data$Subject_ID
)


# ---- Save final analytic dataset ---------------------------------------------

write_csv(
  final_data,
  file.path(
    output_dir,
    "final_analysis_data.csv"
  )
)
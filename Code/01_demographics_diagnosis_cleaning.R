# ==============================================================================
# HCC / Liver Disease Cohort Construction
# UF Health Cancer Institute Biostatistics
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# 1. Load source data
# ------------------------------------------------------------------------------

demo <- read_csv(
  here("Data", "Raw OneFL Data", "CaseCon_demographic.csv"),
  show_col_types = FALSE
)

diagnosis <- read_csv(
  here("Data", "Raw OneFL Data", "CaseCon_diagnosis_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    ADMIT_DATE = as.Date(ADMIT_DATE),
    DX_DATE = as.Date(DX_DATE)
  )

# ------------------------------------------------------------------------------
# 2. Clean demographic data
# ------------------------------------------------------------------------------

demo_clean <- demo %>%
  mutate(
    HISPANIC = na_if(HISPANIC, "R"),
    HISPANIC = na_if(HISPANIC, "NI"),
    HISPANIC = na_if(HISPANIC, "UN"),
    HISPANIC = na_if(HISPANIC, "OT"),
    RACE = case_when(
      RACE == "01" ~ "Other",
      RACE == "02" ~ "Asian",
      RACE == "03" ~ "Black or African American",
      RACE == "04" ~ "Other",
      RACE == "05" ~ "White",
      RACE == "06" ~ "Other",
      RACE == "OT" ~ "Other",
      RACE %in% c("07", "NI", "UN") ~ NA_character_,
      TRUE ~ as.character(RACE)
    )
  ) %>%
  drop_na(RACE, HISPANIC) %>%
  filter(AGE != 200)

# ------------------------------------------------------------------------------
# 3. Standardize diagnosis codes
# ------------------------------------------------------------------------------

# Named vector: names = original undotted codes; values = standardized codes.
dx_corrections <- c(
  # HCC
  "C220" = "C22.0",
  "1550" = "155.0",

  # Alcohol-related liver disease
  "K700" = "K70.0",
  "K7010" = "K70.10",
  "K7011" = "K70.11",
  "K702" = "K70.2",
  "K7030" = "K70.30",
  "K7031" = "K70.31",
  "K7040" = "K70.40",
  "K7041" = "K70.41",
  "K709" = "K70.9",

  # Cirrhosis, fibrosis, fatty liver, and other liver disease
  "5716" = "571.6",
  "K743" = "K74.3",
  "K744" = "K74.4",
  "K745" = "K74.5",
  "5713" = "571.3",
  "5718" = "571.8",
  "5719" = "571.9",
  "K760" = "K76.0",
  "K7581" = "K75.81",
  "K740" = "K74.0",
  "K741" = "K74.1",
  "K742" = "K74.2",
  "K7460" = "K74.60",
  "K7469" = "K74.69",
  "K7689" = "K76.89",
  "K769" = "K76.9",
  "5710" = "571.0",
  "5711" = "571.1",
  "5712" = "571.2",
  "5715" = "571.5",

  # Hepatitis B
  "07020" = "070.20",
  "07022" = "070.22",
  "07023" = "070.23",
  "07030" = "070.30",
  "07032" = "070.32",
  "07033" = "070.33",
  "B180" = "B18.0",
  "B181" = "B18.1",
  "B1910" = "B19.10",
  "B1911" = "B19.11",

  # Hepatitis C
  "07044" = "070.44",
  "07054" = "070.54",
  "07070" = "070.70",
  "07071" = "070.71",
  "B182" = "B18.2",
  "B1920" = "B19.20",
  "B1921" = "B19.21"
)

diagnosis_clean <- diagnosis %>%
  mutate(
    DX_original = as.character(DX),
    DX_was_fixed = DX_original %in% names(dx_corrections),
    DX_standardized = coalesce(
      unname(dx_corrections[DX_original]),
      DX_original
    )
  )

# ------------------------------------------------------------------------------
# 4. Define HCC and liver disease code sets
# ------------------------------------------------------------------------------

hcc_codes <- c("C22.0", "155.0")

liver_codes <- c(
  # Alcohol-related liver disease
  "K70.0", "K70.10", "K70.11", "K70.2", "K70.30", "K70.31",
  "K70.40", "K70.41", "K70.9",

  # Cirrhosis, fibrosis, fatty liver, and other liver disease
  "571.6", "K74.3", "K74.4", "K74.5", "571.3", "571.8", "571.9",
  "K76.0", "K75.81", "K74.0", "K74.1", "K74.2", "K74.60", "K74.69",
  "K76.89", "K76.9", "571.0", "571.1", "571.2", "571.5",

  # Hepatitis B
  "070.20", "070.22", "070.23", "070.30", "070.32", "070.33",
  "B18.0", "B18.1", "B19.10", "B19.11",

  # Hepatitis C
  "070.44", "070.54", "070.70", "070.71", "B18.2", "B19.20", "B19.21"
)

# ------------------------------------------------------------------------------
# 5. Identify first diagnosis per patient
# ------------------------------------------------------------------------------

get_first_diagnosis <- function(data, codes, date_name, dx_name, prefix) {
  data %>%
    filter(
      DX_standardized %in% codes,
      !is.na(ADMIT_DATE)
    ) %>%
    arrange(Subject_ID, ADMIT_DATE, DX_standardized, DX_original) %>%
    group_by(Subject_ID) %>%
    slice_head(n = 1) %>%
    transmute(
      Subject_ID,
      "{date_name}" := ADMIT_DATE,
      "{dx_name}" := DX_standardized,
      "{prefix}_DX_original" := DX_original,
      "{prefix}_DX_fixed" := DX_was_fixed
    ) %>%
    ungroup()
}

first_hcc <- get_first_diagnosis(
  data = diagnosis_clean,
  codes = hcc_codes,
  date_name = "First_HCC_DX_DATE",
  dx_name = "HCC_DX",
  prefix = "First_HCC"
)

first_ld <- get_first_diagnosis(
  data = diagnosis_clean,
  codes = liver_codes,
  date_name = "First_LD_DATE",
  dx_name = "First_LD_DX",
  prefix = "First_LD"
)

# ------------------------------------------------------------------------------
# 6. Construct eligible liver disease cohort
# ------------------------------------------------------------------------------

analytic_cohort <- demo %>%
  semi_join(demo_clean, by = "Subject_ID") %>%
  left_join(first_ld, by = "Subject_ID") %>%
  left_join(first_hcc, by = "Subject_ID") %>%
  filter(!is.na(First_LD_DATE)) %>%
  mutate(
    HCC_group = if_else(is.na(First_HCC_DX_DATE), "Non-HCC", "HCC"),
    LD_YEAR = lubridate::year(First_LD_DATE),
    HCC_YEAR = lubridate::year(First_HCC_DX_DATE)
  ) %>%
  filter(
    !LD_YEAR %in% c(2013, 2014),
    HCC_group == "Non-HCC" | First_HCC_DX_DATE >= First_LD_DATE
  ) %>%
  mutate(
    First_LD_or_HCC_DX_fixed = coalesce(First_LD_DX_fixed, FALSE) |
      coalesce(First_HCC_DX_fixed, FALSE),
    Diagnosis_fix_category = case_when(
      coalesce(First_LD_DX_fixed, FALSE) & coalesce(First_HCC_DX_fixed, FALSE) ~
        "First LD and first HCC codes fixed",
      coalesce(First_LD_DX_fixed, FALSE) ~ "First LD code fixed only",
      coalesce(First_HCC_DX_fixed, FALSE) ~ "First HCC code fixed only",
      TRUE ~ "Neither first diagnosis code fixed"
    )
  )

# ------------------------------------------------------------------------------
# 7. Quality-control summaries
# ------------------------------------------------------------------------------

# Cohort counts
cohort_summary <- analytic_cohort %>%
  summarise(
    total_patients = n_distinct(Subject_ID),
    hcc_patients = n_distinct(Subject_ID[HCC_group == "HCC"]),
    non_hcc_patients = n_distinct(Subject_ID[HCC_group == "Non-HCC"]),
    first_ld_code_fixed = n_distinct(Subject_ID[coalesce(First_LD_DX_fixed, FALSE)]),
    first_hcc_code_fixed = n_distinct(Subject_ID[coalesce(First_HCC_DX_fixed, FALSE)]),
    either_first_code_fixed = n_distinct(Subject_ID[First_LD_or_HCC_DX_fixed])
  )

# Mutually exclusive diagnosis-code correction categories
dx_fix_summary <- analytic_cohort %>%
  distinct(Subject_ID, Diagnosis_fix_category) %>%
  count(Diagnosis_fix_category, name = "n_patients") %>%
  mutate(percent = 100 * n_patients / sum(n_patients)) %>%
  arrange(desc(n_patients))

# Diagnosis record counts per patient
dx_record_counts <- analytic_cohort %>%
  distinct(Subject_ID) %>%
  left_join(
    diagnosis_clean %>%
      filter(DX_standardized %in% liver_codes) %>%
      count(Subject_ID, name = "n_ld_records"),
    by = "Subject_ID"
  ) %>%
  left_join(
    diagnosis_clean %>%
      filter(DX_standardized %in% hcc_codes) %>%
      count(Subject_ID, name = "n_hcc_records"),
    by = "Subject_ID"
  ) %>%
  mutate(
    n_ld_records = replace_na(n_ld_records, 0L),
    n_hcc_records = replace_na(n_hcc_records, 0L),
    LD_status = case_when(
      n_ld_records == 0 ~ "No LD",
      n_ld_records == 1 ~ "Single LD",
      TRUE ~ "Multiple LD"
    ),
    HCC_status = case_when(
      n_hcc_records == 0 ~ "No HCC",
      n_hcc_records == 1 ~ "Single HCC",
      TRUE ~ "Multiple HCC"
    )
  )

record_count_summary <- dx_record_counts %>%
  count(LD_status, HCC_status) %>%
  arrange(LD_status, HCC_status)

# Optional review table for patients whose first LD or HCC code was corrected
patients_with_fixed_first_dx <- analytic_cohort %>%
  filter(First_LD_or_HCC_DX_fixed) %>%
  select(
    Subject_ID,
    HCC_group,
    First_LD_DATE,
    First_LD_DX_original,
    First_LD_DX,
    First_LD_DX_fixed,
    First_HCC_DX_DATE,
    First_HCC_DX_original,
    HCC_DX,
    First_HCC_DX_fixed,
    Diagnosis_fix_category
  ) %>%
  arrange(HCC_group, Subject_ID)

# Print key QC outputs
cohort_summary
dx_fix_summary
record_count_summary

# ------------------------------------------------------------------------------
# 8. Optional export
# ------------------------------------------------------------------------------

# write_csv(
#   analytic_cohort,
#   here("Data", "Processed", "onefl_hcc_liver_disease_cohort.csv")
# )

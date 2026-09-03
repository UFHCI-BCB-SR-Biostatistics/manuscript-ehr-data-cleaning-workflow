library(tidyverse)
library(here)
library(hms)

#------------------------------------------------------------
# 1. Load data
#------------------------------------------------------------

# Restricted patient-level data are not included in this repository.
# Specify the appropriate local paths and file names before running this script.

vital <- read_csv(
  here(
    "SPECIFY_PATH",
    "vital_signs_data.csv"
  ),
  show_col_types = FALSE
)

final_data <- read_csv(
  here(
    "SPECIFY_PATH",
    "analytic_data.csv"
  ),
  show_col_types = FALSE
) %>%
  transmute(
    Subject_ID,
    First_LD_DATE = as.Date(First_LD_DATE),
    HCC_group
  ) %>%
  distinct()


#------------------------------------------------------------
# 2. Restrict vital records to the analytic cohort
#------------------------------------------------------------

bmi_subset <- vital %>%
  semi_join(
    final_data,
    by = "Subject_ID"
  ) %>%
  transmute(
    Subject_ID,
    MEASURE_DATE = as.Date(MEASURE_DATE),
    MEASURE_TIME = na_if(str_trim(MEASURE_TIME), ""),
    HT,
    WT,
    ORIGINAL_BMI
  ) %>%
  left_join(
    final_data,
    by = "Subject_ID"
  ) %>%
  mutate(
    MEASURE_TIME_SORT = suppressWarnings(
      parse_hm(MEASURE_TIME)
    )
  ) %>%
  arrange(
    Subject_ID,
    MEASURE_DATE,
    MEASURE_TIME_SORT
  )


# Optional export of restricted analytic-sample vital records

# write_csv(
#   bmi_subset,
#   here(
#     "SPECIFY_OUTPUT_PATH",
#     "bmi_analytic_sample.csv"
#   )
# )


#------------------------------------------------------------
# 3. Identify the measurement date closest to first LD diagnosis
#------------------------------------------------------------

closest_first <- bmi_subset %>%
  filter(
    !is.na(MEASURE_DATE),
    !is.na(First_LD_DATE),
    MEASURE_DATE <= First_LD_DATE
  ) %>%
  mutate(
    days_before_ld = as.integer(
      First_LD_DATE - MEASURE_DATE
    )
  ) %>%
  group_by(
    Subject_ID
  ) %>%
  filter(
    days_before_ld == min(days_before_ld)
  ) %>%
  ungroup()


#------------------------------------------------------------
# 4. Select BMI, height, and weight from the closest date
#------------------------------------------------------------

closest_bmi_ht_wt <- closest_first %>%
  group_by(
    Subject_ID,
    First_LD_DATE,
    MEASURE_DATE,
    MEASURE_TIME
  ) %>%
  summarise(
    HT_final = first(
      HT[!is.na(HT) & HT > 0],
      default = NA_real_
    ),
    
    WT_final = first(
      WT[!is.na(WT) & WT > 0],
      default = NA_real_
    ),
    
    BMI_original = first(
      ORIGINAL_BMI[!is.na(ORIGINAL_BMI)],
      default = NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  
  # Calculate BMI when valid height and weight are available
  mutate(
    BMI_calculated = if_else(
      !is.na(HT_final) &
        !is.na(WT_final) &
        HT_final > 0 &
        WT_final > 0,
      round(
        (WT_final * 703) / HT_final^2,
        2
      ),
      NA_real_
    ),
    
    original_bmi_valid =
      !is.na(BMI_original) &
      BMI_original > 10 &
      BMI_original < 100,
    
    calculated_bmi_valid =
      !is.na(BMI_calculated) &
      BMI_calculated > 10 &
      BMI_calculated < 100,
    
    has_all_three =
      !is.na(HT_final) &
      !is.na(WT_final) &
      !is.na(BMI_original),
    
    has_original_bmi =
      !is.na(BMI_original),
    
    has_ht_wt =
      !is.na(HT_final) &
      !is.na(WT_final)
  ) %>%
  
  # When multiple records occur on the selected date, prioritize:
  # 1. Records containing BMI, height, and weight
  # 2. Records containing original BMI
  # 3. Records containing height and weight
  # 4. Earliest available measurement time
  arrange(
    Subject_ID,
    MEASURE_DATE,
    desc(has_all_three),
    desc(has_original_bmi),
    desc(has_ht_wt),
    desc(!is.na(MEASURE_TIME)),
    MEASURE_TIME
  ) %>%
  
  group_by(
    Subject_ID,
    First_LD_DATE,
    MEASURE_DATE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  
  # Prefer a valid original BMI.
  # If unavailable or implausible, use calculated BMI.
  mutate(
    BMI_final = case_when(
      original_bmi_valid ~ BMI_original,
      calculated_bmi_valid ~ BMI_calculated,
      TRUE ~ NA_real_
    ),
    
    BMI_source = case_when(
      original_bmi_valid ~
        "Original BMI",
      
      is.na(BMI_original) &
        calculated_bmi_valid ~
        "Calculated BMI - original missing",
      
      !is.na(BMI_original) &
        !original_bmi_valid &
        calculated_bmi_valid ~
        "Calculated BMI - original implausible",
      
      TRUE ~
        "BMI unavailable"
    )
  ) %>%
  
  select(
    Subject_ID,
    First_LD_DATE,
    MEASURE_DATE,
    BMI_final
  )


#------------------------------------------------------------
# 5. Merge BMI with the analytic cohort
#------------------------------------------------------------

bmi_final_data <- final_data %>%
  left_join(
    closest_bmi_ht_wt,
    by = c(
      "Subject_ID",
      "First_LD_DATE"
    )
  )


#------------------------------------------------------------
# 6. Save cleaned BMI data
#------------------------------------------------------------

write_rds(
  bmi_final_data,
  here(
    "SPECIFY_OUTPUT_PATH",
    "bmi_analysis_data.rds"
  )
)

write_csv(
  bmi_final_data,
  here(
    "SPECIFY_OUTPUT_PATH",
    "bmi_analysis_data.csv"
  )
)


#------------------------------------------------------------
# 7. Quality control: patients without an eligible pre-LD BMI
#------------------------------------------------------------

patients_missing_bmi <- final_data %>%
  select(
    Subject_ID,
    First_LD_DATE
  ) %>%
  anti_join(
    closest_bmi_ht_wt %>%
      distinct(Subject_ID),
    by = "Subject_ID"
  )


# Review vital records for patients without a pre-LD measurement

patients_missing_bmi_records <- bmi_subset %>%
  semi_join(
    patients_missing_bmi,
    by = "Subject_ID"
  ) %>%
  mutate(
    measurement_timing = case_when(
      is.na(MEASURE_DATE) |
        is.na(First_LD_DATE) ~
        "Missing date",
      
      MEASURE_DATE > First_LD_DATE ~
        "After first LD diagnosis",
      
      MEASURE_DATE <= First_LD_DATE ~
        "On or before first LD diagnosis"
    )
  ) %>%
  select(
    Subject_ID,
    First_LD_DATE,
    MEASURE_DATE,
    measurement_timing
  ) %>%
  arrange(
    Subject_ID,
    MEASURE_DATE
  )


# Summarize why patients did not receive a baseline BMI

missing_bmi_summary <- patients_missing_bmi_records %>%
  distinct(
    Subject_ID,
    measurement_timing
  ) %>%
  count(
    measurement_timing,
    name = "n_patients"
  )

missing_bmi_summary
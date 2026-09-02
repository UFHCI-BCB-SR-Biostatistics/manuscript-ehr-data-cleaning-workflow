# ==============================================================================
# OneFlorida+ Biomarker Condensing and Monthly Aggregation
# ==============================================================================

# Load libraries
library(tidyverse)
library(here)
library(lubridate)
library(rlang)

# Project paths
# Restricted patient-level data are not included in the public repository.
# Replace placeholder input file names as needed for the local data source.

raw_biomarker_dir <- here(
  "Data",
  "Biomarker Data",
  "UPDATED Raw Biomarker Data"
)

mean_biomarker_dir <- here(
  "Data",
  "Biomarker Data",
  "UPDATED Mean Biomarker Data"
)

exclusion_dir <- here(
  "Data",
  "Exclusion Summary Data"
)

demographic_file <- here(
  "Data",
  "Clean Final Analysis Data",
  "demographics_diagnosis.csv"
)

# Create output folders if they do not already exist
dir.create(mean_biomarker_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exclusion_dir, recursive = TRUE, showWarnings = FALSE)

#### Condensing Lab Results Function ####
process <- function(lab_df,
                    result_col,
                    result_date_col,
                    result_threshold,
                    lab_shortname,
                    dem) {
  
  result_col      <- ensym(result_col)
  result_date_col <- ensym(result_date_col)
  
  result_month_col <- paste0(lab_shortname, "_RESULT_MONTH")
  mean_result_col  <- paste0("mean_", lab_shortname, "_RESULT_NUM")
  
  lab_df0 <- lab_df %>%
    filter(Subject_ID %in% dem$Subject_ID) %>%
    right_join(dem, by = "Subject_ID") %>%
    mutate(
      lab_flag = 1,   # <-- all rows kept; no lab type filtering
      !!result_date_col := as.Date(!!result_date_col),
      First_LD_DATE     = as.Date(First_LD_DATE),
      First_HCC_DX_DATE    = as.Date(First_HCC_DX_DATE)
    )
  
  # Keep only allowed lab names for both inclusion & exclusion accounting
  lab_df1 <- lab_df0 %>% filter(lab_flag == 1)
  
  # ---- Build inclusion flags per row (group-specific) ----
  # HCC: within LD..HCC (inclusive)
  hcc_flagged <- lab_df1 %>%
    filter(HCC_group == "HCC") %>%
    mutate(
      within_period = !is.na(First_LD_DATE) & !is.na(First_HCC_DX_DATE) &
        (!!result_date_col >= First_LD_DATE) & (!!result_date_col <= First_HCC_DX_DATE),
      below_cutoff  = (!!result_col < result_threshold) | is.na(!!result_col),
      include_row   = within_period & below_cutoff,
      period_type   = "LD_to_HCC"
    )
  
  # Non-HCC: date >= LD
  nonhcc_flagged <- lab_df1 %>%
    filter(HCC_group == "Non-HCC") %>%
    mutate(
      within_period = !is.na(First_LD_DATE) & (!!result_date_col >= First_LD_DATE),
      below_cutoff  = (!!result_col < result_threshold) | is.na(!!result_col),
      include_row   = within_period & below_cutoff,
      period_type   = "LD_to_LastObs"
    )
  
  flagged_all <- bind_rows(hcc_flagged, nonhcc_flagged)
  
  # ---- Exclusions log (based on flags above) ----
  exclusions <- flagged_all %>%
    filter(!include_row) %>%
    mutate(
      First_LD_DATE  = as.Date(First_LD_DATE),
      First_HCC_DX_DATE = as.Date(First_HCC_DX_DATE),
      study_period = dplyr::if_else(
        HCC_group == "HCC",
        paste0(as.character(First_LD_DATE), " — ", as.character(First_HCC_DX_DATE)),
        paste0(as.character(First_LD_DATE), " — NA")
      )
    ) %>%
    rowwise() %>%
    mutate(
      reason = paste(
        c(
          if (!below_cutoff)  "Out of Range Lab Value"  else NULL,
          if (!within_period) "Not Within Study Period" else NULL
        ),
        collapse = "; "
      )
    ) %>%
    ungroup() %>%
    transmute(
      Subject_ID,
      HCC_group,
      ld_date  = First_LD_DATE,        # <-- add LD date
      hcc_date = First_HCC_DX_DATE,       # <-- add HCC date (NA for non-HCC)
      study_period,
      result_date = !!result_date_col,
      value       = !!result_col,
      threshold   = result_threshold,
      period_type,
      reason
    ) %>%
    arrange(Subject_ID, result_date)
  
  # ---- Subject-level exclusions: patient in dem but never appears in lab_df ----
  subjects_missing_in_lab <- dem %>%
    dplyr::anti_join(lab_df %>% dplyr::distinct(Subject_ID), by = "Subject_ID") %>%
    dplyr::mutate(
      First_LD_DATE  = as.Date(First_LD_DATE),
      First_HCC_DX_DATE = as.Date(First_HCC_DX_DATE)
    ) %>%
    dplyr::transmute(
      Subject_ID,
      HCC_group,
      result_date = as.Date(NA),
      value       = NA_real_,
      threshold   = result_threshold,
      period_type = "N/A",
      reason      = "No Lab Records in Source Data"
    )
  
  # Append to exclusions (avoid duplicates if already present for any reason)
  exclusions <- exclusions %>%
    dplyr::bind_rows(
      dplyr::anti_join(subjects_missing_in_lab, dplyr::select(., Subject_ID), by = "Subject_ID")
    ) %>%
    dplyr::arrange(Subject_ID, result_date)
  
  # ---- Included rows only (for analysis outputs) ----
  included_all <- flagged_all %>% filter(include_row)
  
  # ----------------- HCC -----------------
  hcc_raw <- included_all %>%
    filter(HCC_group == "HCC") %>%
    arrange(Subject_ID, !!result_date_col) %>%
    mutate(
      !!result_month_col := floor_date(!!result_date_col, "month"),
      First_HCC_MONTH    = floor_date(First_HCC_DX_DATE, "month"),
      First_HCC_YEAR     = year(First_HCC_DX_DATE)
    )
  
  hcc_raw1 <- hcc_raw %>%
    group_by(Subject_ID, .data[[result_month_col]]) %>%
    mutate(!!mean_result_col := round(mean(!!result_col, na.rm = TRUE), 2)) %>%
    ungroup()
  
  hcc_final <- hcc_raw1 %>%
    group_by(Subject_ID, .data[[result_month_col]]) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  hcc_last <- hcc_final %>%
    group_by(Subject_ID) %>%
    filter(!!result_date_col == max(!!result_date_col)) %>%
    mutate(SELECT_FLAG = 1, last_type = "HCC_pre_HCC") %>%
    ungroup()
  
  # ----------------- Non-HCC -----------------
  nonhcc_raw <- included_all %>%
    filter(HCC_group == "Non-HCC") %>%
    arrange(Subject_ID, !!result_date_col) %>%
    mutate(!!result_month_col := floor_date(!!result_date_col, "month"))
  
  nonhcc_raw1 <- nonhcc_raw %>%
    group_by(Subject_ID, .data[[result_month_col]]) %>%
    mutate(!!mean_result_col := round(mean(!!result_col, na.rm = TRUE), 2)) %>%
    ungroup()
  
  nonhcc_final <- nonhcc_raw1 %>%
    group_by(Subject_ID, .data[[result_month_col]]) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(First_HCC_MONTH = as.Date(NA), First_HCC_YEAR = as.integer(NA))
  
  nonhcc_last <- nonhcc_final %>%
    group_by(Subject_ID) %>%
    filter(!!result_date_col == max(!!result_date_col)) %>%
    mutate(SELECT_FLAG = 1, last_type = "NonHCC_latest") %>%
    ungroup()
  
  # ----------------- Combine Outputs -----------------
  combined_final <- bind_rows(hcc_final, nonhcc_final)
  combined_raw   <- bind_rows(hcc_raw1, nonhcc_raw1)
  last_all       <- bind_rows(hcc_last, nonhcc_last)
  
  list(
    raw_all    = combined_raw,    # raw rows with monthly means for included data
    monthly    = combined_final,  # 1 row per Subject_ID per month
    last_all   = last_all,        # last HCC pre-HCC + last Non-HCC
    exclusions = exclusions       # excluded rows with reasons
  )
}

# excl_by_id <- out$exclusions %>%
#   group_by(Subject_ID, HCC_group) %>%
#   summarize(n_excluded = dplyr::n(),
#             reasons = paste(sort(unique(reason)), collapse = ";"),
#             .groups = "drop")

#### Demographics Data ####
dem <- read_csv(
  demographic_file,
  show_col_types = FALSE
)

#### Result Condensing ####

#### Albumin ####

#Read Data
alb <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Albumin Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "ALB")

#Process
alb_out <- process(
  lab_df          = alb,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  4.9,
  lab_shortname   = "ALB",
  dem             = dem
)

#Check Raw Albumin Results Subjects 
alb_out$raw_all %>%
  dplyr::select(Subject_ID) %>%
  group_by(Subject_ID) %>%
  count() #4371 unique patients 

#Save raw albumin results 
alb_r <- alb_out$raw_all

#Monthly 
alb_m <- alb_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = ALB_RESULT_MONTH, 
    mean_RESULT_NUM = mean_ALB_RESULT_NUM)

summary(alb_m$mean_RESULT_NUM)
names(alb_m)
# write_csv(alb_m, file.path(mean_biomarker_dir, "alb_mean_data_110525.csv"))

#Exclusions 
alb_ex <- alb_out$exclusions 
alb_ex$Lab_Name <- "Albumin"
alb_excl_summary <- alb_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(alb_ex, file.path(exclusion_dir, "Albumin Patient Exclusion 12.06.25.csv"))

#Latest lab results 
alb_latest <- alb_out$last_all

#### ALT ####
alt <- read_csv(
  file.path(raw_biomarker_dir, "All Raw ALT Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "ALT")
alt_out <- process(
  lab_df          = alt,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  553,
  lab_shortname   = "ALT",
  dem             = dem
)

#Raw ALT Results 
alt_r <- alt_out$raw_all
alt_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
alt_m <- alt_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = ALT_RESULT_MONTH, 
    mean_RESULT_NUM = mean_ALT_RESULT_NUM)

# write_csv(alt_m, file.path(mean_biomarker_dir, "alt_mean_data_110525.csv"))

#Exclusions 
alt_ex <- alt_out$exclusions
alt_ex$Lab_Name <- "ALT"
alt_excl_summary <- alt_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(alt_ex, file.path(exclusion_dir, "ALT Patient Exclusion 12.06.25.csv"))

#Latest 
alt_latest <- alt_out$last_all

#### AST ####
ast <- read_csv(
  file.path(raw_biomarker_dir, "All Raw AST Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "AST")
ast_out <- process(
  lab_df          = ast,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  676.12,
  lab_shortname   = "AST",
  dem             = dem
)
#Raw AST Results 
ast_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()
ast_r <- ast_out$raw_all

#Monthly 
ast_m <- ast_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = AST_RESULT_MONTH, 
    mean_RESULT_NUM = mean_AST_RESULT_NUM)

# write_csv(ast_m, file.path(mean_biomarker_dir, "ast_mean_data_110525.csv"))

#Exclusions 
ast_ex <- ast_out$exclusions
ast_ex$Lab_Name <- "AST"
ast_excl_summary <- ast_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(ast_ex, file.path(exclusion_dir, "AST Patient Exclusion 12.06.25.csv"))

#Latest 
ast_l <- ast_out$last_all

#### Bilirubin ####
bili <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Bilirubin Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "BILI")
bili_out <- process(
  lab_df          = bili,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  21.9,
  lab_shortname   = "BILI",
  dem             = dem
)
#Raw Bilirubin Results 
bili_r <- bili_out$raw_all 
bili_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
bili_m <- bili_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = BILI_RESULT_MONTH, 
    mean_RESULT_NUM = mean_BILI_RESULT_NUM)

# write_csv(bili_m, file.path(mean_biomarker_dir, "bili_mean_data_110525.csv"))

#Exclusions 
bili_ex <- bili_out$exclusions
bili_ex$Lab_Name <- "Bilirubin"
bili_excl_summary <- bili_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(bili_ex, file.path(exclusion_dir, "Bilirubin Patient Exclusion 12.06.25.csv"))

#Latest 
bili_l <- bili_out$last_all

#### Total Cholesterol ####
tchol <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Total Cholesterol Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "TC")
tchol_out <- process(
  lab_df          = tchol,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  316,
  lab_shortname   = "TCHOL",
  dem             = dem
)

#Raw Total Cholesterol Results 
tchol_r <- tchol_out$raw_all
tchol_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

names(bili_m)

#Monthly 
tchol_m <- tchol_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = TCHOL_RESULT_MONTH, 
    mean_RESULT_NUM = mean_TCHOL_RESULT_NUM)

#Latest 
tchol_l <- tchol_out$last_all

#Exclusions 
tchol_ex <- tchol_out$exclusions
tchol_ex$Lab_Name <- "Total Cholesterol"
tchol_excl_summary <- tchol_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(tchol_ex, file.path(exclusion_dir, "Total Cholesterol Patient Exclusion 12.06.25.csv"))

#### LDL Cholesterol ####
ldl <- read_csv(
  file.path(raw_biomarker_dir, "All Raw LDL Cholesterol Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "LDL")

ldl_out <- process(
  lab_df          = ldl,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  197,
  lab_shortname   = "LDLCHOL",
  dem             = dem
)

#Raw LDL Cholesterol Results 
ldl_r <- ldl_out$raw_all 
ldl_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
ldl_m <- ldl_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = LDLCHOL_RESULT_MONTH, 
    mean_RESULT_NUM = mean_LDLCHOL_RESULT_NUM)

ldl_m1 <- ldl_m %>%
  dplyr::select(Subject_ID:First_LD_DX, HCC_DX:mean_RESULT_NUM)

#Exclusions 
ldl_ex <- ldl_out$exclusions
ldl_ex$Lab_Name <- "LDL Cholesterol"
ldl_excl_summary <- ldl_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(ldl_ex, file.path(exclusion_dir, "LDL Cholesterol Patient Exclusion 12.06.25.csv"))

#Latest
ldl_l <- ldl_out$last_all

#### HDL Cholesterol ####
hdl <- read_csv(
  file.path(raw_biomarker_dir, "All Raw HDL Cholesterol Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "HDL")

hdl_out <- process(
  lab_df          = hdl,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  112,
  lab_shortname   = "HDLCHOL",
  dem             = dem
)

#Raw HDL Cholesterol Results 
hdl_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()
hdl_r <- hdl_out$raw_all

#Monthly 
hdl_m <- hdl_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = HDLCHOL_RESULT_MONTH, 
    mean_RESULT_NUM = mean_HDLCHOL_RESULT_NUM)

#Exclusions 
hdl_ex <- hdl_out$exclusions
hdl_ex$Lab_Name <- "HDL Cholesterol"
hdl_excl_summary <- hdl_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(hdl_ex, file.path(exclusion_dir, "HDL Cholesterol Patient Exclusion 12.06.25.csv"))

#### Triglycerides ####
trig <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Triglycerides Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "TG")
trig_out <- process(
  lab_df          = trig,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  1657,
  lab_shortname   = "TRIG",
  dem             = dem
)

#Raw Triglycerides Results 
trig_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()
trig_r <- trig_out$raw_all

#Monthly 
trig_m <- trig_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = TRIG_RESULT_MONTH, 
    mean_RESULT_NUM = mean_TRIG_RESULT_NUM)

#Latest 
trig_l <- trig_out$last_all

#Exclusions 
trig_ex <- trig_out$exclusions
trig_ex$Lab_Name <- "Triglycerides"
trig_excl_summary <- trig_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(trig_ex, file.path(exclusion_dir, "Triglycerides Patient Exclusion 12.06.25.csv"))


#### Hemoglobin ####
hemo <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Hemoglobin Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "HGB")
hemo_out <- process(
  lab_df          = hemo,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  17,
  lab_shortname   = "HGB",
  dem             = dem
)
#Raw Hemoglobin Results 
hemo_r <- hemo_out$raw_all
hemo_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
hemo_m <- hemo_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = HGB_RESULT_MONTH, 
    mean_RESULT_NUM = mean_HGB_RESULT_NUM)

# write_csv(hemo_m, file.path(mean_biomarker_dir, "hemo_mean_data_110525.csv"))

#Latest 
hemo_l <- hemo_out$last_all

#Exclusions 
hemo_ex <- hemo_out$exclusions
hemo_ex$Lab_Name <- "Hemoglobin"
hemo_excl_summary <- hemo_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(hemo_ex, file.path(exclusion_dir, "Hemoglobin Patient Exclusion 12.06.25.csv"))

#### INR ####
inr <- read_csv(
  file.path(raw_biomarker_dir, "All Raw INR Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "INR")
inr_out <- process(
  lab_df          = inr,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  4.8,
  lab_shortname   = "INR",
  dem             = dem
)
#Raw INR Results 
inr_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()
inr_r <- inr_out$raw_all

#Monthly 
inr_m <- inr_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = INR_RESULT_MONTH, 
    mean_RESULT_NUM = mean_INR_RESULT_NUM)

# write_csv(inr_m, file.path(mean_biomarker_dir, "inr_mean_data_110525.csv"))

#Latest 
inr_l <- inr_out$last_all

#Exclusions 
inr_ex <- inr_out$exclusions
inr_ex$Lab_Name <- "INR"
inr_excl_summary <- inr_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(inr_ex, file.path(exclusion_dir, "INR Patient Exclusion 12.06.25.csv"))

#### Platelet ####
plat <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Platelets Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "PLAT")
plat_out <- process(
  lab_df          = plat,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  513,
  lab_shortname   = "PLAT",
  dem             = dem
)
#Raw Platelets Results 
plat_r <- plat_out$raw_all 
plat_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
plat_m <- plat_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = PLAT_RESULT_MONTH, 
    mean_RESULT_NUM = mean_PLAT_RESULT_NUM)

# write_csv(plat_m, file.path(mean_biomarker_dir, "plat_mean_data_110525.csv"))

#Latest 
plat_l <- plat_out$last_all

#Exclusions 
plat_ex <- plat_out$exclusions
plat_ex$Lab_Name <- "Platelets"
plat_excl_summary <- plat_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(plat_ex, file.path(exclusion_dir, "Platelet Patient Exclusion 12.06.25.csv"))

#### Alpha-Feto Protein ####
afp <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Alpha-Feto Protein Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "AFP")
afp_out <- process(
  lab_df          = afp,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  11619,
  lab_shortname   = "AFP",
  dem             = dem
)
#Raw Alpha-Feto Protein Results 
afp_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
afp_m <- afp_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = AFP_RESULT_MONTH, 
    mean_RESULT_NUM = mean_AFP_RESULT_NUM)

# write_csv(afp_m, file.path(mean_biomarker_dir, "afp_mean_data_110525.csv"))

#Latest 
afp_l <- afp_out$last_all

#Exclusions 
afp_ex <- afp_out$exclusions
afp_ex$Lab_Name <- "Alpha-Fetoprotein"
afp_excl_summary <- afp_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(afp_ex, file.path(exclusion_dir, "Alpha-Feto Protein Patient Exclusion 12.06.25.csv"))

#### Neutrophils ####
neut <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Neutrophils Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "NEUT")
neut_out <- process(
  lab_df          = neut,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  21.16,
  lab_shortname   = "NEUT",
  dem             = dem
)
#Raw Neutrophils Results 
neut_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
neut_m <- neut_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = NEUT_RESULT_MONTH, 
    mean_RESULT_NUM = mean_NEUT_RESULT_NUM)

# write_csv(neut_m, file.path(mean_biomarker_dir, "neut_mean_data_110525.csv"))

#Latest 
neut_l <- neut_out$last_all

#Exclusions 
neut_ex <- neut_out$exclusions
neut_ex$Lab_Name <- "Neutrophils"
neut_excl_summary <- neut_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(neut_ex, file.path(exclusion_dir, "Neutrophil Patient Exclusion 12.06.25.csv"))

#### Lymphocytes ####
lymph <- read_csv(
  file.path(raw_biomarker_dir, "All Raw Lymphocytes Data.csv"),
  show_col_types = FALSE
) %>%
  select(-any_of(c("X", "...1"))) %>%
  mutate(LAB_TYPE = "LYMPH")
lymph_out <- process(
  lab_df          = lymph,
  result_col      = RESULT_NUM,
  result_date_col = RESULT_DATE,
  result_threshold =  5.04,
  lab_shortname   = "LYMPH",
  dem             = dem
)
#Raw Lymphocytes Results 
lymph_out$raw_all %>% dplyr::select(Subject_ID) %>% group_by(Subject_ID) %>% count()

#Monthly 
lymph_m <- lymph_out$monthly %>%
  dplyr::select(-c(8:9, 22, 24:26)) %>%
  dplyr::select(-c(RAW_LAB_NAME, LAB_ORDER_DATE, RESULT_NUM, HCC_YEAR)) %>%
  rename(
    RESULT_MONTH = LYMPH_RESULT_MONTH, 
    mean_RESULT_NUM = mean_LYMPH_RESULT_NUM)

# write_csv(lymph_m, file.path(mean_biomarker_dir, "lymph_mean_data_110525.csv"))

#Latest 
lymph_l <- lymph_out$last_all

#Exclusions 
lymph_ex <- lymph_out$exclusions
lymph_ex$Lab_Name <- "Lymphocytes"
lymph_excl_summary <- lymph_ex %>% group_by(reason) %>% summarise(n = dplyr::n(), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pct = round(100*n/sum(n), 1))

write_csv(lymph_ex, file.path(exclusion_dir, "Lymphocyte Patient Exclusion 12.06.25.csv"))

#### MEAN - Cholesterol Filling Formula ####
tchol_m2 <- tchol_m %>% transmute(Subject_ID, RESULT_MONTH, mean_RESULT_NUM, RESULT_UNIT = "mg/dL", LAB_TYPE = "TC")
ldl_m2   <- ldl_m   %>% transmute(Subject_ID, RESULT_MONTH, mean_RESULT_NUM, RESULT_UNIT = "mg/dL", LAB_TYPE = "LDL")
hdl_m2   <- hdl_m   %>% transmute(Subject_ID, RESULT_MONTH, mean_RESULT_NUM, RESULT_UNIT = "mg/dL", LAB_TYPE = "HDL")
trig_m2  <- trig_m  %>% transmute(Subject_ID, RESULT_MONTH, mean_RESULT_NUM,RESULT_UNIT = "mg/dL", LAB_TYPE = "TG")

all_lipids_m <- bind_rows(tchol_m2, ldl_m2, hdl_m2, trig_m2)

lipids_wide_m <- all_lipids_m %>%
  pivot_wider(
    id_cols = c(Subject_ID, RESULT_MONTH),
    names_from = LAB_TYPE,
    values_from = mean_RESULT_NUM, 
    values_fn = mean
  )


lipids_filled_m <- lipids_wide_m %>%
  mutate(
    value_count = rowSums(!is.na(select(., TC, LDL, HDL, TG))),
    original_na_TC = is.na(TC),
    original_na_LDL = is.na(LDL),
    original_na_HDL = is.na(HDL),
    original_na_TG = is.na(TG),
    TC = ifelse(value_count == 3 & is.na(TC), LDL + HDL + (TG / 5), TC),
    LDL = ifelse(value_count == 3 & is.na(LDL), TC - HDL - (TG / 5), LDL),
    HDL = ifelse(value_count == 3 & is.na(HDL), TC - LDL - (TG / 5), HDL),
    TG = ifelse(value_count == 3 & is.na(TG), (TC - LDL - HDL) * 5, TG)
  )

lipids_long_final_m <- lipids_filled_m %>%
  group_by(Subject_ID, RESULT_MONTH) %>%
  pivot_longer(
    cols = c(TC, LDL, HDL, TG),
    names_to = "LAB_TYPE",
    values_to = "RESULT_NUM",
    values_drop_na = TRUE
  )

imputed_rows_m <- lipids_long_final_m %>%
  filter(
    (LAB_TYPE == "TC" & original_na_TC) |
      (LAB_TYPE == "LDL" & original_na_LDL) |
      (LAB_TYPE == "HDL" & original_na_HDL) |
      (LAB_TYPE == "TG" & original_na_TG)
  )

if (nrow(imputed_rows_m) > 0) {
  imputed_rows_m <- imputed_rows_m %>%
    mutate(
      RAW_LAB_NAME = case_when(
        LAB_TYPE == "TC" ~ "Filled TC",
        LAB_TYPE == "LDL" ~ "Filled LDL",
        LAB_TYPE == "HDL" ~ "Filled HDL",
        LAB_TYPE == "TG" ~ "Filled TRIGLYCERIDES",
        TRUE ~ as.character(LAB_TYPE)
      ),
      RESULT_UNIT = "mg/dL"
    ) %>%
    select(Subject_ID, RAW_LAB_NAME, RESULT_MONTH, RESULT_NUM, RESULT_UNIT, LAB_TYPE)
  
  final_combined_data_m <- bind_rows(all_lipids_m, imputed_rows_m)
} else {
  final_combined_data_m <- all_lipids_m
}

fc1 <- final_combined_data_m %>%
  group_by(Subject_ID, RESULT_MONTH, LAB_TYPE) %>%
  mutate(
    mean_RESULT_NUM = case_when(
      is.na(mean_RESULT_NUM) & str_detect(RAW_LAB_NAME, regex("Filled", ignore_case = TRUE)) ~ RESULT_NUM,
      TRUE ~ mean_RESULT_NUM
    )
  ) %>%
  dplyr::select(-c(RAW_LAB_NAME, RESULT_NUM))

fc2 <- final_combined_data_m %>%
  group_by(Subject_ID, RESULT_MONTH, LAB_TYPE) %>%
  mutate(
    # 1) Use RESULT_NUM when RAW_LAB_NAME contains "Filled"
    .from_filled = is.na(mean_RESULT_NUM) &
      str_detect(RAW_LAB_NAME, regex("Filled", ignore_case = TRUE)) &
      !is.na(RESULT_NUM),
    mean_RESULT_NUM = if_else(.from_filled, RESULT_NUM, mean_RESULT_NUM),
    
    # 2) If still NA, take any non-missing value from the same (Subject, Month, Lab)
    .grp_value = suppressWarnings(dplyr::first(na.omit(mean_RESULT_NUM))),
    .from_group = is.na(mean_RESULT_NUM) & !is.na(.grp_value),
    mean_RESULT_NUM = if_else(is.na(mean_RESULT_NUM), .grp_value, mean_RESULT_NUM),
    
    # Optional: coerce to numeric if mixed types slipped in
    mean_RESULT_NUM = suppressWarnings(as.numeric(mean_RESULT_NUM))
  ) %>%
  ungroup() %>%
  # Keep flags if you want; remove helper column
  select(-c(RAW_LAB_NAME, RESULT_NUM, .grp_value)) %>%
  distinct(Subject_ID, RESULT_MONTH, LAB_TYPE, .keep_all = TRUE)


#05MAR202120200127200201150

### New Total Cholesterol ####
new_tchol_m <- fc2 %>% 
  filter(LAB_TYPE == "TC")

tchol_dem <- tchol_m %>%
  dplyr::select(Subject_ID, BIRTH_DATE:within_period, First_HCC_DX_DATE, First_HCC_MONTH, First_HCC_YEAR) %>%
  group_by(Subject_ID) %>%
  distinct(Subject_ID, .keep_all = TRUE) %>%
  dplyr::select(Subject_ID, within_period, First_HCC_MONTH, First_HCC_YEAR)

tchol_m1 <- tchol_m %>%
  dplyr::select(Subject_ID, RESULT_DATE, RESULT_MONTH) %>%
  group_by(Subject_ID) %>%
  full_join(new_tchol_m, by = c("Subject_ID", "RESULT_MONTH")) %>%
  arrange(Subject_ID, RESULT_MONTH) %>%
  left_join(dem, by = "Subject_ID") %>%
  dplyr::select(-c(7:10, HCC_YEAR)) %>%
  left_join(tchol_dem, by = "Subject_ID")

names(lymph_m)
names(tchol_m1)

write_csv(tchol_m1, file.path(mean_biomarker_dir, "tchol_mean_data_110525.csv"))

### New LDL Cholesterol ####
new_ldl_m <- fc2 %>% 
  filter(LAB_TYPE == "LDL")

ldl_dem <- ldl_m %>%
  dplyr::select(Subject_ID, BIRTH_DATE:within_period, First_HCC_DX_DATE, First_HCC_MONTH, First_HCC_YEAR) %>%
  group_by(Subject_ID) %>%
  distinct(Subject_ID, .keep_all = TRUE) %>%
  dplyr::select(Subject_ID, within_period, First_HCC_MONTH, First_HCC_YEAR)

ldl_m1 <- ldl_m %>%
  dplyr::select(Subject_ID, RESULT_DATE, RESULT_MONTH) %>%
  group_by(Subject_ID) %>%
  full_join(new_ldl_m, by = c("Subject_ID", "RESULT_MONTH")) %>%
  arrange(Subject_ID, RESULT_MONTH) %>%
  left_join(dem, by = "Subject_ID") %>%
  dplyr::select(-c(7:10, HCC_YEAR)) %>%
  left_join(tchol_dem, by = "Subject_ID") %>%
  mutate(mean_RESULT_NUM = ifelse(mean_RESULT_NUM<0, NA, mean_RESULT_NUM))

write_csv(ldl_m1, file.path(mean_biomarker_dir, "ldlchol_mean_data_110525.csv"))

### New HDL Cholesterol ####
new_hdl_m <- fc2 %>% 
  filter(LAB_TYPE == "HDL")

hdl_dem <- hdl_m %>%
  dplyr::select(Subject_ID, BIRTH_DATE:within_period, First_HCC_DX_DATE, First_HCC_MONTH, First_HCC_YEAR) %>%
  group_by(Subject_ID) %>%
  distinct(Subject_ID, .keep_all = TRUE) %>%
  dplyr::select(Subject_ID, within_period, First_HCC_MONTH, First_HCC_YEAR)

hdl_m1 <- hdl_m %>%
  dplyr::select(Subject_ID, RESULT_DATE, RESULT_MONTH) %>%
  group_by(Subject_ID) %>%
  full_join(new_hdl_m, by = c("Subject_ID", "RESULT_MONTH")) %>%
  arrange(Subject_ID, RESULT_MONTH) %>%
  left_join(dem, by = "Subject_ID") %>%
  dplyr::select(-c(7:10, HCC_YEAR)) %>%
  left_join(tchol_dem, by = "Subject_ID")

write_csv(hdl_m1, file.path(mean_biomarker_dir, "hdlchol_mean_data_110525.csv"))


### New Triglycerides ####
new_trig_m <- fc2 %>% 
  filter(LAB_TYPE == "TG")

trig_dem <- trig_m %>%
  dplyr::select(Subject_ID, BIRTH_DATE:within_period, First_HCC_DX_DATE, First_HCC_MONTH, First_HCC_YEAR) %>%
  group_by(Subject_ID) %>%
  distinct(Subject_ID, .keep_all = TRUE) %>%
  dplyr::select(Subject_ID, within_period, First_HCC_MONTH, First_HCC_YEAR)

trig_m1 <- trig_m %>%
  dplyr::select(Subject_ID, RESULT_DATE, RESULT_MONTH) %>%
  group_by(Subject_ID) %>%
  full_join(new_trig_m, by = c("Subject_ID", "RESULT_MONTH")) %>%
  arrange(Subject_ID, RESULT_MONTH) %>%
  left_join(dem, by = "Subject_ID") %>%
  dplyr::select(-c(7:10, HCC_YEAR)) %>%
  left_join(tchol_dem, by = "Subject_ID") %>%
  mutate(mean_RESULT_NUM = ifelse(mean_RESULT_NUM<0, NA, mean_RESULT_NUM))

write_csv(trig_m1, file.path(mean_biomarker_dir, "trig_mean_data_110525.csv"))


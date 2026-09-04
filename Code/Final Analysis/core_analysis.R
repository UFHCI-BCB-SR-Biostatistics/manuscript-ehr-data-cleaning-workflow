suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(survival)
  library(survivalROC)
  library(ggplot2)
  library(purrr)
  library(patchwork)
  library(scales)
})

# Biomarker definitions

biomarker_catalog <- tibble::tribble(
  ~input_name, ~label,                             ~required_cols,
  "ALT",       "ALT (U/L)",                       c("mean_ALT_RESULT_NUM"),
  "AST",       "AST (U/L)",                       c("mean_AST_RESULT_NUM"),
  "ALB",       "Albumin (g/dL)",                  c("mean_ALB_RESULT_NUM"),
  "BILI",      "Total Bilirubin (mg/dL)",         c("mean_BILI_RESULT_NUM"),
  "INR",       "INR",                             c("mean_INR_RESULT_NUM"),
  "PLAT",      "Platelet count (10^3/uL)",        c("mean_PLAT_RESULT_NUM"),
  "HGB",       "Hemoglobin (g/dL)",               c("mean_HGB_RESULT_NUM"),
  "NEUT",      "Neutrophil count (10^3/uL)",      c("mean_NEUT_RESULT_NUM"),
  "LYMPH",     "Lymphocyte count (10^3/uL)",      c("mean_LYMPH_RESULT_NUM"),
  "AFP",       "AFP (ng/mL)",                     c("mean_AFP_RESULT_NUM"),
  "TC",        "Total Cholesterol (mg/dL)",       c("mean_TC_RESULT_NUM"),
  "HDL",       "HDL Cholesterol (mg/dL)",         c("mean_HDL_RESULT_NUM"),
  "LDL",       "LDL Cholesterol (mg/dL)",         c("mean_LDL_RESULT_NUM"),
  "TG",        "Triglycerides (mg/dL)",           c("mean_TG_RESULT_NUM")
)

format_model_spec_label <- function(input_biomarker,
                                     primary_longitudinal_biomarkers) {
  model_key <- paste(input_biomarker, collapse = "||")
  primary_key <- paste(primary_longitudinal_biomarkers, collapse = "||")
  model_labels <- c(
    setNames("all biomarker", primary_key),
    setNames(
      "all biomarker (AST/ALT ratio replaces AST + ALT)",
      paste(c(setdiff(primary_longitudinal_biomarkers, c("ALT", "AST")), "AST_ALT_RATIO"), collapse = "||")
    ),
    setNames(
      "all biomarker (NLR replaces NEUT + LYMPH)",
      paste(c(setdiff(primary_longitudinal_biomarkers, c("NEUT", "LYMPH")), "NLR"), collapse = "||")
    ),
    setNames(
      "all biomarker (PLR replaces PLAT + LYMPH)",
      paste(c(setdiff(primary_longitudinal_biomarkers, c("PLAT", "LYMPH")), "PLR"), collapse = "||")
    )
  )
  dplyr::coalesce(unname(model_labels[model_key]), paste(input_biomarker, collapse = " + "))
}

# Discrimination at the 36-month prediction horizon

compute_performance_metrics <- function(eval_df,
                                        risk_score,
                                        predict_time) {
  perf_df <- eval_df %>% mutate(risk_score = as.numeric(risk_score))
  c_index <- survival::concordance(
    Surv(duration, hcc_status) ~ risk_score,
    data = perf_df,
    reverse = TRUE
  )$concordance
  time_auc <- survivalROC(
    Stime = perf_df$duration,
    status = perf_df$hcc_status,
    marker = perf_df$risk_score,
    predict.time = predict_time,
    method = "NNE",
    span = 0.25 * nrow(perf_df)^(-0.20)
  )$AUC
  eval_horizon_df <- perf_df %>%
    filter(duration >= predict_time | (hcc_status == 1L & duration < predict_time))

  list(
    c_index = as.numeric(c_index),
    auc = time_auc,
    evaluable_n = nrow(eval_horizon_df),
    horizon_events = sum(eval_horizon_df$hcc_status == 1L & eval_horizon_df$duration <= predict_time),
    horizon_nonevents = sum(eval_horizon_df$hcc_status == 0L | eval_horizon_df$duration > predict_time)
  )
}

# Subject-level bootstrap

build_subject_bootstrap_weights <- function(subject_ids) {
  
  sampled_subject_ids <- sample(
    subject_ids,
    size = length(subject_ids),
    replace = TRUE
  )
  
  sampled_table <- table(sampled_subject_ids)
  
  tibble::tibble(
    Subject_ID = names(sampled_table),
    boot_weight = as.integer(sampled_table)
  )
}

expand_subject_level_dataset_by_weights <- function(eval_df,
                                                    bootstrap_weights) {
  
  eval_df %>%
    left_join(bootstrap_weights, by = "Subject_ID") %>%
    mutate(
      boot_weight = dplyr::coalesce(boot_weight, 0L)
    ) %>%
    filter(boot_weight > 0) %>%
    slice(rep(seq_len(n()), boot_weight)) %>%
    select(-boot_weight)
}

compute_corrected_performance <- function(apparent_performance,
                                          optimism_summary,
                                          prefix = "optimism_") {
  list(
    c_index = apparent_performance$c_index - optimism_summary[[paste0(prefix, "c_index")]],
    auc = apparent_performance$auc - optimism_summary[[paste0(prefix, "auc")]],
    evaluable_n = apparent_performance$evaluable_n,
    horizon_events = apparent_performance$horizon_events,
    horizon_nonevents = apparent_performance$horizon_nonevents
  )
}

# Data ingestion and preprocessing

# Laboratory results are represented by one record per patient and calendar
# month. 
collapse_subject_month_records <- function(df) {
  df %>%
    select(-within_period) %>%
    group_by(Subject_ID, RESULT_MONTH) %>%
    summarise(
      across(
        everything(),
        ~ dplyr::first(.x[!is.na(.x)], default = .x[NA_integer_])
      ),
      .groups = "drop"
    )
}

read_and_prepare_data <- function(data_path, bmi_data_path, min_followup_months) {
  
  raw <- read.csv(
    data_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )

  raw <- raw[, -1, drop = FALSE] %>%
    collapse_subject_month_records()
  
  bmi_data <- read.csv(
    file = bmi_data_path,
    stringsAsFactors = FALSE,
    check.names = TRUE,
    na.strings = c("", "NA", "NaN")
  )
  bmi_data$Subject_ID <- as.character(bmi_data$Subject_ID)
  
  raw <- raw %>%
    left_join(bmi_data %>% dplyr::select(Subject_ID, BMI_final), by = "Subject_ID")
  
  df <- raw %>%
    mutate(
      Subject_ID      = as.character(Subject_ID),
      BIRTH_DATE      = ymd(BIRTH_DATE),
      First_LD_DATE   = lubridate::floor_date(ymd(First_LD_DATE), "month"),
      First_HCC_YEAR  = suppressWarnings(as.integer(First_HCC_YEAR)),
      First_HCC_MONTH = suppressWarnings(as.integer(First_HCC_MONTH)),
      First_HCC_DATE  = dplyr::coalesce(
        ymd(First_HCC_DATE),
        make_date(First_HCC_YEAR, First_HCC_MONTH, 1L)
      ),
      RESULT_MONTH    = ymd(RESULT_MONTH),
      SEX             = factor(SEX),
      RACE            = factor(RACE),
      HISPANIC        = factor(HISPANIC),
      HCC_group       = as.character(HCC_group)
    )
  
  numeric_candidates <- c(
    "AGE", "BMI_final",
    "mean_ALT_RESULT_NUM",  "mean_AST_RESULT_NUM",  "mean_ALB_RESULT_NUM",
    "mean_BILI_RESULT_NUM", "mean_INR_RESULT_NUM",  "mean_PLAT_RESULT_NUM",
    "mean_HGB_RESULT_NUM",  "mean_NEUT_RESULT_NUM", "mean_LYMPH_RESULT_NUM",
    "mean_AFP_RESULT_NUM",  "mean_TC_RESULT_NUM",   "mean_HDL_RESULT_NUM",
    "mean_LDL_RESULT_NUM",  "mean_TG_RESULT_NUM"
  )
  
  df <- df %>%
    mutate(across(all_of(numeric_candidates), ~ suppressWarnings(as.numeric(.))))
  
  subject_raw_cohort <- df %>%
    group_by(Subject_ID) %>%
    summarise(
      First_LD_DATE  = first(First_LD_DATE[!is.na(First_LD_DATE)]),
      HCC_group      = first(HCC_group[!is.na(HCC_group)]),
      First_HCC_DATE = first(First_HCC_DATE[!is.na(First_HCC_DATE)]),
      last_lab_date  = max(RESULT_MONTH, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      last_lab_date = if_else(is.infinite(last_lab_date), as.Date(NA), last_lab_date),
      eligibility_months = 12L * (year(last_lab_date) - year(First_LD_DATE)) +
        (month(last_lab_date) - month(First_LD_DATE)),
      hcc_status = if_else(HCC_group == "HCC", 1L, 0L),
      followup_end_date = case_when(
        hcc_status == 1L ~ First_HCC_DATE,
        hcc_status == 0L ~ last_lab_date
      ),
      duration = as.numeric(
        time_length(interval(First_LD_DATE, followup_end_date), "months")
      )
    ) %>%
    filter(
      !is.na(First_LD_DATE),
      !is.na(followup_end_date),
      !is.na(duration),
      duration > 0
    )
  
  subject_followup <- subject_raw_cohort %>%
    filter(eligibility_months >= min_followup_months)
  
  baseline_raw_data <- df %>%
    inner_join(
      subject_raw_cohort %>%
        select(Subject_ID, First_LD_DATE, followup_end_date, duration, hcc_status),
      by = "Subject_ID"
    ) %>%
    arrange(Subject_ID, RESULT_MONTH) %>%
    group_by(Subject_ID) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  baseline_data <- df %>%
    inner_join(
      subject_followup %>%
        select(Subject_ID, First_LD_DATE, followup_end_date, duration, hcc_status),
      by = "Subject_ID"
    ) %>%
    arrange(Subject_ID, RESULT_MONTH) %>%
    group_by(Subject_ID) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  list(
    analysis_data = df,
    subject_raw_cohort = subject_raw_cohort,
    subject_followup = subject_followup,
    baseline_raw_data = baseline_raw_data,
    baseline_data = baseline_data,
    biomarker_catalog = biomarker_catalog
  )
}

# Follow-up distribution
# The reverse Kaplan-Meier estimator treats HCC diagnoses as censored and
# follow-up termination among patients without HCC as the event of interest.

build_followup_summary_table <- function(subject_followup, analysis_data) {
  followup_data <- subject_followup %>%
    transmute(
      Subject_ID,
      hcc_status,
      followup_months = duration,
      followup_years = duration / 12,
      reverse_status = 1L - hcc_status
    )

  reverse_followup_fit <- survival::survfit(
    survival::Surv(followup_months, reverse_status) ~ 1,
    data = followup_data
  )

  followup_quantiles <- as.numeric(
    stats::quantile(
      reverse_followup_fit,
      probs = c(0.25, 0.50, 0.75)
    )$quantile
  )

  laboratory_months <- analysis_data %>%
    select(Subject_ID, RESULT_MONTH) %>%
    inner_join(
      subject_followup %>%
        select(Subject_ID, First_LD_DATE, followup_end_date),
      by = "Subject_ID"
    ) %>%
    filter(
      !is.na(RESULT_MONTH),
      RESULT_MONTH >= First_LD_DATE,
      RESULT_MONTH <= followup_end_date
    ) %>%
    count(Subject_ID, name = "months_with_measurements")

  laboratory_month_quantiles <- quantile(
    laboratory_months$months_with_measurements,
    probs = c(0.25, 0.50, 0.75),
    names = FALSE
  )
  hcc_time_quantiles <- quantile(
    followup_data$followup_months[followup_data$hcc_status == 1L],
    probs = c(0.25, 0.50, 0.75),
    names = FALSE
  )
  person_years <- sum(followup_data$followup_years)
  hcc_events <- sum(followup_data$hcc_status == 1L)

  tibble::tibble(
    `Estimation method` = "Reverse Kaplan-Meier",
    `Subjects` = nrow(subject_followup),
    `Person-years` = person_years,
    `HCC events` = hcc_events,
    `HCC incidence (%)` = 100 * hcc_events / nrow(subject_followup),
    `Incidence rate per 1,000 person-years` = 1000 * hcc_events / person_years,
    `First quartile (months)` = followup_quantiles[1],
    `Median follow-up (months)` = followup_quantiles[2],
    `Third quartile (months)` = followup_quantiles[3],
    `First quartile of months with measurements` = laboratory_month_quantiles[1],
    `Median months with measurements` = laboratory_month_quantiles[2],
    `Third quartile of months with measurements` = laboratory_month_quantiles[3],
    `First quartile of time to HCC (months)` = hcc_time_quantiles[1],
    `Median time to HCC (months)` = hcc_time_quantiles[2],
    `Third quartile of time to HCC (months)` = hcc_time_quantiles[3]
  )
}

# Analysis biomarker construction

append_model_biomarkers <- function(df) {
  df %>%
    mutate(
      bm_ALT = mean_ALT_RESULT_NUM,
      bm_AST = mean_AST_RESULT_NUM,
      bm_ALB = mean_ALB_RESULT_NUM,
      bm_BILI = mean_BILI_RESULT_NUM,
      bm_PLAT = mean_PLAT_RESULT_NUM,
      bm_NEUT = mean_NEUT_RESULT_NUM,
      bm_LYMPH = mean_LYMPH_RESULT_NUM,
      bm_HGB = mean_HGB_RESULT_NUM,
      bm_AST_ALT_RATIO = if_else(
        mean_ALT_RESULT_NUM > 0,
        mean_AST_RESULT_NUM / mean_ALT_RESULT_NUM,
        NA_real_
      ),
      bm_NLR = if_else(
        mean_LYMPH_RESULT_NUM > 0,
        mean_NEUT_RESULT_NUM / mean_LYMPH_RESULT_NUM,
        NA_real_
      ),
      bm_PLR = if_else(
        mean_LYMPH_RESULT_NUM > 0,
        mean_PLAT_RESULT_NUM / mean_LYMPH_RESULT_NUM,
        NA_real_
      )
    )
}

# Baseline biomarker missingness table

build_baseline_missingness_table <- function(prepared_obj) {
  
  reportable_biomarkers <- prepared_obj$biomarker_catalog %>%
    mutate(raw_col = purrr::map_chr(required_cols, 1))
  
  baseline_raw_data <- prepared_obj$baseline_raw_data
  total_subjects <- nrow(baseline_raw_data)
  
  purrr::map_dfr(seq_len(nrow(reportable_biomarkers)), function(i) {
    biomarker_name <- reportable_biomarkers$input_name[i]
    raw_col <- reportable_biomarkers$raw_col[i]
    available_n <- sum(!is.na(baseline_raw_data[[raw_col]]))
    missing_n <- sum(is.na(baseline_raw_data[[raw_col]]))
    
    tibble::tibble(
      Biomarker = biomarker_name,
      `Baseline available n` = available_n,
      `Baseline missing n` = missing_n,
      `Baseline missing %` = 100 * missing_n / total_subjects
    )
  })
}

# Longitudinal biomarker missingness table

build_longitudinal_missingness_table <- function(prepared_obj,
                                                 primary_longitudinal_biomarkers) {
  
  reportable_biomarkers <- prepared_obj$biomarker_catalog %>%
    mutate(raw_col = purrr::map_chr(required_cols, 1))
  
  eligible_monthly_rows <- prepared_obj$analysis_data %>%
    inner_join(
      prepared_obj$subject_raw_cohort %>%
        transmute(
          Subject_ID,
          first_ld_date_ref = First_LD_DATE,
          followup_end_date = followup_end_date
        ),
      by = "Subject_ID"
    ) %>%
    filter(
      !is.na(RESULT_MONTH),
      RESULT_MONTH >= first_ld_date_ref,
      RESULT_MONTH <= followup_end_date
    )
  
  total_rows <- nrow(eligible_monthly_rows)
  
  purrr::map_dfr(seq_len(nrow(reportable_biomarkers)), function(i) {
    biomarker_name <- reportable_biomarkers$input_name[i]
    raw_col <- reportable_biomarkers$raw_col[i]
    observed_n <- sum(!is.na(eligible_monthly_rows[[raw_col]]))
    missing_n <- sum(is.na(eligible_monthly_rows[[raw_col]]))
    
    tibble::tibble(
      Biomarker = biomarker_name,
      `Eligible monthly rows` = total_rows,
      `Observed rows` = observed_n,
      `Missing rows` = missing_n,
      `Missing %` = 100 * missing_n / total_rows,
      `Retained for primary longitudinal model` = if_else(
        biomarker_name %in% primary_longitudinal_biomarkers,
        "Yes",
        "No"
      )
    )
  })
}

# Time-varying start-stop dataset construction
# Complete-case eligibility is based on the primary eight-biomarker panel;
# derived ratios are retained for sensitivity analyses.

build_timevarying_dataset <- function(analysis_df,
                                      subject_followup,
                                      complete_case_cols,
                                      retain_analysis_cols = complete_case_cols) {
  
  analysis_df %>%
    inner_join(subject_followup, by = "Subject_ID", suffix = c("", "_subj")) %>%
    filter(
      !is.na(RESULT_MONTH),
      !is.na(First_LD_DATE_subj),
      RESULT_MONTH >= First_LD_DATE_subj,
      RESULT_MONTH <= followup_end_date
    ) %>%
    mutate(
      t_mo = as.numeric(time_length(interval(First_LD_DATE_subj, RESULT_MONTH), "months")),
      ev_time_mo = if_else(
        hcc_status == 1L,
        as.numeric(time_length(interval(First_LD_DATE_subj, First_HCC_DATE), "months")),
        NA_real_
      )
    ) %>%
    filter(
      !is.na(AGE),
      !is.na(SEX),
      !is.na(RACE),
      !is.na(HISPANIC),
      !is.na(BMI_final),
      if_all(all_of(complete_case_cols), ~ !is.na(.))
    ) %>%
    arrange(Subject_ID, t_mo) %>%
    group_by(Subject_ID) %>%
    mutate(
      duration_mo = first(duration),
      next_t = lead(t_mo, default = first(duration)),
      tstart = t_mo,
      tstop_raw = pmin(next_t, duration_mo),
      event = if_else(
        first(hcc_status) == 1L &
          first(ev_time_mo) > tstart &
          first(ev_time_mo) <= tstop_raw,
        1L,
        0L
      ),
      tstop = if_else(event == 1L, first(ev_time_mo), tstop_raw)
    ) %>%
    ungroup() %>%
    filter(tstop > tstart) %>%
    select(
      Subject_ID,
      AGE, SEX, RACE, HISPANIC, BMI_final,
      duration, hcc_status,
      RESULT_MONTH, First_LD_DATE_subj, First_HCC_DATE,
      tstart, tstop, event,
      all_of(retain_analysis_cols)
    )
}

# Index dataset at the prediction horizon
# Complete-case eligibility is based on the primary eight-biomarker panel;
# derived ratios are retained for sensitivity analyses.

build_index_dataset_at_horizon <- function(analysis_df,
                                           subject_followup,
                                           complete_case_cols,
                                           predict_time,
                                           retain_analysis_cols = complete_case_cols) {
  
  subject_followup %>%
    mutate(horizon_date = First_LD_DATE %m+% months(predict_time)) %>%
    transmute(
      Subject_ID,
      first_ld_date_ref = First_LD_DATE,
      horizon_date = horizon_date,
      duration = duration,
      hcc_status = hcc_status
    ) %>%
    inner_join(analysis_df, by = "Subject_ID") %>%
    filter(
      !is.na(RESULT_MONTH),
      RESULT_MONTH >= first_ld_date_ref,
      RESULT_MONTH <= horizon_date,
      !is.na(AGE),
      !is.na(SEX),
      !is.na(RACE),
      !is.na(HISPANIC),
      !is.na(BMI_final),
      if_all(all_of(complete_case_cols), ~ !is.na(.))
    ) %>%
    arrange(Subject_ID, RESULT_MONTH) %>%
    group_by(Subject_ID) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(
      Subject_ID,
      duration,
      hcc_status,
      AGE, SEX, RACE, HISPANIC, BMI_final,
      all_of(retain_analysis_cols)
    )
}

# Cox model estimation

fit_final_tvcox_model <- function(tv_df,
                                  selected_analysis_cols) {
  
  rhs_terms <- c("AGE", "SEX", "RACE", "HISPANIC", "BMI_final", selected_analysis_cols, "cluster(Subject_ID)")
  
  cox_formula <- as.formula(
    paste0("Surv(tstart, tstop, event) ~ ", paste(rhs_terms, collapse = " + "))
  )
  
  final_model <- coxph(
    formula = cox_formula,
    data = tv_df,
    ties = "efron",
    x = TRUE,
    robust = TRUE
  )
  
  final_summary <- summary(final_model)
  
  coefficient_table <- tibble::tibble(
    term = rownames(final_summary$coefficients),
    beta = final_summary$coefficients[, "coef"],
    hr = final_summary$coefficients[, "exp(coef)"],
    robust_se = final_summary$coefficients[, "robust se"],
    z = final_summary$coefficients[, "z"],
    p_value = final_summary$coefficients[, "Pr(>|z|)"]
  )
  
  list(
    final_model = final_model,
    coefficient_table = coefficient_table
  )
}

fit_final_baseline_cox_model <- function(baseline_df,
                                         selected_analysis_cols) {
  
  rhs_terms <- c("AGE", "SEX", "RACE", "HISPANIC", "BMI_final", selected_analysis_cols)
  
  cox_formula <- as.formula(
    paste0("Surv(duration, hcc_status) ~ ", paste(rhs_terms, collapse = " + "))
  )
  
  final_model <- coxph(
    formula = cox_formula,
    data = baseline_df,
    ties = "efron",
    x = TRUE
  )
  
  final_summary <- summary(final_model)
  
  coefficient_table <- tibble::tibble(
    term = rownames(final_summary$coefficients),
    beta = final_summary$coefficients[, "coef"],
    hr = final_summary$coefficients[, "exp(coef)"],
    se = final_summary$coefficients[, "se(coef)"],
    z = final_summary$coefficients[, "z"],
    p_value = final_summary$coefficients[, "Pr(>|z|)"]
  )
  
  list(
    final_model = final_model,
    coefficient_table = coefficient_table
  )
}

# Bootstrap model estimation

fit_weighted_tvcox_bootstrap_model <- function(tv_df,
                                               selected_analysis_cols,
                                               bootstrap_weights) {
  
  weighted_tv_df <- tv_df %>%
    left_join(bootstrap_weights, by = "Subject_ID") %>%
    mutate(
      boot_weight = dplyr::coalesce(boot_weight, 0L)
    ) %>%
    filter(boot_weight > 0)
  
  rhs_terms <- c("AGE", "SEX", "RACE", "HISPANIC", "BMI_final", selected_analysis_cols)
  
  cox_formula <- as.formula(
    paste0("Surv(tstart, tstop, event) ~ ", paste(rhs_terms, collapse = " + "))
  )
  
  coxph(
    formula = cox_formula,
    data = weighted_tv_df,
    ties = "efron",
    weights = boot_weight,
    x = TRUE
  )
}

fit_weighted_baseline_cox_bootstrap_model <- function(baseline_df,
                                                      selected_analysis_cols,
                                                      bootstrap_weights) {
  
  weighted_baseline_df <- baseline_df %>%
    left_join(bootstrap_weights, by = "Subject_ID") %>%
    mutate(
      boot_weight = dplyr::coalesce(boot_weight, 0L)
    ) %>%
    filter(boot_weight > 0)
  
  rhs_terms <- c("AGE", "SEX", "RACE", "HISPANIC", "BMI_final", selected_analysis_cols)
  
  cox_formula <- as.formula(
    paste0("Surv(duration, hcc_status) ~ ", paste(rhs_terms, collapse = " + "))
  )
  
  coxph(
    formula = cox_formula,
    data = weighted_baseline_df,
    ties = "efron",
    weights = boot_weight,
    x = TRUE
  )
}

# Analysis cohorts

prepare_common_modeling_objects <- function(prepared_obj,
                                            primary_longitudinal_biomarkers,
                                            predict_time) {
  
  sensitivity_biomarkers <- c("AST_ALT_RATIO", "NLR", "PLR")
  all_model_biomarkers <- c(primary_longitudinal_biomarkers, sensitivity_biomarkers)
  
  analysis_data <- append_model_biomarkers(prepared_obj$analysis_data)
  
  primary_analysis_cols <- paste0("bm_", primary_longitudinal_biomarkers)
  retain_analysis_cols <- paste0("bm_", all_model_biomarkers)
  
  tv_df_raw <- build_timevarying_dataset(
    analysis_df = analysis_data,
    subject_followup = prepared_obj$subject_followup,
    complete_case_cols = primary_analysis_cols,
    retain_analysis_cols = retain_analysis_cols
  )
  
  index_df <- build_index_dataset_at_horizon(
    analysis_df = analysis_data,
    subject_followup = prepared_obj$subject_followup,
    complete_case_cols = primary_analysis_cols,
    predict_time = predict_time,
    retain_analysis_cols = retain_analysis_cols
  )
  
  analytic_subject_ids <- unique(index_df$Subject_ID)
  
  tv_df <- tv_df_raw %>%
    filter(Subject_ID %in% analytic_subject_ids)

  baseline_augmented_data <- append_model_biomarkers(prepared_obj$baseline_data)
  baseline_complete_data <- baseline_augmented_data %>%
    filter(
      !is.na(AGE),
      !is.na(SEX),
      !is.na(RACE),
      !is.na(HISPANIC),
      !is.na(BMI_final),
      if_all(all_of(primary_analysis_cols), ~ !is.na(.))
    ) %>%
    select(
      Subject_ID,
      duration,
      hcc_status,
      AGE, SEX, RACE, HISPANIC, BMI_final,
      all_of(retain_analysis_cols)
    )
  
  list(
    primary_longitudinal_biomarkers = primary_longitudinal_biomarkers,
    selected_analysis_cols = primary_analysis_cols,
    tv_data = tv_df,
    index_data = index_df,
    baseline_complete_data = baseline_complete_data
  )
}

# Joint bootstrap validation in the time-varying Cox cohort

run_joint_tvcox_harrell_bootstrap_validation <- function(tv_df,
                                                         index_df,
                                                         model_spec_list,
                                                         primary_longitudinal_biomarkers,
                                                         predict_time,
                                                         n_bootstrap,
                                                         bootstrap_seed) {
  subject_ids <- unique(index_df$Subject_ID)
  
  bootstrap_detail <- future.apply::future_lapply(
    seq_len(n_bootstrap),
    function(b) {
      bootstrap_weights <- build_subject_bootstrap_weights(subject_ids)
      
      purrr::map_dfr(model_spec_list, function(model_spec) {
      selected_analysis_cols <- paste0("bm_", model_spec)
      model_key <- paste(model_spec, collapse = "||")
      model_label <- format_model_spec_label(
        input_biomarker = model_spec,
        primary_longitudinal_biomarkers = primary_longitudinal_biomarkers
      )
      
      model_index_df <- index_df %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      model_subject_ids <- unique(model_index_df$Subject_ID)
      
      model_tv_df <- tv_df %>%
        filter(Subject_ID %in% model_subject_ids) %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      apparent_eval_df <- expand_subject_level_dataset_by_weights(
        eval_df = model_index_df,
        bootstrap_weights = bootstrap_weights
      )
      
      boot_model <- fit_weighted_tvcox_bootstrap_model(
        tv_df = model_tv_df,
        selected_analysis_cols = selected_analysis_cols,
        bootstrap_weights = bootstrap_weights
      )
      
      apparent_lp <- as.numeric(
        predict(boot_model, newdata = apparent_eval_df, type = "lp")
      )
      
      test_lp <- as.numeric(
        predict(boot_model, newdata = model_index_df, type = "lp")
      )
      
      apparent_perf <- compute_performance_metrics(
        eval_df = apparent_eval_df,
        risk_score = apparent_lp,
        predict_time = predict_time
      )
      
      test_perf <- compute_performance_metrics(
        eval_df = model_index_df,
        risk_score = test_lp,
        predict_time = predict_time
      )
      
      tibble::tibble(
        bootstrap_id = b,
        model_key = model_key,
        model_label = model_label,
        optimism_c_index = apparent_perf$c_index - test_perf$c_index,
        optimism_auc = apparent_perf$auc - test_perf$auc
      )
      })
    },
    future.seed = bootstrap_seed
  ) %>%
    bind_rows()
  
  optimism_summary_table <- bootstrap_detail %>%
    group_by(model_key, model_label) %>%
    summarise(
      optimism_c_index = mean(optimism_c_index),
      optimism_auc = mean(optimism_auc),
      .groups = "drop"
    )
  
  list(
    n_bootstrap = n_bootstrap,
    bootstrap_detail = bootstrap_detail,
    optimism_summary_table = optimism_summary_table
  )
}

# Bootstrap validation in the baseline Cox cohort

run_joint_baseline_cox_harrell_bootstrap_validation <- function(baseline_df,
                                                                model_spec_list,
                                                                primary_longitudinal_biomarkers,
                                                                predict_time,
                                                                n_bootstrap,
                                                                bootstrap_seed) {
  subject_ids <- unique(baseline_df$Subject_ID)
  
  bootstrap_detail <- future.apply::future_lapply(
    seq_len(n_bootstrap),
    function(b) {
      bootstrap_weights <- build_subject_bootstrap_weights(subject_ids)
      
      purrr::map_dfr(model_spec_list, function(model_spec) {
      selected_analysis_cols <- paste0("bm_", model_spec)
      model_key <- paste(model_spec, collapse = "||")
      model_label <- format_model_spec_label(
        input_biomarker = model_spec,
        primary_longitudinal_biomarkers = primary_longitudinal_biomarkers
      )
      
      model_baseline_df <- baseline_df %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      apparent_eval_df <- expand_subject_level_dataset_by_weights(
        eval_df = model_baseline_df,
        bootstrap_weights = bootstrap_weights
      )
      
      boot_model <- fit_weighted_baseline_cox_bootstrap_model(
        baseline_df = model_baseline_df,
        selected_analysis_cols = selected_analysis_cols,
        bootstrap_weights = bootstrap_weights
      )
      
      apparent_lp <- as.numeric(
        predict(boot_model, newdata = apparent_eval_df, type = "lp")
      )
      
      test_lp <- as.numeric(
        predict(boot_model, newdata = model_baseline_df, type = "lp")
      )
      
      apparent_perf <- compute_performance_metrics(
        eval_df = apparent_eval_df,
        risk_score = apparent_lp,
        predict_time = predict_time
      )
      
      test_perf <- compute_performance_metrics(
        eval_df = model_baseline_df,
        risk_score = test_lp,
        predict_time = predict_time
      )
      
      tibble::tibble(
        bootstrap_id = b,
        model_key = model_key,
        model_label = model_label,
        optimism_c_index = apparent_perf$c_index - test_perf$c_index,
        optimism_auc = apparent_perf$auc - test_perf$auc
      )
      })
    },
    future.seed = bootstrap_seed
  ) %>%
    bind_rows()
  
  optimism_summary_table <- bootstrap_detail %>%
    group_by(model_key, model_label) %>%
    summarise(
      optimism_c_index = mean(optimism_c_index),
      optimism_auc = mean(optimism_auc),
      .groups = "drop"
    )
  
  list(
    n_bootstrap = n_bootstrap,
    bootstrap_detail = bootstrap_detail,
    optimism_summary_table = optimism_summary_table
  )
}

# Bootstrap comparison in the aligned cohort

run_joint_matched_harrell_bootstrap_comparison <- function(matched_baseline_data,
                                                           matched_tv_data,
                                                           matched_index_data,
                                                           model_spec_list,
                                                           primary_longitudinal_biomarkers,
                                                           predict_time,
                                                           n_bootstrap,
                                                           bootstrap_seed) {
  subject_ids <- unique(matched_index_data$Subject_ID)
  
  bootstrap_detail <- future.apply::future_lapply(
    seq_len(n_bootstrap),
    function(b) {
      bootstrap_weights <- build_subject_bootstrap_weights(subject_ids)
      
      purrr::map_dfr(model_spec_list, function(model_spec) {
      selected_analysis_cols <- paste0("bm_", model_spec)
      model_key <- paste(model_spec, collapse = "||")
      model_label <- format_model_spec_label(
        input_biomarker = model_spec,
        primary_longitudinal_biomarkers = primary_longitudinal_biomarkers
      )
      
      model_baseline_df <- matched_baseline_data %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      model_index_df <- matched_index_data %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      model_subject_ids <- intersect(
        model_baseline_df %>% distinct(Subject_ID) %>% pull(Subject_ID),
        model_index_df %>% distinct(Subject_ID) %>% pull(Subject_ID)
      )
      
      model_baseline_df <- model_baseline_df %>%
        filter(Subject_ID %in% model_subject_ids)
      
      model_index_df <- model_index_df %>%
        filter(Subject_ID %in% model_subject_ids)
      
      model_tv_df <- matched_tv_data %>%
        filter(Subject_ID %in% model_subject_ids) %>%
        filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
      
      baseline_apparent_eval_df <- expand_subject_level_dataset_by_weights(
        eval_df = model_baseline_df,
        bootstrap_weights = bootstrap_weights
      )
      
      tv_apparent_eval_df <- expand_subject_level_dataset_by_weights(
        eval_df = model_index_df,
        bootstrap_weights = bootstrap_weights
      )
      
      baseline_boot_model <- fit_weighted_baseline_cox_bootstrap_model(
        baseline_df = model_baseline_df,
        selected_analysis_cols = selected_analysis_cols,
        bootstrap_weights = bootstrap_weights
      )
      
      tv_boot_model <- fit_weighted_tvcox_bootstrap_model(
        tv_df = model_tv_df,
        selected_analysis_cols = selected_analysis_cols,
        bootstrap_weights = bootstrap_weights
      )
      
      baseline_apparent_lp <- as.numeric(
        predict(baseline_boot_model, newdata = baseline_apparent_eval_df, type = "lp")
      )
      
      baseline_test_lp <- as.numeric(
        predict(baseline_boot_model, newdata = model_baseline_df, type = "lp")
      )
      
      tv_apparent_lp <- as.numeric(
        predict(tv_boot_model, newdata = tv_apparent_eval_df, type = "lp")
      )
      
      tv_test_lp <- as.numeric(
        predict(tv_boot_model, newdata = model_index_df, type = "lp")
      )
      
      baseline_apparent_perf <- compute_performance_metrics(
        eval_df = baseline_apparent_eval_df,
        risk_score = baseline_apparent_lp,
        predict_time = predict_time
      )
      
      baseline_test_perf <- compute_performance_metrics(
        eval_df = model_baseline_df,
        risk_score = baseline_test_lp,
        predict_time = predict_time
      )
      
      tv_apparent_perf <- compute_performance_metrics(
        eval_df = tv_apparent_eval_df,
        risk_score = tv_apparent_lp,
        predict_time = predict_time
      )
      
      tv_test_perf <- compute_performance_metrics(
        eval_df = model_index_df,
        risk_score = tv_test_lp,
        predict_time = predict_time
      )
      
      tibble::tibble(
        bootstrap_id = b,
        model_key = model_key,
        model_label = model_label,
        baseline_optimism_c_index = baseline_apparent_perf$c_index - baseline_test_perf$c_index,
        baseline_optimism_auc = baseline_apparent_perf$auc - baseline_test_perf$auc,
        tv_optimism_c_index = tv_apparent_perf$c_index - tv_test_perf$c_index,
        tv_optimism_auc = tv_apparent_perf$auc - tv_test_perf$auc
      )
      })
    },
    future.seed = bootstrap_seed
  ) %>%
    bind_rows()
  
  optimism_summary_table <- bootstrap_detail %>%
    group_by(model_key, model_label) %>%
    summarise(
      baseline_optimism_c_index = mean(baseline_optimism_c_index),
      baseline_optimism_auc = mean(baseline_optimism_auc),
      tv_optimism_c_index = mean(tv_optimism_c_index),
      tv_optimism_auc = mean(tv_optimism_auc),
      .groups = "drop"
    )
  
  list(
    n_bootstrap = n_bootstrap,
    bootstrap_detail = bootstrap_detail,
    optimism_summary_table = optimism_summary_table
  )
}

# Time-varying Cox model estimation
# Derived-ratio models exclude rows for which the ratio is undefined.

run_tvcox_model_suite <- function(common_obj,
                                  model_spec_list,
                                  predict_time,
                                  n_bootstrap,
                                  bootstrap_seed) {
  
  result_list <- purrr::map(model_spec_list, function(model_spec) {
    selected_analysis_cols <- paste0("bm_", model_spec)
    model_label <- format_model_spec_label(
      input_biomarker = model_spec,
      primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers
    )
    
    model_index_df <- common_obj$index_data %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    analytic_subject_ids <- unique(model_index_df$Subject_ID)
    
    model_tv_df <- common_obj$tv_data %>%
      filter(Subject_ID %in% analytic_subject_ids) %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    final_fit <- fit_final_tvcox_model(
      tv_df = model_tv_df,
      selected_analysis_cols = selected_analysis_cols
    )
    
    apparent_lp <- as.numeric(
      predict(final_fit$final_model, newdata = model_index_df, type = "lp")
    )
    
    apparent_info <- compute_performance_metrics(
      eval_df = model_index_df,
      risk_score = apparent_lp,
      predict_time = predict_time
    )
    
    list(
      model_label = model_label,
      model_spec = model_spec,
      selected_analysis_cols = selected_analysis_cols,
      tv_data = model_tv_df,
      index_data = model_index_df,
      model_result = list(
        final_model = final_fit$final_model,
        coefficient_table = final_fit$coefficient_table,
        apparent_performance = apparent_info,
        corrected_performance = NULL,
        bootstrap_validation = NULL
      )
    )
  })
  
  joint_bootstrap <- run_joint_tvcox_harrell_bootstrap_validation(
    tv_df = common_obj$tv_data,
    index_df = common_obj$index_data,
    model_spec_list = model_spec_list,
    primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers,
    predict_time = predict_time,
    n_bootstrap = n_bootstrap,
    bootstrap_seed = bootstrap_seed
  )
  
  result_list <- purrr::map(result_list, function(result_obj) {
    model_key <- paste(result_obj$model_spec, collapse = "||")
    
    optimism_row <- joint_bootstrap$optimism_summary_table %>%
      filter(model_key == !!model_key)
    
    optimism_summary <- as.list(
      optimism_row %>%
        select(
          optimism_c_index,
          optimism_auc
        ) %>%
        slice(1)
    )
    
    result_obj$model_result$corrected_performance <- compute_corrected_performance(
      apparent_performance = result_obj$model_result$apparent_performance,
      optimism_summary = optimism_summary
    )
    
    result_obj$model_result$bootstrap_validation <- list(
      n_bootstrap = joint_bootstrap$n_bootstrap,
      bootstrap_detail = joint_bootstrap$bootstrap_detail %>%
        filter(model_key == !!model_key) %>%
        select(-model_key, -model_label),
      optimism_summary = optimism_summary
    )
    
    result_obj
  })
  
  result_list
}

# Baseline Cox model estimation
# Derived-ratio models exclude rows for which the ratio is undefined.

run_baseline_cox_model_suite <- function(common_obj,
                                         model_spec_list,
                                         predict_time,
                                         n_bootstrap,
                                         bootstrap_seed) {
  
  result_list <- purrr::map(model_spec_list, function(model_spec) {
    selected_analysis_cols <- paste0("bm_", model_spec)
    model_label <- format_model_spec_label(
      input_biomarker = model_spec,
      primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers
    )
    
    model_baseline_df <- common_obj$baseline_complete_data %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    final_fit <- fit_final_baseline_cox_model(
      baseline_df = model_baseline_df,
      selected_analysis_cols = selected_analysis_cols
    )
    
    apparent_lp <- as.numeric(
      predict(final_fit$final_model, newdata = model_baseline_df, type = "lp")
    )
    
    apparent_info <- compute_performance_metrics(
      eval_df = model_baseline_df,
      risk_score = apparent_lp,
      predict_time = predict_time
    )
    
    list(
      model_label = model_label,
      model_spec = model_spec,
      selected_analysis_cols = selected_analysis_cols,
      baseline_complete_data = model_baseline_df,
      model_result = list(
        final_model = final_fit$final_model,
        coefficient_table = final_fit$coefficient_table,
        apparent_performance = apparent_info,
        corrected_performance = NULL,
        bootstrap_validation = NULL
      )
    )
  })
  
  joint_bootstrap <- run_joint_baseline_cox_harrell_bootstrap_validation(
    baseline_df = common_obj$baseline_complete_data,
    model_spec_list = model_spec_list,
    primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers,
    predict_time = predict_time,
    n_bootstrap = n_bootstrap,
    bootstrap_seed = bootstrap_seed
  )
  
  result_list <- purrr::map(result_list, function(result_obj) {
    model_key <- paste(result_obj$model_spec, collapse = "||")
    
    optimism_row <- joint_bootstrap$optimism_summary_table %>%
      filter(model_key == !!model_key)
    
    optimism_summary <- as.list(
      optimism_row %>%
        select(
          optimism_c_index,
          optimism_auc
        ) %>%
        slice(1)
    )
    
    result_obj$model_result$corrected_performance <- compute_corrected_performance(
      apparent_performance = result_obj$model_result$apparent_performance,
      optimism_summary = optimism_summary
    )
    
    result_obj$model_result$bootstrap_validation <- list(
      n_bootstrap = joint_bootstrap$n_bootstrap,
      bootstrap_detail = joint_bootstrap$bootstrap_detail %>%
        filter(model_key == !!model_key) %>%
        select(-model_key, -model_label),
      optimism_summary = optimism_summary
    )
    
    result_obj
  })
  
  result_list
}

# Baseline and time-varying Cox comparisons in the aligned cohort

run_matched_model_comparison_suite <- function(common_obj,
                                               model_spec_list,
                                               predict_time,
                                               n_bootstrap,
                                               bootstrap_seed) {
  
  common_matched_subject_ids <- intersect(
    common_obj$baseline_complete_data %>% distinct(Subject_ID) %>% pull(Subject_ID),
    common_obj$index_data %>% distinct(Subject_ID) %>% pull(Subject_ID)
  )
  
  matched_baseline_data_common <- common_obj$baseline_complete_data %>%
    filter(Subject_ID %in% common_matched_subject_ids)
  
  matched_tv_data_common <- common_obj$tv_data %>%
    filter(Subject_ID %in% common_matched_subject_ids)
  
  matched_index_data_common <- common_obj$index_data %>%
    filter(Subject_ID %in% common_matched_subject_ids)
  
  result_list <- purrr::map(model_spec_list, function(model_spec) {
    selected_analysis_cols <- paste0("bm_", model_spec)
    model_label <- format_model_spec_label(
      input_biomarker = model_spec,
      primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers
    )
    
    matched_baseline_data <- matched_baseline_data_common %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    matched_index_data <- matched_index_data_common %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    matched_subject_ids <- intersect(
      matched_baseline_data %>% distinct(Subject_ID) %>% pull(Subject_ID),
      matched_index_data %>% distinct(Subject_ID) %>% pull(Subject_ID)
    )
    
    matched_baseline_data <- matched_baseline_data %>%
      filter(Subject_ID %in% matched_subject_ids)
    
    matched_index_data <- matched_index_data %>%
      filter(Subject_ID %in% matched_subject_ids) %>%
      mutate(
        horizon_evaluable = duration >= predict_time | (hcc_status == 1L & duration < predict_time),
        incident_hcc_by_horizon = hcc_status == 1L & duration <= predict_time
      )
    
    matched_tv_data <- matched_tv_data_common %>%
      filter(Subject_ID %in% matched_subject_ids) %>%
      filter(if_all(all_of(selected_analysis_cols), ~ !is.na(.)))
    
    matched_summary <- list(
      matched_subject_n = length(matched_subject_ids),
      matched_event_n = sum(matched_index_data$hcc_status == 1L),
      matched_horizon_event_n = sum(matched_index_data$incident_hcc_by_horizon),
      matched_evaluable_n = sum(matched_index_data$horizon_evaluable)
    )
    
    baseline_final_fit <- fit_final_baseline_cox_model(
      baseline_df = matched_baseline_data,
      selected_analysis_cols = selected_analysis_cols
    )
    
    baseline_apparent_lp <- as.numeric(
      predict(baseline_final_fit$final_model, newdata = matched_baseline_data, type = "lp")
    )
    
    baseline_apparent_info <- compute_performance_metrics(
      eval_df = matched_baseline_data,
      risk_score = baseline_apparent_lp,
      predict_time = predict_time
    )
    
    matched_tv_final_fit <- fit_final_tvcox_model(
      tv_df = matched_tv_data,
      selected_analysis_cols = selected_analysis_cols
    )
    
    matched_tv_apparent_lp <- as.numeric(
      predict(
        matched_tv_final_fit$final_model,
        newdata = matched_index_data %>%
          select(-horizon_evaluable, -incident_hcc_by_horizon),
        type = "lp"
      )
    )
    
    matched_tv_apparent_info <- compute_performance_metrics(
      eval_df = matched_index_data %>% select(-horizon_evaluable, -incident_hcc_by_horizon),
      risk_score = matched_tv_apparent_lp,
      predict_time = predict_time
    )
    
    list(
      model_label = model_label,
      model_spec = model_spec,
      matched_subject_ids = matched_subject_ids,
      matched_summary = matched_summary,
      matched_baseline_data = matched_baseline_data,
      matched_tv_data = matched_tv_data,
      matched_index_data = matched_index_data,
      baseline_cox_result = list(
        final_model = baseline_final_fit$final_model,
        coefficient_table = baseline_final_fit$coefficient_table,
        apparent_performance = baseline_apparent_info,
        corrected_performance = NULL
      ),
      matched_tvcox_result = list(
        final_model = matched_tv_final_fit$final_model,
        coefficient_table = matched_tv_final_fit$coefficient_table,
        apparent_performance = matched_tv_apparent_info,
        corrected_performance = NULL
      ),
      bootstrap_validation = NULL
    )
  })
  
  joint_bootstrap <- run_joint_matched_harrell_bootstrap_comparison(
    matched_baseline_data = matched_baseline_data_common,
    matched_tv_data = matched_tv_data_common,
    matched_index_data = matched_index_data_common,
    model_spec_list = model_spec_list,
    primary_longitudinal_biomarkers = common_obj$primary_longitudinal_biomarkers,
    predict_time = predict_time,
    n_bootstrap = n_bootstrap,
    bootstrap_seed = bootstrap_seed
  )
  
  result_list <- purrr::map(result_list, function(result_obj) {
    model_key <- paste(result_obj$model_spec, collapse = "||")
    
    optimism_row <- joint_bootstrap$optimism_summary_table %>%
      filter(model_key == !!model_key)
    
    baseline_optimism_summary <- as.list(
      optimism_row %>%
        select(
          baseline_optimism_c_index,
          baseline_optimism_auc
        ) %>%
        slice(1)
    )
    
    tv_optimism_summary <- as.list(
      optimism_row %>%
        select(
          tv_optimism_c_index,
          tv_optimism_auc
        ) %>%
        slice(1)
    )
    
    result_obj$baseline_cox_result$corrected_performance <- compute_corrected_performance(
      apparent_performance = result_obj$baseline_cox_result$apparent_performance,
      optimism_summary = baseline_optimism_summary,
      prefix = "baseline_optimism_"
    )
    
    result_obj$matched_tvcox_result$corrected_performance <- compute_corrected_performance(
      apparent_performance = result_obj$matched_tvcox_result$apparent_performance,
      optimism_summary = tv_optimism_summary,
      prefix = "tv_optimism_"
    )
    
    result_obj$bootstrap_validation <- list(
      n_bootstrap = joint_bootstrap$n_bootstrap,
      bootstrap_detail = joint_bootstrap$bootstrap_detail %>%
        filter(model_key == !!model_key) %>%
        select(-model_key, -model_label),
      optimism_summary = as.list(optimism_row[1, ])
    )
    
    result_obj
  })
  
  result_list
}

# Summary tables

build_bootstrap_validation_summary_table <- function(result_list, model_type) {
  purrr::map_dfr(result_list, function(result_obj) {
    tibble::tibble(
      `Model type` = model_type,
      `Biomarker set` = result_obj$model_label,
      `Bootstrap replicates` = result_obj$model_result$bootstrap_validation$n_bootstrap,
      `Apparent C-index` = result_obj$model_result$apparent_performance$c_index,
      `Optimism-corrected C-index` = result_obj$model_result$corrected_performance$c_index,
      `Apparent AUC at horizon` = result_obj$model_result$apparent_performance$auc,
      `Optimism-corrected AUC at horizon` = result_obj$model_result$corrected_performance$auc,
      `Evaluable-at-horizon subjects` = result_obj$model_result$apparent_performance$evaluable_n,
      `Horizon events` = result_obj$model_result$apparent_performance$horizon_events
    )
  })
}

build_matched_performance_comparison_table <- function(comparison_result_list) {
  purrr::map_dfr(comparison_result_list, function(comparison_result) {
    tibble::tibble(
      `Model type` = c("Baseline Cox", "Time-varying Cox"),
      `Biomarker set` = c(comparison_result$model_label, comparison_result$model_label),
      `Bootstrap replicates` = c(
        comparison_result$bootstrap_validation$n_bootstrap,
        comparison_result$bootstrap_validation$n_bootstrap
      ),
      `Apparent C-index` = c(
        comparison_result$baseline_cox_result$apparent_performance$c_index,
        comparison_result$matched_tvcox_result$apparent_performance$c_index
      ),
      `Optimism-corrected C-index` = c(
        comparison_result$baseline_cox_result$corrected_performance$c_index,
        comparison_result$matched_tvcox_result$corrected_performance$c_index
      ),
      `Apparent AUC at horizon` = c(
        comparison_result$baseline_cox_result$apparent_performance$auc,
        comparison_result$matched_tvcox_result$apparent_performance$auc
      ),
      `Optimism-corrected AUC at horizon` = c(
        comparison_result$baseline_cox_result$corrected_performance$auc,
        comparison_result$matched_tvcox_result$corrected_performance$auc
      ),
      `Evaluable-at-horizon subjects` = c(
        comparison_result$baseline_cox_result$apparent_performance$evaluable_n,
        comparison_result$matched_tvcox_result$apparent_performance$evaluable_n
      ),
      `Horizon events` = c(
        comparison_result$baseline_cox_result$apparent_performance$horizon_events,
        comparison_result$matched_tvcox_result$apparent_performance$horizon_events
      )
    )
  })
}

# Hazard ratios for the primary joint models in the aligned cohort

build_model_coefficients_table <- function(matched_joint_result) {
  variable_labels <- c(
    AGE = "Age",
    SEXMale = "Male sex",
    `RACEBlack or African American` = "Race: Black or African American",
    RACEOther = "Race: Other",
    RACEWhite = "Race: White",
    HISPANICYes = "Hispanic ethnicity",
    BMI_final = "BMI",
    bm_ALT = "ALT",
    bm_AST = "AST",
    bm_ALB = "Albumin",
    bm_BILI = "Total bilirubin",
    bm_PLAT = "Platelet count",
    bm_NEUT = "Neutrophil count",
    bm_LYMPH = "Lymphocyte count",
    bm_HGB = "Hemoglobin"
  )

  extract_coefs <- function(model_fit, model_type) {
    model_summary <- summary(model_fit)
    coefficients <- model_summary$coefficients
    confidence_intervals <- model_summary$conf.int
    p_values <- as.numeric(coefficients[, "Pr(>|z|)"])

    tibble::tibble(
      `Model type` = model_type,
      Variable = recode(rownames(coefficients), !!!variable_labels),
      `HR (95% CI)` = sprintf(
        "%.2f (%.2f-%.2f)",
        exp(coefficients[, "coef"]),
        confidence_intervals[, "lower .95"],
        confidence_intervals[, "upper .95"]
      ),
      `P value` = if_else(
        p_values < 0.001,
        "<.001",
        sub("^0", "", sprintf("%.3f", p_values))
      )
    )
  }

  bind_rows(
    extract_coefs(matched_joint_result$baseline_cox_result$final_model, "Baseline Cox"),
    extract_coefs(matched_joint_result$matched_tvcox_result$final_model, "Time-varying Cox")
  ) %>%
    pivot_wider(
      names_from = `Model type`,
      values_from = c(`HR (95% CI)`, `P value`),
      names_glue = "{`Model type`} {.value}"
    ) %>%
    select(
      Variable,
      `Baseline Cox HR (95% CI)`,
      `Baseline Cox P value`,
      `Time-varying Cox HR (95% CI)`,
      `Time-varying Cox P value`
    )
}

# Apparent calibration of the aligned baseline model

build_calibration_results <- function(matched_joint_result, predict_time) {
  calibration_model <- matched_joint_result$baseline_cox_result$final_model
  calibration_data <- matched_joint_result$matched_baseline_data
  baseline_hazard <- basehaz(calibration_model, centered = FALSE)
  hazard_at_horizon <- tail(
    baseline_hazard$hazard[baseline_hazard$time <= predict_time],
    1L
  )

  calibration_data <- calibration_data %>%
    mutate(
      linear_predictor = as.numeric(
        predict(calibration_model, newdata = calibration_data, type = "lp", reference = "zero")
      ),
      predicted_risk = 1 - exp(-hazard_at_horizon * exp(linear_predictor)),
      calibration_group = ntile(predicted_risk, 4L)
    )

  calibration_table <- calibration_data %>%
    group_by(calibration_group) %>%
    group_modify(~ {
      km <- summary(
        survfit(Surv(duration, hcc_status) ~ 1, data = .x),
        times = predict_time,
        extend = TRUE
      )
      tibble(
        n = nrow(.x),
        events_by_horizon = sum(.x$hcc_status == 1L & .x$duration <= predict_time),
        mean_predicted_risk = mean(.x$predicted_risk),
        observed_risk = 1 - km$surv,
        observed_lower_95 = 1 - km$upper,
        observed_upper_95 = 1 - km$lower
      )
    }) %>%
    ungroup()

  axis_limit <- min(
    1,
    max(
      0.10,
      ceiling(1.05 * max(calibration_table$mean_predicted_risk, calibration_table$observed_risk) / 0.05) * 0.05
    )
  )

  calibration_plot <- ggplot(
    calibration_table,
    aes(x = mean_predicted_risk, y = observed_risk)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey55") +
    geom_point(shape = 21, size = 2.8, stroke = 0.7, color = "black", fill = "white") +
    scale_x_continuous(
      limits = c(0, axis_limit),
      labels = label_percent(accuracy = 1),
      breaks = breaks_pretty(n = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, axis_limit),
      labels = label_percent(accuracy = 1),
      breaks = breaks_pretty(n = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_equal() +
    labs(
      x = paste0("Predicted ", predict_time, "-month risk"),
      y = paste0("Observed ", predict_time, "-month risk")
    ) +
    theme_classic(base_size = 11)

  list(table = calibration_table, plot = calibration_plot)
}

# Proportional hazards tests for the primary joint models

build_ph_test_table <- function(tv_joint_result, baseline_joint_result) {
  term_labels <- c(
    AGE = "Age",
    SEX = "Sex",
    RACE = "Race",
    HISPANIC = "Ethnicity",
    BMI_final = "BMI",
    bm_ALT = "ALT",
    bm_AST = "AST",
    bm_ALB = "Albumin",
    bm_BILI = "Total bilirubin",
    bm_PLAT = "Platelet count",
    bm_NEUT = "Neutrophil count",
    bm_LYMPH = "Lymphocyte count",
    bm_HGB = "Hemoglobin",
    GLOBAL = "Overall test"
  )

  extract_test <- function(model, model_label) {
    as_tibble(
      cox.zph(model, transform = "km", terms = TRUE, singledf = FALSE, global = TRUE)$table,
      rownames = "term"
    ) %>%
      transmute(
        Model = model_label,
        `Test term` = recode(term, !!!term_labels, .default = term),
        `Chi-square` = chisq,
        df,
        `P value` = p
      )
  }

  bind_rows(
    extract_test(baseline_joint_result$model_result$final_model, "Baseline Cox"),
    extract_test(tv_joint_result$model_result$final_model, "Time-varying Cox")
  )
}

# Illustration of counting-process record construction

build_timevarying_record_figure <- function(tv_data) {
  example_id <- tv_data %>%
    group_by(Subject_ID) %>%
    summarise(
      intervals = n(),
      event = any(event == 1L),
      end_time = max(tstop),
      .groups = "drop"
    ) %>%
    filter(intervals == 4L, event) %>%
    arrange(end_time, Subject_ID) %>%
    slice(1L) %>%
    pull(Subject_ID)

  example_data <- tv_data %>%
    filter(Subject_ID == example_id) %>%
    arrange(tstart) %>%
    mutate(
      Interval = row_number(),
      `Visit month` = format(RESULT_MONTH, "%Y.%m"),
      `Start month` = round(tstart),
      `Stop month` = round(tstop),
      `Platelet count (10^3/uL)` = round(bm_PLAT),
      `AST (U/L)` = round(bm_AST),
      `HCC event` = event,
      visit_label = paste0(
        "Lab visit ", Interval,
        "\n", `Visit month`,
        "\nPlatelets = ", `Platelet count (10^3/uL)`,
        "\nAST = ", `AST (U/L)`
      )
    )

  record_table <- example_data %>%
    select(
      Interval,
      `Visit month`,
      `Start month`,
      `Stop month`,
      `Platelet count (10^3/uL)`,
      `AST (U/L)`,
      `HCC event`
    )

  timeline_plot <- ggplot(example_data, aes(x = tstart, y = 0)) +
    annotate(
      "segment",
      x = min(example_data$tstart),
      xend = max(example_data$tstop),
      y = 0,
      yend = 0,
      linewidth = 0.8,
      color = "grey35"
    ) +
    geom_point(size = 2.9, color = "grey20") +
    geom_text(aes(y = 0.24, label = visit_label), size = 3, lineheight = 0.95, vjust = 0) +
    geom_point(
      data = filter(example_data, event == 1L),
      aes(x = tstop),
      shape = 23,
      size = 3.5,
      fill = "firebrick",
      color = "firebrick"
    ) +
    geom_text(
      data = filter(example_data, event == 1L),
      aes(x = tstop, y = -0.22, label = "HCC diagnosis"),
      size = 3.1,
      fontface = "bold",
      color = "firebrick",
      vjust = 1
    ) +
    coord_cartesian(ylim = c(-0.45, 0.9), clip = "off") +
    labs(title = "Clinical timeline for one patient", x = "Months since liver disease index date", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0)
    )

  table_long <- record_table %>%
    mutate(across(everything(), as.character), row = row_number()) %>%
    pivot_longer(-row, names_to = "column", values_to = "value") %>%
    mutate(column = factor(column, levels = names(record_table)))
  header <- tibble(
    row = 0,
    column = factor(names(record_table), levels = names(record_table)),
    value = names(record_table)
  )

  record_plot <- ggplot() +
    geom_text(data = header, aes(column, -row, label = value), fontface = "bold", size = 3) +
    geom_text(data = table_long, aes(column, -row, label = value), size = 3) +
    geom_hline(yintercept = -seq(0.5, nrow(record_table) + 0.5), linewidth = 0.25, color = "grey75") +
    scale_y_continuous(limits = c(-(nrow(record_table) + 0.6), 0.5)) +
    labs(title = "Corresponding counting-process records") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0))

  list(
    table = record_table,
    plot = timeline_plot / record_plot +
      plot_layout(heights = c(1.15, 1)) +
      plot_annotation(
        title = "Illustration of time-varying laboratory record construction",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
  )
}

# Analysis settings

project_dir <- normalizePath(".", winslash = "/")
data_path <- file.path(
  project_dir,
  "data",
  "Final Clean Analysis Data with Biomarkers_v2.csv"
)
bmi_data_path <- file.path(project_dir, "data", "OneFL BMI Data - 05.08.25.csv")
output_dir <- file.path(project_dir, "code", "result_final")

predict_time <- 36
min_followup_months <- 6
n_bootstrap <- 500
bootstrap_seed <- 20260315
bootstrap_workers <- 5L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
future::plan(future::multisession, workers = bootstrap_workers)

primary_longitudinal_biomarkers <- c(
  "ALT", "AST", "ALB", "BILI",  "PLAT", "NEUT", "LYMPH", "HGB"
)

joint_model_specs <- list(
  primary_longitudinal_biomarkers,
  c(setdiff(primary_longitudinal_biomarkers, c("ALT", "AST")), "AST_ALT_RATIO"),
  c(setdiff(primary_longitudinal_biomarkers, c("NEUT", "LYMPH")), "NLR"),
  c(setdiff(primary_longitudinal_biomarkers, c("PLAT", "LYMPH")), "PLR")
)
tv_model_specs <- c(as.list(primary_longitudinal_biomarkers), joint_model_specs)
baseline_model_specs <- list(primary_longitudinal_biomarkers)
matched_model_specs <- tv_model_specs

# Data ingestion and preprocessing

prepared_obj <- read_and_prepare_data(
  data_path = data_path,
  bmi_data_path = bmi_data_path,
  min_followup_months = min_followup_months
)

# Earlier cohort-selection counts were obtained during data extraction.
cohort_selection_table <- tibble::tibble(
  Stage = c(
    "Patients with a liver disease diagnosis",
    "After excluding diagnoses confined to 2013-2014",
    "After excluding incomplete demographic information",
    "After excluding HCC before the index date",
    "After excluding patients without eligible laboratory data",
    "After excluding laboratory records outside the study period",
    "After excluding biologically implausible laboratory values",
    "After requiring at least 6 months of follow-up"
  ),
  `Excluded at stage` = c(
    NA, 3632, 275, 6, 289, 241, 42,
    nrow(prepared_obj$subject_raw_cohort) - nrow(prepared_obj$subject_followup)
  ),
  Remaining = c(
    9290, 5658, 5383, 5377, 5088, 4847,
    nrow(prepared_obj$subject_raw_cohort),
    nrow(prepared_obj$subject_followup)
  )
)

common_obj <- prepare_common_modeling_objects(
  prepared_obj = prepared_obj,
  primary_longitudinal_biomarkers = primary_longitudinal_biomarkers,
  predict_time = predict_time
)

followup_summary_table <- build_followup_summary_table(
  prepared_obj$subject_followup,
  prepared_obj$analysis_data
)

baseline_missingness_table <- build_baseline_missingness_table(prepared_obj)

longitudinal_missingness_table <- build_longitudinal_missingness_table(
  prepared_obj = prepared_obj,
  primary_longitudinal_biomarkers = primary_longitudinal_biomarkers
)

tvcox_result_list <- run_tvcox_model_suite(
  common_obj = common_obj,
  model_spec_list = tv_model_specs,
  predict_time = predict_time,
  n_bootstrap = n_bootstrap,
  bootstrap_seed = bootstrap_seed
)

evaluable_tv_subjects <- common_obj$index_data %>%
  filter(duration >= predict_time | (hcc_status == 1L & duration < predict_time))

tvcox_joint_model_summary_table <- tibble::tibble(
  Subjects = n_distinct(common_obj$tv_data$Subject_ID),
  Intervals = nrow(common_obj$tv_data),
  Events = sum(common_obj$tv_data$event),
  `Evaluable-at-horizon subjects` = nrow(evaluable_tv_subjects),
  `Horizon events` = sum(
    evaluable_tv_subjects$hcc_status == 1L &
      evaluable_tv_subjects$duration <= predict_time
  )
)

tvcox_bootstrap_validation_summary_table <- build_bootstrap_validation_summary_table(
  tvcox_result_list,
  "Time-varying Cox"
)

baseline_cox_result_list <- run_baseline_cox_model_suite(
  common_obj = common_obj,
  model_spec_list = baseline_model_specs,
  predict_time = predict_time,
  n_bootstrap = n_bootstrap,
  bootstrap_seed = bootstrap_seed + 1000L
)

baseline_joint_result <- detect(
  baseline_cox_result_list,
  ~ .x$model_label == "all biomarker"
)

baseline_cox_joint_model_summary_table <- tibble::tibble(
  Subjects = nrow(baseline_joint_result$baseline_complete_data),
  Events = sum(baseline_joint_result$baseline_complete_data$hcc_status),
  `Evaluable-at-horizon subjects` = baseline_joint_result$model_result$apparent_performance$evaluable_n,
  `Horizon events` = baseline_joint_result$model_result$apparent_performance$horizon_events
)

baseline_cox_bootstrap_validation_summary_table <- build_bootstrap_validation_summary_table(
  baseline_cox_result_list,
  "Baseline Cox"
)

matched_comparison_result_list <- run_matched_model_comparison_suite(
  common_obj = common_obj,
  model_spec_list = matched_model_specs,
  predict_time = predict_time,
  n_bootstrap = n_bootstrap,
  bootstrap_seed = bootstrap_seed + 2000L
)

matched_performance_comparison_table <- build_matched_performance_comparison_table(
  comparison_result_list = matched_comparison_result_list
)

tv_joint_result <- detect(tvcox_result_list, ~ .x$model_label == "all biomarker")
matched_joint_result <- detect(matched_comparison_result_list, ~ .x$model_label == "all biomarker")

model_display_labels <- c(
  ALT = "ALT",
  AST = "AST",
  ALB = "Albumin",
  BILI = "Total bilirubin",
  PLAT = "Platelet count",
  NEUT = "Neutrophil count",
  LYMPH = "Lymphocyte count",
  HGB = "Hemoglobin",
  `all biomarker` = "Joint eight-biomarker model",
  `all biomarker (AST/ALT ratio replaces AST + ALT)` = "AST/ALT ratio replaces AST and ALT",
  `all biomarker (NLR replaces NEUT + LYMPH)` = "NLR replaces neutrophil and lymphocyte counts",
  `all biomarker (PLR replaces PLAT + LYMPH)` = "PLR replaces platelet and lymphocyte counts"
)

table3_model_performance <- tvcox_bootstrap_validation_summary_table %>%
  filter(`Biomarker set` %in% c(primary_longitudinal_biomarkers, "all biomarker")) %>%
  mutate(`Biomarker set` = recode(`Biomarker set`, !!!model_display_labels))

table4_joint_model_performance <- bind_rows(
  baseline_cox_bootstrap_validation_summary_table %>%
    mutate(
      Cohort = "Baseline-complete cohort",
      Subjects = nrow(baseline_joint_result$baseline_complete_data),
      Events = sum(baseline_joint_result$baseline_complete_data$hcc_status)
    ),
  matched_performance_comparison_table %>%
    filter(`Biomarker set` == "all biomarker") %>%
    mutate(
      Cohort = "Aligned comparison cohort",
      Subjects = matched_joint_result$matched_summary$matched_subject_n,
      Events = matched_joint_result$matched_summary$matched_event_n
    )
) %>%
  mutate(`Biomarker set` = recode(`Biomarker set`, !!!model_display_labels)) %>%
  select(Cohort, `Model type`, `Biomarker set`, Subjects, Events, everything())

table4_single_biomarker_performance <- matched_performance_comparison_table %>%
  filter(`Biomarker set` %in% primary_longitudinal_biomarkers) %>%
  mutate(`Biomarker set` = recode(`Biomarker set`, !!!model_display_labels))

supplementary_table2_coefficients <- build_model_coefficients_table(matched_joint_result)

supplementary_table3_sensitivity <- bind_rows(
  tvcox_bootstrap_validation_summary_table %>%
    filter(grepl("all biomarker", `Biomarker set`, fixed = TRUE)) %>%
    mutate(Cohort = "Time-varying Cox cohort"),
  matched_performance_comparison_table %>%
    filter(grepl("all biomarker", `Biomarker set`, fixed = TRUE)) %>%
    mutate(Cohort = "Aligned comparison cohort")
) %>%
  mutate(`Biomarker set` = recode(`Biomarker set`, !!!model_display_labels)) %>%
  select(Cohort, `Model type`, `Biomarker set`, everything())

supplementary_table4_ph_tests <- build_ph_test_table(
  tv_joint_result = tv_joint_result,
  baseline_joint_result = baseline_joint_result
)

calibration_results <- build_calibration_results(
  matched_joint_result = matched_joint_result,
  predict_time = predict_time
)

timevarying_record_figure <- build_timevarying_record_figure(common_obj$tv_data)

output_tables <- list(
  cohort_selection = cohort_selection_table,
  cohort_followup_summary = followup_summary_table,
  baseline_biomarker_missingness = baseline_missingness_table,
  longitudinal_biomarker_missingness = longitudinal_missingness_table,
  tvcox_joint_model_summary = tvcox_joint_model_summary_table,
  baseline_cox_joint_model_summary = baseline_cox_joint_model_summary_table,
  table3_model_performance = table3_model_performance,
  table4_joint_model_performance = table4_joint_model_performance,
  table4_single_biomarker_performance = table4_single_biomarker_performance,
  supplementary_table2_coefficients = supplementary_table2_coefficients,
  supplementary_table3_sensitivity = supplementary_table3_sensitivity,
  supplementary_table4_ph_tests = supplementary_table4_ph_tests,
  supplementary_figure2_calibration_data = calibration_results$table,
  figure2_timevarying_records = timevarying_record_figure$table
)

purrr::iwalk(
  output_tables,
  ~ write.csv(
    .x,
    file.path(output_dir, paste0(.y, ".csv")),
    row.names = FALSE,
    na = ""
  )
)

sink(file.path(output_dir, "core_analysis_console.txt"))
iwalk(output_tables, function(table, table_name) {
  cat("\n", table_name, "\n", sep = "")
  print(table, n = Inf, width = Inf)
})
sink()

ggsave(
  file.path(output_dir, "Figure2_Time_Varying_Record_Construction.png"),
  timevarying_record_figure$plot,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Supplementary_Figure2_Calibration.png"),
  calibration_results$plot,
  width = 4.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

saveRDS(
  list(
    settings = list(
      predict_time = predict_time,
      min_followup_months = min_followup_months,
      n_bootstrap = n_bootstrap,
      bootstrap_seed = bootstrap_seed,
      bootstrap_workers = bootstrap_workers,
      primary_longitudinal_biomarkers = primary_longitudinal_biomarkers,
      tv_model_specs = tv_model_specs,
      baseline_model_specs = baseline_model_specs,
      matched_model_specs = matched_model_specs
    ),
    tables = output_tables,
    tvcox_results = tvcox_result_list,
    baseline_cox_results = baseline_cox_result_list,
    matched_comparison_results = matched_comparison_result_list
  ),
  file.path(output_dir, "core_analysis_results.rds")
)

future::plan(future::sequential)

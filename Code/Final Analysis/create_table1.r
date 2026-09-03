# Descriptive statistics and baseline comparisons
# Continuous variables: Median [IQR], Wilcoxon rank-sum test
# Categorical variables: n (%), Pearson chi-square test

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(lubridate)
})

biomarker_spec <- tibble::tibble(
  analysis_var = c(
    "ALT", "AST", "Bilirubin", "Albumin", "INR",
    "Platelets", "Neutrophils", "Lymphocytes", "Hemoglobin",
    "TC", "HDL", "LDL", "TG", "AFP"
  ),
  source_var = c(
    "mean_ALT_RESULT_NUM", "mean_AST_RESULT_NUM", "mean_BILI_RESULT_NUM",
    "mean_ALB_RESULT_NUM", "mean_INR_RESULT_NUM", "mean_PLAT_RESULT_NUM",
    "mean_NEUT_RESULT_NUM", "mean_LYMPH_RESULT_NUM", "mean_HGB_RESULT_NUM",
    "mean_TC_RESULT_NUM", "mean_HDL_RESULT_NUM", "mean_LDL_RESULT_NUM",
    "mean_TG_RESULT_NUM", "mean_AFP_RESULT_NUM"
  ),
  label = c(
    "ALT (U/L)", "AST (U/L)", "Total Bilirubin (mg/dL)",
    "Albumin (g/dL)", "INR", "Platelet count (10^3/uL)",
    "Neutrophil count (10^3/uL)", "Lymphocyte count (10^3/uL)",
    "Hemoglobin (g/dL)", "Total Cholesterol (mg/dL)",
    "HDL Cholesterol (mg/dL)", "LDL Cholesterol (mg/dL)",
    "Triglycerides (mg/dL)", "AFP (ng/mL)"
  ),
  table_section = c(
    rep("Liver Enzymes and Function", 5),
    rep("Haematologic Indices", 4),
    rep("Lipid Profile and Tumor Markers", 5)
  ),
  digits = c(0, 0, 2, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1)
)

longitudinal_plot_order <- c(
  "ALT", "AST", "Bilirubin", "Albumin", "INR", "Platelets",
  "Hemoglobin", "Neutrophils", "Lymphocytes", "TC", "HDL", "LDL", "TG", "AFP"
)

# Each patient contributes at most one record per calendar month. Results from
# complementary laboratory panels are combined without altering the monthly
# summary values reported in the source data.
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

format_count_pct <- function(n, denom, pct_digits = 1) {
  pct <- 100 * n / denom
  sprintf(
    "%s (%.*f%%)",
    format(n, big.mark = ",", trim = TRUE),
    pct_digits,
    pct
  )
}

format_median_iqr <- function(x, digits = 1) {
  q <- quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
  q <- formatC(
    q,
    format = "f",
    digits = digits,
    big.mark = ",",
    drop0trailing = TRUE
  )
  sprintf(
    "%s [%s, %s]",
    q[2],
    q[1],
    q[3]
  )
}

format_pvalue <- function(p) {
  if (p < 0.001) {
    return("< .001")
  }
  sub("^0", "", sprintf("%.3f", p))
}

prepare_analysis_data <- function(data_path, bmi_data_path) {
  df <- read.csv(
    file = data_path,
    stringsAsFactors = FALSE,
    check.names = TRUE,
    na.strings = c("", "NA", "NaN")
  )

  df <- df[, -1, drop = FALSE] %>%
    collapse_subject_month_records()
  
  bmi_data <- read.csv(
    file = bmi_data_path,
    stringsAsFactors = FALSE,
    check.names = TRUE,
    na.strings = c("", "NA", "NaN")
  )
  bmi_data$Subject_ID <- as.character(bmi_data$Subject_ID)
  
  df <- df %>%
    left_join(bmi_data %>% dplyr::select(Subject_ID, BMI_final), by = "Subject_ID")
  
  df <- df %>%
    mutate(
      Subject_ID     = as.character(Subject_ID),
      RESULT_MONTH   = suppressWarnings(ymd(as.character(RESULT_MONTH))),
      First_LD_DATE  = suppressWarnings(ymd(as.character(First_LD_DATE))),
      First_HCC_DATE = suppressWarnings(ymd(as.character(First_HCC_DATE)))
    )
  
  # Eligibility requires at least six elapsed calendar months between the
  # index month and the last observed laboratory month.
  eligible_patients <- df %>%
    group_by(Subject_ID) %>%
    summarise(
      index_month = floor_date(min(First_LD_DATE, na.rm = TRUE), "month"),
      last_lab_month = floor_date(max(RESULT_MONTH, na.rm = TRUE), "month"),
      .groups = "drop"
    ) %>%
    mutate(
      eligibility_months = 12 * (year(last_lab_month) - year(index_month)) +
        (month(last_lab_month) - month(index_month))
    ) %>%
    filter(eligibility_months >= 6)
    
  df <- df %>% filter(Subject_ID %in% eligible_patients$Subject_ID)

  df_baseline <- df %>%
    filter(!is.na(Subject_ID), !is.na(RESULT_MONTH)) %>%
    group_by(Subject_ID) %>%
    arrange(RESULT_MONTH, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  df_baseline <- df_baseline %>%
    mutate(
      Age = AGE,
      Sex = factor(SEX, levels = c("Female", "Male")),
      Race = factor(RACE, levels = c("White", "Black or African American", "Asian", "Other")),
      Ethnicity = factor(
        HISPANIC,
        levels = c("Yes", "No"),
        labels = c("Hispanic", "Non-Hispanic")
      ),
      Group = factor(
        HCC_group,
        levels = c("Non-HCC", "HCC"),
        labels = c("Non-HCC", "Incident HCC")
      ),
      BMI = as.numeric(BMI_final)
    )
  
  for (i in seq_len(nrow(biomarker_spec))) {
    df_baseline[[biomarker_spec$analysis_var[i]]] <- df_baseline[[biomarker_spec$source_var[i]]]
  }
  
  list(
    df = df,
    df_baseline = df_baseline,
    biomarker_spec = biomarker_spec
  )
}

append_row <- function(rows, characteristic, total = "", non_hcc = "", incident_hcc = "", p_value = "") {
  rows[[length(rows) + 1]] <- data.frame(
    Characteristic = characteristic,
    Total = total,
    Non_HCC = non_hcc,
    Incident_HCC = incident_hcc,
    P_value = p_value,
    stringsAsFactors = FALSE
  )
  rows
}

append_continuous_rows <- function(rows, data, variable, display_label, digits) {
  total_vec <- data[[variable]]
  non_hcc_vec <- data[[variable]][data$Group == "Non-HCC"]
  incident_hcc_vec <- data[[variable]][data$Group == "Incident HCC"]
  complete_rows <- !is.na(data[[variable]]) & !is.na(data$Group)
  observed_values <- data[[variable]][complete_rows]
  observed_groups <- droplevels(data$Group[complete_rows])
  p_val <- wilcox.test(observed_values ~ observed_groups, exact = FALSE)$p.value
  
  rows <- append_row(
    rows = rows,
    characteristic = paste0(display_label, ", Median [IQR]"),
    total = format_median_iqr(total_vec, digits = digits),
    non_hcc = format_median_iqr(non_hcc_vec, digits = digits),
    incident_hcc = format_median_iqr(incident_hcc_vec, digits = digits),
    p_value = format_pvalue(p_val)
  )
  
  rows <- append_row(
    rows = rows,
    characteristic = "   Missing, n (%)",
    total = format_count_pct(sum(is.na(total_vec)), length(total_vec)),
    non_hcc = format_count_pct(sum(is.na(non_hcc_vec)), length(non_hcc_vec)),
    incident_hcc = format_count_pct(sum(is.na(incident_hcc_vec)), length(incident_hcc_vec)),
    p_value = ""
  )
  
  rows
}

append_categorical_rows <- function(rows, data, variable, display_label) {
  x <- data[[variable]]
  g <- data$Group
  complete_rows <- !is.na(x) & !is.na(g)
  p_val <- suppressWarnings(
    chisq.test(table(x[complete_rows], g[complete_rows]), correct = FALSE)$p.value
  )
  
  rows <- append_row(
    rows = rows,
    characteristic = paste0(display_label, ", n (%)"),
    total = "",
    non_hcc = "",
    incident_hcc = "",
    p_value = format_pvalue(p_val)
  )
  
  level_values <- levels(droplevels(x))
  
  total_denom <- sum(!is.na(x))
  non_hcc_denom <- sum(!is.na(x[g == "Non-HCC"]))
  incident_hcc_denom <- sum(!is.na(x[g == "Incident HCC"]))
  
  for (lev in level_values) {
    total_n <- sum(x == lev, na.rm = TRUE)
    non_hcc_n <- sum(x[g == "Non-HCC"] == lev, na.rm = TRUE)
    incident_hcc_n <- sum(x[g == "Incident HCC"] == lev, na.rm = TRUE)
    
    rows <- append_row(
      rows = rows,
      characteristic = paste0("   ", lev),
      total = format_count_pct(total_n, total_denom),
      non_hcc = format_count_pct(non_hcc_n, non_hcc_denom),
      incident_hcc = format_count_pct(incident_hcc_n, incident_hcc_denom),
      p_value = ""
    )
  }
  
  rows
}

create_baseline_summary_table <- function(df_baseline, biomarker_spec) {
  rows <- list()
  
  rows <- append_row(rows, "Demographics")
  rows <- append_continuous_rows(rows, df_baseline, "Age", "Age (years)", digits = 0)
  rows <- append_continuous_rows(rows, df_baseline, "BMI", "BMI (kg/m^2)", digits = 1)
  rows <- append_categorical_rows(rows, df_baseline, "Sex", "Sex")
  rows <- append_categorical_rows(rows, df_baseline, "Race", "Race")
  rows <- append_categorical_rows(rows, df_baseline, "Ethnicity", "Ethnicity")
  
  section_order <- c(
    "Liver Enzymes and Function",
    "Haematologic Indices",
    "Lipid Profile and Tumor Markers"
  )
  
  for (section_name in section_order) {
    rows <- append_row(rows, section_name)
    
    section_spec <- biomarker_spec %>%
      filter(table_section == section_name)
    
    for (i in seq_len(nrow(section_spec))) {
      rows <- append_continuous_rows(
        rows = rows,
        data = df_baseline,
        variable = section_spec$analysis_var[i],
        display_label = section_spec$label[i],
        digits = section_spec$digits[i]
      )
    }
  }
  
  table_df <- bind_rows(rows)
  
  n_total <- nrow(df_baseline)
  n_non_hcc <- sum(df_baseline$Group == "Non-HCC", na.rm = TRUE)
  n_incident_hcc <- sum(df_baseline$Group == "Incident HCC", na.rm = TRUE)
  
  colnames(table_df) <- c(
    "Characteristic",
    paste0("Total (N = ", format(n_total, big.mark = ","), ")"),
    paste0("Non-HCC (N = ", format(n_non_hcc, big.mark = ","), ")"),
    paste0("Incident HCC (N = ", format(n_incident_hcc, big.mark = ","), ")"),
    "P-value"
  )
  
  table_df
}

create_baseline_distribution_plot <- function(df_baseline, biomarker_spec) {
  panel_design <- "
ABCD
EFGH
IJKL
#MN#
"
  
  plot_order <- intersect(longitudinal_plot_order, biomarker_spec$analysis_var)
  label_map <- stats::setNames(biomarker_spec$label, biomarker_spec$analysis_var)
  
  plot_data_baseline <- df_baseline %>%
    mutate(
      Plot_Group = factor(
        Group,
        levels = c("Incident HCC", "Non-HCC"),
        labels = c("HCC", "Non-HCC")
      )
    ) %>%
    dplyr::select(Plot_Group, all_of(plot_order)) %>%
    pivot_longer(
      cols = -Plot_Group,
      names_to = "Biomarker",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(
      Biomarker = factor(
        Biomarker,
        levels = plot_order,
        labels = label_map[plot_order]
      )
    )
  
  ggplot(plot_data_baseline, aes(x = Plot_Group, y = Value, fill = Plot_Group)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(
      values = c("HCC" = "#E15759", "Non-HCC" = "#868686"),
      name = "HCC Status"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#D9D9D9", color = "#666666")
    ) +
    labs(
      y = "Biomarker level",
      title = "Baseline Distributions of All Biomarkers by HCC Status"
    ) +
    ggh4x::facet_manual(
      vars(Biomarker),
      design = panel_design,
      scales = "free_y"
    )
}

create_correlation_heatmap <- function(df_baseline, biomarker_spec) {
  heatmap_vars <- biomarker_spec$analysis_var
  
  short_label_map <- c(
    ALT = "ALT",
    AST = "AST",
    Bilirubin = "Bilirubin",
    Albumin = "Albumin",
    INR = "INR",
    Platelets = "Platelets",
    Neutrophils = "Neutrophils",
    Lymphocytes = "Lymphocytes",
    Hemoglobin = "Hemoglobin",
    TC = "Total Cholesterol",
    HDL = "HDL Cholesterol",
    LDL = "LDL Cholesterol",
    TG = "Triglycerides",
    AFP = "Alpha-fetoprotein"
  )
  
  heatmap_data <- df_baseline %>%
    dplyr::select(all_of(heatmap_vars))
  
  colnames(heatmap_data) <- short_label_map[heatmap_vars]
  
  cormat <- round(
    cor(
      heatmap_data,
      use = "pairwise.complete.obs",
      method = "spearman"
    ),
    2
  )
  
  cormat[upper.tri(cormat, diag = TRUE)] <- NA_real_
  
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
  colnames(melted_cormat) <- c("Row", "Col", "value")
  
  melted_cormat$Col <- factor(
    as.character(melted_cormat$Col),
    levels = colnames(cormat)
  )
  melted_cormat$Row <- factor(
    as.character(melted_cormat$Row),
    levels = rev(rownames(cormat))
  )
  
  ggplot(melted_cormat, aes(x = Col, y = Row, fill = value)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
    scale_fill_gradient2(
      low = "#075AFF",
      mid = "white",
      high = "#FF0000",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Spearman Correlation Heatmap"
    )
}
prepare_longitudinal_data <- function(df, biomarker_spec, max_months = 60) {
  longitudinal_raw <- df %>%
    filter(
      !is.na(Subject_ID),
      !is.na(HCC_group),
      !is.na(RESULT_MONTH),
      !is.na(First_LD_DATE)
    ) %>%
    mutate(
      First_LD_Calendar_Month = floor_date(First_LD_DATE, unit = "month"),
      Result_Calendar_Month = floor_date(RESULT_MONTH, unit = "month"),
      Time_Months =
        12 * (year(Result_Calendar_Month) - year(First_LD_Calendar_Month)) +
        (month(Result_Calendar_Month) - month(First_LD_Calendar_Month)),
      Time_Years = Time_Months / 12
    ) %>%
    filter(Time_Months >= 0, Time_Months <= max_months)
  
  longitudinal_df <- longitudinal_raw %>%
    dplyr::select(
      Subject_ID,
      HCC_group,
      Time_Months,
      Time_Years,
      all_of(biomarker_spec$source_var)
    )
  
  names(longitudinal_df)[5:ncol(longitudinal_df)] <- biomarker_spec$analysis_var
  
  plot_order <- intersect(longitudinal_plot_order, biomarker_spec$analysis_var)
  label_map <- stats::setNames(biomarker_spec$label, biomarker_spec$analysis_var)
  
  longitudinal_df %>%
    pivot_longer(
      cols = all_of(biomarker_spec$analysis_var),
      names_to = "Biomarker",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    group_by(Subject_ID, HCC_group, Time_Months, Biomarker) %>%
    summarise(
      Value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Time_Years = Time_Months / 12,
      HCC_status = factor(HCC_group, levels = c("HCC", "Non-HCC")),
      Biomarker = factor(
        Biomarker,
        levels = plot_order,
        labels = label_map[plot_order]
      )
    ) %>%
    filter(!is.na(Biomarker), !is.na(HCC_status))
}

create_longitudinal_plot <- function(trajectory_summary, afp_display_upper = 200) {
  panel_design <- "
ABCD
EFGH
IJKL
#MN#
"
  
  plot_df <- trajectory_summary %>%
    mutate(
      median_value_plot = ifelse(
        Biomarker == "AFP (ng/mL)",
        pmin(median_value, afp_display_upper),
        median_value
      ),
      q1_plot = ifelse(
        Biomarker == "AFP (ng/mL)",
        pmin(q1, afp_display_upper),
        q1
      ),
      q3_plot = ifelse(
        Biomarker == "AFP (ng/mL)",
        pmin(q3, afp_display_upper),
        q3
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = Time_Years,
      y = median_value_plot,
      color = HCC_status,
      fill = HCC_status,
      group = HCC_status
    )
  ) +
    geom_ribbon(
      aes(ymin = q1_plot, ymax = q3_plot),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1.05) +
    scale_x_continuous(
      limits = c(0, 3),
      breaks = 0:3,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_color_manual(
      values = c("HCC" = "#E15759", "Non-HCC" = "#8A8A8A"),
      breaks = c("HCC", "Non-HCC"),
      name = "HCC Status"
    ) +
    scale_fill_manual(
      values = c("HCC" = "#E15759", "Non-HCC" = "#8A8A8A"),
      breaks = c("HCC", "Non-HCC"),
      name = "HCC Status"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#D9D9D9", color = "#666666"),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.3),
      panel.grid.minor = element_line(color = "#ECECEC", linewidth = 0.2),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      axis.title = element_text(size = 12)
    ) +
    labs(
      title = "Longitudinal Trajectories of Biomarkers by HCC Status",
      x = "Years since first LD",
      y = "Biomarker Level"
    ) +
    guides(
      fill = "none",
      color = guide_legend(
        override.aes = list(
          linewidth = 1.2,
          alpha = 1
        )
      )
    ) +
    ggh4x::facet_manual(
      vars(Biomarker),
      design = panel_design,
      scales = "free_y"
    )
}

project_dir <- normalizePath(".", winslash = "/")
data_path <- file.path(
  project_dir,
  "data",
  "Final Clean Analysis Data with Biomarkers_v2.csv"
)
bmi_data_path <- file.path(project_dir, "data", "OneFL BMI Data - 05.08.25.csv")
output_dir <- file.path(project_dir, "code", "result_final")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- prepare_analysis_data(data_path, bmi_data_path)
baseline_summary <- create_baseline_summary_table(
  prepared$df_baseline,
  prepared$biomarker_spec
)
section_column <- names(baseline_summary)[1]
total_column <- names(baseline_summary)[2]
laboratory_section_row <- which(
  baseline_summary[[section_column]] == "Liver Enzymes and Function"
)
table1_baseline_characteristics <- baseline_summary[2:(laboratory_section_row - 1), ] %>%
  filter(!(Characteristic == "   Missing, n (%)" & .data[[total_column]] == "0 (0.0%)"))
table2_baseline_biomarkers <- baseline_summary[laboratory_section_row:nrow(baseline_summary), ]

figure3 <- create_baseline_distribution_plot(
  prepared$df_baseline,
  prepared$biomarker_spec
)

supplementary_figure1 <- create_correlation_heatmap(
  prepared$df_baseline,
  prepared$biomarker_spec
)

longitudinal_data <- prepare_longitudinal_data(
  prepared$df,
  prepared$biomarker_spec,
  max_months = 36
)
trajectory_summary <- longitudinal_data %>%
  group_by(HCC_status, Biomarker, Time_Months, Time_Years) %>%
  summarise(
    n_subjects = n_distinct(Subject_ID),
    median_value = median(Value),
    q1 = quantile(Value, 0.25, names = FALSE),
    q3 = quantile(Value, 0.75, names = FALSE),
    .groups = "drop"
  )
figure4 <- create_longitudinal_plot(trajectory_summary)

write.csv(
  table1_baseline_characteristics,
  file.path(output_dir, "table1_baseline_characteristics.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  table2_baseline_biomarkers,
  file.path(output_dir, "table2_baseline_biomarkers.csv"),
  row.names = FALSE,
  na = ""
)

sink(file.path(output_dir, "table1_console.txt"))
print(as.data.frame(table1_baseline_characteristics), row.names = FALSE)
cat("\n")
print(as.data.frame(table2_baseline_biomarkers), row.names = FALSE)
sink()

ggsave(
  file.path(output_dir, "Figure3_Baseline_Biomarker_Distributions.png"),
  figure3,
  width = 10,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Supplementary_Figure1_Biomarker_Correlations.png"),
  supplementary_figure1,
  width = 12,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure4_Longitudinal_Biomarker_Trajectories.png"),
  figure4,
  width = 15,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)

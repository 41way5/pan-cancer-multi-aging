rm(list = ls())
# Lib===================================================
{
  library(lubridate)
  library(survival)
  library(patchwork)
  library(gridExtra)
  library(cowplot)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(scales)
  library(magick)
  library(tidyverse)
  library(data.table)
  library(Hmisc)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(missRanger)

  setwd("E:\\Tutor\\UK Biobank\\")
}

load("./Project/PAC and cancer/Data/UKB/ukb_aging_ca_Dec.RData")
load("./Project/PAC and cancer/Data/UKB/ukb_aging_ca_reg_Dec.RData") # Added dataset

# Add cancer registry data----
ca_diseases <- hesin_all_clean %>%
  filter(Category == "Cancers") %>%
  distinct(Disease) %>%
  pull() %>%
  setdiff(c("Liver Mets.", "Lymph Node Mets.")) %>%
  union("Brain Cancer")
ca_diseases

date_cols_reg_death <- paste0(ca_diseases, "_reg_death_date") # Hesin only (from hallmark_all_final)
date_cols_reg <- paste0(ca_diseases, "_reg_date") # Cancer registry only (from ukb_ca)
date_cols_hesin <- paste0(ca_diseases, "_date") # Hesin only (from ukb_ca_death)

ukb_ca <- ukb_ca %>%
  mutate(ca_reg_date = ca_date) %>%
  mutate(across(any_of(c(date_cols_reg, "ca_reg_date", "date_of_death")), ~ ifelse(is.na(.), "", as.character(.))))
ukb_ca_death <- ukb_ca_death %>%
  mutate(ca_reg_death_date = ca_date) %>%
  mutate(across(any_of(c(date_cols_reg_death, "ca_reg_death_date", "date_of_death")), ~ ifelse(is.na(.), "", as.character(.))))

# Pan-cancer associations with multilayer and multiorgan biological aging across diverse populations

#----Figure 1: study design and overview----
# Panel a: overview plot-----
# multilayer and multiorgan aging clocks and disease of cellular aging hallmarks#

# All disease/phenotype category
# AIC, Altered intercellular communication, (+)
# CS, Cellular senescence (?)
# DNS, Deregulated nutrient sensing (?)
# EA, Epigenetic alterations (=)
# GI, Genomic instability (=)
# LP, Loss of proteostasis (+)
# MD, Mitochondrial dysfunction (+)
# SCE, Stem cell exhaustion (-)
# TA, Telomere attrition (-)

hallmark <- c("GI", "EA", "LP", "TA", "DNS", "CS", "MD", "AIC", "SCE")
clock <- c("Proteomics", "Metabolomics", "Clinical", "Disease")
ca <- c("=", "=", "+", "-", "±", "±", "+", "+", "-")

mapping <- list(
  Proteomics        = c("LP", "AIC", "CS"),
  Metabolomics      = c("DNS", "MD"),
  Clinical          = c("DNS", "CS", "SCE", "MD"),
  Disease           = c("GI", "EA", "LP", "TA", "DNS", "CS", "MD", "AIC", "SCE")
)

mat <- matrix(0,
  nrow = length(hallmark), ncol = length(clock),
  dimnames = list(hallmark, clock)
)

for (cl in names(mapping)) {
  mat[mapping[[cl]], cl] <- 1
}

mat <- matrix(0,
  nrow = length(hallmark), ncol = length(clock),
  dimnames = list(hallmark, clock)
)
for (cl in names(mapping)) {
  mat[mapping[[cl]], cl] <- 1
}

mat2 <- cbind(mat, Cancer = NA)

library("scales")
show_col(pal_npg("nrc")(10))
# Parl
Heatmap(
  t(mat2),
  col = c("0" = "white", "1" = "#00A087FF"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  column_title = "Hallmarks of Aging",
  row_title = "Aging clocks",
  height = unit(1, "inch"),
  width = unit(nrow(mat2) / ncol(mat2) * 1, "inch"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6, fontface = "bold"),
  column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  column_names_rot = 30,
  column_names_side = "top",
  show_heatmap_legend = FALSE,
  row_split = factor(c(rep("Clock", ncol(mat)), "Cancer"),
    levels = c("Clock", "Cancer")
  ),
  column_split = factor(c(rep("Causes", 4), rep("Responses", 2), rep("Culprits", 3)),
    levels = c("Causes", "Responses", "Culprits")
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (rownames(t(mat2))[i] == "Cancer") {
      grid.rect(x, y, w, h, gp = gpar(col = "black", fill = "white", lwd = 0.5))
      grid.text(ca[j], x, y, gp = gpar(fontsize = 7))
    } else {
      grid.rect(x, y, w, h, gp = gpar(col = "black", fill = fill, lwd = 0.5))
    }
  }
)

# Panel b -----
organ_list <- c(
  "Organismal",
  "Immune", "Artery", "Adipose", "Muscle", "Skin",
  "Brain", "Esophagus",
  "Lung", "Heart",
  "Liver", "Stomach", "Pancreas", "Intestine", "Kidney", "Male"
)
length(organ_list) # 16
organ_pac <- paste0(organ_list, "_res")

cox <- function(data, exposures, outcomes, covar, condition = T, early = F) {
  result_list <- data.frame()
  p_values_all <- list()
  select <- dplyr::select
  data <- data %>%
    mutate(across(where(is.character), ~ ifelse(. == "", NA, .))) %>%
    mutate(death = ifelse(is.na(death), 0, 1)) %>%
    mutate(death_date = date_of_death, death_s = NA)

  # Add count and running time
  start_time <- Sys.time()
  total_tasks <- length(outcomes) * length(exposures)
  counter <- 0

  # Calculated the incident case and follow-up years
  for (var in outcomes) {
    if (early) {
      data <- data %>%
        mutate(
          !!paste0(var, "_base") := ifelse(as.Date(!!sym(paste0(var, "_date"))) <= as.Date(recruit_date), 1, NA),
          !!paste0(var, "_inci") := ifelse(as.Date(!!sym(paste0(var, "_date"))) > as.Date(recruit_date) & !!paste0(var, "_s") != 0, 1, NA)
        ) %>%
        mutate(
          !!paste0(var, "_date") := ifelse(!is.na(!!sym(paste0(var, "_date"))), !!sym(paste0(var, "_date")), ifelse(death == 1, date_of_death, max(!!sym(paste0(var, "_date")), na.rm = T))),
          !!paste0(var, "_fu") := as.numeric(as.Date(!!sym(paste0(var, "_date"))) - as.Date(recruit_date)) / 365.25,
          !!paste0(var, "_inci") := ifelse(is.na(!!sym(paste0(var, "_inci"))), 0, !!sym(paste0(var, "_inci"))),
          !!paste0(var, "_base2") := ifelse(is.na(!!sym(paste0(var, "_base"))), 0, !!sym(paste0(var, "_base")))
        ) %>%
        mutate(
          !!paste0(var, "_diag_age") := as.numeric(!!sym(paste0(var, "_fu"))) + age,
          !!paste0(var, "_inci") := ifelse(!!sym(paste0(var, "_inci")) == 1 & !!sym(paste0(var, "_diag_age")) < 50, 1, 0),
          !!paste0(var) := ifelse(!!sym(paste0(var, "_inci")) == 1, 1, 0)
        ) %>%
        select(-!!sym(paste0(var, "_diag_age")))
    } else {
      data <- data %>%
        mutate(
          !!paste0(var, "_base") := ifelse(as.Date(!!sym(paste0(var, "_date"))) <= as.Date(recruit_date), 1, NA),
          !!paste0(var, "_inci") := ifelse(as.Date(!!sym(paste0(var, "_date"))) > as.Date(recruit_date) & !!paste0(var, "_s") != 0, 1, NA)
        ) %>%
        mutate(
          !!paste0(var, "_date") := ifelse(!is.na(!!sym(paste0(var, "_date"))), !!sym(paste0(var, "_date")), ifelse(death == 1, date_of_death, max(!!sym(paste0(var, "_date")), na.rm = T))),
          !!paste0(var, "_fu") := as.numeric(as.Date(!!sym(paste0(var, "_date"))) - as.Date(recruit_date)) / 365.25,
          !!paste0(var, "_inci") := ifelse(is.na(!!sym(paste0(var, "_inci"))), 0, !!sym(paste0(var, "_inci"))),
          !!paste0(var, "_base2") := ifelse(is.na(!!sym(paste0(var, "_base"))), 0, !!sym(paste0(var, "_base"))),
          !!paste0(var) := ifelse(!!sym(paste0(var, "_inci")) == 1, 1, 0)
        )
    }
  }

  for (outcome in outcomes) {
    for (exposure in exposures) {
      # Print count##################################################################
      counter <- counter + 1
      percent_done <- round(counter / total_tasks * 100, 1)
      message(
        "Running: outcome = ", outcome,
        ", exposure = ", exposure,
        " ... (", counter, "/", total_tasks, ", ", percent_done, "% done)"
      )
      ##############################################################################

      formula_input <- as.formula(paste0(paste0("survival::Surv( `", outcome, "_fu`,`", outcome, "_inci`", ") ~"), paste(c(paste0("`", exposure, "`"), covar), collapse = "+")))

      data2 <- data %>%
        filter(is.na(!!sym(paste0(outcome, "_base")))) %>%
        filter(!is.na(!!sym(paste0(exposure))))

      if (grepl("Male", exposure)) {
        data2 <- data2 %>% filter(sex == 1)
      }

      if (grepl("Male", exposure) && outcome %in% c(
        "breast", "uterine", "ovarian",
        "Breast Cancer", "Uterine Cancer", "Ovarian Cancer"
      )) {
        next
      }

      if (outcome %in% c(
        "breast", "uterine", "ovarian",
        "Breast Cancer", "Uterine Cancer", "Ovarian Cancer"
      )) {
        data2 <- data2 %>% filter(sex == 0)
      }

      if (outcome %in% c("prostate", "Prostate Cancer")) {
        data2 <- data2 %>% filter(sex == 1)
      }

      cox_model <- survival::coxph(formula = formula_input, data = data2)
      cox_summary <- summary(cox_model)
      coef_data <- as.data.frame(cox_summary$coefficients[1, "coef"])
      se_data <- as.data.frame(cox_summary$coefficients[1, "se(coef)"])
      z_data <- as.data.frame(cox_summary$coefficients[1, "z"])
      p_values <- cox_summary$coefficients[1, "Pr(>|z|)"]
      p_values_all[[paste(outcome, exposure)]] <- p_values

      n_cases_inci <- sum(data2[, paste0(outcome, "_inci")] == 1)
      n_cases_base <- sum(data[!is.na(data[, exposure]), paste0(outcome, "_base2")] == 1)

      N <- nrow(data2)

      mean_diff <- data %>%
        summarise(
          mean_diff = mean(ifelse(!!sym(outcome) == 1, !!sym(exposure), NA), na.rm = TRUE) -
            mean(ifelse(!!sym(outcome) == 0, !!sym(exposure), NA), na.rm = TRUE)
        ) %>%
        pull(mean_diff)

      median_age <- data2 %>%
        filter(!!sym(paste0(outcome, "_inci")) == 1) %>%
        summarise(median_fu = median(as.numeric(!!sym(paste0(outcome, "_fu"))) + age, na.rm = TRUE)) %>%
        pull(median_fu)

      median_follow_up <- median(as.numeric(data2[[paste0(outcome, "_fu")]]), na.rm = TRUE)
      hr_data <- exp(coef_data)
      hr_ci_data <- exp(confint(cox_model))

      # Store results in a list
      result <- data.frame(
        Outcome = outcome,
        Exposure = exposure,
        N_Cases_inci = n_cases_inci,
        N_Cases_base = n_cases_base,
        N = nrow(data2),
        Median_age = median_age,
        Mean_diff = mean_diff,
        Median_Follow_Up = median_follow_up,
        Z = z_data[1, ],
        Beta = coef_data[1, ],
        SE = se_data[1, ],
        HR = hr_data[1, ],
        CI_Lower = hr_ci_data[1, "2.5 %"],
        CI_Upper = hr_ci_data[1, "97.5 %"],
        P_Value = p_values
      )

      result_list <- rbind(result_list, result)
    }
  }

  all_p_values <- unlist(p_values_all)
  fdr_p_values <- p.adjust(all_p_values, method = "fdr")

  result_list$FDR_P_Value <- fdr_p_values[paste(result_list$Outcome, result_list$Exposure)]

  end_time <- Sys.time()
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (run_time < 60) {
    msg <- paste0(round(run_time, 2), " seconds")
  } else if (run_time < 3600) {
    msg <- paste0(round(run_time / 60, 2), " minutes")
  } else {
    msg <- paste0(round(run_time / 3600, 2), " hours")
  }

  message("\nTotal runtime: ", msg)

  if (condition) {
    return(list(result_list = result_list, data = data))
  } else {
    return(list(result_list = result_list))
  }
}

# Check cases overlapping between Hesin and Can reg-----

cancer_overlap_summary <- do.call(rbind, lapply(ca_diseases, function(ca) {
  dat <- dt_aging_2_cancer

  hcol <- paste0(ca, "_date")
  rcol <- paste0(ca, "_reg_date")

  hflag <- dat[[hcol]] != ""
  rflag <- dat[[rcol]] != ""

  hesin_case <- n_distinct(dat$eid[hflag])
  reg_case <- n_distinct(dat$eid[rflag])
  both_case <- n_distinct(dat$eid[hflag & rflag])
  hesin_only <- n_distinct(dat$eid[hflag & !rflag])
  reg_only <- n_distinct(dat$eid[!hflag & rflag])

  data.frame(
    cancer = ca,
    reg_case = reg_case,
    hesin_case = hesin_case,
    both_case = both_case,
    reg_only_case = reg_only,
    hesin_only_case = hesin_only,
    stringsAsFactors = FALSE
  )
})) %>% as_tibble()

cancer_overlap_summary <- cancer_overlap_summary %>%
  mutate(
    both_over_reg_pct   = 100 * both_case / reg_case,
    both_over_hesin_pct = 100 * both_case / hesin_case
  )

# HeatPlot view function
ht_ca <- function(data) {
  mat <- xtabs(Beta ~ Outcome + Exposure, data)
  mat_HR <- xtabs(HR ~ Outcome + Exposure, data)

  mat_p <- xtabs(P_Value ~ Outcome + Exposure, data)
  mat_q <- xtabs(FDR_P_Value ~ Outcome + Exposure, data)
  col_fun = colorRamp2(c(-1, 0, 1), c("#007F7F", "white", "#D85F4D"))

  Heatmap(mat,
    col = col_fun, na_col = "white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    cluster_rows = F, cluster_columns = F,
    row_names_side = "left", row_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_names_side = "bottom", show_column_names = T, column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 6),
    column_title_side = "bottom",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (mat_q[i, j] < 0.05) {
        grid.text(sprintf("%.1f", mat_HR[i, j]), x, y, gp = gpar(fontsize = 5))
      } else if ((mat_p)[i, j] < 0.05) {
        if ((mat)[i, j] > 1) {
          grid.text("+", x, y, gp = gpar(fontsize = 4))
        }
      }
    },
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      title = NULL,
      legend_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "mm"),
      labels_gp = gpar(col = "black", fontsize = 5)
    ),

    row_title = NULL
  ) %>% draw()
}

# Association test
## *****Section A: aging 2 cancer***** ----
### A.1 (hallmark diseases) ----
# defined dataset including multiomic aging clocks, covar and cancer registry data

dt_aging_2_cancer <- hallmark_all_final %>%
  select(any_of(c(
    "eid", "non_ca_base",
    c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi"), # covar
    c(paste0("base_", categories), "any_base"), # exposure
    date_cols_hesin
  ))) %>%
  left_join(ukb_ca_death)

cox_disease_ca <- cox(
  data = dt_aging_2_cancer %>% filter(ca_base_any == 0), # Without any malignant cancer at baseline
  exposures = c(paste0("base_", categories), "any_base"),
  outcomes = c(
    "ca_reg_death",
    paste0((ca_diseases), "_reg_death"),
    "death"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

disease_ca <- cox_disease_ca$result_list
ht_ca(disease_ca)

### A.2 (Proteomics, clinical, and metabolomic aging clocks) ----

dt_aging_clock_2_cancer <- hallmark_all_final %>%
  select(any_of(c(
    "eid", "t1", "t2",
    c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi", "non_ca_base"), # covar
    c(organ_pac, "phenoage_res", "kdm_res", "met_res"), # exposure
    date_cols_hesin
  ))) %>%
  left_join(ukb_ca_death) %>%
  mutate(`Liver Cancer_reg_death_date` = pmax(`Liver Cancer_date`, `Liver Cancer_reg_death_date`))

table(dt_aging_clock_2_cancer$t1) # 46785
table(is.na(dt_aging_clock_2_cancer$Organismal_res)) # 53014

cox_pac_ca <- cox(
  data = dt_aging_clock_2_cancer %>% filter(ca_base_any == 0)
  ,
  exposures = c(
    organ_pac, "phenoage_res", "kdm_res", "met_res"
  ),
  outcomes = c(
    "ca_reg_death",
    paste0(ca_diseases, "_reg_death"),
    "death"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok")
)

pac_ca <- cox_pac_ca$result_list
ht_ca(pac_ca)

# Define levels for exposure and outcome
# Cancer and death
ca_levels <- c(
  # Overall
  "Any cancer",
  # Digestive system cancers
  "Oesophageal Cancer", "Stomach Cancer", "Bowel Cancer",
  "Liver Cancer", "Biliary Cancer", "Pancreatic Cancer",
  # Respiratory system cancers
  "Lung Cancer",
  # Skin cancers
  "Skin Cancer", "Melanoma", "Mesothelioma",
  # Female-specific cancers
  "Breast Cancer", "Uterine Cancer", "Ovarian Cancer",
  # Male-specific cancers
  "Prostate Cancer",
  # Genitourinary system cancers
  "Kidney Cancer", "Bladder Cancer",
  # CNS system cancers
  "Brain Cancer",
  # Endocrine cancers
  "Thyroid Cancer",
  # Hematological malignancies
  "Non-Hodgkins lymphoma", "Plasma Cell Cancer", "Leukaemia",
  "Polycythaemia Vera", "Myelodysplasia", "MGUS",
  # Outcome
  "All-cause mortality"
)
# Aging clocks
pac_levels <- c(
  "Organismal", "Immune", "Artery", "Adipose", "Muscle", "Skin",
  "Brain", "Esophagus", "Lung", "Heart",
  "Liver", "Stomach", "Pancreas", "Intestine", "Kidney", "Male",
  "MetAge",
  "PhenoAge", "KDMAge"
)

# Aging diseases
disease_levels <- c("Any hallmarks", "GI", "TA", "EA", "LP", "DNS", "CS", "MD", "SCE", "AIC")

# Plot view
disease_ca_fomrat <- function(result_list) {
  result_disease <- result_list$result_list %>%
    mutate(
      Outcome_label = case_when(
        Outcome == "ca" ~ "Any cancer",
        Outcome == "Cancers_all" ~ "Any cancer",
        Outcome == "brain" ~ "Brain cancer",
        Outcome == "esophagus" ~ "Oesophageal cancer",
        Outcome == "lung" ~ "Lung cancer",
        Outcome == "stomach" ~ "Stomach cancer",
        Outcome == "liver" ~ "Liver cancer",
        Outcome == "pancreas" ~ "Pancreatic cancer",
        Outcome == "kidney" ~ "Kidney cancer",
        Outcome == "bowel" ~ "Bowel cancer",
        Outcome == "melanoma" ~ "Melanoma",
        Outcome == "skin" ~ "Skin cancer",
        Outcome == "nhl" ~ "Non-Hodgkins lymphoma",
        Outcome == "NHL" ~ "Non-Hodgkins lymphoma",
        Outcome == "leukaemia" ~ "Leukaemia",
        Outcome == "pcc" ~ "Plasma cell cancer",
        Outcome == "breast" ~ "Breast cancer",
        Outcome == "uterine" ~ "Uterine cancer",
        Outcome == "ovarian" ~ "Ovarian cancer",
        Outcome == "prostate" ~ "Prostate cancer",
        Outcome == "death" ~ "All-cause mortality",
        TRUE ~ as.character(Outcome)
      ),
      Exposure_label = sub("_res", "", Exposure),
      Exposure_label = sub("_super$", "", Exposure_label),
      Exposure_label = sub("_extra$", "", Exposure_label),
      Exposure_label = case_when(
        Exposure == "phenoage_res" ~ "PhenoAge",
        Exposure == "kdm_res" ~ "KDMAge",
        Exposure == "met_res" ~ "MetAge",
        Exposure == "Telomere_length" ~ "Telomere length",
        TRUE ~ as.character(Exposure_label)
      ),
    ) %>%
    mutate(
      Outcome_label = factor(Outcome_label, levels = ca_levels),
      Exposure_label = factor(Exposure_label, levels = pac_levels)
    )

  return(result_disease)
}

# 1a, BA to CA plotting------
result_cancer_fomrat <- function(result_list) {
  result_disease <- result_list %>%
    mutate(
      Outcome_label = case_when(
        Outcome == "ca" ~ "Any cancer",
        Outcome == "cancer" ~ "Any cancer",
        Outcome == "Cancers_all" ~ "Any cancer",
        Outcome == "brain" ~ "Brain cancer",
        Outcome == "esophagus" ~ "Oesophageal cancer",
        Outcome == "lung" ~ "Lung cancer",
        Outcome == "stomach" ~ "Stomach cancer",
        Outcome == "liver" ~ "Liver cancer",
        Outcome == "pancreas" ~ "Pancreatic cancer",
        Outcome == "kidney" ~ "Kidney cancer",
        Outcome == "bowel" ~ "Bowel cancer",
        Outcome == "melanoma" ~ "Melanoma",
        Outcome == "skin" ~ "Skin cancer",
        Outcome == "nhl" ~ "Non-Hodgkins lymphoma",
        Outcome == "NHL" ~ "Non-Hodgkins lymphoma",
        Outcome == "leukaemia" ~ "Leukaemia",
        Outcome == "pcc" ~ "Plasma cell cancer",
        Outcome == "breast" ~ "Breast cancer",
        Outcome == "uterine" ~ "Uterine cancer",
        Outcome == "ovarian" ~ "Ovarian cancer",
        Outcome == "prostate" ~ "Prostate cancer",
        Outcome == "death" ~ "All-cause mortality",
        TRUE ~ as.character(Outcome)
      ),
      Exposure_label = sub("_res", "", Exposure),
      Exposure_label = sub("_super$", "", Exposure_label),
      Exposure_label = sub("_extra$", "", Exposure_label),
      Exposure_label = case_when(
        Exposure == "phenoage_res" ~ "PhenoAge",
        Exposure == "kdm_res" ~ "KDMAge",
        Exposure == "met_res" ~ "MetAge",
        Exposure == "Telomere_length" ~ "Telomere length",
        TRUE ~ as.character(Exposure_label)
      ),
      Exposure_label = case_when(
        Exposure == "any_base" ~ "Any hallmarks",
        TRUE ~ as.character(Exposure_label)
      ),
      Exposure_label = sub("base_", "", Exposure_label),
    ) %>%
    mutate(
      Outcome_label = factor(Outcome_label, levels = ca_levels),
      Exposure_label = factor(Exposure_label, levels = c(disease_levels, pac_levels))
    )

  return(result_disease)
}

aging_cancer_result <- rbind(cox_disease_ca$result_list, cox_pac_ca$result_list) %>%
  mutate(Outcome = gsub("_reg(_death)?$", "", Outcome)) %>%
  result_cancer_fomrat() %>%
  mutate(N_Cases_base = 0)

cancer_inci <- cox_disease_ca$data %>% select(-any_of(c(date_cols_hesin, "ca_date", "ca_base", "ca_inci")))

colnames(cancer_inci) <- gsub("_reg_death", "", colnames(cancer_inci))

disease_prefixes <- c("ca", "death", ca_diseases)

for (prefix in disease_prefixes) {
  age_col <- paste0(prefix, "_age")
  cancer_inci <- cancer_inci %>%
    mutate(!!age_col := if_else(!!sym(paste0(prefix, "_inci")) == 1,
      as.numeric(!!sym(paste0(prefix, "_fu")) + age),
      NA
    ))
}

aging_cancer_ht <- function(data) {
  ht_opt$ROW_ANNO_PADDING = unit(0.3, "mm")
  disease_prefixes <- c("ca", "death", ca_diseases)

  mat <- xtabs(Beta ~ Outcome_label + Exposure_label, data)
  mat[abs(mat) == 0] <- NA
  mat_p_fdr <- xtabs(FDR_P_Value ~ Outcome_label + Exposure_label, data)
  mat_p <- xtabs(P_Value ~ Outcome_label + Exposure_label, data)
  mat[(mat_p) > 0.05] <- 0
  mat[(mat_p) < 0.05 & (mat_p_fdr) >= 0.05 & (mat) < 0] <- 0

  mat_hr <- xtabs(HR ~ Outcome_label + Exposure_label, data)

  mat_case <- data %>%
    dplyr::select(Outcome_label, N_Cases_base, N_Cases_inci) %>%
    distinct(Outcome_label, .keep_all = TRUE) %>%
    dplyr::select(2:3)
  mat_disease_gap <- xtabs(Mean_diff ~ Outcome_label + Exposure_label, data)
  mat_disease_gap <- mat_disease_gap[, c(11:21)]

  mat_distri <- cancer_inci %>%
    dplyr::select(paste0(disease_prefixes, "_age")) %>%
    rename_with(~ gsub("_age", "", .), contains("_age")) %>%
    as.data.frame() %>%
    rename(
      ca = "Any cancer", NHL = "Non-Hodgkins lymphoma", death = "All-cause mortality"
    ) %>%
    t()

  mat_distri_order <- mat_distri[match(rownames(mat), rownames(mat_distri)), ]

  sig_count <- rowSums(mat_p_fdr < 0.05, na.rm = TRUE)

  sig_count_mat <- cbind(
    group1 = rowSums(mat_p_fdr[, 1:9, drop = FALSE] < 0.05, na.rm = TRUE),
    group2 = rowSums(mat_p_fdr[, 10:21, drop = FALSE] < 0.05, na.rm = TRUE),
    group3 = rowSums(mat_p_fdr[, 22, drop = FALSE] < 0.05, na.rm = TRUE),
    group4 = rowSums(mat_p_fdr[, 23:24, drop = FALSE] < 0.05, na.rm = TRUE)
  )

  col_fun = colorRamp2(c(-1, 0, 1), c("#007F7F", "white", "#D85F4D"))
  color_panel <- c("Black", c("#A52A2A", "#4B0082", "#FF7F50", "#FFFF00", "#808080", "#808000", "#D2B48C", "#00FFFF", "#0000FF", "#008000"))

  cancer_class <- c(
    "Any cancer" = "Any",
    "All-cause mortality" = "Any",
    "Oesophageal Cancer" = "Digestive",
    "Stomach Cancer" = "Digestive",
    "Bowel Cancer" = "Digestive",
    "Liver Cancer" = "Digestive",
    "Biliary Cancer" = "Digestive",
    "Pancreatic Cancer" = "Digestive",
    "Lung Cancer" = "Lung",
    "Skin Cancer" = "Skin",
    "Melanoma" = "Skin",
    "Mesothelioma" = "Mesothelial",
    "Breast Cancer" = "Breast",
    "Uterine Cancer" = "Genital",
    "Ovarian Cancer" = "Genital",
    "Prostate Cancer" = "Genital",
    "Kidney Cancer" = "Urinary tract",
    "Bladder Cancer" = "Urinary tract",
    "Brain Cancer" = "Brain",
    "Thyroid Cancer" = "Thyroid",
    "Non-Hodgkins lymphoma" = "Haematological/\nImmunological",
    "Plasma Cell Cancer" = "Haematological/\nImmunological",
    "Leukaemia" = "Haematological/\nImmunological",
    "Polycythaemia Vera" = "Haematological/\nImmunological",
    "Myelodysplasia" = "Haematological/\nImmunological",
    "MGUS" = "Haematological/\nImmunological"
  )

  class_col <- c(
    "Digestive" = "#D2B48C",
    "Lung" = "#00FFFF",
    "Skin" = "#008000",
    "Mesothelial" = "#0000FF",
    "Breast" = "#F4A3B4",
    "Genital" = "#FF7F50",
    "Urinary tract" = "#808000",
    "Brain" = "#4B0082",
    "Thyroid" = "#808080",
    "Haematological/\nImmunological" = "#FFFF00",
    "Any" = "white"
  )

  cat_levels <- c(
    "Any",
    "Digestive",
    "Lung",
    "Skin",
    "Mesothelial",
    "Breast",
    "Genital",
    "Urinary tract",
    "Brain",
    "Thyroid",
    "Haematological/\nImmunological"
  )

  Heatmap(mat,
    col = col_fun, na_col = "lightgrey", show_column_dend = F, show_row_dend = F,
    border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    cluster_rows = F, cluster_columns = F,
    width = ncol(mat) * unit(3, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_title = "Multi-omics, multi-organ aging hallmarks and incident cancer and mortality",
    column_title_gp = gpar(fontsize = 6, fontface = "bold"),
    column_names_side = "bottom", show_column_names = T, column_names_gp = gpar(fontsize = 6),
    row_split = factor(c(rep("Any", 1), rep("Cancer subtypes", 24), rep("All-cause mortality", 1)), levels = c("Any", "Cancer subtypes", "All-cause mortality")),
    column_split = factor(c(rep("Aging hallmark disease", each = 10), rep("Organ proteomic aging", each = 16), rep("Metabolic aging", each = 1), rep("Phenotypic aging", each = 2)), levels = c("Aging hallmark disease", "Organ proteomic aging", "Metabolic aging", "Phenotypic aging")),
    column_title_side = "top",

    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(mat[i, j]) && !is.na(mat_p_fdr[i, j]) && mat_p_fdr[i, j] < 0.05) {
        grid.text(sprintf("%.1f", mat_hr[i, j]), x, y, gp = gpar(fontsize = 5))
      } else if (!is.na(mat[i, j]) && !is.na(mat_p[i, j]) && mat_p[i, j] < 0.05 && mat[i, j] > 0) {
        grid.text("+", x, y, gp = gpar(fontsize = 4))
      }
    },
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      title = "Log (HR)",
      legend_height = unit(0.6, "cm"),
      grid_width = unit(0.5, "mm"),
      labels_gp = gpar(col = "black", fontsize = 5),
      title_gp = gpar(col = "black", fontsize = 5)
    ),

    left_annotation = rowAnnotation(
      "FDR hits" = anno_lines(
        sig_count,
        gp = gpar(col = "#808000", lwd = 0.4),
        add_points = TRUE,
        pt_gp = gpar(col = "black", pch = 16),
        size = unit(0.5, "mm"),
        border = F,
        axis = TRUE,
        axis_param = list(
          side = "top",
          gp = gpar(fontsize = 5, lwd = 0.2), labels_rot = 0
        ),
        width = unit(0.7, "cm")
      ),
      " " = factor(cancer_class[rownames(mat)], levels = cat_levels),
      col = list(" " = class_col),
      show_legend = c("FDR hits" = FALSE, " " = TRUE),
      annotation_name_gp = gpar(fontsize = 5),
      annotation_name_rot = c("FDR hits" = 0),
      simple_anno_size = unit(1, "mm"),
      gap = unit(0.5, "mm"),
      annotation_legend_param = list(
        " " = list(
          title_gp = gpar(fontsize = 5),
          labels_gp = gpar(fontsize = 5),
          grid_width = unit(1, "mm"),
          grid_height = unit(1, "mm")
        )
      )
    ),
    right_annotation = rowAnnotation(
      "Incident\ncases" = anno_barplot(mat_case,
        border = F, ylim = c(0, 70000), gp = gpar(fill = c("#00a087ff", "#3c5488ff"), col = c("#00a087ff", "#3c5488ff")),
        axis_param = list(gp = gpar(fontsize = 4, col = "black", lwd = 0.2), at = c(1000, 20000, 40000, 75000), labels = c("1k", "20k", "40k", "75k"), labels_rot = 45),
      ),
      "Age\nat diagnosis" = anno_density(mat_distri_order,
        type = "violin", joyplot_scale = 2,
        gp = gpar(fill = "#CCCCCC80", lwd = 0.2, col = c("black", rep("black", 10))), border = F,
        axis_param = list(gp = gpar(fontsize = 5, col = "black", lwd = 0.2), labels_rot = 0)
      ),

      width = unit(2, "cm"),
      annotation_name_gp = gpar(fontsize = 5),
      annotation_width = unit(c(1, 1), c("null", "null")),
      gap = unit(2.5, "mm")
    ),

    row_title = NULL
  ) %>% return()
}

aging_cancer_ht(aging_cancer_result) %>% draw()

cairo_pdf("./Project/PAC and cancer/Figures/fig.1/reg/BA_2_CA.pdf",
  width = 7, height = 3
)

ht <- aging_cancer_ht(aging_cancer_result) %>% draw()
for (i in seq_len(3)) {
  decorate_annotation("FDR hits",
    {
      grid.rect(gp = gpar(col = "black", lwd = 0.2, fill = NA))
    },
    slice = i
  )
}

dev.off()

## *****Section B: CA to BA***** ----
### B.1 CA to BA (Any hallmark diseases) ----
#### Overall ----

table(hallmark_all_final$cancer_base)
table(ukb_ca_death$ca_base)
table(ukb_ca_death$ca_base_any)

dt_cancer_2_aging <- hallmark_all_final %>%
  select(any_of(c(
    "eid", "non_ca_base", "cancer_base",
    c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi"), # covar
    paste0(c(
      "AIC", "CS", "DNS", "EA", "GI", "LP", "MD", "SCE", "TA",
      labels,
      "any_1", "any_2", "any_3", "any_4"
    ), "_date"), # outcome
    c("ca_base_any") # exposure
  ))) %>%
  left_join(ukb_ca_death)

cox_ca_disease <- cox(
  data = dt_cancer_2_aging %>% filter(non_ca_base == 0), # Without any non-ca age-related diseases at baseline
  exposures = c(
    "ca_base_any"
  ),
  outcomes = c(
    "AIC", "CS", "DNS", "EA", "GI", "LP", "MD", "SCE", "TA",
    labels,
    "any_1", "any_2", "any_3", "any_4"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

describe(dt_cancer_2_aging %>% filter(non_ca_base == 0) %>% select(ca_base_any)) # 35819
nrow(dt_cancer_2_aging %>% filter(non_ca_base == 0)) # 427810
describe(hallmark_all$sex %>% as.factor())
#### By sex ----

cox_ca_disease_male <- cox(
  data = dt_cancer_2_aging %>% filter(non_ca_base == 0) %>% filter(sex == 1), # Without any non-ca age-related diseases at baseline
  exposures = c(
    "ca_base_any"
  ),
  outcomes = c(
    "AIC", "CS", "DNS", "EA", "GI", "LP", "MD", "SCE", "TA",
    labels,
    "any_1", "any_2", "any_3", "any_4"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

cox_ca_disease_female <- cox(
  data = dt_cancer_2_aging %>% filter(non_ca_base == 0) %>% filter(sex == 0), # Without any non-ca age-related diseases at baseline
  exposures = c(
    "ca_base_any"
  ),
  outcomes = c(
    "AIC", "CS", "DNS", "EA", "GI", "LP", "MD", "SCE", "TA",
    labels,
    "any_1", "any_2", "any_3", "any_4"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

# Formatting Cox results
process_res <- function(res_list) {
  res_list %>%
    filter(Outcome != "any") %>%
    mutate(
      label = case_when(
        grepl("_1$", Outcome) ~ "First disease",
        grepl("_2$", Outcome) ~ "Second disease",
        grepl("_3$", Outcome) ~ "Third disease",
        grepl("_4$", Outcome) ~ "Fourth disease",
        TRUE ~ "Any disease"
      ),
      group = gsub("_[0-9]+$", "", Outcome),
      group = case_when(
        group == "GI" ~ "Genomic instability",
        group == "TA" ~ "Telomere attrition",
        group == "EA" ~ "Epigenetic alterations",
        group == "LP" ~ "Loss of proteostasis",
        group == "DNS" ~ "Deregulated nutrient sensing",
        group == "CS" ~ "Cellular senescence",
        group == "MD" ~ "Mitochondrial dysfunction",
        group == "SCE" ~ "Stem cell exhaustion",
        group == "AIC" ~ "Altered intercellular communication",
        group == "any" ~ "Any aging hallmarks",
        TRUE ~ NA_character_
      ),
      subgroup = Outcome,
      diamond = case_when(label == "Any disease" ~ TRUE, TRUE ~ FALSE),
      fill = case_when(diamond ~ "white", TRUE ~ "black")
    ) %>%
    group_by(group) %>%
    mutate(Avg_HR = mean(HR, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(-Avg_HR, group) %>%
    select(-Avg_HR)
}
ca_disease <- process_res(cox_ca_disease$result_list)
ca_disease_male <- process_res(cox_ca_disease_male$result_list)
ca_disease_female <- process_res(cox_ca_disease_female$result_list)

# Add test heterogeneity for sex
het_test_one <- function(outcome, df_male, df_female) {
  m <- df_male %>% filter(Outcome == outcome)
  f <- df_female %>% filter(Outcome == outcome)

  if (nrow(m) == 0 | nrow(f) == 0) {
    return(data.frame(Outcome = outcome, Q = NA, p_het = NA, Beta_m = NA, Beta_f = NA))
  }

  Q <- (m$Beta - f$Beta)^2 / (m$SE^2 + f$SE^2)
  p <- 1 - pchisq(Q, df = 1)

  data.frame(Outcome = outcome, Q = Q, p_het = p, Beta_m = m$Beta, Beta_f = f$Beta)
}

outcomes_of_interest <- unique(ca_disease$subgroup)

phet_results <- do.call(rbind, lapply(outcomes_of_interest, het_test_one,
  df_male = ca_disease_male,
  df_female = ca_disease_female
)) %>%
  mutate(
    color = case_when(
      !is.na(p_het) & p_het < 0.05 & Beta_f > Beta_m ~ " F",
      !is.na(p_het) & p_het < 0.05 & Beta_f < Beta_m ~ " M",
      TRUE ~ ""
    )
  )

ca_disease_female <- ca_disease_female %>% left_join(phet_results %>% select(Outcome, color),
  by = c("subgroup" = "Outcome")
)

CA_2_BA_disease_plot <- ckbplotr::forest_plot(
  panels = list(
    Overall = process_res(cox_ca_disease$result_list) %>% mutate(color = ""),
    Male = process_res(cox_ca_disease_male$result_list) %>% mutate(color = ""),
    Female = process_res(cox_ca_disease_female$result_list) %>% left_join(phet_results %>% select(Outcome, color),
      by = c("subgroup" = "Outcome")
    )
  ),
  col.key = "subgroup",
  row.labels = process_res(cox_ca_disease$result_list) %>% select(subgroup, group, label),
  row.labels.space = c(0, 0.5, 0, 0),
  col.left = "N_Cases_inci",
  col.left.heading = "No.\nof events",
  col.estimate = "Beta",
  col.stderr = "SE",
  digits = 2,
  bold.labels = TRUE,
  logscale = FALSE,
  nullval = 1,
  ci.delim = "-",
  base_size = 7,
  base_line_size = 0.1,
  pointsize = 1.5,
  shape = 15,
  col.diamond = "diamond",
  fill = "fill",
  xlim = c(0.5, 7),
  xticks = c(1, 2, 4, 6),
  left.space = unit(5, "mm"),
  right.space = unit(15, "mm"),
  col.right = "color",
  col.right.heading = ""
)

cairo_pdf("./Project/PAC and cancer/Figures/fig.1/reg/CA_2_BA (disease)2.pdf",
  width = 7.5, height = 6
)
print(CA_2_BA_disease_plot$plot)
dev.off()

### B.2, CA to BA (pan cancer and hallmark diseases) ----
# Define the data for baseline disease
{
  disease_prefixes <- c("ca", ca_diseases)

  disease_prefixes_all <- c(
    "ca",
    paste0(ca_diseases, "_reg_death"),
    hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Disease) %>% pull(),
    hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Category) %>% pull(),
    "death"
  )

  hallmark_base <- hallmark_all_final %>%
    select(any_of(c(
      "eid", "non_ca_base",
      c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi"), # covar
      c(paste0("base_", categories), "any_base"), # exposure
      paste0(c(
        labels_nca[str_detect(labels_nca, "_1")],
        "any_1", "any_2", "any_3", "any_4"
      ), "_date"),
      paste0(c(hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Disease) %>% pull()), "_date"),
      paste0(c(hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Category) %>% pull()), "_date")
    ))) %>%
    left_join(ukb_ca_death) %>%
    mutate(death_date = date_of_death)

  for (prefix in disease_prefixes_all) {
    interval_col <- paste0(prefix, "_interval22")
    base_col <- paste0(prefix, "_base22")
    inci_col <- paste0(prefix, "_inci22")

    hallmark_base <- hallmark_base %>%
      mutate(
        !!interval_col := if_else((!!sym(paste0(prefix, "_date"))) != "",
          as.numeric(difftime(as.Date(.$recruit_date), as.Date(!!sym(paste0(prefix, "_date"))), units = "days")), 0
        ),
        !!base_col := if_else(!!sym(interval_col) > 0, 1, 0) %>% as.numeric(),
        !!inci_col := if_else(!!sym(interval_col) < 0, 1, 0) %>% as.numeric(),
      )
  }
}

cox_panca_disease <- cox(
  data = hallmark_base %>% filter(non_ca_base == 0), # Without any non-ca age-related diseases at baseline
  outcomes = c(
    labels_nca[str_detect(labels_nca, "_1")],
    "any_1", "any_2", "any_3", "any_4",
    "death"
  ),
  exposures = c(
    "ca_base_any", ca_diseases %>% paste0(., "_reg_death_base22")
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

panca_disease <- cox_panca_disease$result_list %>%
  mutate(
    Outcome = gsub("_1_nca", "", Outcome),
    Exposure = gsub("_reg_death", "", Exposure),
    Exposure = gsub("_base22", "", Exposure),
    Exposure = gsub("_base", "", Exposure)
  ) %>%
  mutate(
    Exposure = ifelse(Exposure == "ca_any", "Any cancer", Exposure),
    Exposure = ifelse(Exposure == "NHL", "Non-Hodgkins lymphoma", Exposure)
  ) %>%
  filter(Exposure != "Mesothelioma") %>% # Exclude Mesothelioma with <10 baseline cases
  mutate(Outcome = case_when(
    Outcome == "any_1" ~ "Any (1st)",
    Outcome == "any_2" ~ "Any (2nd)",
    Outcome == "any_3" ~ "Any (3rd)",
    Outcome == "any_4" ~ "Any (4th)",
    Outcome == "death" ~ "Death",
    T ~ Outcome
  )) %>%
  mutate(
    Exposure = factor(Exposure, levels = setdiff(ca_levels, c("All-cause mortality", "Mesothelioma"))),
    Outcome = factor(Outcome, levels = c("Any (1st)", "Any (2nd)", "Any (3rd)", "Any (4th)", "GI", "TA", "EA", "LP", "DNS", "CS", "MD", "SCE", "AIC", "Death"))
  )

# Baseline cancer cases
mat <- xtabs(HR ~ Outcome + Exposure, panca_disease) %>% t()
mat_base <- hallmark_base %>%
  filter(non_ca_base == 0) %>%
  select(all_of(c("ca_base_any", ca_diseases %>% paste0(., "_reg_death_base22")))) %>%
  rename_with(~ gsub("_reg_death", "", .), contains("_reg_death_base22")) %>%
  rename_with(~ gsub("_base22", "", .), contains("_base22")) %>%
  rename_with(~ gsub("_base", "", .), contains("_base")) %>%
  rename(ca_any = "Any cancer", NHL = "Non-Hodgkins lymphoma") %>%
  select(row.names(mat)) %>%
  summarise_all(~ sum(. == 1)) %>%
  as.matrix() %>%
  t()

ht_panca_disease <- function(data) {
  mat <- xtabs(Beta ~ Outcome + Exposure, data) %>% t()
  mat_hr <- xtabs(HR ~ Outcome + Exposure, data) %>% t()
  mat_p_fdr <- xtabs(FDR_P_Value ~ Outcome + Exposure, data) %>% t()
  mat_p <- xtabs(P_Value ~ Outcome + Exposure, data) %>% t()

  # Highlight sig cell only
  mat[(mat_p) > 0.05] <- 0
  mat[(mat_p) < 0.05 & (mat_p_fdr) >= 0.05 & (mat) < 0] <- 0

  col_fun = colorRamp2(c(-1, 0, 1, 2), c("#007F7F", "white", "#ef6548", "#b30000"))

  Heatmap(mat,
    col = col_fun, na_col = "white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    cluster_rows = F, cluster_columns = F,
    width = ncol(mat) * unit(3.3, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_names_side = "bottom", show_column_names = T, column_names_gp = gpar(fontsize = 6),
    column_split = factor(c(rep("Overall", each = 4), rep(" ", each = 9), rep("", each = 1)), levels = c("Overall", " ", "")),
    column_title_gp = gpar(fontsize = 6),
    column_title_side = "bottom",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (mat_p_fdr[i, j] < 0.05)
        grid.text(sprintf("%.1f", mat_hr[i, j]), x, y, gp = gpar(fontsize = 5))
    },
    heatmap_legend_param = list(
      title = "Log (HR)",
      legend_height = unit(0.5, "cm"),
      grid_width = unit(0.5, "mm"),
      labels_gp = gpar(col = "black", fontsize = 5),
      title_gp = gpar(col = "black", fontsize = 5)
    ),

    row_title = NULL,
    column_title = NULL,
  ) %>% return()
}
h1 <- ht_panca_disease(panca_disease)

### B.3, CA to BA (aging clocks) ----

ca_progression <- hallmark_base %>%
  left_join(hallmark_all_final %>% select(c("eid", c(organ_pac, "met_res", "phenoage_res", "kdm_res")))) %>%
  mutate(across(ends_with("_interval22"), ~ ifelse(. < 0, NA, .))) %>%
  mutate(across(contains("_interval22"), ~ (. / 365.25)))

prevalent_cases <- hallmark_base %>%
  mutate(across(ends_with("_interval22"), ~ ifelse(. < 0, NA, .) / 365.25)) %>%
  {
    interval_cols <- grep("_base22$", names(.), value = TRUE)
    case_mat <- mutate(select(., all_of(interval_cols)), across(everything(), ~ ifelse(!is.na(.) & . == 1, 1, 0)))
    data.frame(
      Outcome = interval_cols,
      Overall_Cases = colSums(case_mat, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

# Define function
calculate_beta_matrix <- function(data, var_list1, var_list2, covariates = NULL) {
  beta_matrix <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  fdr_p_values <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  n_matrix <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))

  for (i in seq_along(var_list1)) {
    for (j in seq_along(var_list2)) {

      formula <- as.formula(paste0(var_list1[i], " ~ `", var_list2[j], "`"))
      if (!is.null(covariates)) {
        covariate_formula <- paste(covariates, collapse = " + ")
        formula <- update(formula, paste(". ~ . +", covariate_formula))
      }
      lm_model <- lm(formula, data = data)
      beta_matrix[i, j] <- coef(lm_model)[2]

      p_value <- summary(lm_model)$coefficients[2, "Pr(>|t|)"]
      fdr_p_values[i, j] <- p.adjust(p_value, method = "fdr")

      base_var <- sub("_interval22$", "_base22", var_list2[j])
      if (base_var %in% names(data)) {
        n_matrix[i, j] <- sum(data[[base_var]] == 1 & !is.na(data[[var_list1[i]]]), na.rm = TRUE)
      } else {
        n_matrix[i, j] <- NA
      }
    }
  }

  colnames(beta_matrix) <- var_list2
  rownames(beta_matrix) <- var_list1
  colnames(fdr_p_values) <- var_list2
  rownames(fdr_p_values) <- var_list1
  colnames(n_matrix) <- var_list2
  rownames(n_matrix) <- var_list1

  return(list(
    beta_matrix = beta_matrix,
    fdr_p_values = fdr_p_values,
    n_matrix = n_matrix
  ))
}

disease_prefixes <- c("ca", ca_diseases)

table(ca_progression$ca_base_any) # 39470
table(ca_progression$ca_base) # 30775
table(ca_progression$non_ca_base) # 74370

result_ca_pac <- calculate_beta_matrix(ca_progression %>% filter(ca_base_any == 1),
  var_list1 = c(organ_pac, "met_res", "phenoage_res", "kdm_res"),
  var_list2 = c("ca_interval22", paste0(setdiff(ca_diseases, "Myelodysplasia"), "_reg_death_interval22")),
  covariates = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok")
)

result_ca_pac <- calculate_beta_matrix(ca_progression %>% filter(ca_base == 1),
  var_list1 = c(organ_pac, "met_res", "phenoage_res", "kdm_res"),
  var_list2 = c("ca_interval22", paste0(setdiff(ca_diseases, "Myelodysplasia"), "_reg_death_interval22")),
  covariates = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi", "non_ca_base")
)

result_ca_pac <- calculate_beta_matrix(ca_progression,
  var_list1 = c(organ_pac, "met_res", "phenoage_res", "kdm_res"),
  var_list2 = c("ca_interval22", paste0(setdiff(ca_diseases, "Myelodysplasia"), "_reg_death_interval22")),
  covariates = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "ca_D")
)

# Plot view
{
  disease_prefixes <- c("ca", ca_diseases)
  disease_prefixes2 <- gsub("ca", "Any cancer", disease_prefixes)
  disease_prefixes2 <- gsub("NHL", "Non-Hodgkins lymphoma", disease_prefixes2)
  disease_prefixes2 <- setdiff(disease_prefixes2, "Mesothelioma")

  rownames(result_ca_pac$beta_matrix) <- c(organ_list, "MetAge", "PhenoAge", "KDMAge")
  colnames(result_ca_pac$beta_matrix) <- disease_prefixes2

  rownames(result_ca_pac$fdr_p_values) <- c(organ_list, "MetAge", "PhenoAge", "KDMAge")
  colnames(result_ca_pac$fdr_p_values) <- disease_prefixes2

  rownames(result_ca_pac$n_matrix) <- c(organ_list, "MetAge", "PhenoAge", "KDMAge")
  colnames(result_ca_pac$n_matrix) <- disease_prefixes2

  mat_beta <- (result_ca_pac$beta_matrix) %>% t()
  mat_beta <- mat_beta[match(rownames(xtabs(HR ~ Outcome + Exposure, panca_disease) %>% t()), rownames(mat_beta)), ]
  mat_p <- result_ca_pac$fdr_p_values %>% t()
  mat_p <- mat_p[match(rownames(xtabs(HR ~ Outcome + Exposure, panca_disease) %>% t()), rownames(mat_p)), ]
  mat_n <- result_ca_pac$n_matrix %>% t()
  mat_n <- mat_n[match(rownames(xtabs(HR ~ Outcome + Exposure, panca_disease) %>% t()), rownames(mat_n)), ]

  mat_beta[(mat_p) > 0.05] <- 0
  mat_beta[abs(mat_beta) < 0.01] <- 0
  mat_p[abs(mat_beta) < 0.01] <- NA

  mat_beta[mat_n <= 3] <- 0
  mat_beta[mat_n == 0] <- NA
  mat_p[mat_n <= 3] <- NA

  h2 <- Heatmap(mat_beta,
    na_col = "lightgrey",
    name = "beta", border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    col = colorRamp2(c(-0.1, 0, 0.1), c("#007F7F", "white", "#D85F4D")),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = F, cluster_columns = F,
    row_dend_width = unit(0.3, "cm"),
    row_dend_gp = gpar(lwd = 0.2),
    column_dend_height = unit(0.3, "cm"),
    column_dend_gp = gpar(lwd = 0.2),
    show_column_dend = F, show_row_dend = F,
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 5),
    column_split = factor(c(rep("Organ-specific aging", each = 16), rep("Met aging", each = 1), rep("Pheno aging", each = 2)), levels = c("Organ-specific aging", "Met aging", "Pheno aging")),
    column_title_side = "top",
    row_names_side = "left",
    row_title = NULL,
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 4))
      } else if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 4))
      } else if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 4))
      }
    },
    heatmap_legend_param = list(
      at = c(-0.1, 0, 0.1),
      title = "beta",
      title_gp = gpar(fontsize = 5),
      legend_height = unit(1, "cm"),
      grid_width = unit(0.5, "mm"),
      legend_width = unit(1, "mm"),
      labels_gp = gpar(col = "black", fontsize = 5)
    ),
    right_annotation = rowAnnotation(
      "Prevalent cases" = anno_barplot(mat_base,
        border = F, ylim = c(0, 40000), gp = gpar(fill = c("#3c5488ff"), col = c("#3c5488ff")),
        axis_param = list(gp = gpar(fontsize = 4, col = "black", lwd = 0.2), at = c(500, 10000, 20000, 40000), labels = c("500", "10k", "20k", "40k"), labels_rot = 90),
      ),

      width = unit(1, "cm"),
      annotation_name_gp = gpar(fontsize = 6),
      gap = unit(2, "mm")
    ),
    width = ncol(mat_beta) * unit(3, "mm"),
    height = nrow(mat_beta) * unit(2, "mm")
  )
}

grid.grabExpr(
  draw(h1 + h2,
    ht_gap = unit(1, "mm"), auto_adjust = T,
    column_title = "Multi-omics, multi-organ aging hallmarks and incident cancer and mortality",
    column_title_gp = gpar(fontsize = 6, fontface = "bold")
  )
)

## *****Section C: BA to disease in CA***** ----
### C.1 BA to disease in all participants (validation of aging clocks)----

organ_mri_list <- c(
  "Adipose", "Brain", "Heart",
  "Liver", "Pancreas", "Kidney", "Spleen"
) # 7

organ_mri_pac <- paste0(organ_mri_list, "_mri_res")

# Non ca aging-related diseases
non_ca_date_cols <- setdiff(
  hesin_all_clean %>%
    filter(Category != "Cancers") %>%
    distinct(Disease) %>%
    pull(Disease),
  c("High Cholesterol, High Triglycerides, Low HDL, Raised LDL", "Myasthenia Gravis")
) %>% paste0("_date") # 56

dt_aging_2_disease_imaging <- hallmark_all_final %>%
  select(
    eid, "any_base", "non_ca_base", "date_visit_imaging1",
    c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "age_at_imaging"), # covar
    c(
      organ_pac, "phenoage_res", "kdm_res", "met_res",
      organ_mri_pac
    ), # exposure
    c(
      setdiff(
        hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Disease) %>% pull(),
        "High Cholesterol, High Triglycerides, Low HDL, Raised LDL"
      )
    ) %>% paste0(., "_date")
  ) %>%
  left_join(ukb_ca_death) %>%
  mutate(recruit_date = as.character(date_visit_imaging1)) %>%
  filter(!is.na(date_visit_imaging1))

tmp <- dt_aging_2_disease_imaging %>%
  filter(!is.na(Brain_res) & !is.na(Brain_mri_res) & !is.na(met_res) & !is.na(phenoage_res))

library(ComplexHeatmap)
library(circlize)

# Correlation patterns of diverse aging clocks
vars_all <- c(
  organ_pac,
  "phenoage_res",
  "kdm_res",
  "met_res",
  organ_mri_pac
)

cor_mat <- cor(
  dt_aging_2_disease_imaging[, vars_all],
  use = "complete.obs"
)

cor_tri <- cor_mat
cor_tri[upper.tri(cor_tri, diag = TRUE)] <- NA

col_fun <- colorRamp2(
  c(-0.3, 0, 0.3),
  c("#4575b4", "white", "#d73027")
)

Heatmap(
  cor_tri,
  name = "r",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  na_col = "white",
  rect_gp = gpar(col = "grey90"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- cor_tri[i, j]
    if (!is.na(val) && val > 0.1) {
      grid.text(sprintf("%.2f", val), x, y,
        gp = gpar(fontsize = 8)
      )
    }
  }
)

dt_aging_2_disease_imaging$non_ca_date <- do.call(
  pmin,
  c(
    lapply(
      dt_aging_2_disease_imaging[non_ca_date_cols],
      function(x) ifelse(x == "", NA, x)
    ),
    na.rm = TRUE
  )
)

dt_aging_2_disease_imaging <- dt_aging_2_disease_imaging %>%
  mutate(across(any_of(c("non_ca_date")), ~ ifelse(is.na(.), "", as.character(.)))) %>%
  mutate(
    non_ca_imaging = ifelse(non_ca_date <= date_visit_imaging1, 1, 0),
    ca_imaging = ifelse(ca_date <= date_visit_imaging1, 1, 0),
    non_ca_imaging = ifelse(is.na(non_ca_imaging), 0, non_ca_imaging),
    ca_imaging = ifelse(is.na(ca_imaging), 0, ca_imaging)
  )

nrow(dt_aging_2_disease_imaging) # 102730
table(dt_aging_2_disease_imaging$non_ca_imaging)
table(dt_aging_2_disease_imaging$ca_imaging)

# A, Imaging aging clocks are lack of power
cox_pac_disease_imaging <- cox(
  data = dt_aging_2_disease_imaging %>% filter(ca_imaging != 1),
  outcomes = c(

    "ca",
    paste0(ca_diseases %>% setdiff(., "MGUS"), "_reg_death"),
    "death"
  ),
  exposures = c(
    organ_mri_pac
  ),
  covar = c(
    "age",
    "age_at_imaging",
    "sex", "ethnic", "edu", "tdi", "center", "pa", "smok"
  )
)

ht_ca(cox_pac_disease_imaging$result_list)

# B, Multiomic Multiorgan aging clocks and aging hallmark disease

hallmark_base2 <- hallmark_base %>%
  left_join(hallmark_all_final %>% select("eid", organ_pac, "phenoage_res", "kdm_res", "met_res", "t1", "t2"))

cox_pac_disease <- cox(
  data = hallmark_base2 %>% filter(ca_base != 1 & non_ca_base == 0) %>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .))),
  outcomes = c(
    setdiff(
      hesin_all_clean %>%
        filter(Category != "Cancers") %>%
        distinct(Disease) %>%
        pull(Disease),
      c("High Cholesterol, High Triglycerides, Low HDL, Raised LDL", "Myasthenia Gravis")
    )
  ),
  exposures = c(
    organ_pac, "phenoage_res", "kdm_res", "met_res"
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

ht_pac_disease <- function(data) {
  tmp <- disease_map %>%
    select(Disease, Category) %>%
    filter(Category != "Cancers") %>%
    filter(Disease != "High Cholesterol, High Triglycerides, Low HDL, Raised LDL") %>%
    filter(Disease != "Myasthenia Gravis") %>%
    mutate(Category = ifelse(Category == "Haematological/ Immunological", "Haematological/\nImmunological", Category))

  cat_col <- setNames(c(pal_aaas("default")(10), c("#FFDC91ff"))[1:length(unique(tmp$Category))], as.character(unique(tmp$Category)))

  tmp$Category <- factor(
    tmp$Category,
    levels = c(
      "Infections",
      "Haematological/\nImmunological",
      "Endocrine",
      "Cardiovascular",
      "Respiratory",
      "Digestive",
      "Neurological",
      "Eye",
      "Skin",
      "Musculoskeletal",
      "Genitourinary"
    )
  )
  tmp$Disease <- factor(tmp$Disease, levels = unique(tmp$Disease))

  incident_table <- hallmark_base %>%
    summarise(across(
      all_of(
        paste0(
          hesin_all_clean %>%
            filter(Category != "Cancers") %>%
            filter(Disease != "High Cholesterol, High Triglycerides, Low HDL, Raised LDL") %>%
            filter(Disease != "Myasthenia Gravis") %>%
            distinct(Disease) %>%
            pull(),
          "_inci22"
        )
      ),
      ~ sum(.x, na.rm = TRUE)
    )) %>%
    tidyr::pivot_longer(everything(),
      names_to = "Disease",
      values_to = "Incident_cases"
    ) %>%
    mutate(Disease = gsub("_inci22$", "", Disease)) %>%
    mutate(Disease <- factor(Disease, levels = unique(tmp$Disease))) %>%
    column_to_rownames("Disease") %>%
    mutate(log2(Incident_cases))

  data <- data %>%
    mutate(Outcome = factor(Outcome, levels = tmp$Disease)) %>%
    mutate(
      Exposure = sub("_res", "", Exposure),
      Exposure = sub("_super$", "", Exposure),
      Exposure = sub("_extra$", "", Exposure),
      Exposure = case_when(
        Exposure == "phenoage" ~ "PhenoAge",
        Exposure == "kdm" ~ "KDMAge",
        Exposure == "met" ~ "MetAge",
        Exposure == "Telomere_length" ~ "Telomere length",
        TRUE ~ as.character(Exposure)
      )
    ) %>%
    mutate(Exposure = factor(Exposure, levels = pac_levels))

  mat <- xtabs(Beta ~ Outcome + Exposure, data)
  mat_HR <- xtabs(HR ~ Outcome + Exposure, data)
  mat[abs(mat) == 0] <- NA
  mat_p <- xtabs(P_Value ~ Outcome + Exposure, data)
  mat_q <- xtabs(FDR_P_Value ~ Outcome + Exposure, data)

  mat[(mat_p) >= 0.05] <- 0
  mat[(mat_p) < 0.05 & (mat_q) >= 0.05 & (mat) < 0] <- 0

  col_fun = colorRamp2(c(-1, 0, 1), c("#007F7F", "white", "#D85F4D"))

  ht_opt$ROW_ANNO_PADDING = unit(0.6, "mm")

  Heatmap(mat,
    col = col_fun, na_col = "white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    cluster_rows = F, cluster_columns = F,
    row_names_side = "left", row_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_names_side = "bottom", show_column_names = T, column_names_gp = gpar(fontsize = 6),
    row_split = factor(tmp$Category),
    row_title_gp = gpar(fontsize = 5, fontface = "bold"),
    column_split = factor(c(rep("Proteomic", each = 16), rep("Metabolomic", each = 1), rep("Clinical", each = 2)),
      levels = c("Proteomic", "Metabolomic", "Clinical")
    ),
    column_title_gp = gpar(fontsize = 5, fontface = "bold"),
    column_title_side = "bottom",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (mat_q[i, j] < 0.05) {
        grid.text(sprintf("%.1f", mat_HR[i, j]), x, y, gp = gpar(fontsize = 5))
      } else if ((mat_p)[i, j] < 0.05) {
        if ((mat_HR)[i, j] > 1) {
          grid.text("+", x, y, gp = gpar(fontsize = 4))
        }
      }
    },
    heatmap_legend_param = list(
      title = "Log (HR)",
      legend_height = unit(1, "cm"),
      grid_width = unit(0.5, "mm"),
      labels_gp = gpar(col = "black", fontsize = 5),
      title_gp = gpar(col = "black", fontsize = 5)
    ),
    left_annotation = rowAnnotation(
      " " = anno_simple(
        as.character(tmp$Category),
        col = cat_col,
        width = unit(0.03, "cm")
      )
    ),
    right_annotation = rowAnnotation(
      "Incident\ncases" = anno_barplot(
        incident_table$Incident_cases,
        border = FALSE,
        ylim = c(0, 60000),
        gp = gpar(fill = "#3c5488ff", col = "#3c5488ff"),
        axis_param = list(
          gp = gpar(fontsize = 5, col = "black", lwd = 0.2),
          at = c(0, 30000, 60000),
          labels = c("0", "30k", "60k"),
          labels_rot = 0
        )
      ), width = unit(1, "cm"),
      annotation_name_gp = gpar(fontsize = 6)
    ),

    column_title = NULL
  ) %>% draw()
}

ht_pac_disease(cox_pac_disease$result_list)

cairo_pdf("./Project/PAC and cancer/Figures/fig.1/reg/BA_2_disease.pdf",
  width = 4.5, height = 6
)
ht_pac_disease(cox_pac_disease$result_list)
dev.off()

### C.2 BA to disease in CA patients----
# Note: too few cases
table(hallmark_base2$ca_base == 1 & hallmark_base2$non_ca_base == 0) # 28043
table(hallmark_base2$ca_base_any == 1 & hallmark_base2$non_ca_base == 0) # 28043

cox_pac_disease_ca <- cox(
  data = hallmark_base2 %>% filter(ca_base_any == 1 & non_ca_base == 0) %>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .))),
  outcomes = c(setdiff(
    hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Disease) %>% pull(),
    "High Cholesterol, High Triglycerides, Low HDL, Raised LDL"
  )),
  exposures = c(organ_pac, "phenoage_res", "kdm_res", "met_res"),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)
ht_pac_disease(cox_pac_disease_ca$result_list)

aging_disease_inca_results <- cox_pac_disease_ca$result_list %>%
  filter(Exposure == "phenoage_res") %>%
  left_join(hesin_all_clean %>% filter(Category != "Cancers") %>% distinct(Disease, .keep_all = T) %>% select(Disease, Category), by = c("Outcome" = "Disease")) %>%
  mutate(sig = ifelse(FDR_P_Value >= 0.05, "F", "T")) %>%
  group_by(Category) %>%
  arrange(desc(-log10(FDR_P_Value))) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(
    label = ifelse(rank <= 3 & sig == "T", Outcome, ""),
    Category = ifelse(Category == "Haematological/ Immunological", "Haematological/\nImmunological", Category)
  )

cat_levels <- levels(factor(aging_disease_inca_results$Category))

bar_df <- data.frame(
  Category = cat_levels,
  group = cat_levels,
  x = seq_along(cat_levels)
)

fig2d <- ggplot(aging_disease_inca_results, aes(x = Category, y = pmin(-log10(FDR_P_Value), 100))) +
  geom_point(aes(size = abs(Beta), color = Category, shape = sig),
    alpha = 0.7, stroke = 0,
    position = position_jitter(width = 0.3)
  ) +
  scale_shape_manual(values = c("F" = 1, "T" = 16), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", lwd = 0.2) +
  geom_text_repel(aes(label = label),
    size = 1.5,
    box.padding = 0.1,
    point.padding = 0.1, max.overlaps = 50
  ) +
  scale_color_manual(values = c(
    "Cardiovascular" = "#A52A2A",
    "Neurological" = "#4B0082",
    "Respiratory" = "#FF7F50",
    "Haematological/\nImmunological" = "#FFFF00",
    "Musculoskeletal" = "#808080",
    "Genitourinary" = "#808000",
    "Digestive" = "#D2B48C",
    "Eye" = "#00FFFF",
    "Infections" = "#0000FF",
    "Endocrine" = "#008000",
    "Skin" = "black"
  )) +
  labs(x = "", y = expression(-log[10] * "(q)"), title = "") +
  scale_size_continuous(range = c(0.5, 2), name = "Log (HR)", ) +
  theme_bw(base_size = 6, base_line_size = 0.1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(color = "black", size = 0.1),
    legend.key.width = unit(0.1, "cm"),
    axis.text.x = element_text(color = "black", size = 5, angle = 45, vjust = 0.5),
    axis.text.y = element_text(color = "black", size = 5),
    legend.key.size = unit(0.1, "cm"), # Adjust key size to reduce spacing
    legend.spacing = unit(0, "cm"), # Adjust spacing between legend items
    legend.key.spacing.y = unit(0, "pt"),
    legend.box.spacing = margin(0),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    strip.background = element_blank(),
    legend.title = element_text(size = 5),
  ) +
  scale_y_continuous(
    limits = c(
      -4, 101
    ), expand = expansion(add = c(0, 0.5))
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = bar_df,
    aes(x = x, y = -3.5, fill = group),
    height = 1,
    width = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Cardiovascular" = "#A52A2A",
      "Neurological" = "#4B0082",
      "Respiratory" = "#FF7F50",
      "Haematological/\nImmunological" = "#FFFF00",
      "Musculoskeletal" = "#808080",
      "Genitourinary" = "#808000",
      "Digestive" = "#D2B48C",
      "Eye" = "#00FFFF",
      "Infections" = "#0000FF",
      "Endocrine" = "#008000",
      "Skin" = "black"
    ),
    guide = "none"
  ) +
  coord_cartesian(clip = "off")

### C.3 BA to disease (all) in CA patients----
cox_pac_disease_ca_all <- cox(
  data = hallmark_base2 %>% filter(ca_base_any == 1 & non_ca_base == 0),
  outcomes = c(
    labels_nca[str_detect(labels_nca, "_1")],
    "any_1", "any_2", "any_3", "any_4", "death"
  ),
  exposures = c(organ_pac, "phenoage_res", "kdm_res", "met_res"),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")
)

aging_disease_inca_all_results <- cox_pac_disease_ca_all$result_list %>%
  mutate(
    Outcome = gsub("_1_nca", "", Outcome),
    Exposure = gsub("_res", "", Exposure)
  ) %>%
  mutate(
    Exposure = ifelse(Exposure == "kdm", "KDMAge", Exposure),
    Exposure = ifelse(Exposure == "phenoage", "PhenoAge", Exposure),
    Exposure = ifelse(Exposure == "met", "MetAge", Exposure)
  ) %>%
  mutate(Outcome = case_when(
    Outcome == "any_1" ~ "Any (1st)",
    Outcome == "any_2" ~ "Any (2nd)",
    Outcome == "any_3" ~ "Any (3rd)",
    Outcome == "any_4" ~ "Any (4th)",
    Outcome == "death" ~ "Death",
    TRUE ~ Outcome
  )) %>%
  mutate(
    Exposure = factor(Exposure, levels = pac_levels),
    Outcome = factor(Outcome, levels = c("Any (1st)", "Any (2nd)", "Any (3rd)", "Any (4th)", "GI", "TA", "EA", "LP", "DNS", "CS", "MD", "SCE", "AIC", "Death"))
  )

ht_ca(aging_disease_inca_all_results)

# Add interaction between aging and death/disease

result_list_death <- cox(
  data = hallmark_base2,
  exposures = "bmi",
  outcomes = c(
    "any_1", "any_2", "any_3", "any_4",
    "death",
    labels_nca[str_detect(labels_nca, "_1")]
  ),
  covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok")
)

table(result_list_death$data$ca_base_any) # 39470
table(result_list_death$data$ca_base) # 30775

cox_intearction <- function(var_list) {
  library(glue)
  library(survival)

  outcomes <- list(
    death = c("death_fu", "death"),
    any1  = c("any_1_fu", "any_1_inci"),
    any2  = c("any_2_fu", "any_2_inci"),
    any3  = c("any_3_fu", "any_3_inci"),
    any4  = c("any_4_fu", "any_4_inci"),
    AIC   = c("AIC_1_nca_fu", "AIC_1_nca_inci"),
    GI    = c("GI_1_nca_fu", "GI_1_nca_inci"),
    TA    = c("TA_1_nca_fu", "TA_1_nca_inci"),
    EA    = c("EA_1_nca_fu", "EA_1_nca_inci"),
    LP    = c("LP_1_nca_fu", "LP_1_nca_inci"),
    DNS   = c("DNS_1_nca_fu", "DNS_1_nca_inci"),
    CS    = c("CS_1_nca_fu", "CS_1_nca_inci"),
    MD    = c("MD_1_nca_fu", "MD_1_nca_inci"),
    SCE   = c("SCE_1_nca_fu", "SCE_1_nca_inci")
  )

  p_ls <- data.frame()

  for (var in var_list) {
    data <- result_list_death$data %>%
      dplyr::filter(!is.na(.data[[var]]), non_ca_base == 0)

    row_res <- list(var = var)

    for (nm in names(outcomes)) {
      fu <- outcomes[[nm]][1]
      inci <- outcomes[[nm]][2]
      formula_str <- glue("Surv({fu}, {inci}) ~ ca_base_any + {var} + ca_base_any * {var} + age + sex + ethnic + edu + smok + tdi + center")
      fit <- coxph(as.formula(formula_str), data = data)
      pval <- summary(fit)$coefficients[paste("ca_base_any", var, sep = ":"), "Pr(>|z|)"]

      row_res[[paste0("p_", nm)]] <- pval
    }

    p_ls <- rbind(p_ls, as.data.frame(row_res))
  }

  return(p_ls)
}

p_interaction <- cox_intearction(c(organ_pac, "met_res", "kdm_res", "phenoage_res"))
p_interaction_ordered <- p_interaction %>%
  select(c(paste0("p_", c("any1", "any2", "any3", "any4", "GI", "TA", "EA", "LP", "DNS", "CS", "MD", "SCE", "AIC", "death"))))

p_interaction_result <- p_interaction %>%
  mutate(var = gsub("_res", "", var)) %>%
  mutate(
    var = ifelse(var == "kdm", "KDMAge", var),
    var = ifelse(var == "phenoage", "PhenoAge", var),
    var = ifelse(var == "met", "MetAge", var)
  ) %>%
  mutate(var = factor(var, levels = pac_levels)) %>%
  mutate(
    death = ifelse(p_death < 0.05, 1, 0),
    any = ifelse(
      rowSums(across(3:15, ~ .x < 0.05), na.rm = TRUE) > 0,
      1, 0
    )
  )

ht_ca_ungroup <- function(data) {
  ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
  mat <- xtabs(Beta ~ Outcome + Exposure, data) %>% t()
  mat_HR <- xtabs(HR ~ Outcome + Exposure, data) %>% t()
  mat_p <- xtabs(FDR_P_Value ~ Outcome + Exposure, data) %>% t()
  interaction_mat <- p_interaction_ordered
  col_fun = colorRamp2(c(-1, 0, 1), c("#4DBBD5B2", "white", "#DC0000B2"))

  mat[(mat_p) > 0.05] <- 0

  Heatmap(mat,
    col = col_fun, na_col = "white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black", lwd = 0.2), rect_gp = gpar(col = "black", lwd = 0.2),
    cluster_rows = F, cluster_columns = F,
    width = ncol(mat) * unit(3, "mm"),
    height = nrow(mat) * unit(2, "mm"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_names_side = "bottom", show_column_names = T, column_names_gp = gpar(fontsize = 6),
    column_split = factor(c(rep("Overall", each = 4), rep(" ", each = 9), rep("", each = 1)), levels = c("Overall", " ", "")),
    column_title = NULL,
    column_title_gp = gpar(fontsize = 6, fontface = "bold"),
    column_title_side = "top",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (mat_p[i, j] < 0.05) {
        grid.text(sprintf("%.1f", mat_HR[i, j]), x, y, gp = gpar(fontsize = 5))
      }

      if (mat_p[i, j] < 0.05 & interaction_mat[i, j] < 0.05) {
        grid.rect(x, y, width, height,
          gp = gpar(col = "black", lwd = 0.8, fill = NA)
        )
      }
    },
    heatmap_legend_param = list(
      at = c(-1, 0, 1),
      title = "Log (HR)",
      legend_height = unit(0.6, "cm"),
      grid_width = unit(0.5, "mm"),
      labels_gp = gpar(col = "black", fontsize = 4.5),
      title_gp = gpar(col = "black", fontsize = 4.5)
    ),
    right_annotation = rowAnnotation(
      `Interaction (death)` = (p_interaction_result$death),
      `Interaction (any)` = (p_interaction_result$any),

      gp = gpar(col = "black", lwd = 0.3),
      gap = unit(1, "mm"),
      col = list(
        `Interaction (death)` = colorRamp2(c(0, 1), c("white", "#4575b4")),
        `Interaction (any)` = colorRamp2(c(0, 1), c("white", "#4575b4"))
      ),
      width = unit(0.3, "cm"),
      simple_anno_size_adjust = TRUE,
      show_legend = c(
        `Interaction (death)` = FALSE, `Interaction (any)` = FALSE
      ),
      annotation_name_gp = gpar(fontsize = 4.5),
      annotation_name_rot = 45
    ),

    row_title = NULL,
  ) %>% return()
}

ht_ca_ungroup(aging_disease_inca_all_results)

# Combined <Figure 2>---------------
fig2a <- grid.grabExpr({
  ht <- aging_cancer_ht(aging_cancer_result) %>% draw()
  for (i in seq_len(3)) {
    decorate_annotation("FDR hits",
      {
        grid.rect(gp = gpar(col = "black", lwd = 0.2, fill = NA))
      },
      slice = i
    )
  }
})

fig2b <- grid.grabExpr(
  draw(h1 + h2,
    ht_gap = unit(1, "mm"), auto_adjust = T,
    column_title = "Prevalent cancer and multi-omics, multi-organ aging hallmarks and mortality",
    column_title_gp = gpar(fontsize = 6, fontface = "bold")
  )
)

fig2c <- grid.grabExpr(
  {
    ht <- ht_ca_ungroup(aging_disease_inca_all_results)
    draw(ht)
  },
  wrap = TRUE
)

cairo_pdf("./Project/PAC and cancer/Figures/fig.1/reg/Fig.2.pdf",
  width = 7.2, height = 7.8
)

plot_grid(fig2a, fig2b,
  plot_grid(fig2c, fig2d, nrow = 1, rel_widths = c(1, 1.7), labels = c("c", "d"), label_size = 8),
  ncol = 1, nrow = 3, labels = c("a", "b", ""), label_size = 8, rel_heights = c(1, 1, 1)
) +
  theme(plot.margin = margin(0.3, 0.1, 0, 0.1, "cm"))

dev.off()

#### Traj analysis----
# Aging trajectory analysis by aging clocks and cancer types

traj <- hallmark_all_final %>%
  select(any_of(c(
    "eid", "t1", "t2",
    c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi", "non_ca_base"), # covar
    c(organ_pac, "phenoage_res", "kdm_res", "met_res"), # exposure
    date_cols_hesin
  ))) %>%
  left_join(ukb_ca_death)

# Get disease interval
# (prevalent diseases: + days)
# (incident disease: -days)
traj_disease_prefixes <- c(
  "ca",
  hesin_all_clean %>% filter(Category == "Cancers") %>% distinct(Disease) %>% pull() %>% setdiff(c("Liver Mets.", "Lymph Node Mets."))
)

for (prefix in traj_disease_prefixes) {
  interval_col <- paste0(prefix, "_interval")
  disease_col <- paste0(prefix, "_con")
  base_col <- paste0(prefix, "_base")
  inci_col <- paste0(prefix, "_inci")
  date_col <- paste0(prefix, "_reg_death_date")

  traj <- traj %>%
    mutate(
      !!interval_col := if_else(
        !!sym(date_col) != "",
        as.numeric(difftime(as.Date(recruit_date), as.Date(!!sym(date_col)), units = "days")),
        0
      ),
      !!disease_col := if_else(!!sym(date_col) != "", 1, 0) %>% as.factor(),
      !!base_col := if_else(!!sym(interval_col) > 0, 1, 0) %>% as.factor(),
      !!inci_col := if_else(!!sym(interval_col) < 0, 1, 0) %>% as.factor()
    )
}

# Truncated to 40 years
traj <- traj %>% mutate(across(contains("_interval"), ~ ifelse(. / 365.25 >= 40, 40, . / 365.25)))

# Match people for those with incident disease
case_control_match <- function(data, disease_status, disease_interval) {
  # Subset cases (individuals where disease_status_var == 1)
  cases <- data[data[[disease_status]] == 1, ]

  # Initialize an empty list to store matched controls
  controls <- list()

  # Loop through each case
  for (i in 1:nrow(cases)) {
    # Extract the case information
    case <- cases[i, ]
    age_range <- case$age + c(-2, 2)
    sex <- case$sex

    # Subset potential controls based on age and sex
    potential_controls <- data[data[[disease_status]] == 0 &
      data$age >= age_range[1] &
      data$age <= age_range[2] &
      data$sex == sex, ]

    # Randomly select 5 controls from potential_controls
    matched_controls <- potential_controls[sample(nrow(potential_controls), 5), ]

    # Assign the time_to_disease of cases to matched controls

    matched_controls[[disease_interval]] <- case[[disease_interval]]
    # Add matched controls to the list
    controls[[i]] <- matched_controls
  }

  # Combine all matched controls into a single data frame
  matched_controls_df <- do.call(rbind, controls)

  # Return the matched controls data frame
  return(matched_controls_df)
}

control_cancer <- case_control_match(traj, "Cancers_inci", "Cancers_interval")
matched_cancer <- rbind(traj %>% filter(Cancers_inci == 1 | Cancers_base == 1), control_cancer) %>%
  mutate(Cancers_con = factor(Cancers_con, levels = c(0, 1), labels = c("Control", "Cancers")))

cancer_vars <- paste0(traj_disease_prefixes, "_interval")
traj_long <- traj %>%
  pivot_longer(
    cols = all_of(cancer_vars),
    names_to = "Cancer_type",
    values_to = "Cancer_interval"
  ) %>%
  mutate(Cancer_type = sub("_interval$", "", Cancer_type))

table(traj_long$Cancer_type)

# Any cancer type and phenoAge

traj$ca_interval

ggplot(
  traj,
  aes(x = `ca_interval`, y = phenoage_res),
) +
  geom_smooth(method = "loess", se = F, size = 0.2) + # Linear regression lines with confidence intervals
  labs(
    x = "Time since diagnosis",
    y = "Age gap",
    color = NULL, # Remove legend title
    fill = NULL # Remove fill legend title
  ) +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = alpha(c("#4292c6", "#fc9272"), 0.3)) +
  theme_classic(base_size = 6, base_line_size = 0.1) +
  theme(
    legend.position = c(0.3, 0.92),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2))

### Saving working image--------
save.image("./Project/PAC and cancer/Data/UKB/ukb_aging_ca_fig1.RData")

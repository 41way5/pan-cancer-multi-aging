#Contribution of age+sex, PRS, and aging exposome to pan cancers-------------------
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
  library(CoxR2)
  
  linear <- function(data, a, b,covar=NULL) {
    
    betas <- c()
    ses <- c()
    ps <- c()
    var_a_list <- c()
    var_b_list <- c()
    
    # Iterate over all combinations of variables in lists 'a' and 'b'
    for (var_a_name in a) {
      for (var_b_name in b) {
        
        formula_str <- paste(var_a_name, "~", var_b_name)
        if (!is.null(covar)) {
          formula_str <- paste(formula_str, paste(covar, collapse = " + "), sep = " + ")
        }
        formula <- as.formula(formula_str)
        model <- lm(formula, data = data)
        
        # Extract coefficient, standard error, and p-value
        beta <- coef(model)[2]
        se <- sqrt(diag(vcov(model)))[2]
        p <- summary(model)$coefficients[2, 4]
        
        # Append results to lists
        betas <- c(betas, beta)
        ses <- c(ses, se)
        ps <- c(ps, p)
        var_a_list <- c(var_a_list, var_a_name)
        var_b_list <- c(var_b_list, var_b_name)
      }
    }
    
    # Create a data frame with the results
    result <- data.frame(
      var_a = var_a_list,
      var_b = var_b_list,
      beta = betas,
      se = ses,
      p = ps
    )
    
    return(result)
  }
  select <- dplyr::select
}
# Set the FIG4_* environment variables to local paths before running.
# FIG4_FIG2_WORKSPACE should match FIG2_OUTPUT_WORKSPACE from fig2.R.
# FIG4_PRS_SCORE_DIR contains the phenotype subdirectories with .sscore files.
# Set both output file paths in existing directories outside the repository.
fig4_path <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path)) {
    stop(sprintf("Set %s to a local path before running fig4.R.", variable),
         call. = FALSE)
  }
  path.expand(path)
}

load(fig4_path("FIG4_FIG2_WORKSPACE"))


#Contribution of age+sex, PRS, and aging exposome to pan cancers-------------------
# Get PRS data
#From UKB standard set
#OR
#Construct manually

#Get PRS for diseases and trais----- 
# a, standard set using external SNPs
{
  ukb_prs <- read.csv(file = fig4_path("FIG4_PRS_CSV"),check.names = FALSE)
  ukb_pc <- read.csv(file = fig4_path("FIG4_PC_CSV"),check.names = FALSE)
  
  #rename PRS colname
  labs <- colnames(ukb_prs)
  prs_idx <- grepl("^Standard PRS for ", labs)
  
  colnames(ukb_prs)[prs_idx] <- paste0(
    "PRS_",
    sub(".*\\(([^()]*)\\)$", "\\1", labs[prs_idx])
  )
  
  for (i in which(prs_idx)) {
    attr(ukb_prs[[i]], "label") <- sub("^Standard PRS for (.*) \\([^()]*\\)$", "\\1", labs[i])
  }
  
  #rename PC colname
  colnames(ukb_pc) <- sub("^p22009_a([0-9]+)$", "pc\\1", colnames(ukb_pc))
  
  
  
  ukb_prs <- ukb_prs %>% dplyr::rename(eid=`Participant ID`,gen_batch=`Genotype measurement batch`) %>% mutate(prs_in=1) %>% 
    mutate(gen_batch = ifelse(str_detect(gen_batch, "UKBiLEVEAX"), "BiLEVE", "Axiom")) %>% 
    left_join(ukb_pc %>% select(eid, pc1:pc40), by = "eid")
  
}
# b, manually constructed based on PGC

{
  prs_list <- bind_rows(fread(file = fig4_path("FIG4_PRS_LIST"), header = FALSE, col.names = c("pgs_id", "pheno")), 
                        fread(file = fig4_path("FIG4_PRS_LIST2"), header = FALSE, col.names = c("pgs_id", "pheno"))
  )
  
  prs_add <- lapply(seq_len(nrow(prs_list)), function(i) {
    pheno <- prs_list$pheno[i]
    f <- file.path(fig4_path("FIG4_PRS_SCORE_DIR"), pheno, paste0("prs_", pheno, ".sscore"))
    fread(file = f) %>%
      select(IID, SCORE1_AVG) %>%
      dplyr::rename(eid = IID) %>%
      dplyr::rename(!!paste0("PRS_", pheno) := SCORE1_AVG)
  }) %>%
    Reduce(function(x, y) full_join(x, y, by = "eid"), .)
  
  #Combined all PRS 
  ukb_prs_all <- ukb_prs %>% left_join(prs_add) %>% map_eid %>% filter(!is.na(eid))
  
}
#Figure 4 part A-------------
#For cancers
tmp <- cox_pac_ca$data %>% 
  left_join(ukbb_add_outcome2 %>% select(eid,smoke_cat,alcohol_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat)) %>% 
  mutate(pa_cat=ifelse(pa==0,1,0),obe_cat=ifelse(bmi>=30,1,0)) %>% 
  mutate(lifestyle=rowSums(dplyr::select(., smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,obe_cat,fruit_cat,fish_cat,red_cat,process_cat)),
         life_cat=ifelse(lifestyle<=2,"Favorable",ifelse(lifestyle<=4,"Intermediate","Unfavorable"))
  ) %>% 
  left_join(ukb_prs_all %>% select(-(pc11:pc40)), by = "eid") %>% 
  mutate(across(starts_with("PRS_"), ~ as.numeric(scale(.))))



#P3
get_block_wald_all <- function(data, diseases, rhs = NULL) {
  
  res_all <- data.frame()
  
  lifestyle_var <- c("smoke_cat", "alcohol_cat", "pa_cat", "obe_cat", "tv_cat",
                     "sleep_cat", "fruit_cat", "fish_cat", "red_cat", "process_cat")
  
  prs_map <- c(
    "Breast Cancer" = "PRS_BC",
    "Bowel Cancer"  = "PRS_CRC",
    "Ovarian Cancer"  = "PRS_EOC",
    "Prostate Cancer"  = "PRS_PC",
    "NHL"="PRS_NHL",
    "Uterine Cancer"="PRS_Uterine_Cancer",
    "Skin Cancer" ="PRS_Skin_Cancer",
    "Oesophageal Cancer"="PRS_Esophageal_Cancer",
    "Pancreatic Cancer"="PRS_Pancreatic_Cancer",
    "Melanoma"="PRS_Melanoma",
    "Bladder Cancer"="PRS_Bladder_Cancer",
    "Lung Cancer"="PRS_Lung_Cancer",
    "Thyroid Cancer"="PRS_Thyroid_Cancer",
    "Kidney Cancer"="PRS_Kidney_Cancer",
    "Leukaemia" ="PRS_Leukemia",
    "Polycythaemia Vera"="PRS_P_vera",
    "Stomach Cancer"="PRS_Stomach_Cancer",
    "Plasma Cell Cancer"="PRS_PCC",
    "Brain Cancer"="PRS_Brain_Cancer"
  )
  
  
  for (d in diseases) {
    
    time_col   <- paste0(d, "_fu")
    status_col <- paste0(d, "_inci")
    base_col <- paste0(d, "_base2")
    
    if (d == "death") {
      prs_var <- colnames(ukb_prs)[2:40]
    } else {
      prs_var <- unname(prs_map[sapply(names(prs_map), grepl, x = d)])[1]
      if (length(prs_var) == 0 || is.na(prs_var)) prs_var <- NULL
    }
    
    if (!is.null(prs_var)) {
      prs_var <- c(prs_var, colnames(ukb_prs)[c(41,43:52)])
    }
    
    vars_all <- na.omit(c("age", "sex", prs_var, lifestyle_var,
                          "phenoage_res", "kdm_res", "met_res","non_ca_base"))
    
    data2 <- data %>% select(all_of(vars_all),time_col,status_col,base_col) %>%  
      filter(if_all(all_of(vars_all), ~ !is.na(.)))
    
    if (d %in% paste0(c("Uterine Cancer", "Ovarian Cancer", "Breast Cancer"), "_reg_death")) {
      data2 <- data2 %>% filter(sex == 0)
    }
    
    if (d %in% paste0("Prostate Cancer", "_reg_death")) {
      data2 <- data2 %>% filter(sex == 1)
    }
    
    rhs <- paste(
      na.omit(c("age", "sex", prs_var, lifestyle_var,
                "phenoage_res", "kdm_res", "met_res","non_ca_base")),
      collapse = " + "
    )
    
    form_txt <- paste0(
      "Surv(`", time_col, "`, `", status_col, "` == 1) ~ ", rhs
    )
    
    fit <- coxph(
      as.formula(form_txt),
      data = data2
    )
    
    aov_tab <- as.data.frame(anova(fit, test = "Chisq"))
    aov_tab$term <- rownames(aov_tab)
    
    get_term_chisq <- function(term_names) {
      sum(aov_tab$Chisq[aov_tab$term %in% term_names], na.rm = TRUE)
    }
    
    block_chisq <- c(
      age_sex   = get_term_chisq(c("age", "sex")),
      prs       = if (is.null(prs_var)) 0 else get_term_chisq(prs_var),
      lifestyle = get_term_chisq(lifestyle_var),
      aging     = get_term_chisq(c("phenoage_res", "kdm_res", "met_res", "non_ca_base"))
    )
    
    total_chisq <- sum(block_chisq, na.rm = TRUE)
    
    res <- data.frame(
      block = names(block_chisq),
      wald_chisq = as.numeric(block_chisq),
      total_chisq = total_chisq,
      prop = block_chisq / total_chisq,
      disease = d
    )
    
    res_all <- rbind(res_all, res)
  }
  
  rownames(res_all) <- NULL
  res_all
}


wald_all <- get_block_wald_all(
  data = tmp,
  diseases = c(paste0(c(ca_diseases),"_reg_death"),"death")
)

p3 <- {
  plot_dat <- wald_all %>%
  mutate(
    disease = gsub("_reg_death", "", disease),
    disease = ifelse(disease == "death", "All-cause mortality", disease),
    block = factor(block, levels = c("aging", "lifestyle", "prs", "age_sex"))
  )

disease_order <- plot_dat %>%
  group_by(disease) %>%
  summarise(age_related = sum(prop[block == "age_sex"])) %>%
  arrange(age_related) %>%
  pull(disease)

plot_dat <- plot_dat %>%
  mutate(disease = factor(disease, levels = disease_order))


ggplot(plot_dat, aes(x = disease, y = prop, fill = block)) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c(
      "age_sex" = "#122356",
      "prs" = "#EEDC92",
      "lifestyle" = "#75BABB",
      "aging" = "#882C23"
    ),
    breaks = c("age_sex", "prs", "lifestyle", "aging"),
    labels = c(
      "age_sex" = "Age and sex",
      "prs" = "PRS",
      "lifestyle" = "Lifestyle",
      "aging" = "Aging"
    )
  )+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  labs(
    x = NULL,
    y = expression("Relative proportion of total " * chi^2),
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )
}



get_r2_tab <- function(data, diseases, fit_list) {
  
  r2_tab <- data.frame()
  lifestyle_var = c("smoke_cat", "alcohol_cat", "pa_cat", "obe_cat", "tv_cat", "sleep_cat", "fruit_cat", "fish_cat", "red_cat", "process_cat")
  
  for (d in diseases) {
    
    fu_col   <- paste0(d, "_fu")
    inci_col <- paste0(d, "_inci")
    
    prs_map <- c(
      "Breast Cancer" = "PRS_BC",
      "Bowel Cancer"  = "PRS_CRC",
      "Ovarian Cancer"  = "PRS_EOC",
      "Prostate Cancer"  = "PRS_PC",
      "NHL"="PRS_NHL",
      "Uterine Cancer"="PRS_Uterine_Cancer",
      "Skin Cancer" ="PRS_Skin_Cancer",
      "Oesophageal Cancer"="PRS_Esophageal_Cancer",
      "Pancreatic Cancer"="PRS_Pancreatic_Cancer",
      "Melanoma"="PRS_Melanoma",
      "Bladder Cancer"="PRS_Bladder_Cancer",
      "Lung Cancer"="PRS_Lung_Cancer",
      "Thyroid Cancer"="PRS_Thyroid_Cancer",
      "Kidney Cancer"="PRS_Kidney_Cancer",
      "Leukaemia" ="PRS_Leukemia",
      "Polycythaemia Vera"="PRS_P_vera",
      "Stomach Cancer"="PRS_Stomach_Cancer",
      "Plasma Cell Cancer"="PRS_PCC",
      "Brain Cancer"="PRS_Brain_Cancer"
    )
    
    
    if (d == "death") {
      prs_var <- colnames(ukb_prs)[2:40]
    } else {
      prs_var <- unname(prs_map[sapply(names(prs_map), grepl, x = d)])[1]
      if (length(prs_var) == 0 || is.na(prs_var)) prs_var <- NULL
    }
    
    if (!is.null(prs_var)) {
      prs_var <- c(prs_var, colnames(ukb_prs)[c(41,43:52)])
    }
    
    # fit_list <- list(
    #   model1 = paste(na.omit(c("sex", prs_var)), collapse = " + "), 
    #   # Model 1: sex + PRS
    #   model2 = paste(na.omit(c("sex", prs_var, lifestyle_var)), collapse = " + "), 
    #   # Model 2: Model 1 + lifestyle
    #   model3 = paste(na.omit(c("sex", prs_var, lifestyle_var, "age")), collapse = " + "), 
    #   # Model 3: Model 2 + age
    #   model4 = paste(na.omit(c("age", prs_var, "sex", lifestyle_var, "phenoage_res", "kdm_res", "met_res", "non_ca_base")), collapse = " + ") 
    #   # Model 4: Model 3 + phenotypic/metabolic aging + baseline non-cancer disease
    # )
    
    
    fit_list <- list(
      model1 = paste(na.omit(c("age","sex")), collapse = " + "), 
      # Model 1: age + sex
      model2 = paste(na.omit(c("age","sex",prs_var)), collapse = " + "), 
      # Model 2: Model 1 + PRS
      model3= paste(na.omit(c("age","sex", prs_var, lifestyle_var)), collapse = " + "), 
      # Model 3: Model 2 + lifestyle
      model4 = paste(na.omit(c("age", "sex", prs_var, lifestyle_var, "phenoage_res", "kdm_res", "met_res", "non_ca_base")), collapse = " + ") 
      # Model 4: Model 3 + phenotypic/metabolic aging + baseline non-cancer disease
    )
    
    vars_all <- na.omit(c("sex", prs_var, lifestyle_var, "age", "phenoage_res", "kdm_res", "met_res", "non_ca_base"))
    data2 <- data %>% filter(if_all(all_of(vars_all), ~ !is.na(.)))
    
    if (d %in% paste0(c("Uterine Cancer", "Ovarian Cancer", "Breast Cancer"),"_reg_death")) {
      data2 <- data2 %>% filter(sex == 0)
    }
    
    if (d %in% paste0(c("Prostate Cancer"),"_reg_death")) {
      data2 <- data2 %>% filter(sex == 1)
    }
    
    for (m in names(fit_list)) {
      
      form_txt <- paste0(
        "Surv(`", fu_col, "`, `", inci_col, "` == 1) ~ ",
        fit_list[[m]]
      )
      
      fit <- coxph(
        as.formula(form_txt),
        data = data2,
        na.action = na.omit
      )
      
      x <- coxr2(fit)
      
      r2_tab <- rbind(
        r2_tab,
        data.frame(
          disease = d,
          model   = m,
          N = nrow(model.frame(fit)),
          nevent  = x$nevent,
          R2      = x$rsq
        )
      )
    }
  }
  
  rownames(r2_tab) <- NULL
  r2_tab
}

r2_tab <- get_r2_tab(
  data = tmp,
  diseases = c(paste0(c(ca_diseases),"_reg_death"),"death")
)

p1 <- {
r2_tab_plot <- r2_tab %>%
  mutate(
    disease = ifelse(disease == "death", "All-cause mortality", disease),
    disease = factor(gsub("_reg_death", "", disease), levels = disease_order),
    model = factor(model, levels = c("model1", "model2", "model3", "model4"))
  )

ggplot(r2_tab_plot, aes(x = disease, y = R2, group = model, color = model)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = c(
      model1 = "#122356",
      model2 = "#EEDC92",
      model3 = "#75BABB",
      model4 = "#882C23"
    ),
    labels = c(
      model1 = "Age and sex",
      model2 = "Age, sex, and PRS",
      model3 = "Age, sex, PRS, and lifestyle",
      model4 = "Age, sex, PRS, lifestyle, and aging"
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))
  ) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  labs(x = NULL, y = expression("Model " * italic(R)^2), color = NULL)+
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.margin = margin(b = -5),
    legend.box.margin = margin(b = -5),
    legend.key.width = unit(0.2, "cm"),
  )

}

p2 <- {p2_dat <- plot_dat %>%
  filter(block %in% c("aging", "lifestyle")) %>%
  select(disease, block, prop) %>%
  tidyr::pivot_wider(names_from = block, values_from = prop) %>%
  mutate(prop = aging / (aging + lifestyle))

ggplot(p2_dat, aes(x = disease, y = prop)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.2, color = "grey40") +
  geom_segment(aes(xend = disease, y = 0, yend = prop), linewidth = 0.2) +
  geom_point(size = 0.7, color = "#882C23") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1), expand = c(0, 0),
                     labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))) +
  theme_classic(base_size = 6, base_line_size = 0.2) +
  labs(
    x = NULL,
    y = "Aging/\nmodifiable risk"
  ) +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x = element_blank(),axis.ticks.x = element_blank()
  )
}

f3a <- plot_grid(p1, p2, p3, ncol = 1, align = "v",rel_heights = c(1,0.5,1.5))


#For non-cancer chronic diseases
tmp2 <-cox_pac_disease$data %>% 
  left_join(ukbb_add_outcome2 %>% select(eid,smoke_cat,alcohol_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat)) %>% 
  mutate(pa_cat=ifelse(pa==0,1,0),obe_cat=ifelse(bmi>=30,1,0)) %>% 
  mutate(lifestyle=rowSums(dplyr::select(., smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,obe_cat,fruit_cat,fish_cat,red_cat,process_cat)),
         life_cat=ifelse(lifestyle<=2,"Favorable",ifelse(lifestyle<=4,"Intermediate","Unfavorable"))
  ) %>% 
  left_join(ukb_prs_all %>% select(-(pc11:pc40)), by = "eid") %>% 
  mutate(across(starts_with("PRS_"), ~ as.numeric(scale(.))))


non_ca_diseases <- c("LRTI",
                     "AMD","Diabetic Eye Dz","Anaemia",
                     "Dementia", "Parkinson's", 
                     "AF","Ischaemic Stroke","PAD","CHF","CHD","MI","Cardiomyopathy","Hypertension",
                     "T2DM",
                     "COPD","Pulmonary Fibrosis",
                     "Cirrhosis","Fatty Liver",
                     "RA","Gout","Osteoarthritis","Osteoporosis",
                     "CKD","Glomerulonephritis"
)

#p2&3
get_block_wald_all2 <- function(data, diseases, rhs = NULL) {
  
  res_all <- data.frame()
  
  lifestyle_var <- c("smoke_cat", "alcohol_cat", "pa_cat", "obe_cat", "tv_cat",
                     "sleep_cat", "fruit_cat", "fish_cat", "red_cat", "process_cat")
  
  prs_map <- c(
    "AMD" = "PRS_AMD","Diabetic Eye Dz"="PRS_Diabetic_Eye_Dz","Anaemia"="PRS_Anaemia","LRTI"="PRS_LRTI","Osteoarthritis"="PRS_Osteoarthritis",
    "Dementia" = "PRS_AD","Parkinson's" = "PRS_PD",
    "AF" = "PRS_AF",
    "Ischaemic Stroke" = "PRS_ISS",
    "CHD" = "PRS_CAD", "MI" = "PRS_CAD", "Cardiomyopathy"="PRS_Cardiomyopathy","CHF" = "PRS_CHF",
    "Hypertension" = "PRS_HT",
    "PAD" = "PRS_CVD",
    "T2DM" = "PRS_T2D","COPD"="PRS_COPD","Pulmonary Fibrosis"="PRS_Pulmonary_Fibrosis","Cirrhosis"="PRS_Cirrhosis","Fatty Liver"="PRS_Fatty_Liver",
    "Gout"="PRS_Gout","CKD"="PRS_CKD","Glomerulonephritis"="PRS_Glomerulonephritis",
    "Osteoporosis" = "PRS_OP",
    "RA" = "PRS_RA"
  )
  
  
  for (d in diseases) {
    
    time_col   <- paste0(d, "_fu")
    status_col <- paste0(d, "_inci")
    base_col <- paste0(d, "_base22")
    
    if (d == "death") {
      prs_var <- colnames(ukb_prs)[2:40]
    } else {
      prs_var <- unname(prs_map[sapply(names(prs_map), grepl, x = d)])[1]
      if (length(prs_var) == 0 || is.na(prs_var)) prs_var <- NULL
    }
    
    if (!is.null(prs_var)) {
      prs_var <- c(prs_var, colnames(ukb_prs)[c(41,43:52)])
    }
    
    vars_all <- na.omit(c("age", "sex", prs_var, lifestyle_var,
                          "phenoage_res", "kdm_res", "met_res"))
    
    data2 <- data %>% select(all_of(vars_all),time_col,status_col,base_col) %>%  
      filter(if_all(all_of(vars_all), ~ !is.na(.)))%>% 
      filter(.data[[base_col]] == 0)
    
    rhs <- paste(
      na.omit(c("age", "sex", prs_var, lifestyle_var,
                "phenoage_res", "kdm_res", "met_res")),
      collapse = " + "
    )
    
    form_txt <- paste0(
      "Surv(`", time_col, "`, `", status_col, "` == 1) ~ ", rhs
    )
    
    fit <- coxph(
      as.formula(form_txt),
      data = data2
    )
    
    aov_tab <- as.data.frame(anova(fit, test = "Chisq"))
    aov_tab$term <- rownames(aov_tab)
    
    get_term_chisq <- function(term_names) {
      sum(aov_tab$Chisq[aov_tab$term %in% term_names], na.rm = TRUE)
    }
    
    block_chisq <- c(
      age_sex   = get_term_chisq(c("age", "sex")),
      prs       = if (is.null(prs_var)) 0 else get_term_chisq(prs_var),
      lifestyle = get_term_chisq(lifestyle_var),
      aging     = get_term_chisq(c("phenoage_res", "kdm_res", "met_res"))
    )
    
    total_chisq <- sum(block_chisq, na.rm = TRUE)
    
    res <- data.frame(
      block = names(block_chisq),
      wald_chisq = as.numeric(block_chisq),
      total_chisq = total_chisq,
      prop = block_chisq / total_chisq,
      disease = d
    )
    
    res_all <- rbind(res_all, res)
  }
  
  rownames(res_all) <- NULL
  res_all
}


wald_all2 <- get_block_wald_all2(
  data = tmp2,
  diseases = non_ca_diseases
)
p3b  <- {
plot_dat <- wald_all2 %>%
  mutate(
    disease = gsub("_reg_death", "", disease),
    disease = ifelse(disease == "death", "All-cause mortality", disease),
    block = factor(block, levels = c("aging", "lifestyle", "prs", "age_sex"))
  )

disease_order <- plot_dat %>%
  group_by(disease) %>%
  summarise(age_related = sum(prop[block == "age_sex"])) %>%
  arrange(age_related) %>%
  pull(disease)

plot_dat <- plot_dat %>%
  mutate(disease = factor(disease, levels = disease_order))

ggplot(plot_dat, aes(x = disease, y = prop, fill = block)) +
  geom_col(width = 0.9) +
  scale_fill_manual(
    values = c(
      "age_sex" = "#122356",
      "prs" = "#EEDC92",
      "lifestyle" = "#75BABB",
      "aging" = "#882C23"
    ),
    breaks = c("age_sex", "prs", "lifestyle", "aging"),
    labels = c(
      "age_sex" = "Age and sex",
      "prs" = "PRS",
      "lifestyle" = "Lifestyle",
      "aging" = "Aging"
    )
  )+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  labs(
    x = NULL,
    y = expression("Relative proportion of total " * chi^2),
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )
}
#p2
p2b <- {
  p2_dat <- plot_dat %>%
  filter(block %in% c("aging", "lifestyle")) %>%
  select(disease, block, prop) %>%
  tidyr::pivot_wider(names_from = block, values_from = prop) %>%
  mutate(prop = aging / (aging + lifestyle))

  ggplot(p2_dat, aes(x = disease, y = prop)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.2, color = "grey40") +
  geom_segment(aes(xend = disease, y = 0, yend = prop), linewidth = 0.2) +
  geom_point(size = 0.7, color = "#882C23") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0,0.5,1), expand = c(0, 0),
                     labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))) +
  theme_classic(base_size = 6, base_line_size = 0.2) +
  labs(
    x = NULL,
    y = "Aging/\nmodifiable risk"
  ) +
  theme(
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x = element_blank(),axis.ticks.x = element_blank()
  )
}

#p1
get_r2_tab2 <- function(data, diseases, fit_list) {
  
  r2_tab <- data.frame()
  lifestyle_var = c("smoke_cat", "alcohol_cat", "pa_cat", "obe_cat", "tv_cat", "sleep_cat", "fruit_cat", "fish_cat", "red_cat", "process_cat")
  
  for (d in diseases) {
    
    fu_col   <- paste0(d, "_fu")
    inci_col <- paste0(d, "_inci")
    base_col <- paste0(d, "_base22")
    
    
    prs_map <- c(
      "AMD" = "PRS_AMD",
      "Dementia" = "PRS_AD","Parkinson's" = "PRS_PD",
      "AF" = "PRS_AF",
      "Ischaemic Stroke" = "PRS_ISS",
      "CHD" = "PRS_CAD", "MI" = "PRS_CAD", "Cardiomyopathy"="PRS_Cardiomyopathy","CHF" = "PRS_CHF",
      "Hypertension" = "PRS_HT",
      "PAD" = "PRS_CVD",
      "T2DM" = "PRS_T2D","COPD"="PRS_COPD","Pulmonary Fibrosis"="PRS_Pulmonary_Fibrosis","Cirrhosis"="PRS_Cirrhosis","Fatty Liver"="PRS_Fatty_Liver",
      "Gout"="PRS_Gout","CKD"="PRS_CKD","Glomerulonephritis"="PRS_Glomerulonephritis",
      "Osteoporosis" = "PRS_OP",
      "RA" = "PRS_RA"
    )
    
    if (d == "death") {
      prs_var <- colnames(ukb_prs)[2:40]
    } else {
      prs_var <- unname(prs_map[sapply(names(prs_map), grepl, x = d)])[1]
      if (length(prs_var) == 0 || is.na(prs_var)) prs_var <- NULL
    }
    
    if (!is.null(prs_var)) {
      prs_var <- c(prs_var, colnames(ukb_prs)[c(41,43:52)])
    }
    
    
    fit_list <- list(
      model1 = paste(na.omit(c("age","sex")), collapse = " + "), 
      # Model 1: age + sex
      model2 = paste(na.omit(c("age","sex",prs_var)), collapse = " + "), 
      # Model 2: Model 1 + PRS
      model3= paste(na.omit(c("age","sex", prs_var, lifestyle_var)), collapse = " + "), 
      # Model 3: Model 2 + lifestyle
      model4 = paste(na.omit(c("age", "sex", prs_var, lifestyle_var, "phenoage_res", "kdm_res", "met_res")), collapse = " + ") 
      # Model 4: Model 3 + phenotypic/metabolic aging + baseline non-cancer disease
    )
    
    vars_all <- na.omit(c("sex", prs_var, lifestyle_var, "age", "phenoage_res", "kdm_res", "met_res"))
    data2 <- data %>% filter(if_all(all_of(vars_all), ~ !is.na(.))) %>% 
      filter(.data[[base_col]] == 0)
    
    if (d %in% paste0(c("Uterine Cancer", "Ovarian Cancer", "Breast Cancer"),"_reg_death")) {
      data2 <- data2 %>% filter(sex == 0)
    }
    
    if (d %in% paste0(c("Prostate Cancer"),"_reg_death")) {
      data2 <- data2 %>% filter(sex == 1)
    }
    
    for (m in names(fit_list)) {
      
      form_txt <- paste0(
        "Surv(`", fu_col, "`, `", inci_col, "` == 1) ~ ",
        fit_list[[m]]
      )
      
      fit <- coxph(
        as.formula(form_txt),
        data = data2,
        na.action = na.omit
      )
      
      x <- coxr2(fit)
      
      r2_tab <- rbind(
        r2_tab,
        data.frame(
          disease = d,
          model   = m,
          N = nrow(model.frame(fit)),
          nevent  = x$nevent,
          R2      = x$rsq
        )
      )
    }
  }
  
  rownames(r2_tab) <- NULL
  r2_tab
}

r2_tab2 <- get_r2_tab2(
  data = tmp2,
  diseases = c(non_ca_diseases)
)

p1b <- {
r2_tab_plot2 <- r2_tab2 %>%
  mutate(
    # disease = ifelse(disease == "death", "All-cause mortality", disease),
    disease = factor(gsub("_reg_death", "", disease), levels = disease_order),
    model = factor(model, levels = c("model1", "model2", "model3", "model4"))
  )


ggplot(r2_tab_plot2, aes(x = disease, y = R2, group = model, color = model)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = c(
      model1 = "#122356",
      model2 = "#EEDC92",
      model3 = "#75BABB",
      model4 = "#882C23"
    ),
    labels = c(
      # model1 = "Age and sex",
      # model2 = "Age, sex, and PRS",
      # model3 = "Age, sex, PRS, and lifestyle",
      # model4 = "Age, sex, PRS, lifestyle, and aging",
      model1 = "",
      model2 = "",
      model3 = "",
      model4 = ""
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),labels = function(x) ifelse(x %in% c(0, 1), x, sprintf("%.2f", x))
  ) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  labs(x = NULL, y = expression("Model " * italic(R)^2), color = NULL)+
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(b = -5),
    legend.box.margin = margin(b = -5),
    legend.key.width = unit(0.2, "cm"),
  )
}

f3b <- plot_grid(p1b, p2b, p3b, ncol = 1, align = "v",rel_heights = c(1,0.5,1.5))

save.image(fig4_path("FIG4_OUTPUT_WORKSPACE"))

#Get stats of contribution
#Total variance 
r2_tab %>%
  filter(model == "model4") %>%
  summarise(
    R2_summary = sprintf(
      "%.2f (%.2f–%.2f)",
      median(R2, na.rm = TRUE),
      quantile(R2, 0.25, na.rm = TRUE),
      quantile(R2, 0.75, na.rm = TRUE)
    )
  )

r2_tab2 %>%
  filter(model == "model4") %>%
  summarise(
    R2_summary = sprintf(
      "%.2f (%.2f–%.2f)",
      median(R2, na.rm = TRUE),
      quantile(R2, 0.25, na.rm = TRUE),
      quantile(R2, 0.75, na.rm = TRUE)
    )
  )

wald_all %>%
  filter(block == "aging") %>%
  summarise(
    R2_summary = sprintf(
      "%.2f (%.2f–%.2f)",
      median(prop, na.rm = TRUE),
      quantile(prop, 0.25, na.rm = TRUE),
      quantile(prop, 0.75, na.rm = TRUE)
    )
  )

wald_all2 %>%
  filter(block == "aging") %>%
  summarise(
    R2_summary = sprintf(
      "%.2f (%.2f–%.2f)",
      median(prop, na.rm = TRUE),
      quantile(prop, 0.25, na.rm = TRUE),
      quantile(prop, 0.75, na.rm = TRUE)
    )
  )

#Figure 4 part B-------------
#Combined effect of lifestyle, PRS, and biological aging for cancers------
#Main data
tmp_interplay <- cox_pac_ca$data %>% #data ex ca_base_any 
  left_join(ukbb_add_outcome2 %>% select(eid,smoke_cat,alcohol_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat)) %>% 
  mutate(pa_cat=ifelse(pa==0,1,0),obe_cat=ifelse(bmi>=30,1,0)) %>% 
  mutate(lifestyle=rowSums(dplyr::select(., smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,obe_cat,fruit_cat,fish_cat,red_cat,process_cat)),
         life_cat=ifelse(lifestyle<=2,"Favorable",ifelse(lifestyle<=4,"Intermediate","Unfavorable")),
         life_cat2=ifelse(lifestyle<=2,"Low",ifelse(lifestyle<=4,"Mid","High")),
         
  ) %>% 
  left_join(ukb_prs_all %>% select(-(pc11:pc40)), by = "eid") %>% 
  mutate(across(starts_with("PRS_"), ~ as.numeric(scale(.)))) %>% 
  mutate(
    lung_prs_3 = case_when(
      PRS_Lung_Cancer <= quantile(PRS_Lung_Cancer, 0.2,na.rm=T) ~ "Low",
      PRS_Lung_Cancer > quantile(PRS_Lung_Cancer, 0.2,na.rm=T) & PRS_Lung_Cancer <= quantile(PRS_Lung_Cancer, 0.8,na.rm=T) ~ "Mid",
      PRS_Lung_Cancer > quantile(PRS_Lung_Cancer, 0.8,na.rm=T) ~ "High"
    ),
    lung_prs_cat=ifelse(PRS_Lung_Cancer<median(PRS_Lung_Cancer,na.rm=T),0,1)
  )


#KM plot
forest_km <- function(
    data = tmp_interplay,
    surv_formula = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)",
    var_forest = "lung_prs_3",
    var_km = "lung_prs_3"
){
  
  dat <- data %>%
    dplyr::filter(!is.na(.data[[var_km]])) %>%
    dplyr::mutate(
      !!var_km := as.character(.data[[var_km]])
    )
  
  if (stringr::str_detect(var_km, "prs")) {
    wanted_levels <- c("Low", "Mid", "High")
    ref_km <- "Low"
    legend.title <- "Genetic risk"
    palette_vals <- c("Low" = "#5C9E3E", "Mid" = "#8FB6D8", "High" = "#B22222")
  } else if (stringr::str_detect(var_km, "life")) {
    wanted_levels <- c("Low", "Mid", "High")
    ref_km <- "Low"
    legend.title <- "Lifestyle risk"
    palette_vals <- c("Low" = "#5C9E3E", "Mid" = "#8FB6D8", "High" = "#B22222")
  } else {
    wanted_levels <- c("Low", "Mid", "High")
    ref_km <- "Mid"
    legend.title <- "Aging risk"
    palette_vals <- c("Low" = "#5C9E3E", "Mid" = "#8FB6D8", "High" = "#B22222")
  }
  
  dat[[var_km]] <- factor(dat[[var_km]], levels = wanted_levels)
  
  f2_txt <- paste0(
    surv_formula,
    " ~ relevel(as.factor(", var_km, "), ref = '", ref_km,
    "') + age + sex + ethnic + edu + tdi + center"
  )
  f2 <- stats::as.formula(f2_txt)
  
  cox <- survival::coxph(f2, data = dat)
  
  coef_tab <- summary(cox)$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Variable") %>%
    dplyr::mutate(
      HR = `exp(coef)`,
      CI_Lower = exp(stats::confint(cox)[, 1]),
      CI_Upper = exp(stats::confint(cox)[, 2])
    ) %>%
    dplyr::filter(stringr::str_detect(Variable, var_km)) %>%
    dplyr::mutate(
      level = gsub(".*\\)", "", Variable),
      level = factor(level, levels = wanted_levels)
    ) %>%
    dplyr::arrange(level)
  
  hr_data <- data.frame(
    variable = wanted_levels,
    HR = NA_real_,
    Lower_CI = NA_real_,
    Upper_CI = NA_real_,
    P_value = NA_real_,
    stringsAsFactors = FALSE
  )
  
  hr_data$HR[match(as.character(coef_tab$level), hr_data$variable)] <- coef_tab$HR
  hr_data$Lower_CI[match(as.character(coef_tab$level), hr_data$variable)] <- coef_tab$CI_Lower
  hr_data$Upper_CI[match(as.character(coef_tab$level), hr_data$variable)] <- coef_tab$CI_Upper
  hr_data$P_value[match(as.character(coef_tab$level), hr_data$variable)] <- coef_tab$`Pr(>|z|)`
  
  hr_data$HR[hr_data$variable == ref_km] <- 1
  hr_data$Lower_CI[hr_data$variable == ref_km] <- 1
  hr_data$Upper_CI[hr_data$variable == ref_km] <- 1
  
  legend_labels <- c(
    paste0(
      wanted_levels[1], ": ",
      ifelse(
        wanted_levels[1] == ref_km,
        "HR 1 (ref)",
        paste0(
          round(hr_data$HR[hr_data$variable == wanted_levels[1]], 2),
          " (",
          round(hr_data$Lower_CI[hr_data$variable == wanted_levels[1]], 2),
          "–",
          round(hr_data$Upper_CI[hr_data$variable == wanted_levels[1]], 2),
          ")"
        )
      )
    ),
    paste0(
      wanted_levels[2], ": ",
      ifelse(
        wanted_levels[2] == ref_km,
        "HR 1 (ref)",
        paste0(
          round(hr_data$HR[hr_data$variable == wanted_levels[2]], 2),
          " (",
          round(hr_data$Lower_CI[hr_data$variable == wanted_levels[2]], 2),
          "–",
          round(hr_data$Upper_CI[hr_data$variable == wanted_levels[2]], 2),
          ")"
        )
      )
    ),
    paste0(
      wanted_levels[3], ": ",
      ifelse(
        wanted_levels[3] == ref_km,
        "HR 1 (ref)",
        paste0(
          round(hr_data$HR[hr_data$variable == wanted_levels[3]], 2),
          " (",
          round(hr_data$Lower_CI[hr_data$variable == wanted_levels[3]], 2),
          "–",
          round(hr_data$Upper_CI[hr_data$variable == wanted_levels[3]], 2),
          ")"
        )
      )
    )
  )
  
  # km_txt <- paste0(surv_formula, " ~ ", var_km)
  # km_fit <- survival::survfit(
  #   stats::as.formula(km_txt),
  #   data = dat
  # )
  
  km_formula <- stats::as.formula(
    paste0(surv_formula, " ~ `", var_km, "`"),
    env = environment()
  )
  
  km_fit <- survival::survfit(
    formula = km_formula,
    data = dat
  )
  
  km_fit$call$formula <- km_formula
  
  strata_names <- names(km_fit$strata)
  strata_levels <- gsub(".*=", "", strata_names)
  palette_use <- unname(palette_vals[strata_levels])
  
  legend_labels_use <- legend_labels[match(strata_levels, wanted_levels)]
  
  km_plot <- survminer::ggsurvplot(
    km_fit,
    data = dat,
    conf.int = TRUE,
    censor = FALSE,
    risk.table = FALSE,
    fun = "event",
    palette = palette_use,
    xlab = "Follow-up time (years)",
    ylab = "Cumulative event",
    size = 0.2,xlim=c(0,15.5),
    legend.title = legend.title,
    legend.labs = legend_labels_use
  )
  
  km_plot$plot <- km_plot$plot +
    ggplot2::theme_classic(base_line_size = 0.2, base_size = 6) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 5, color = "black"),
      legend.title = ggplot2::element_text(size = 5, hjust = 0),
      legend.key.width = grid::unit(0.15, "cm"),
      axis.text.x = ggplot2::element_text(color = "black", size = 6),
      axis.text.y = ggplot2::element_text(color = "black", size = 6),
      legend.position = c(0.01, 0.95),
      legend.justification = c("left", "top"),
      legend.spacing.y = grid::unit(-0.1, "cm"),
      legend.key.size = grid::unit(0.15, "cm"),
      legend.title.align = 1,
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 5)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3))
  
  return(km_plot$plot)
}


prs_km <- forest_km()
life_km <- forest_km(var_forest="life_cat2",var_km="life_cat2")


interplay_forest <- function(data=tmp_interplay, var1, var2, surv="Surv(ca_reg_death_fu, ca_reg_death_inci)", title=F){
  
  data_interplay <- data %>% 
    mutate(
      var1_cat = case_when(
        var1 != "non_ca_base" ~ ifelse(!!sym(var1) < median(!!sym(var1), na.rm = TRUE), 0, 1),
        TRUE ~ !!sym(var1)
      ),
      interplay_var = paste0(var1_cat,!!sym(var2))
    ) %>% filter(!is.na(!!sym(var1))&!is.na(!!sym(var2)))
  
  
  ref <- case_when(
    str_detect(var2,"prs")~"1High",
    TRUE ~ "1Unfavorable"
  )
  
  
  f1 <- as.formula(
    glue::glue("{surv} ~ relevel(as.factor(interplay_var), ref = '{ref}') + age + sex + ethnic + edu + tdi + center")
  )
  
  trait <- case_when(var1 =="phenoage_res"~"PhenoAge",
                     var1 =="non_ca_base"~"Aging hallmarks",
                     T~paste0(gsub("_res","",var1)," Aging"))
  
  
  cox_format <- function(cox) {summary(cox)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("Variable") %>%
      mutate(
        HR = `exp(coef)`,
        CI_Lower = exp(confint(cox)[, 1]),
        CI_Upper = exp(confint(cox)[, 2])
      ) %>% filter(str_detect(Variable, "interplay_var")) %>% 
      mutate(Variable=gsub("relevel\\(.*?\\), ref = \".*?\"\\)", "", Variable)) %>% 
      bind_rows(.,data.frame(Variable=ref,HR=1,stringsAsFactors = FALSE)) %>% 
      mutate(group=ifelse(str_detect(Variable, "0"),"Accelerated aging (-)","Accelerated aging (+)"))
    # mutate(Variable= sub("_res", "", Variable) %>% str_to_title)
  }
  
  calculate_interaction_pvalue <- function(data, surv, var1, var2, covariates) {
    # Evaluate the survival formula correctly
    # surv_formula <- eval(parse(text = surv))
    
    # Create formula strings
    formula_with_interaction <- as.formula(
      glue::glue("{surv} ~ {var1} * {var2} + {covariates}")
    )
    formula_without_interaction <- as.formula(
      glue::glue("{surv} ~ {var1} + {var2} + {covariates}")
    )
    
    # Fit models
    model_with_interaction <- coxph(formula_with_interaction, data = data)
    model_without_interaction <- coxph(formula_without_interaction, data = data)
    
    # Perform ANOVA
    anova_result <- anova(model_without_interaction, model_with_interaction)
    # Extract the p-value for the interaction term
    interaction_pvalue <-anova_result[2, "Pr(>|Chi|)"] %>% sprintf("%.3f", .) 
    
    return(interaction_pvalue)
  }
  
  p_interaction <- calculate_interaction_pvalue(
    data = data_interplay %>% filter(!is.na(var1) & !is.na(var2)),
    surv = surv,
    var1 = "var1_cat",
    var2 = var2, 
    covariates = "age + sex + ethnic + edu + tdi + center"
  )
  
  # y_limits  <- case_when(
  #   str_detect(var2,"prs")~c("Low", "Inter.", "High"),
  #   TRUE ~ c("Favorable", "Inter.", "Unfavorable")
  # )
  
  y_limits <- c("Low", "Mid", "High")
  
  
  plot <- coxph(f1, data = data_interplay) %>% 
    cox_format() %>% mutate(var=gsub("^.", "", Variable)) %>% 
    mutate(var = case_when(
      str_detect(var2, "prs") ~ factor(var, levels = c("Low", "Mid", "High")),
      TRUE ~ var
    )) %>%
    mutate(
      var = recode(var,
                   "Favorable" = "Low",
                   "Intermediate" = "Mid",
                   "Unfavorable" = "High"
      )
    ) %>% 
    ggplot(.,aes(x = HR, y=(var))) +
    geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
    geom_point(aes(color = "black"),color = "black",size = 1,position=position_dodge(0.7),shape=15)+
    scale_y_discrete(limits=y_limits)+
    labs(x = "HR (95% CI)", y = ifelse(str_detect(var2, "prs"),"Genetic risk","Lifestyle risk"), title = if (title) trait else NULL) +
    theme_classic(base_size = 6,base_line_size = 0.2)+
    facet_grid(.~group,scales = "free")+
    theme(            legend.position = "right",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                      axis.text.x = element_text(color="black",size = 5),
                      axis.text.y = element_text(color="black",size = 5),
                      
                      legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                      legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                      legend.key.spacing.y = unit(0, "pt"),
                      legend.box.spacing = margin(0),
                      # legend.margin = margin(t = 1, b =1, unit = "pt"),
                      # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                      plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"),
                      strip.background = element_blank(),
                      strip.text = if (title) element_text() else element_blank(),
                      plot.title = if (title) element_text(hjust = 0.5, size = 6, margin = margin(b = -2)) else element_blank(),
                      # aspect.ratio = 2
                      # panel.spacing = unit(0.01, "cm")
                      
    )+coord_flip()+
    geom_text(
      data = data.frame(
        label = c(paste0("italic(P)[interaction]~'='~", p_interaction), ""),
        group = c("Accelerated aging (-)", "Accelerated aging (+)")
      ),
      mapping = aes(x = max(plot$data$HR, na.rm = TRUE), 
                    y = -Inf * 0.95, label = label),
      parse = TRUE,
      hjust = -0.05,
      vjust = 2,
      size = 2,
      inherit.aes = FALSE
    )
  return(plot)
}
#PAC or aging hallmarks x lifestyle ---- 
#Lung cancer
#Aging x lifestyle
f4b_1.1 <- interplay_forest(var1="non_ca_base",var2 = "life_cat",title = T)
f4b_1.2 <-interplay_forest(var1="Organismal_res",var2 = "life_cat",title = T)
f4b_1.3 <-interplay_forest(var1="phenoage_res",var2 = "life_cat",title = T)
#Aging x PRS
f4b_1.4 <- interplay_forest(var1="non_ca_base",var2 = "lung_prs_3",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")
f4b_1.5 <-interplay_forest(var1="Organismal_res",var2 = "lung_prs_3",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")
f4b_1.6 <-interplay_forest(var1="phenoage_res",var2 = "lung_prs_3",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")

#Plotting----

plot_grid(life_km,prs_km,ncol=1,labels = c("c","e"),label_size = 8) %>% plot_grid(.,
                                                                                  plot_grid(f4b_1.1,f4b_1.2,f4b_1.3,  
                                                                                            f4b_1.4,f4b_1.5,f4b_1.6,  
                                                                                            nrow=2,ncol=3,align = "hv",labels = c("d","","","f","",""),label_size = 8),
                                                                                  nrow=1,rel_widths = c(1,3)
)


#Aging x PRS X lifestyle
f4b_1.7 <- interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==0),var1="non_ca_base",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)",title = T)
f4b_1.8 <-interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==0),var1="Organismal_res",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)",title = T)
f4b_1.9 <-interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==0),var1="phenoage_res",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)",title = T)
f4b_1.10 <- interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==1),var1="non_ca_base",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")
f4b_1.11 <-interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==1),var1="Organismal_res",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")
f4b_1.12 <-interplay_forest(data=tmp_interplay %>% filter(lung_prs_cat==1),var1="phenoage_res",var2 = "life_cat",surv = "Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci`)")


plot_grid(prs_km, f4b_1.7,f4b_1.8,f4b_1.9,  
          life_km, f4b_1.10,f4b_1.11,f4b_1.12,  
          nrow=2,ncol=4,align = "hv",labels = c("c","d","","","","","",""),label_size = 8
)


#Any cancer
#Aging x lifestyle
f4b_2.1 <- interplay_forest(var1="non_ca_base",var2 = "life_cat",title = T)
f4b_2.2 <-interplay_forest(var1="Organismal_res",var2 = "life_cat",title = T)
f4b_2.3 <-interplay_forest(var1="phenoage_res",var2 = "life_cat",title = T)


#Plotting output-----

cairo_pdf(fig4_path("FIG4_OUTPUT_PDF"), 
          width = 7.19, height = 7)

(plot_grid(f3a,f3b,nrow = 1,labels = c("a","b"),label_size = 8)+
  theme(plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"))) %>% 
  plot_grid(.,plot_grid(life_km,prs_km,ncol=1,labels = c("c","e"),label_size = 8) %>% plot_grid(.,
                                                                                                plot_grid(f4b_1.1,f4b_1.2,f4b_1.3,  
                                                                                                          f4b_1.4,f4b_1.5,f4b_1.6,  
                                                                                                          nrow=2,ncol=3,align = "hv",labels = c("d","","","f","",""),label_size = 8),
                                                                                                nrow=1,rel_widths = c(1,3)
  )+theme(plot.margin = margin(0, 0.2, 0, 0.2, "cm")),
            ncol=1,rel_heights=c(4.2,3))

dev.off()    

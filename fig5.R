#Fig.6, CHIP, aging, and hemato cancer

#Get previous data
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
# Set the FIG5_* environment variables to local paths before running.
# FIG5_FIG2_WORKSPACE should match FIG2_OUTPUT_WORKSPACE from fig2.R.
# FIG5_PRS_SCORE_DIR contains the phenotype subdirectories with .sscore files.
# Set FIG5_OUTPUT_PDF to a file in an existing directory outside the repository.
fig5_path <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path)) {
    stop(sprintf("Set %s to a local path before running fig5.R.", variable),
         call. = FALSE)
  }
  path.expand(path)
}

load(fig5_path("FIG5_FIG2_WORKSPACE"))

#Get PRS for diseases and trais----- 
# a, standard set using external SNPs
{
  ukb_prs <- read.csv(file = fig5_path("FIG5_PRS_CSV"),check.names = FALSE)
  ukb_pc <- read.csv(file = fig5_path("FIG5_PC_CSV"),check.names = FALSE)
  
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
  prs_list <- bind_rows(fread(file = fig5_path("FIG5_PRS_LIST"), header = FALSE, col.names = c("pgs_id", "pheno")), 
                        fread(file = fig5_path("FIG5_PRS_LIST2"), header = FALSE, col.names = c("pgs_id", "pheno"))
  )
  
  prs_add <- lapply(seq_len(nrow(prs_list)), function(i) {
    pheno <- prs_list$pheno[i]
    f <- file.path(fig5_path("FIG5_PRS_SCORE_DIR"), pheno, paste0("prs_", pheno, ".sscore"))
    fread(file = f) %>%
      select(IID, SCORE1_AVG) %>%
      dplyr::rename(eid = IID) %>%
      dplyr::rename(!!paste0("PRS_", pheno) := SCORE1_AVG)
  }) %>%
    Reduce(function(x, y) full_join(x, y, by = "eid"), .)
  
  #Combined all PRS 
  ukb_prs_all <- ukb_prs %>% left_join(prs_add) %>% map_eid %>% filter(!is.na(eid))
  
}


#Get CHIP data
ukb_chip <- read.csv(file = fig5_path("FIG5_CHIP_CSV"),check.names = FALSE)

#Define number of CHIP variants
ukb_chip2 <- ukb_chip %>%
  mutate(CHIP_n = as.numeric(`30105-0.0`)) %>%
  map_eid() %>%
  filter(!is.na(eid), !is.na(CHIP_n)) %>%
  mutate(
    chip = ifelse(CHIP_n >= 1, 1, 0),
    CHIP = factor(chip, levels = c(0, 1), labels = c("No CHIP", "CHIP"))
  )
table(ukb_chip2$chip)
nrow(ukb_chip2)

mut_cols <- grep("^30106", colnames(ukb_chip2), value = TRUE)

ukb_chip3 <- ukb_chip2 %>%
  unite("mut_all", all_of(mut_cols), sep = ";", remove = FALSE, na.rm = TRUE) %>% 
  mutate(gene = sub(":.*", "", mut_all))


tmp_chip <- dt_aging_clock_2_cancer %>% left_join(ukb_chip3)

#Check number 
ca_diseases_hemato <- c(
  "NHL", "Plasma Cell Cancer", "Leukaemia",
  "Polycythaemia Vera", "Myelodysplasia", "MGUS"
)

hemato_date_cols <- paste0(ca_diseases_hemato, "_reg_death_date")

tmp_chip %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c("recruit_date", hemato_date_cols)),
      ~ as.Date(dplyr::na_if(as.character(.x), ""))
    )
  ) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(hemato_date_cols),
      ~ is.na(.x) | .x > recruit_date
    )
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_CHIP = sum(CHIP == "CHIP", na.rm = TRUE),
    CHIP_prevalence = round(100 * n_CHIP / n, 1)
  ) #15616



tmp_chip %>% filter(ca_base==0&CHIP=="CHIP") %>%  filter(!is.na(Organismal_res)) %>% nrow() #1609

nrow(tmp_chip %>% filter(ca_base==0&!is.na(CHIP))) #425180


#Any CHIP
cox_chip <- cox(
  data=tmp_chip %>% filter(ca_base==0) #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c(
    # organ_pac,"phenoage_res","kdm_res","met_res"
    "chip"
  ),
  outcomes=c(
    "ca_reg_death",
    paste0(ca_diseases,"_reg_death"),
    "death"
  ), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)

#DNMT3A
cox_chip_DNMT3A <- cox(
  data=tmp_chip %>% filter(ca_base==0) %>%  filter((gene == "DNMT3A" & CHIP == "CHIP") | CHIP == "No CHIP") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c(
    # organ_pac,"phenoage_res","kdm_res","met_res"
    "chip"
  ),
  outcomes=c(
    cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)
  ), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

#TET2
cox_chip_TET2 <- cox(
  data=tmp_chip %>% filter(ca_base==0) %>%  filter((gene == "TET2" & CHIP == "CHIP") | CHIP == "No CHIP") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c(
    # organ_pac,"phenoage_res","kdm_res","met_res"
    "chip"
  ),
  outcomes=c(
    cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)
  ), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

#ASXL1
cox_chip_ASXL1 <- cox(
  data=tmp_chip %>% filter(ca_base==0) %>%  filter((gene == "ASXL1" & CHIP == "CHIP") | CHIP == "No CHIP") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c(
    # organ_pac,"phenoage_res","kdm_res","met_res"
    "chip"
  ),
  outcomes=c(
    cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)
  ), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))




cox_chip_1 <- cox(
  data=tmp_chip %>% filter(ca_base==0&CHIP=="CHIP") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c("phenoage_res","kdm_res","met_res","non_ca_base"),
  outcomes=c(cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

cox_chip_1_DNMT3A <- cox(
  data=tmp_chip %>% filter(ca_base==0&gene=="DNMT3A") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c("phenoage_res","kdm_res","met_res","non_ca_base"),
  outcomes=c(cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

cox_chip_1_TET2 <- cox(
  data=tmp_chip %>% filter(ca_base==0&gene=="TET2") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c("phenoage_res","kdm_res","met_res","non_ca_base"),
  outcomes=c(cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

cox_chip_1_ASXL1 <- cox(
  data=tmp_chip %>% filter(ca_base==0&gene=="ASXL1") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c("phenoage_res","kdm_res","met_res","non_ca_base"),
  outcomes=c(cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))

cox_chip_0 <- cox(
  data=tmp_chip %>% filter(ca_base==0&CHIP=="No CHIP") #%>% mutate(across(all_of(organ_pac), ~ ifelse(t1 == 0, NA, .)))
  ,
  exposures=c("phenoage_res","kdm_res","met_res","non_ca_base"),
  outcomes=c(cox_chip$result_list %>% filter(FDR_P_Value<0.05) %>% pull(Outcome)), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))



# CHIP, aging, and hemato cancer------
outcome_levels <- c(
  "Any cancer","Death","Leukaemia","Myelodysplasia",
  "Polycythaemia Vera","NHL","Lung Cancer","Oesophageal Cancer"
)

sig_outcomes <- cox_chip$result_list %>%
  filter(FDR_P_Value < 0.05) %>%
  pull(Outcome) %>%
  unique()

plot_dat <- bind_rows(
  cox_chip$result_list %>%
    filter(Outcome %in% sig_outcomes) %>%
    mutate(group = "Any CHIP"),
  
  cox_chip_DNMT3A$result_list %>%
    filter(Outcome %in% sig_outcomes) %>%
    mutate(group = "DNMT3A"),
  
  cox_chip_TET2$result_list %>%
    filter(Outcome %in% sig_outcomes) %>%
    mutate(group = "TET2"),
  
  cox_chip_ASXL1$result_list %>%
    filter(Outcome %in% sig_outcomes) %>%
    mutate(group = "ASXL1")
) %>%
  filter(Exposure == "chip") %>%
  mutate(
    Outcome = case_when(
      Outcome == "ca_reg_death" ~ "Any cancer",
      Outcome == "death" ~ "Death",
      grepl("leuk|leuc|leukaemia|leukemia", Outcome, ignore.case = TRUE) ~ "Leukaemia",
      grepl("myelodys", Outcome, ignore.case = TRUE) ~ "Myelodysplasia",
      grepl("polycythaemia|polycythemia|vera", Outcome, ignore.case = TRUE) ~ "Polycythaemia Vera",
      grepl("nhl|non_hodgkin|non-hodgkin", Outcome, ignore.case = TRUE) ~ "NHL",
      grepl("lung", Outcome, ignore.case = TRUE) ~ "Lung Cancer",
      grepl("oesoph|esoph", Outcome, ignore.case = TRUE) ~ "Oesophageal Cancer",
      TRUE ~ Outcome
    ),
    group = factor(group, levels = c("Any CHIP", "DNMT3A", "TET2", "ASXL1"))
  ) %>%
  filter(Outcome %in% outcome_levels)


plot_dat2 <- plot_dat %>%
  mutate(
    Outcome = factor(Outcome, levels = rev(outcome_levels)),
    group = factor(group, levels = c("Any CHIP", "DNMT3A", "TET2", "ASXL1"))
  )

p_left <- ggplot(plot_dat2, aes(x = HR, y = Outcome, color = group, shape = group)) +
  geom_vline(xintercept = 1, color = "grey75", linewidth = 0.1,linetype = "dashed") +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper),
    height = 0,
    linewidth = 0.2,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 1.2, position = position_dodge(width = 0.55)) +
  scale_x_log10(
    breaks = c(0.3, 1, 3, 10, 30),
    labels = c("0.3", "1.0", "3.0", "10.0", "30.0")
  ) +
  scale_y_discrete(
    limits = rev(outcome_levels),
    drop = FALSE,
    expand = expansion(add = c(0.6, 0.6))
  ) +
  scale_color_manual(values = c(
    "Any CHIP" = "black",
    "DNMT3A"   = "#77AFC8",
    "TET2"     = "#D79A7E",
    "ASXL1"    = "#9E9E9E"
  )) +
  scale_shape_manual(values = c(
    "Any CHIP" = 15,
    "DNMT3A"   = 17,
    "TET2"     = 17,
    "ASXL1"    = 17
  )) +
  labs(x = "HR (95% CI)", y = NULL, color = NULL, shape = NULL) +
  theme_classic(base_size = 6, base_line_size = 0.2) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.width = unit(0.1, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    axis.text.x = element_text(color = "black", size = 6),
    axis.text.y = element_text(color = "black", size = 6, margin = margin(r = 2)),
    plot.margin = margin(6, 2, 6, 6)
  )

make_heat_dat <- function(res){
  
  res %>%  mutate(
    Outcome = gsub("_reg_death", "", Outcome),
    Outcome = recode(
      Outcome,
      ca = "Any cancer",
      death = "Death"
    ),
    Exposure = recode(
      Exposure,
      met_res = "MetAge",
      non_ca_base = "Any hallmarks",
      kdm_res = "KDMAge",
      phenoage_res = "PhenoAge"
    )
  ) %>%
    filter(
      Outcome %in% outcome_levels,
      Exposure %in% c("Any hallmarks", "MetAge", "KDMAge", "PhenoAge")
    ) %>%
    mutate(
      Outcome = factor(Outcome, levels = rev(outcome_levels)),
      Exposure = factor(Exposure, levels = c("Any hallmarks", "MetAge", "KDMAge", "PhenoAge"))
    )
}

make_heatmap <- function(dat,title,legend=F){
  
  ggplot(make_heat_dat(dat), aes(x = Exposure, y = Outcome, fill = Beta)) +
    geom_tile(color = "black", linewidth = 0.1, width = 1, height = 1) +
    geom_text(
      aes(label = ifelse(FDR_P_Value < 0.05, sprintf("%.1f", HR), "")),
      size = 2
    ) +
    scale_fill_gradient2(
      low = "#5AA38E",
      mid = "grey90",
      high = "#CF563D",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      name = NULL
    ) +
    scale_y_discrete(
      limits = rev(outcome_levels),
      drop = FALSE,
      expand = expansion(add = c(0.6, 0.6))
    ) +
    labs(x = NULL, y = NULL,title=title) +
    theme_classic(base_size = 6, base_line_size = 0.25) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 6),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 6, hjust = 0.5, margin = margin(b = 2)),
      plot.margin = margin(8, 6, 6, 0),
      legend.position = ifelse(legend, "right", "none"),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.1, "cm"),
      legend.text = element_text(size = 5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
}

p_right1 <- make_heatmap(dat=(cox_chip_1$result_list),title="Any CHIP")
p_right2 <- make_heatmap(dat=(cox_chip_1_DNMT3A$result_list),title="DNMT3A")
p_right3 <- make_heatmap(dat=(cox_chip_1_TET2$result_list),title="TET2")
p_right4 <- make_heatmap(dat=(cox_chip_1_ASXL1$result_list),title="ASXL1")
p_right5 <- make_heatmap(dat=(cox_chip_0$result_list),title="No CHIP",legend = T)


p_gene <- {

gene_df <- data.frame(table(ukb_chip3$gene), stringsAsFactors = FALSE)
colnames(gene_df) <- c("gene", "n")

gene_df <- gene_df %>%
  mutate(
    gene = as.character(gene),
    n = as.numeric(as.character(n))
  ) %>%
  filter(!gene %in% c(";;;;;;", "No CHIP", "")) %>%
  arrange(desc(n))

gene_top <- head(gene_df, 10)
gene_top$gene[10] <- "Other"
gene_top$n[10] <- sum(gene_df$n[10:nrow(gene_df)])

gene_top <- gene_top[1:10, ] %>%
  mutate(
    pct = 100 * n / sum(n),
    ymax = cumsum(pct),
    ymin = c(0, head(ymax, -1)),
    mid = (ymax + ymin) / 2,
    label = ""
  )

gene_top$label[1:3] <- paste0(round(gene_top$pct[1:3], 1), "%")
gene_top$gene <- factor(gene_top$gene, levels = gene_top$gene)

chip_cols <- c(
  "#77AFC8",  # 1
           "#D79A7E",  # 2
           "#9E9E9E",  # 3
           "#CF563D",  # 4
           "#E6C0C9",  # 5
           "#5AA38E",  # 6
           "#5B6DAF",  # 7
           "#E3B58C",  # 8
           "#BFDCE8",  # 9
           "#D9D9D9"   # 10
)

ggplot(gene_top, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = gene)) +
  geom_rect(color = "black", linewidth = 0.1) +
  coord_polar(theta = "y") +
  xlim(c(2.2, 4)) +
  geom_text(
    data = gene_top[1:3, ],
    aes(x = 3.5, y = mid, label = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_fill_manual(values = chip_cols) +
  theme_void(base_size = 6, base_line_size = 0.2) +
  labs(title = NULL) +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.25, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    # plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "right",
    legend.box.margin = margin(0, -6, 0, -6),
    legend.margin = margin(0, -4, 0, -4),
    plot.margin = margin(0, 0.3, 0, 0, "cm")
  )
}

p_final <- p_gene+p_left + p_right1+ p_right2 + p_right3+ p_right4+ p_right5+plot_layout(widths = c(2,2.5, 1,1,1,1,1))+
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0, 0.9)
  )

p_final

#Combined effect of aging and CHIP (KM)

tmp_chip_km <- cox_pac_ca$data %>% 
  left_join(ukbb_add_outcome2 %>% select(eid,smoke_cat,alcohol_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat)) %>% 
  mutate(pa_cat=ifelse(pa==0,1,0),obe_cat=ifelse(bmi>=30,1,0)) %>% 
  mutate(lifestyle=rowSums(dplyr::select(., smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,obe_cat,fruit_cat,fish_cat,red_cat,process_cat)),
         life_cat=ifelse(lifestyle<=2,"Favorable",ifelse(lifestyle<=4,"Intermediate","Unfavorable"))
  ) %>% 
  left_join(ukb_prs_all %>% select(-(pc11:pc40)), by = "eid") %>% 
  mutate(across(starts_with("PRS_"), ~ as.numeric(scale(.)))) %>% 
  left_join(ukb_chip3)

tmp_chip_km$phenoage_cat <- as.integer(tmp_chip_km$phenoage_res > median(tmp_chip_km$phenoage_res, na.rm = TRUE))

#Only P for CHIP+ Aging
make_chip_km <- function(data=tmp_chip_km,
                         disease,
                         aging_var = "non_ca_base",
                         title = NULL,legend=F,
                         covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")) {
  
  if (disease == "death") {
    time_var  <- paste0(disease, "_fu")
    event_var <- paste0(disease, "_inci")
  } else {
    time_var  <- paste0(disease, "_reg_death_fu")
    event_var <- paste0(disease, "_reg_death_inci")
  }
  
  dat <- data %>%
    mutate(
      grp = factor(
        case_when(
          chip == 0 & .data[[aging_var]] == 0 ~ "No CHIP, accelerated aging-",
          chip == 0 & .data[[aging_var]] == 1 ~ "No CHIP, accelerated aging+",
          chip == 1 & .data[[aging_var]] == 0 ~ "CHIP, accelerated aging-",
          chip == 1 & .data[[aging_var]] == 1 ~ "CHIP, accelerated aging+"
        ),
        levels = c(
          "No CHIP, accelerated aging-",
          "No CHIP, accelerated aging+",
          "CHIP, accelerated aging-",
          "CHIP, accelerated aging+"
        )
      )
    )
  
  s <- survfit(
    as.formula(paste0("Surv(`", time_var, "`, `", event_var, "` == 1) ~ grp")),
    data = dat
  )
  
  fit <- coxph(
    as.formula(
      paste0(
        "Surv(`", time_var, "`, `", event_var, "` == 1) ~ grp + ",
        paste(covar, collapse = " + ")
      )
    ),
    data = dat
  )
  
  ci <- summary(fit)$conf.int
  hr_lab <- sprintf(
    "HR, CHIP + aging:\n%.2f (%.2f–%.2f)",
    ci["grpCHIP, accelerated aging+", "exp(coef)"],
    ci["grpCHIP, accelerated aging+", "lower .95"],
    ci["grpCHIP, accelerated aging+", "upper .95"]
  )
  
  hr_tab <- ci %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(stringr::str_detect(term, "^grp")) %>%
    dplyr::mutate(
      grp = sub("^grp", "", term),
      HR = `exp(coef)`,
      Lower_CI = `lower .95`,
      Upper_CI = `upper .95`,
      HR_CI = sprintf("%.2f (%.2f–%.2f)", HR, Lower_CI, Upper_CI)
    ) %>%
    dplyr::select(grp, HR, Lower_CI, Upper_CI, HR_CI)
  
  hr_tab <- tibble::tibble(
    grp = levels(dat$grp),
    HR = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    Lower_CI = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    Upper_CI = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    HR_CI = c("1.00 (ref)", rep(NA_character_, length(levels(dat$grp)) - 1))
  ) %>%
    dplyr::rows_update(hr_tab, by = "grp")
  
  print(hr_tab)
  
  df <- data.frame(
    time = summary(s)$time,
    surv = summary(s)$surv,
    grp = sub(".*=", "", summary(s)$strata)
  )
  
  ggplot(df, aes(x = time, y = 1 - surv, color = grp)) +
    geom_step(linewidth = 0.3) +
    scale_color_manual(values = c(
      "No CHIP, accelerated aging-" = "#8FB6D8",
      "No CHIP, accelerated aging+" = "#E6B8C3",
      "CHIP, accelerated aging-" = "#D79A7E",
      "CHIP, accelerated aging+" = "#B22222"
    )) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Time (years)",
      y = "Cumulative event",
      color = NULL,
      title = title
    ) +
    annotate(
      "text",
      x = min(df$time, na.rm = TRUE) + 0.05 * diff(range(df$time, na.rm = TRUE)),
      y = max(1 - df$surv, na.rm = TRUE) - 0.1 * diff(range(1 - df$surv, na.rm = TRUE)),
      label = hr_lab,
      hjust = 0,
      size = 2
    ) +
    theme_classic(base_size = 6, base_line_size = 0.2) +
    theme(
      legend.position = ifelse(legend, "top", "none"),
      legend.direction = "horizontal",
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.1, "cm"),
      legend.spacing.y = unit(-0.3, "cm"),
      legend.text = element_text(size = 5),
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5, size = 6),
      plot.margin = margin(0.2, 0.1, 0.1, 0.1, "cm")
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

#Add p for interaction
make_chip_km <- function(data = tmp_chip_km,
                         disease,
                         aging_var = "non_ca_base",
                         title = NULL,
                         legend = FALSE,
                         covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi")) {
  
  if (disease == "death") {
    time_var  <- paste0(disease, "_fu")
    event_var <- paste0(disease, "_inci")
  } else {
    time_var  <- paste0(disease, "_reg_death_fu")
    event_var <- paste0(disease, "_reg_death_inci")
  }
  
  bt <- function(x) {
    ifelse(
      grepl("^[A-Za-z.][A-Za-z0-9._]*$", x),
      x,
      paste0("`", x, "`")
    )
  }
  
  covar_terms <- bt(covar)
  
  dat <- data %>%
    dplyr::mutate(
      time_i = as.numeric(.data[[time_var]]),
      event_i = as.integer(.data[[event_var]] == 1),
      
      chip_i = dplyr::case_when(
        as.character(chip) %in% c("0", "FALSE", "false") ~ 0L,
        as.character(chip) %in% c("1", "TRUE", "true") ~ 1L,
        TRUE ~ NA_integer_
      ),
      
      aging_i = dplyr::case_when(
        as.character(.data[[aging_var]]) %in% c("0", "FALSE", "false") ~ 0L,
        as.character(.data[[aging_var]]) %in% c("1", "TRUE", "true") ~ 1L,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::filter(
      !is.na(time_i),
      !is.na(event_i),
      !is.na(chip_i),
      !is.na(aging_i),
      dplyr::if_all(dplyr::all_of(covar), ~ !is.na(.x))
    ) %>%
    dplyr::mutate(
      grp = factor(
        dplyr::case_when(
          chip_i == 0 & aging_i == 0 ~ "No CHIP, accelerated aging-",
          chip_i == 0 & aging_i == 1 ~ "No CHIP, accelerated aging+",
          chip_i == 1 & aging_i == 0 ~ "CHIP, accelerated aging-",
          chip_i == 1 & aging_i == 1 ~ "CHIP, accelerated aging+"
        ),
        levels = c(
          "No CHIP, accelerated aging-",
          "No CHIP, accelerated aging+",
          "CHIP, accelerated aging-",
          "CHIP, accelerated aging+"
        )
      )
    ) %>%
    dplyr::filter(!is.na(grp)) %>%
    droplevels()
  
  surv_form <- "Surv(time_i, event_i)"
  
  rhs_grp <- paste(
    c("grp", covar_terms),
    collapse = " + "
  )
  
  rhs_int <- paste(
    c("chip_i", "aging_i", covar_terms, "chip_i:aging_i"),
    collapse = " + "
  )
  
  ## KM curve
  s <- survfit(
    as.formula(paste0(surv_form, " ~ grp")),
    data = dat
  )
  
  ## 4-group HR model
  fit_grp <- coxph(
    as.formula(paste0(surv_form, " ~ ", rhs_grp)),
    data = dat,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE,
    na.action = na.fail
  )
  
  ci <- summary(fit_grp)$conf.int
  
  ## Single full Cox model for ANOVA interaction P
  ## interaction term is deliberately placed LAST
  fit_int <- coxph(
    as.formula(paste0(surv_form, " ~ ", rhs_int)),
    data = dat,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE,
    na.action = na.fail
  )
  
  ## ANOVA on ONE model, not anova(fit_main, fit_int)
  anova_int <- anova(fit_int, test = "Chisq")
  
  int_row <- grep(
    "^chip_i:aging_i$|^aging_i:chip_i$",
    rownames(anova_int),
    value = TRUE
  )
  
  p_col <- grep("^Pr", colnames(anova_int), value = TRUE)[1]
  p_interaction <- anova_int[int_row, p_col]
  
  p_int_lab <- dplyr::case_when(
    length(p_interaction) == 0 ~ "P interaction = NA",
    is.na(p_interaction) ~ "P interaction = NA",
    p_interaction < 0.001 ~ "P interaction < 0.001",
    TRUE ~ sprintf("P interaction = %.3f", p_interaction)
  )
  
  hr_lab <- sprintf(
    "HR, CHIP + aging:\n%.2f (%.2f–%.2f)\n%s",
    ci["grpCHIP, accelerated aging+", "exp(coef)"],
    ci["grpCHIP, accelerated aging+", "lower .95"],
    ci["grpCHIP, accelerated aging+", "upper .95"],
    p_int_lab
  )
  
  hr_tab <- ci %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(stringr::str_detect(term, "^grp")) %>%
    dplyr::mutate(
      grp = sub("^grp", "", term),
      HR = `exp(coef)`,
      Lower_CI = `lower .95`,
      Upper_CI = `upper .95`,
      HR_CI = sprintf("%.2f (%.2f–%.2f)", HR, Lower_CI, Upper_CI)
    ) %>%
    dplyr::select(grp, HR, Lower_CI, Upper_CI, HR_CI)
  
  hr_tab <- tibble::tibble(
    grp = levels(dat$grp),
    HR = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    Lower_CI = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    Upper_CI = c(1, rep(NA_real_, length(levels(dat$grp)) - 1)),
    HR_CI = c("1.00 (ref)", rep(NA_character_, length(levels(dat$grp)) - 1))
  ) %>%
    dplyr::rows_update(hr_tab, by = "grp")
  
  print(hr_tab)
  print(anova_int)
  print(paste0("N used = ", nrow(dat), "; events = ", sum(dat$event_i == 1)))
  
  ss <- summary(s)
  
  df <- data.frame(
    time = ss$time,
    surv = ss$surv,
    grp = sub(".*=", "", ss$strata)
  )
  
  ggplot(df, aes(x = time, y = 1 - surv, color = grp)) +
    geom_step(linewidth = 0.3) +
    scale_color_manual(values = c(
      "No CHIP, accelerated aging-" = "#8FB6D8",
      "No CHIP, accelerated aging+" = "#E6B8C3",
      "CHIP, accelerated aging-" = "#D79A7E",
      "CHIP, accelerated aging+" = "#B22222"
    )) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Time (years)",
      y = "Cumulative event",
      color = NULL,
      title = title
    ) +
    annotate(
      "text",
      x = min(df$time, na.rm = TRUE) + 0.05 * diff(range(df$time, na.rm = TRUE)),
      y = max(1 - df$surv, na.rm = TRUE) - 0.1 * diff(range(1 - df$surv, na.rm = TRUE)),
      label = hr_lab,
      hjust = 0,
      size = 2
    ) +
    theme_classic(base_size = 6, base_line_size = 0.2) +
    theme(
      legend.position = ifelse(legend, "top", "none"),
      legend.direction = "horizontal",
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.1, "cm"),
      legend.spacing.y = unit(-0.3, "cm"),
      legend.text = element_text(size = 5),
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5, size = 6),
      plot.margin = margin(0.2, 0.1, 0.1, 0.1, "cm")
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

p_km_hallmark_ca <- make_chip_km(disease = "ca",title = "Any cancer")
p_km_hallmark_death <- make_chip_km(disease = "death",title = "Any-cause mortality")
p_km_hallmark_1 <- make_chip_km(disease = "Leukaemia",title = "Leukaemia")
p_km_hallmark_2 <- make_chip_km(disease = "Myelodysplasia",title = "Myelodysplasia")
p_km_hallmark_3 <- make_chip_km(disease = "Polycythaemia Vera",title = "Polycythaemia Vera")

p_km_pheno_ca <- make_chip_km(disease = "ca",aging_var="phenoage_cat")
p_km_pheno_death <- make_chip_km(disease = "death",aging_var="phenoage_cat")
p_km_pheno_1 <- make_chip_km(disease = "Leukaemia",aging_var="phenoage_cat")
p_km_pheno_2 <- make_chip_km(disease = "Myelodysplasia",aging_var="phenoage_cat")
p_km_pheno_3 <- make_chip_km(disease = "Polycythaemia Vera",aging_var="phenoage_cat")


cairo_pdf(fig5_path("FIG5_OUTPUT_PDF"), 
          width = 7.3, height = 5.3)
#Main plot
p_final %>% plot_grid(.,
plot_grid(p_km_hallmark_ca,p_km_hallmark_death,p_km_hallmark_1,p_km_hallmark_2,p_km_hallmark_3,
          p_km_pheno_ca,p_km_pheno_death,p_km_pheno_1,p_km_pheno_2,p_km_pheno_3,align="hv",
          ncol = 5,nrow = 2, 
          labels = c("d", "e", "f", "g", "h", "", "", "", "", ""),  label_size = 8, label_fontface = "bold"
),ncol = 1
)
dev.off()    

#Legend
ggplot(
  data.frame(
    x = 1:4,
    y = 1,
    grp = factor(
      c(
        "No CHIP, accelerated aging -",
        "No CHIP, accelerated aging +",
        "CHIP, accelerated aging -",
        "CHIP, accelerated aging +"
      ),
      levels = c(
        "No CHIP, accelerated aging -",
        "No CHIP, accelerated aging +",
        "CHIP, accelerated aging -",
        "CHIP, accelerated aging +"
      )
    )
  ),
  aes(x, y, color = grp)
) +
  geom_segment(aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y), linewidth = 0.3) +
  scale_color_manual(values = c(
    "No CHIP, accelerated aging -" = "#8FB6D8",
    "No CHIP, accelerated aging +" = "#E6B8C3",
    "CHIP, accelerated aging -" = "#D79A7E",
    "CHIP, accelerated aging +" = "#B22222"
  )) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.direction = "horizontal",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.text = element_text(size = 5,margin = margin(t = -3, b = -3, l = 1, r = 1)),
    legend.box.margin = margin(3, -6, 0, -6),
    legend.margin = margin(2, -4, 0, -4),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    legend.box = "horizontal",
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


#CHIP and PAC
gap_vars <- c(organ_pac,"phenoage_res","kdm_res","met_res")

test_res_chip <- lapply(gap_vars, function(v) {
  dat <- tmp_chip %>%
    filter(!is.na(CHIP), !is.na(age), !is.na(sex), !is.na(.data[[v]])) %>% 
    filter(ca_base_any==0&non_ca_base==0)
  
  if (nrow(dat) == 0) return(NULL)
  if (length(unique(dat$CHIP)) < 2) return(NULL)
  
  fit <- lm(as.formula(paste0("`", v, "` ~ CHIP + age + sex + ca_base_any + non_ca_base")), data = dat)
  sm <- summary(fit)$coefficients
  
  data.frame(
    var = v,
    beta = sm["CHIPCHIP", "Estimate"],
    se   = sm["CHIPCHIP", "Std. Error"],
    p    = sm["CHIPCHIP", "Pr(>|t|)"],
    n    = nrow(dat)
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, method = "fdr"))

plot_list_chip <- lapply(gap_vars, function(v) {
  ggplot(
    tmp_chip %>% filter(!is.na(.data[[v]]), !is.na(CHIP)),
    aes(x = CHIP, y = .data[[v]], fill = CHIP)
  ) +
    geom_boxplot(width = 0.6, outlier.size = 0.3, linewidth = 0.2) +
    # geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
    theme_classic(base_size = 7) +
    labs(
      x = NULL,
      y = "Aging gap",
      title = gsub("_gap|_res", "", v)
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    )
})

plot_grid(plotlist = plot_list_chip, ncol = 4)


#IL1A,DNMT3A, hemato aging to cancer------
tmp_chip_pac <- tmp_chip_km %>% filter(ca_base==0&!is.na(Brain_res)) %>% 
  left_join(olink_full, by = "eid") #48750

tmp_chip_pac$Immune_cat <- as.integer(tmp_chip_pac$Immune_res > median(tmp_chip_pac$Immune_res, na.rm = TRUE))
tmp_chip_pac$Lung_cat <- as.integer(tmp_chip_pac$Lung_res > median(tmp_chip_pac$Lung_res, na.rm = TRUE))


#Immune aging x CHIP
coxph(
  Surv(`Lung Cancer_reg_death_fu`, `Lung Cancer_reg_death_inci` == 1) ~
    IL1RN+ 
    age + sex + ethnic + edu + tdi + center + pa + smok + bmi,
  data = tmp_chip_pac %>% filter(CHIP=="CHIP")
) %>% summary()


#Define markers 
  # IL-1 mRNA program
  IL1_core <- c("IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RAP", "IL1RN")
    IL1_core <- IL1_core[IL1_core %in% colnames(tmp_chip_pac)]
  
  # granulocyte / neutrophil-skewed myelopoiesis
  granulo_core<- c( "CXCL8", "S100A8", "S100A9", "MPO", "LCN2", "MMP8", "MMP9", "OLR1")
    granulo_core <- granulo_core[granulo_core %in% colnames(tmp_chip_pac)]
  
  # monocyte-derived macrophage program
  mono_mac_core <- c("CD14", "LYZ", "LST1", "FCGR3A", "CSF1R", "CD68", "CD163", "SPP1", "MMP12", "TREM2")
    mono_mac_core <- mono_mac_core[mono_mac_core %in% colnames(tmp_chip_pac)]
  # IL-1 related program
  IL1_related_core <- c("IL18", "IL18BP", "IL18R1", "IL33", "IL1RL1", "IL36A", "IL36B", "IL36G", "IL36RN", "IL1RL2", "IL37", "IL1F10", "SIGIRR")
    IL1_related_core <- IL1_related_core[IL1_related_core %in% colnames(tmp_chip_pac)]

  protein_list <- list(
      IL1_core = IL1_core,
      mono_mac_core = mono_mac_core,
      granulo_core = granulo_core,
      IL1_related_core=IL1_related_core
    )
  
  {
    run_protein_cox_by_chip_strata <- function(
    data = tmp_chip_pac,
    disease = "Lung Cancer",
    proteins,
    strat_var = c("Immune_cat"),
    covar = c("age", "sex", "ethnic", "edu", "tdi", "center", "pa", "smok", "bmi"),
    min_n = 50,
    min_event = 5
    ) {
      
      library(dplyr)
      library(purrr)
      library(survival)
      library(broom)
      library(tibble)
      
      time_var  <- paste0(disease, "_reg_death_fu")
      event_var <- paste0(disease, "_reg_death_inci")
      
      out <- purrr::map_dfr(strat_var, function(strat_var_i) {
        
        purrr::imap_dfr(proteins, function(protein_vec, category_name) {
          
          category_clean <- category_name
          
          protein_vec <- protein_vec[protein_vec %in% colnames(data)]
          
          if (length(protein_vec) == 0) {
            return(tibble(
              category = category_clean,
              disease = disease,
              protein = NA_character_,
              strat_var = strat_var_i,
              strat_level = NA_character_,
              CHIP = NA_character_,
              n = NA_integer_,
              n_event = NA_integer_,
              HR = NA_real_,
              lower95 = NA_real_,
              upper95 = NA_real_,
              p = NA_real_,
              HR_CI = NA_character_,
              note = "no proteins found in data"
            ))
          }
          
          purrr::map_dfr(protein_vec, function(protein) {
            
            strat_levels <- sort(unique(as.character(data[[strat_var_i]])))
            strat_levels <- strat_levels[!is.na(strat_levels)]
            
            purrr::map_dfr(strat_levels, function(strat_level_i) {
              
              purrr::map_dfr(c("No CHIP", "CHIP"), function(chip_level) {
                
                dat_use <- data %>%
                  filter(
                    as.character(.data[[strat_var_i]]) == strat_level_i,
                    CHIP == chip_level,
                    !is.na(.data[[time_var]]),
                    !is.na(.data[[event_var]]),
                    !is.na(.data[[protein]])
                  ) %>%
                  mutate(
                    prot_z = as.numeric(scale(.data[[protein]]))
                  )
                
                n_event <- sum(dat_use[[event_var]] == 1, na.rm = TRUE)
                
                if (nrow(dat_use) < min_n || n_event < min_event) {
                  return(tibble(
                    category = category_clean,
                    disease = disease,
                    protein = protein,
                    strat_var = strat_var_i,
                    strat_level = strat_level_i,
                    CHIP = chip_level,
                    n = nrow(dat_use),
                    n_event = n_event,
                    HR = NA_real_,
                    lower95 = NA_real_,
                    upper95 = NA_real_,
                    p = NA_real_,
                    HR_CI = NA_character_,
                    note = "too few samples/events"
                  ))
                }
                
                fml <- as.formula(
                  paste0(
                    "Surv(`", time_var, "`, `", event_var, "` == 1) ~ ",
                    "prot_z + ",
                    paste(covar, collapse = " + ")
                  )
                )
                
                fit <- tryCatch(
                  survival::coxph(fml, data = dat_use),
                  error = function(e) NULL
                )
                
                if (is.null(fit)) {
                  return(tibble(
                    category = category_clean,
                    disease = disease,
                    protein = protein,
                    strat_var = strat_var_i,
                    strat_level = strat_level_i,
                    CHIP = chip_level,
                    n = nrow(dat_use),
                    n_event = n_event,
                    HR = NA_real_,
                    lower95 = NA_real_,
                    upper95 = NA_real_,
                    p = NA_real_,
                    HR_CI = NA_character_,
                    note = "cox failed"
                  ))
                }
                
                res <- broom::tidy(
                  fit,
                  exponentiate = TRUE,
                  conf.int = TRUE
                ) %>%
                  filter(term == "prot_z")
                
                if (nrow(res) == 0) {
                  return(tibble(
                    category = category_clean,
                    disease = disease,
                    protein = protein,
                    strat_var = strat_var_i,
                    strat_level = strat_level_i,
                    CHIP = chip_level,
                    n = nrow(dat_use),
                    n_event = n_event,
                    HR = NA_real_,
                    lower95 = NA_real_,
                    upper95 = NA_real_,
                    p = NA_real_,
                    HR_CI = NA_character_,
                    note = "protein term missing"
                  ))
                }
                
                HR_i    <- res$estimate
                lower_i <- res$conf.low
                upper_i <- res$conf.high
                p_i     <- res$p.value
                
                tibble(
                  category = category_clean,
                  disease = disease,
                  protein = protein,
                  strat_var = strat_var_i,
                  strat_level = strat_level_i,
                  CHIP = chip_level,
                  n = nrow(dat_use),
                  n_event = n_event,
                  HR = HR_i,
                  lower95 = lower_i,
                  upper95 = upper_i,
                  p = p_i,
                  HR_CI = sprintf("%.2f (%.2f–%.2f)", HR_i, lower_i, upper_i),
                  note = NA_character_
                )
              })
            })
          })
        })
      })
      
      out %>%
        group_by(strat_var, category) %>%
        mutate(
          p_fdr = p.adjust(p, method = "BH")
        ) %>%
        ungroup() %>%
        mutate(
          group = paste0(strat_var, "=", strat_level, ", ", CHIP)
        ) %>%
        arrange(strat_var, category, protein, strat_level, CHIP)
    }
    
    res_mono_mac_lung <- run_protein_cox_by_chip_strata(
      data = tmp_chip_pac,
      disease = "Lung Cancer",
      proteins = protein_list,
      strat_var = c("Immune_cat","Lung_cat","phenoage_cat","non_ca_base")
    )
    
    
    category_order <- c(
      "IL1_core",
      "mono_mac_core",
      "granulo_core",
      "IL1_related_core"
    )
    
    category_labs <- c(
      IL1_core         = "IL-1 axis",
      mono_mac_core    = "Monocyte/macrophage",
      granulo_core     = "Granulocytic myeloid",
      IL1_related_core = "Extended IL-1"
    )
    
    strat_labs <- c(
      Immune_cat   = "Immune aging",
      Lung_cat     = "Lung aging",
      phenoage_cat = "Pheno aging",
      non_ca_base  = "Aging hallmark"
    )
    
    plot_dat <- res_mono_mac_lung %>%
      filter(!is.na(HR), !is.na(lower95), !is.na(upper95)) %>%
      mutate(
        category = as.character(category),
        
        category = case_when(
          category %in% c("IL1", "IL1_core") ~ "IL1_core",
          category %in% c("mono_mac", "mono_mac_core") ~ "mono_mac_core",
          category %in% c("granulo", "granulo_core") ~ "granulo_core",
          category %in% c("IL1_related", "IL1_related_core") ~ "IL1_related_core",
          TRUE ~ category
        ),
        
        category = factor(
          category,
          levels = category_order,
          labels = category_labs[category_order]
        ),
        
        strat_var = as.character(strat_var),
        strat_var_lab = ifelse(
          strat_var %in% names(strat_labs),
          strat_labs[strat_var],
          strat_var
        ),
        strat_var_lab = factor(
          strat_var_lab,
          levels = strat_labs[unique(strat_var[strat_var %in% names(strat_labs)])]
        ),
        
        aging_lab = case_when(
          as.character(strat_level) == "0" ~ "accelerated aging -",
          as.character(strat_level) == "1" ~ "accelerated aging +",
          TRUE ~ as.character(strat_level)
        ),
        
        group = case_when(
          aging_lab == "accelerated aging -" & CHIP == "No CHIP" ~ "No CHIP, accelerated aging -",
          aging_lab == "accelerated aging -" & CHIP == "CHIP"    ~ "CHIP, accelerated aging -",
          aging_lab == "accelerated aging +"    & CHIP == "No CHIP" ~ "No CHIP, accelerated aging +",
          aging_lab == "accelerated aging +"    & CHIP == "CHIP"    ~ "CHIP, accelerated aging +",
          TRUE ~ NA_character_
        ),
        
        group = factor(
          group,
          levels = c(
            "No CHIP, accelerated aging -","No CHIP, accelerated aging +",
            "CHIP, accelerated aging -","CHIP, accelerated aging +"
          )
        ),
        
        sig_lab = case_when(
          group == "CHIP, accelerated aging +" & p < 0.001 ~ "***",
          group == "CHIP, accelerated aging +" & p < 0.01  ~ "**",
          group == "CHIP, accelerated aging +" & p < 0.05  ~ "*",
          TRUE ~ ""
        ),
        star_y = upper95 + 0.08
      ) %>%
      filter(!is.na(category), !is.na(group), !is.na(strat_var_lab))
    
    plot_dat <- plot_dat %>%
      group_by(category) %>%
      mutate(
        protein = factor(protein, levels = unique(protein))
      ) %>%
      ungroup()
    
    p <- ggplot(
      plot_dat,
      aes(
        x = protein,
        y = HR,
        ymin = lower95,
        ymax = upper95,
        color = group
      )
    ) +
      geom_hline(
        yintercept = 1,
        linetype = 2,
        linewidth = 0.2,
        color = "black"
      ) +
      geom_errorbar(
        position = position_dodge(width = 0.7),
        width = 0,
        linewidth = 0.2
      ) +
      geom_point(
        position = position_dodge(width = 0.7), shape=16,
        size = 0.7
      ) +
      geom_text(
        aes(
          y = star_y,
          label = sig_lab,
          group = group
        ),color="black",
        position = position_dodge(width = 0.7),
        size = 2.3,
        vjust = 0,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      facet_grid(
        strat_var_lab ~ category,
        scales = "free_x",
        space = "free_x"
      ) +
      scale_y_continuous(
        breaks = c(1, 2, 3, 4, 5),limits = c(0,4.5)
      ) +
      scale_color_manual(
        values = c(
          "No CHIP, accelerated aging -" = "#8FB6D8",
          "No CHIP, accelerated aging +" = "#E6B8C3",
          "CHIP, accelerated aging -" = "#D79A7E",
          "CHIP, accelerated aging +" = "#B22222"
        )
      ) +
      labs(
        x = NULL,
        y = "HR (95% CI)",
        color = NULL
      ) +
      theme_classic(base_size = 6, base_line_size = 0.2) +
      theme(
        legend.position = "top",
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(color = "black",margin = margin(t = 0, r = 0, b = 3, l = 0, unit = "mm")),
        strip.text.y = element_text(color = "black",angle = 270),
        axis.text.x = element_text(color = "black",angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),legend.key.width = unit(0.3, "cm"),
        panel.spacing.x = unit(0.3, "lines"),
        panel.spacing.y = unit(0.3, "lines"),
        plot.margin = margin(0, 0.2, 0.2, 0.2, "cm"),legend.box.margin = margin(0, 0, -6, 0)
      )
    
    p
    
  }
  
  
  #xfig.10-------------------------------------------
  #Extended from lung cancer (Fig. 6) to other cancers.
  plot_protein_cox_by_chip <- function(
    disease,
    data=tmp_chip_pac,
    proteins=protein_list,
    strat_var=c("Immune_cat","Lung_cat","phenoage_cat","non_ca_base")
  ){
    
    res <- run_protein_cox_by_chip_strata(
      data=data,
      disease=disease,
      proteins=proteins,
      strat_var=strat_var
    )
    
    category_order <- c(
      "IL1_core",
      "mono_mac_core",
      "granulo_core",
      "IL1_related_core"
    )
    
    category_labs <- c(
      IL1_core="IL-1 axis",
      mono_mac_core="Monocyte/macrophage",
      granulo_core="Granulocytic myeloid",
      IL1_related_core="Extended IL-1"
    )
    
    # automatically generate labels from strat_var
    strat_labs <- setNames(
      ifelse(
        strat_var=="phenoage_cat","Pheno aging",
        ifelse(
          strat_var=="non_ca_base","Aging hallmark",
          paste0(sub("_cat$","",strat_var)," aging")
        )
      ),
      strat_var
    )
    
    plot_dat <- res %>%
      filter(!is.na(HR),!is.na(lower95),!is.na(upper95)) %>%
      mutate(
        
        category=case_when(
          category %in% c("IL1","IL1_core") ~ "IL1_core",
          category %in% c("mono_mac","mono_mac_core") ~ "mono_mac_core",
          category %in% c("granulo","granulo_core") ~ "granulo_core",
          category %in% c("IL1_related","IL1_related_core") ~ "IL1_related_core",
          TRUE ~ category
        ),
        
        category=factor(
          category,
          levels=category_order,
          labels=category_labs[category_order]
        ),
        
        strat_var_lab=recode(strat_var,!!!strat_labs),
        
        strat_var_lab=factor(
          strat_var_lab,
          levels=unname(strat_labs)
        ),
        
        aging_lab=case_when(
          as.character(strat_level)=="0" ~ "accelerated aging -",
          as.character(strat_level)=="1" ~ "accelerated aging +",
          TRUE ~ as.character(strat_level)
        ),
        
        group=case_when(
          aging_lab=="accelerated aging -" & CHIP=="No CHIP" ~
            "No CHIP, accelerated aging -",
          aging_lab=="accelerated aging +" & CHIP=="No CHIP" ~
            "No CHIP, accelerated aging +",
          aging_lab=="accelerated aging -" & CHIP=="CHIP" ~
            "CHIP, accelerated aging -",
          aging_lab=="accelerated aging +" & CHIP=="CHIP" ~
            "CHIP, accelerated aging +",
          TRUE ~ NA_character_
        ),
        
        group=factor(
          group,
          levels=c(
            "No CHIP, accelerated aging -",
            "No CHIP, accelerated aging +",
            "CHIP, accelerated aging -",
            "CHIP, accelerated aging +"
          )
        ),
        
        sig_lab=case_when(
          group=="CHIP, accelerated aging +" & p<0.001 ~ "***",
          group=="CHIP, accelerated aging +" & p<0.01 ~ "**",
          group=="CHIP, accelerated aging +" & p<0.05 ~ "*",
          TRUE ~ ""
        ),
        
        star_y=upper95+0.08
      ) %>%
      filter(
        !is.na(category),
        !is.na(group),
        !is.na(strat_var_lab)
      ) %>%
      group_by(category) %>%
      mutate(
        protein=factor(protein,levels=unique(protein))
      ) %>%
      ungroup()
    
    p <- ggplot(
      plot_dat,
      aes(
        x=protein,
        y=HR,
        ymin=lower95,
        ymax=upper95,
        color=group
      )
    )+
      geom_hline(
        yintercept=1,
        linetype=2,
        linewidth=0.2,
        color="black"
      )+
      geom_errorbar(
        position=position_dodge(width=0.7),
        width=0,
        linewidth=0.2
      )+
      geom_point(
        position=position_dodge(width=0.7),
        shape=16,
        size=0.7
      )+
      geom_text(
        aes(
          y=star_y,
          label=sig_lab,
          group=group
        ),
        color="black",
        position=position_dodge(width=0.7),
        size=2.3,
        vjust=0,
        show.legend=FALSE,
        na.rm=TRUE
      )+
      facet_grid(
        strat_var_lab~category,
        scales="free_x",
        space="free_x"
      )+
      scale_y_continuous(
        breaks=c(1,2,3,4,5),
        limits=c(0,4.5)
      )+
      scale_color_manual(
        values=c(
          "No CHIP, accelerated aging -"="#8FB6D8",
          "No CHIP, accelerated aging +"="#E6B8C3",
          "CHIP, accelerated aging -"="#D79A7E",
          "CHIP, accelerated aging +"="#B22222"
        )
      )+
      labs(
        x=NULL,
        y="HR (95% CI)",
        color=NULL
      )+
      theme_classic(
        base_size=6,
        base_line_size=0.2
      )+
      theme(
        legend.position="top",
        legend.box="horizontal",
        strip.background=element_blank(),
        strip.text.x=element_text(
          color="black",
          margin=margin(t=0,r=0,b=3,l=0,unit="mm")
        ),
        strip.text.y=element_text(
          color="black",
          angle=270
        ),
        axis.text.x=element_text(
          color="black",
          angle=45,
          hjust=1,
          vjust=1
        ),
        axis.text.y=element_text(color="black"),
        axis.title.y=element_text(color="black"),
        legend.key.width=unit(0.3,"cm"),
        panel.spacing.x=unit(0.3,"lines"),
        panel.spacing.y=unit(0.3,"lines"),
        plot.margin=margin(0,0.2,0.2,0.2,"cm"),
        legend.box.margin=margin(0,0,-6,0)
      )
    
    list(
      result=res,
      plot_dat=plot_dat,
      plot=p
    )
  }


  organ_res <- paste0(organ_list,"_res")
  
  tmp_chip_pac <- tmp_chip_km %>%
    filter(ca_base==0 & !is.na(Brain_res)) %>%
    left_join(olink_full,by="eid") %>%
    mutate(across(
      any_of(organ_res),
      ~as.integer(.x>median(.x,na.rm=TRUE)),
      .names="{.col}_cat"
    )) %>%
    rename_with(~sub("_res_cat$","_cat",.x),ends_with("_res_cat"))
  
  
  lung_chip <- plot_protein_cox_by_chip(
    "Lung Cancer",
    strat_var=c("phenoage_cat","non_ca_base",paste0(organ_list,"_cat"))
  )
  lung_chip$plot
  
  
  
  liver_chip <- plot_protein_cox_by_chip(
    "Liver Cancer",
    strat_var=c("Immune_cat","Artery_cat","phenoage_cat","non_ca_base")
  )
  liver_chip$plot
  
  
  bowel_chip <- plot_protein_cox_by_chip(
    "Bowel Cancer",
    strat_var=c("Immune_cat","Intestine_cat","phenoage_cat","non_ca_base")
  )
  bowel_chip$plot
  
  
  
  pancreas_chip <- plot_protein_cox_by_chip(
    "Pancreatic Cancer",
    strat_var=c("Immune_cat","Pancreas_cat","phenoage_cat","non_ca_base")
  )
  pancreas_chip$plot
  
  
  
  
#Prevalent cancer cluster and aging-related diseases

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
# Set FIG3_WORKSPACE to the local RData file loaded and saved by this script.
# Set FIG3_OUTPUT_DIR to an existing directory outside the repository.
fig3_path <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path)) {
    stop(sprintf("Set %s to a local path before running fig3.R.", variable),
         call. = FALSE)
  }
  path.expand(path)
}

load(fig3_path("FIG3_WORKSPACE"))

#Cluster 1: Digestive cancer 
describe(hallmark_base %>% select(c(paste0(c("Oesophageal Cancer","Stomach Cancer","Bowel Cancer","Liver Cancer"),"_reg_death_base22"))))


result_ca_hallmark_c1 <- cox(
  data=hallmark_base %>% filter(non_ca_base!=1),
  exposures=paste0(c("Oesophageal Cancer","Stomach Cancer","Bowel Cancer","Liver Cancer"),"_reg_death_base22"),
  outcomes=c(
    hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>%
      distinct(Disease) %>% pull()
  ),
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
  condition = F
)


C1_aging <- result_ca_hallmark_c1$result_list %>% filter(FDR_P_Value<0.05) %>%
  # filter(!str_detect(Exposure, "Myelodysplasia")) %>% 
  group_by(Outcome) %>%
  filter(n() >= 2) %>%
  ungroup() %>% left_join(hesin_all_clean %>% 
                            filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% dplyr::select(Disease,Category),by=c("Outcome"="Disease")) %>%
  mutate(Category=ifelse(Category=="Haematological/ Immunological","Haematological/\nImmunological",Category)) %>% 
  mutate(Category = case_when(
    Category %in% c("Cardiovascular", "Haematological/\nImmunological", "Infections") ~ Category,
    TRUE ~ "Others"
  )) %>% 
  arrange(Category, desc(-FDR_P_Value)) %>%
  mutate(Category = factor(Category, levels = unique(Category[order(FDR_P_Value)])))  



ggplot(C1_aging %>% mutate(Exposure=gsub("_reg_death_base22","",Exposure)),
       aes(x = HR, y = reorder(Outcome, -FDR_P_Value), shape = Exposure)) + 
  geom_vline(xintercept = 1, color = "grey", linewidth = 0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                 height = 0, linewidth = 0.1,
                 position = position_dodge(0.35),
                 color = "black") +   
  geom_point(size = 0.7,                
             position = position_dodge(0.35),
             color = "black") +         
  scale_shape_manual(values = c(16,17,15,18)) +  
  scale_x_continuous(trans = "log",
                     limits = c(1, 250),
                     breaks = c(1,5,10,25,50,100,200)) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic(base_size = 5, base_line_size = 0.1) +
  theme(
    legend.position = c(0.96, 0.01),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank(),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),   
    strip.text = element_text(size = 4.5),  
    panel.border = element_blank(),       
    panel.spacing = unit(0.05, "cm")      
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free")







ggplot(C1_aging %>% mutate(Exposure=gsub("_reg_death_base22","",Exposure)),
       aes(x = HR, y=(reorder(Outcome, -FDR_P_Value)),color=Exposure)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.1,position=position_dodge(0.35))+
  geom_point(size = 0.5,shape=19,position=position_dodge(0.35))+
  scale_x_continuous(trans = "log",limits = c(1, 250),breaks = c(1,5,10,25,50,100,200)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 5,base_line_size = 0.1)+
  theme(            legend.position = "top",
                    legend.title =element_blank(),
                    legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 5),
                    axis.text.y = element_text(color="black",size = 5),
                    strip.background = element_blank(),)+
  facet_grid(Category~., scales = "free_y",space = "free")+ theme(legend.text = element_text(size = 4, color = "black"),
                                                                  legend.title =element_blank(),
                                                                  legend.key.width = unit(0.15, "cm"),
                                                                  axis.text.x = element_text(color="black",size = 5),
                                                                  axis.text.y = element_text(color="black",size = 5),
                                                                  legend.position = c(0.95, 0.01), 
                                                                  legend.justification = c("right", "bottom"),
                                                                  legend.spacing.y = unit(0.1, "cm"),
                                                                  legend.key.size = unit(0.15, "cm"),   
                                                                  legend.title.align = 1,
                                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                                  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )




#Cluster 2: Hematological cancer
describe(hallmark_base %>% select(c(paste0(c("NHL","Leukaemia","Plasma Cell Cancer","Polycythaemia Vera"),"_reg_death_base22"))))

result_ca_hallmark_c2 <- cox(
  data=hallmark_base %>% filter(non_ca_base!=1),
  exposures=paste0(c("NHL","Leukaemia","Plasma Cell Cancer","Polycythaemia Vera"),"_reg_death_base22"),
  outcomes=c(
    hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>%
      distinct(Disease) %>% pull()
  ),
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
  condition = F
)


C2_aging <- result_ca_hallmark_c2$result_list %>% filter(FDR_P_Value<0.05) %>%
  group_by(Outcome) %>%
  filter(n() >= 2) %>% 
  ungroup() %>% left_join(hesin_all_clean %>% 
                            filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% dplyr::select(Disease,Category),by=c("Outcome"="Disease")) %>%
  mutate(Category=ifelse(Category=="Haematological/ Immunological","Haematological/\nImmunological",Category)) %>% 
  mutate(Category = case_when(
    Category %in% c("Cardiovascular", "Haematological/\nImmunological", "Infections") ~ Category,
    TRUE ~ "Others"
  )) %>% 
  arrange(Category, desc(-FDR_P_Value)) %>%
  mutate(Category = factor(Category, 
                           levels = unique(Category[order(FDR_P_Value)])))  



ggplot(C2_aging %>% mutate(Exposure=gsub("_reg_death_base22","",Exposure)),
       aes(x = HR, y = reorder(Outcome, -FDR_P_Value), shape = Exposure)) + 
  geom_vline(xintercept = 1, color = "grey", linewidth = 0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                 height = 0, linewidth = 0.1,
                 position = position_dodge(0.5),
                 color = "black") +   
  geom_point(size = 0.7,                
             position = position_dodge(0.5),
             color = "black") +         
  scale_shape_manual(values = c(16,17,15,18)) +  
  scale_x_continuous(trans = "log",
                     limits = c(1, 250),
                     breaks = c(1,5,10,25,50,100,200)) +
  labs(x = "HR (95% CI)", y = NULL) +
  theme_classic(base_size = 5, base_line_size = 0.1) +
  theme(
    legend.position = c(0.96, 0.01),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank(),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),   
    strip.text = element_text(size = 4.5),  
    panel.border = element_blank(),       
    panel.spacing = unit(0.05, "cm")      
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free")



#Cluster 3: Others
describe(hallmark_base %>% select(c(paste0(c("Lung Cancer","Breast Cancer","Kidney Cancer","Brain Cancer"),"_reg_death_base22"))))

result_ca_hallmark_c3 <- cox(
  data=hallmark_base %>% filter(non_ca_base!=1),
  exposures=paste0(c("Lung Cancer","Breast Cancer","Kidney Cancer","Brain Cancer"),"_reg_death_base22"),
  outcomes=c(
    hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>%
      distinct(Disease) %>% pull()
  ),
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
  condition = F
)





#Extract the top 10 prot for each organ
top_proteins_list <- list()
for (organ in organ_list) {
  top_proteins <- extract_top_proteins(organ, 10)
  top_proteins_list[[organ]] <- top_proteins
}
top_proteins_list

prot_in_hemato <- top_proteins_list[c("Organismal","Immune","Artery","Brain","Lung","Heart","Stomach","Pancreas","Intestine","Kidney")] %>% 
  unlist() %>%  as.character() %>% unique() 


hallmark_base_prot <- hallmark_base %>% 
  mutate(across(ends_with("_interval22"),~ ., .names = "{sub('_interval22$', '_interval33', .col)}")) %>%
  mutate(across(contains("_interval33"), ~ ( . / 365.25))) %>% 
  mutate(across(ends_with("_interval22"), ~ ifelse(. < 0, NA, .))) %>% 
  mutate(across(contains("_interval22"), ~ ( . / 365.25))) %>% 
  left_join(hallmark_all_final %>% select(c("eid",c(organ_pac,"met_res","phenoage_res","kdm_res")))) %>%  
  left_join(olink_full %>% dplyr::select(all_of((prot_in_hemato)), "eid"), by = "eid") %>%
  mutate(across(all_of((prot_in_hemato)), ~ scale(as.numeric(.))))

calculate_beta_matrix <- function(data, var_list1, var_list2, covariates = NULL) {
  beta_matrix <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  fdr_p_values <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  
  for (i in seq_along(var_list1)) {
    for (j in seq_along(var_list2)) {
      formula <- as.formula(paste0(var_list1[i], "~", "`",var_list2[j],"`"))
      if (!is.null(covariates)) {
        covariate_formula <- paste(covariates, collapse = " + ")
        formula <- update(formula, paste(". ~ . +", covariate_formula))
      }
      lm_model <- lm(formula, data = data)
      beta_matrix[i, j] <- coef(lm_model)[2]  # Store the beta coefficient
      
      p_value <- summary(lm_model)$coefficients[2, "Pr(>|t|)"]
      fdr_p_values[i, j] <- p.adjust(p_value, method = "fdr")
      
    }
  }
  
  colnames(beta_matrix) <- var_list2
  rownames(beta_matrix) <- var_list1
  
  return(list(beta_matrix = beta_matrix, fdr_p_values = fdr_p_values))
}


result_hemato_pac <- calculate_beta_matrix(hallmark_base_prot %>% filter(!is.na(Organismal_res)) %>% filter(non_ca_base!=1),
                                           var_list1=c(prot_in_hemato), 
                                           var_list2=paste0(c("Polycythaemia Vera","Plasma Cell Cancer","Leukaemia","NHL",#Hemo cancer
                                                              "Oesophageal Cancer","Stomach Cancer","Bowel Cancer","Liver Cancer",
                                                              "Lung Cancer","Breast Cancer","Kidney Cancer","Brain Cancer"
                                                              
                                           ),"_reg_death_interval22"),
                                           covariates = c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"))


mat_hemato <- cbind(result_hemato_pac$beta_matrix,result_hemato_pac$fdr_p_values) %>% as.data.frame() %>% 
  filter(rowSums(across(13:24, ~ . < 0.05)) >= 2) %>% 
  mutate(Cat=case_when(      rownames(.)  %in% top_proteins_list$Immune~"Immune",
                             rownames(.)  %in% top_proteins_list$Artery~"Artery",
                             rownames(.)  %in% top_proteins_list$Brain~"Brain",
                             rownames(.)  %in% top_proteins_list$Lung~"Lung",
                             rownames(.)  %in% top_proteins_list$Heart~"Heart",
                             rownames(.)  %in% top_proteins_list$Stomach~"Stomach",
                             rownames(.)  %in% top_proteins_list$Pancreas~"Pancreas",
                             rownames(.)  %in% top_proteins_list$Intestine~"Intestine",
                             rownames(.)  %in% top_proteins_list$Kidney~"Kidney",
                             rownames(.)  %in% top_proteins_list$Organismal~"Organismal",
                             T~rownames(.))
  )

mat_beta <- (as.matrix(mat_hemato[ , c(1:12)])) %>% t() 
rownames(mat_beta) <- gsub("_reg_death_interval22", "", rownames(mat_beta))
mat_p <- (as.matrix(mat_hemato[ , c(13:24)]))%>% t()


#PAC and survival in cancer-specific patients----------------------------
pac_cancer_survival <- function(ca_list,outcome_list){
  results <- data.frame()
  for (ca in ca_list) {
    ca <- ifelse(ca == "ca",
                 paste0(ca, "_base"),
                 paste0(ca, "_base22"))
    
    result_in_ca_aging_disease_main<- cox(
      data=hallmark_base_prot %>% filter(non_ca_base==0)%>% filter(!!rlang::sym(ca)==1) ,
      outcomes =outcome_list,
      exposures=c(row.names(mat_hemato)),
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
      condition = F)
    
    results <- rbind(results,result_in_ca_aging_disease_main$result_list %>% mutate(cat=ca))
  }
  return(results)
}

pac_nonca_survival <- function(ca_list,outcome_list){
  results <- data.frame()
  for (ca in ca_list) {
    ca <- ifelse(ca == "ca",
                 paste0(ca, "_base"),
                 paste0(ca, "_base22"))
    
    result_in_ca_aging_disease_main<- cox(
      data=hallmark_base_prot %>% filter(non_ca_base==0)%>% filter(!!rlang::sym(ca)==0) ,
      outcomes =outcome_list,
      exposures=c(row.names(mat_hemato)),
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
      condition = F)
    
    results <- rbind(results,result_in_ca_aging_disease_main$result_list )
  }
  return(results)
}

pac_cancer_survival_results <- pac_cancer_survival(ca_list = c("ca"),outcome_list="death")
pac_nonca_survival_results <- pac_nonca_survival(ca_list = c("ca"),outcome_list="death")

pac_cancer_any_results <- pac_cancer_survival(ca_list = c("ca"),outcome_list="any_1")
pac_nonca_any_results <- pac_nonca_survival(ca_list = c("ca"),outcome_list="any_1")

pac_cancer_age_results <-linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                                b=row.names(mat_hemato),
                                a="age",covar = c("sex","ethnic","edu","tdi","center","pa","smok","bmi"))


pac_cancer_chronic_results <- cox( 
  data=hallmark_base_prot %>% filter(non_ca_base==0)%>% filter(ca_base==0),
  outcomes= c(
    "Cardiovascular",
    "Neurological",
    "Respiratory",
    "Haematological/ Immunological",
    "Musculoskeletal",
    "Genitourinary",
    "Digestive",
    "Eye",
    "Infections",
    "Endocrine",
    "Skin"
  ),
  exposures=c(row.names(mat_hemato)),
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi"),
  condition = F
)


chronic_mat <- pac_cancer_chronic_results$result_list %>%
  dplyr::mutate(sig = FDR_P_Value < 0.05) %>%
  dplyr::select(Exposure, Outcome, sig) %>%
  tidyr::pivot_wider(
    names_from = Outcome,
    values_from = sig,
    values_fill = FALSE
  ) %>%
  tibble::column_to_rownames("Exposure") %>%
  as.matrix()


values = c(
  "Cardiovascular" = "#A52A2A",
  "Neurological" = "#4B0082",
  "Respiratory" = "#00FFFF",
  "Haematological/\nImmunological" = "#FFFF00",
  "Musculoskeletal" = "#808080",
  "Genitourinary" = "#808000",
  "Digestive" = "#008000",
  "Eye" = "#FF69B4",
  "Infections" = "#0000FF",
  "Endocrine" = "#D2B48C",
  "Skin"="black"
)


mat_hemato_all <- mat_hemato %>% mutate(Exposure=row.names(.)) %>% 
  left_join(.,pac_cancer_survival_results %>% dplyr::select(Exposure,Beta,FDR_P_Value) %>% dplyr::rename("Death (CA)"=Beta),by="Exposure") %>% 
  left_join(.,pac_nonca_survival_results %>% dplyr::select(Exposure,Beta,FDR_P_Value) %>% dplyr::rename("Death (HC)"=Beta),by="Exposure") %>% 
  left_join(.,pac_cancer_any_results %>% dplyr::select(Exposure,Beta,FDR_P_Value) %>% dplyr::rename("Any Hallmarks (CA)"=Beta),by="Exposure") %>%
  left_join(.,pac_nonca_any_results %>% dplyr::select(Exposure,Beta,FDR_P_Value) %>% dplyr::rename("Any Hallmarks (HC)"=Beta),by="Exposure") %>%
  left_join(.,pac_cancer_age_results %>% dplyr::select(var_b,beta,p) %>% mutate(q=p.adjust(p, method = "fdr")) %>% dplyr::rename("Organ aging"=beta,"Exposure"=var_b),by="Exposure")

row.names(mat_hemato_all) <- mat_hemato_all$Exposure

mat_beta <- mat_hemato_all %>%
  dplyr::select(
    ends_with("_interval22"),
    `Death (CA)`,
    `Death (HC)`,
    `Any Hallmarks (CA)`,
    `Any Hallmarks (HC)`,
    `Organ aging`
  ) %>%
  as.matrix() %>%
  t()

rownames(mat_beta) <- gsub("_reg_death_interval22", "", rownames(mat_beta))

mat_p <- (as.matrix(mat_hemato_all[, c(13:24,28,30,32,34,37)]))%>% t()


## protein -> organ category
protein_cat <- mat_hemato$Cat
names(protein_cat) <- rownames(mat_hemato)

sig_count <- t(apply(mat_p, 1, function(p) {
  sig_protein <- names(p)[p < 0.05]
  sig_cat <- protein_cat[sig_protein]
  
  table(factor(sig_cat, levels = unique(protein_cat)))
}))

sig_count <- as.matrix(sig_count)
rownames(sig_count) <- rownames(mat_p)

## convert count to row proportion
sig_prop <- sig_count / rowSums(sig_count)
sig_prop[is.na(sig_prop)] <- 0

organ_cols <- c(
  Artery     = "#A52A2A",
  Stomach    = "#008000",##0000FF
  Lung       = "#00FFFF",
  Immune     = "#FFFF00",
  Heart      = "#FF7F50",
  Intestine  = "#008000",#"#808080",
  Pancreas   = "#008000",#"#008000",
  Brain      = "#4B0082",
  Kidney     = "#808000",
  Organismal = "black"
)

right_anno <- rowAnnotation(
  `Organ protein\nproportion` = anno_barplot(
    sig_prop,
    gp = gpar(fill = organ_cols[colnames(sig_prop)],lwd = 0.2),border = FALSE, 
    bar_width = 0.8,
    width = unit(1, "cm"),
    axis_param = list(
      at = c(0, 0.5, 1),labels_rot = 45, 
      labels = c("0", "0.5", "1"),
      gp = gpar(fontsize = 4,lwd=0.2)
    )
  ),
  annotation_name_gp = gpar(fontsize = 4)
)

top_anno <- HeatmapAnnotation(
  Cat = mat_hemato$Cat,
  col = list(Cat = organ_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "mm"),
  show_legend = FALSE,
  annotation_legend_param = list(
    Cat = list(
      title = "Organ",                         
      title_gp = gpar(fontsize = 5),          
      labels_gp = gpar(fontsize = 4.5),
      grid_width  = unit(1, "mm"),
      grid_height = unit(1, "mm")
    )
  )
)
colnames(chronic_mat) <- gsub("Haematological/\\s*Immunological", "Immunological", colnames(chronic_mat))
chronic_mat <- chronic_mat[colnames(mat_beta), , drop = FALSE]

draw_chronic_label <- TRUE
bottom_anno <- HeatmapAnnotation(
  Chronic = AnnotationFunction(
    which = "column",
    n = ncol(mat_beta),
    height = unit(ncol(chronic_mat) * 1.8, "mm"),
    fun = function(index) {
      
      m <- chronic_mat[index, , drop = FALSE]
      
      nr <- ncol(m)
      nc <- nrow(m)
      
      pushViewport(viewport(
        xscale = c(0, nc),
        yscale = c(0, nr)
      ))
      
      ## category names: only draw once
      if (draw_chronic_label) {
        for (j in seq_len(nr)) {
          y0 <- nr - j + 0.5
          
          grid.text(
            colnames(chronic_mat)[j],
            x = unit(-0.5, "native"),
            y = unit(y0, "native"),
            just = "right",
            gp = gpar(fontsize = 5)
          )
        }
        draw_chronic_label <<- FALSE
      }
      
      for (j in seq_len(nr)) {
        y0 <- nr - j + 0.5
        
        grid.rect(
          x = unit(nc / 2, "native"),
          y = unit(y0, "native"),
          width = unit(nc, "native"),
          height = unit(1, "native"),
          gp = gpar(fill = "#EAEAEA", col = "white", lwd = 0.2)
        )
      }
      
      ## white dots + significant black dots
      for (i in seq_len(nc)) {
        x0 <- i - 0.5
        y_all <- nr - seq_len(nr) + 0.5
        
        grid.points(
          x = unit(rep(x0, nr), "native"),
          y = unit(y_all, "native"),
          pch = 16,
          size = unit(1, "mm"),
          gp = gpar(col = "white")
        )
        
        sig_y <- which(m[i, ])
        
        if (length(sig_y) > 0) {
          y_sig <- nr - sig_y + 0.5
          
          if (length(y_sig) > 1) {
            grid.lines(
              x = unit(rep(x0, length(y_sig)), "native"),
              y = unit(y_sig, "native"),
              gp = gpar(col = "black", lwd = 0.4)
            )
          }
          
          grid.points(
            x = unit(rep(x0, length(y_sig)), "native"),
            y = unit(y_sig, "native"),
            pch = 16,
            size = unit(1, "mm"),
            gp = gpar(col = "black")
          )
        }
      }
      
      popViewport()
    }
  ),
  show_annotation_name = FALSE
)


p_heat <- grid.grabExpr(draw( 
  Heatmap(mat_beta,
          name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
          col = colorRamp2(c(-0.35,0,0.35), c("#007F7F","white","#D85F4D")),
          show_row_names = TRUE, 
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 5),
          cluster_rows = F, cluster_columns = T,
          row_dend_width = unit(0.3, "cm"),
          row_dend_gp = gpar(lwd = 0.2),
          column_dend_height  = unit(0.3, "cm"),
          column_dend_gp = gpar(lwd = 0.2),
          
          show_column_dend = T, show_row_dend = F,
          
          right_annotation = right_anno,
          top_annotation = top_anno,
          bottom_annotation = bottom_anno,
          
          column_names_rot=45,
          row_split=factor(c(rep("Cancer",12),rep("Survival",4),rep("Age",1)),levels = c("Age","Cancer", "Survival")),
          column_title_gp = gpar(fontsize = 4.5),
          column_split = ((mat_hemato %>% dplyr::select(Cat))),
          column_title_side = "top",
          row_names_side = "left",
          row_title = NULL,
          #column_title = NULL,
          # cell_fun = function(j, i, x, y, w, h, fill) {
          #   if((mat_p)[i, j] < 0.001) {
          #     grid.text("***", x, y,gp = gpar(fontsize = 3.5))
          #   } else if((mat_p)[i, j] < 0.01) {
          #     grid.text("**", x, y,gp = gpar(fontsize = 3.5))
          #   } else if((mat_p)[i, j] < 0.05) {
          #     grid.text("*", x, y,gp = gpar(fontsize = 3.5))
          #   }
          # },
          cell_fun = function(j, i, x, y, w, h, fill) {
            
            if (is.na(mat_p[i, j]) || mat_p[i, j] >= 0.05) {
              grid.rect(
                x, y, w, h,
                gp = gpar(fill = "white", col = "black", lwd = 0.2)
              )
            }
            
            if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.001) {
              grid.text("***", x, y, gp = gpar(fontsize = 3.5))
            } else if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.01) {
              grid.text("**", x, y, gp = gpar(fontsize = 3.5))
            } else if (!is.na(mat_p[i, j]) && mat_p[i, j] < 0.05) {
              grid.text("*", x, y, gp = gpar(fontsize = 3.5))
            }
          },
          heatmap_legend_param = list(
            at = c(-0.35, 0, 0.35),
            # labels = c(-0.1, 0, 0.5),
            title = "Beta",
            legend_height = unit(1, "cm"),
            grid_width = unit(0.5, "mm"),
            legend_width =unit(1, "mm"),
            #title_position = "lefttop-rot",
            labels_gp = gpar(col = "black", fontsize = 5),
            title_gp= gpar(col = "black", fontsize = 5)
            
          ),
          width = ncol(mat_beta)*unit(2.3, "mm"), 
          height =nrow(mat_beta)*unit(1.8, "mm"),
          
          
  ),padding = unit(c(0, 0, 0, 0), "mm"),heatmap_legend_side="right",gap = unit(0, "mm")
)
) %>% plot_grid()+theme(plot.margin = margin(0, 0, 0, 0))



cairo_pdf(file.path(fig3_path("FIG3_OUTPUT_DIR"), "Fig4a.pdf"), 
          width = 7.2, height = 3)
p_heat

dev.off()  


#Time-to-diagnosis plot

plot_interval_hist <- function(data, cancer, xlabel=F){
  
  var <- paste0(cancer, "_reg_death_interval33")
  
  group_levels <- c(
    ">5 years before",
    "3–5 years before",
    "1–3 years before",
    "1 year before/after",
    "1–3 years after",
    "3–5 years after",
    ">5 years after"
  )
  
  df <- data %>%
    filter(!is.na(.data[[var]]),
           .data[[var]] != 0) %>%
    mutate(
      x = .data[[var]],
      x = ifelse(x > 25, 25, x),
      group = cut(
        x,
        breaks = c(-Inf, -5, -3, -1, 1, 3, 5, Inf),
        labels = group_levels,
        right = FALSE
      ),
      group = factor(group, levels = group_levels)
    )
  
  ggplot(df, aes(x = x, fill = group)) +
    geom_histogram(binwidth = 1, boundary = 0, color = NA) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c(
      ">5 years before"      = "#556B4F",
      "3–5 years before"    = "#D8D0A3",
      "1–3 years before"    = "#B9854A",
      "1 year before/after" = "#9B2F2F",
      "1–3 years after"     = "#6B1A14",
      "3–5 years after"     = "#4B0B06",
      ">5 years after"      = "#2B0302"
    ), breaks = group_levels) +
    labs(
      title = paste0(cancer),
      x = if (xlabel) "Years to diagnosis" else NULL,
      y = "Count",
      fill = "Group"
    ) +
    theme_classic(base_size = 5,base_line_size = 0.2) +
    theme(
      plot.title = element_text(hjust = 0.5,size=5,face = "bold"),
      legend.position = "none",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(0.1, "cm"),
      legend.key.height = unit(0.005, "cm"),
      # legend.spacing.x = unit(0.05, "cm")
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

plot_protein_by_interval <- function(data, cancer, protein){
  
  var <- paste0(cancer, "_reg_death_interval33")
  
  group_levels <- c(
    "Control",
    ">5 years before",
    "3–5 years before",
    "1–3 years before",
    "1 year before/after",
    "1–3 years after",
    "3–5 years after",
    ">5 years after"
  )
  
  group_cols <- c(
    "Control"              = "grey80",
    ">5 years before"      = "#556B4F",
    "3–5 years before"    = "#D8D0A3",
    "1–3 years before"    = "#B9854A",
    "1 year before/after" = "#9B2F2F",
    "1–3 years after"     = "#6B1A14",
    "3–5 years after"     = "#4B0B06",
    ">5 years after"      = "#2B0302"
  )
  
  df <- data %>%
    filter(!is.na(.data[[var]]),
           !is.na(.data[[protein]])) %>%
    mutate(
      x = .data[[var]],
      x = ifelse(x > 25, 25, x),
      group = ifelse(
        x == 0,
        "Control",
        as.character(cut(
          x,
          breaks = c(-Inf, -5, -3, -1, 1, 3, 5, Inf),
          labels = group_levels[-1],
          right = FALSE
        ))
      ),
      group = factor(group, levels = group_levels)
    )
  
  protein_organ <- protein_cat[[protein]]
  organ_col <- organ_cols[[protein_organ]]
  
  get_text_col <- function(col){
    rgb <- col2rgb(col) / 255
    lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
    if (lum > 0.6) "black" else "white"
  }
  txt_col <- get_text_col(organ_col)
  
  p <- ggplot(df, aes(x = group, y = .data[[protein]], fill = group)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.7,
      linewidth = 0.2
    ) +
    scale_fill_manual(values = group_cols) +
    labs(
      x = NULL,
      y = paste0(protein, " (SD)"),
      title = protein
    ) +
    theme_classic(base_size = 6, base_line_size = 0.2) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(size = 5),
      axis.text.y  = element_text(size = 5),
      plot.title = element_text(hjust = 0.5, size = 5, margin = margin(b = 0))
    )
  
  protein_organ <- protein_cat[[protein]]
  organ_col <- organ_cols[[protein_organ]]
  
  get_text_col <- function(col){
    rgb <- col2rgb(col)/255
    lum <- 0.299*rgb[1] + 0.587*rgb[2] + 0.114*rgb[3]
    if(lum > 0.6) "black" else "white"
  }
  txt_col <- get_text_col(organ_col)
  
  p +
    annotation_custom(
      grob = grobTree(
        rectGrob(
          x = 0.5,
          y = 0.95,
          width = 1,
          height = 0.1,
          gp = gpar(fill = organ_col, col = NA)
        ),
        textGrob(
          protein_organ,
          x = 0.5,
          y = 0.95,
          gp = gpar(col = txt_col, fontsize = 4.5)
        )
      ),
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    )
}


p_leu <- plot_grid(
  plot_interval_hist(hallmark_base_prot, "Leukaemia"),
  plot_protein_by_interval(hallmark_base_prot, "Leukaemia", "FCRL1"),
  plot_protein_by_interval(hallmark_base_prot, "Leukaemia", "FCAMR"),
  plot_protein_by_interval(hallmark_base_prot, "Leukaemia", "SPINK1"),
  plot_protein_by_interval(hallmark_base_prot, "Leukaemia", "LTBP2"),
  ncol = 5, align = "hv"
)

p_sto <- plot_grid(
  plot_interval_hist(hallmark_base_prot, "Stomach Cancer"),
  plot_protein_by_interval(hallmark_base_prot, "Stomach Cancer", "TFF2"),
  plot_protein_by_interval(hallmark_base_prot, "Stomach Cancer", "GUCA2A"),
  plot_protein_by_interval(hallmark_base_prot, "Stomach Cancer", "SPINK1"),
  plot_protein_by_interval(hallmark_base_prot, "Stomach Cancer", "SPINK4"),
  ncol = 5, align = "hv"
)

p_kid <- plot_grid(
  plot_interval_hist(hallmark_base_prot, "Kidney Cancer",xlabel = T),
  plot_protein_by_interval(hallmark_base_prot, "Kidney Cancer", "TFF2"),
  plot_protein_by_interval(hallmark_base_prot, "Kidney Cancer", "GUCA2A"),
  plot_protein_by_interval(hallmark_base_prot, "Kidney Cancer", "NEFL"),
  plot_protein_by_interval(hallmark_base_prot, "Kidney Cancer", "CALCA"),
  ncol = 5, align = "hv"
)


cairo_pdf(file.path(fig3_path("FIG3_OUTPUT_DIR"), "Fig4b.pdf"), 
          width = 6, height = 2.7)

plot_grid(
  p_leu,p_sto,p_kid,ncol = 1,align = "hv",labels = c("d", "e", "f"),  label_size = 8, label_fontface = "bold"
)

dev.off()  




#Set legend 
group_cols <- c(
  ">5 years before"      = "#556B4F",
  "3–5 years before"    = "#D8D0A3",
  "1–3 years before"    = "#B9854A",
  "1 year before/after" = "#9B2F2F",
  "1–3 years after"     = "#6B1A14",
  "3–5 years after"     = "#4B0B06",
  ">5 years after"      = "#2B0302",
  "Control"             = "grey80"
)

# dummy data
df_leg <- data.frame(
  x = 1,
  y = 1,
  grp = factor(names(group_cols), levels = names(group_cols))
)

p_leg <- ggplot(df_leg, aes(x, y, fill = grp)) +
  geom_tile() +
  scale_fill_manual(values = group_cols) +
  guides(fill = guide_legend(
    nrow = 1,
    byrow = TRUE,label.vjust = 0.5, 
    override.aes = list(
      colour = NA
    )
  )) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.width  = unit(0.18, "cm"),
    legend.key.height = unit(0.12, "cm"),
    legend.text = element_text(size = 5, margin = margin(l = 2)),
    legend.key.spacing.x = unit(0.03, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.text.align = 0.5              
  )


cairo_pdf(file.path(fig3_path("FIG3_OUTPUT_DIR"), "Fig4_lgd.pdf"), 
          width = 5, height = 1)
p_leg

dev.off() 



### Saving working image--------
save.image(fig3_path("FIG3_WORKSPACE"))





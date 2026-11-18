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
}
# Set the FIG2_* environment variables to your local paths before running.
# FIG2_SHAP_DIR contains the shap_<organ>.csv inputs.
# Set both output file paths in an existing directory outside the repository.
fig2_path <- function(variable) {
  path <- Sys.getenv(variable)
  if (!nzchar(path)) {
    stop(sprintf("Set %s to a local path before running fig2.R.", variable),
         call. = FALSE)
  }
  path.expand(path)
}

load(fig2_path("FIG2_FIG1_WORKSPACE"))

# Individual Aging hallmark disease -> incident cancer-------------------------------------------------------------------------------------------  

#Lung cancer----------
hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% 
  filter(if_any(c(GI, TA, EA, DNS, MD, CS, SCE, AIC), ~ . == 1)) %>% distinct(Disease) %>% pull()
hallmark_base$`Autonomic Neuropathy_base22` %>% as.factor() %>% describe()

table(hallmark_base$ca_base)

result_hallmark_disease_lung <- cox(
  data=hallmark_base %>% filter(ca_base==0),
  exposures=c(
    paste0(hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% 
             filter(if_any(c(GI, TA, EA, DNS, MD, CS, SCE, AIC), ~ . == 1)) %>% distinct(Disease) %>% pull(),"_base22")
  ),
  outcomes=c(
    "Lung Cancer_reg_death"
  ), 
  covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi")
)

aging_lung <- result_hallmark_disease_lung$result_list %>% filter(FDR_P_Value<0.05)

#Individual aging disease with lung cancer
  lung_a <- ggplot(aging_lung %>% mutate(Exposure=gsub("_base22","",Exposure)) %>% 
                     dplyr::arrange(desc(-FDR_P_Value)) %>%
                     dplyr::slice(1:10),
                   aes(x = HR, y=(reorder(Exposure, HR)))) +
    geom_vline(xintercept = 1, color = "grey",linewidth=0.2) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
    geom_point(size = 1,shape=15,position=position_dodge(0.7))+
    #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
    scale_y_discrete()+
    labs(x = "HR (95% CI)", y = NULL, title = NULL) +
    theme_classic(base_size = 5,base_line_size = 0.2)+
    theme(            legend.position = "top",
                      legend.title =element_blank(),
                      legend.key.width = unit(0.1, "cm"),
                      axis.text.x = element_text(color="black",size = 5),
                      axis.text.y = element_text(color="black",size = 5))

  #Mapping significant disease to Aging hallmarks
  disease_map %>% 
    dplyr::filter(
      Disease %in% (
        aging_lung %>% 
          dplyr::mutate(Exposure = gsub("_base22", "", Exposure)) %>% 
          dplyr::arrange(FDR_P_Value) %>% 
          dplyr::pull(Exposure)
      )
    ) %>% 
    dplyr::summarise(
      dplyr::across(c(GI, TA, EA, DNS, MD, CS, SCE, AIC), ~ sum(. == 1, na.rm = TRUE))
    )
##Individual association of aging protein with lung cancer    
  ##Extract top n proteins---------------------------------------------------------
  organ_list <- c(
    "Organismal",
    "Immune","Artery","Adipose","Muscle","Skin",
    "Brain","Esophagus",
    "Lung","Heart",
    "Liver","Stomach","Pancreas","Intestine","Kidney","Male"
  )
  
  organ_pac <- paste0(organ_list, "_res")
  organ_extra <- paste0(organ_pac, "_extra")
  organ_super <- paste0(organ_pac, "_super")
  
  extract_top_proteins <- function(organ, num_proteins) {
    file_name <- file.path(fig2_path("FIG2_SHAP_DIR"), paste0("shap_", organ, ".csv"))
    data <- read.csv(file = file_name) %>% 
      dplyr::select(-eid) %>% 
      mutate_all(abs) %>%  
      summarise_all(mean, na.rm = TRUE) %>% 
      gather(key = "protein", value = "mean_value") %>%
      arrange(desc(mean_value))
    
    if (nrow(data) < num_proteins) {
      num_proteins <- nrow(data)
    }
    
    top_proteins <- data$protein[1:num_proteins]
    return(top_proteins)
  }
  #Extract the top 20 prot for each organ
  top_proteins_list <- list()
  for (organ in organ_list) {
    top_proteins <- extract_top_proteins(organ, 20)
    top_proteins_list[[organ]] <- top_proteins
  }
  top_proteins_list
  
  #Extract the top 20 prot for selected organ
  lung_pac <- c("Organismal","Skin", "Brain", "Esophagus", "Lung", "Liver")
  
  prot_in <- top_proteins_list[lung_pac] %>% unlist() %>%  as.character() %>% unique() 
  
  hallmark_base_prot <- hallmark_base %>%
    left_join(hallmark_all_final %>% select(c("eid",c(organ_pac,"met_res","phenoage_res","kdm_res")))) %>%  
    left_join(olink_full %>% dplyr::select(all_of((prot_in)), "eid"), by = "eid") %>%
    mutate(across(all_of((prot_in)), ~ scale(as.numeric(.))))
  
  hallmark_base_prot %>% filter(!is.na(Organismal_res)) %>% filter(ca_base!=1) %>% nrow() #n=49702

  prot_lung_result <- cox(data=hallmark_base_prot %>% filter(!is.na(Organismal_res)) %>% filter(ca_base!=1) ,
                          exposures=prot_in,
                          outcomes=c("Lung Cancer_reg_death"),
                          covar=c("age","sex","ethnic","tdi","center","pa","smok","bmi","non_ca_base"))

  
  add_age_results <-linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                           b=prot_in,
                           a="age", covar = c("sex","ethnic","edu","tdi","center","pa","smok","bmi"))
  
  
  prot_lung <- prot_lung_result$result_list %>% 
    mutate(Cat=case_when(Exposure %in% top_proteins_list$Brain~"Brain",
                         Exposure %in% top_proteins_list$Liver~"Liver",
                         Exposure %in% top_proteins_list$Lung~"Lung",
                         Exposure %in% top_proteins_list$Skin~"Skin",
                         Exposure %in% top_proteins_list$Esophagus~"Esophagus",
                         Exposure %in% top_proteins_list$Organismal~"Organismal"
  )) %>% 
    left_join(
      add_age_results %>% dplyr::select(var_b, beta),
      by = c("Exposure" = "var_b")
    ) %>%
    mutate(
      aging = ifelse(Beta / beta > 0, "Synergistic", "Antagonistic"),
      aging = factor(aging)
    )
  
  table(prot_lung$Cat)
  
  lung_b <- ggplot(prot_lung, aes(x = HR, y = -log10(FDR_P_Value), color = Cat,shape = aging)) +
      geom_point(size = 1) +
      theme_bw() +
      ggtitle("Lung cancer")+
      labs(x = "HR", y = expression(-log[10](italic(q)~value)), color = NULL) +
      theme(legend.position = "right") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black",lwd=0.2) + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",lwd=0.2)+
      geom_text_repel(data = prot_lung %>% filter(FDR_P_Value<0.05) %>%  arrange(FDR_P_Value) %>% head(15), aes(label = Exposure,color=Cat), size =2 ,segment.size = 0.2, max.overlaps = 15, show.legend = FALSE)+
      theme_classic(base_size = 6,base_line_size = 0.2)+
      theme(            legend.position = c(.01, .9), legend.justification = c("left", "top"),
                        plot.title = element_text(size = 7, face = "bold"),
                        legend.title =element_blank(),
                        legend.key.width = unit(0.1, "cm"),
                        axis.text.x = element_text(color="black",size = 5,angle = 0,vjust = 0.5),
                        axis.text.y = element_text(color="black",size = 6),
                        legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                        legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                        legend.key.spacing.y = unit(0, "pt"),
                        legend.box.spacing = margin(0),
                        # legend.margin = margin(t = 1, b =1, unit = "pt"),
                        # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                        strip.background = element_blank(),
                        # panel.spacing = unit(0.01, "cm")
                        
      )+scale_color_manual(values=c("Organismal"="black","Brain"="#4B0082","Skin"="lightgrey","Lung"="#00FFFF","Liver"="#D2B48C","Immune"="#feb24c","Esophagus"="#008000"))+
      scale_shape_manual(values = c("Synergistic" = 17, "Antagonistic" = 16))  
  
  
  prot_lung_top_lung <- prot_lung %>% 
    filter((Cat %in% c("Organismal","Lung","Esophagus","Liver"))&FDR_P_Value<0.05) %>% 
    dplyr::select(Exposure) %>% pull
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  go_enrichment <- enrichGO(gene = prot_lung_top_lung, 
                            OrgDb = org.Hs.eg.db,
                            keyType = "SYMBOL",   
                            ont = "BP",
                            pvalueCutoff = 0.05
  )
  
  enrichment_bp_lung <- as.data.frame(go_enrichment)
  dotplot(go_enrichment, showCategory =5)+theme_bw()
  barplot(go_enrichment, showCategory =5)+theme_bw()
  
  go_enrichment@result <- go_enrichment@result %>%
    dplyr::arrange(qvalue) %>%
    dplyr::mutate(has_IL1RN = grepl("(^|/)IL1RN(/|$)", geneID)) %>%
    dplyr::filter(!has_IL1RN | row_number() == which(has_IL1RN)[1]) %>%
    dplyr::select(-has_IL1RN)
  
  lung_c <- {
    res2 <- go_enrichment@result %>%
      as.data.frame() %>%
      filter(Description!="positive regulation of heterotypic cell-cell adhesion") %>% 
      mutate(
        Description = dplyr::recode(
          Description,
          "respiratory gaseous exchange by respiratory system" = "respiratory gaseous exchange"
        )
      ) %>% 
      dplyr::arrange(qvalue) %>%
      dplyr::slice(1:5) %>%
      dplyr::mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      ) %>%
      dplyr::arrange(GeneRatio_num) %>%
      dplyr::mutate(
        Description = paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description))),
        Description = stringr::str_wrap(Description, width = 28),
        Description = factor(Description, levels = Description)
      )
    
    q_breaks <- signif(seq(min(res2$qvalue), max(res2$qvalue), length.out = 3), 2)
    
    ggplot(res2, aes(x = GeneRatio_num, y = Description)) +
      geom_point(aes(size = Count, color = qvalue)) +
      scale_size_continuous(range = c(0.5, 1.5), breaks = c(2, 4, 6),name = "Gene count") +
      scale_color_gradient(
        low = "#023858", high = "#d0d1e6",expression(italic(q)~value),
        breaks = q_breaks,
        labels = format(q_breaks, scientific = TRUE, digits = 1)
      ) +
      labs(x = "Gene ratio", y = NULL) +
      guides(
        color = guide_colorbar(
          direction = "horizontal",
          barwidth = unit(0.75, "cm"),
          barheight = unit(0.06, "cm"),
          ticks = FALSE
        ),
        size = guide_legend(
          direction = "horizontal",
          override.aes = list(color = "black"),
          keywidth = unit(0.1, "cm"),
          keyheight = unit(0.1, "cm")
        )
      ) +
      theme_bw(base_size = 5, base_line_size = 0.1) +
      theme(
        panel.grid = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.direction = "horizontal",
        legend.title = element_text(size = 3.5),
        legend.text = element_text(size = 3.5),
        axis.text.x = element_text(size = 5, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        legend.spacing.x = unit(0.01, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(0.08, "cm"),
        legend.key.height = unit(0.08, "cm"),
        legend.background = element_rect(fill = "white", colour = NA)
      )
  }
  
  lung_c
  
  
  #sc-RNAseq of lung protein (annotated with lung aging and lung cancer)
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  lung_seq <-  read.csv(file = fig2_path("FIG2_LUNG_RNA_CSV"),skip=9) %>%
    filter(Cell.Type %in% c(
      "lung macrophage",
      # "stem cell",
      # "mucus secreting cell",
      "respiratory goblet cell",
      "lung secretory cell",
      "type II pneumocyte",
      "club cell",
      "multi-ciliated epithelial cell",
      "blood vessel endothelial cell",
      "type I pneumocyte",
      "T cell","B cell","fibroblast"
    )) %>%
    mutate(Cell.Type=case_when(
      # Cell.Type=="aggregated"~"Brain",
      Cell.Type=="microglial cell"~"Microglia",
      Cell.Type=="blood vessel endothelial cell"~"Vascular endothelial cell",
      T~ firstup(Cell.Type))) %>% 
    # mutate(Cell.Type=factor(Cell.Type,levels=c("Brain","Neuron","Oligodendrocyte","Astrocyte","Microglia","Endothelial cell"))) %>%
    mutate(Tissue.Composition=gsub("%", "", Tissue.Composition)) %>% mutate(Tissue.Composition=as.numeric(Tissue.Composition) / 100) %>% 
    mutate(Tissue.Composition=ifelse(Tissue.Composition==0,0.001,Tissue.Composition)) %>%
    mutate(expression_all=ifelse(!is.na(Tissue.Composition)&!is.na(Expression),Expression*(Tissue.Composition),NA)) %>% 
    group_by(Gene.Symbol) %>% mutate(expression_nom=scale(expression_all)) %>% ungroup() %>% 
    # mutate(`Synaptic signaling`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[1]] %>% unlist(),1,0),
    #        `Axon development`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[2]] %>% unlist(),1,0)) %>%
    mutate(`Lung cancer`=ifelse(Gene.Symbol %in% (prot_lung %>% filter(Cat=="Lung"&FDR_P_Value<0.05) %>% dplyr::select(Exposure) %>% pull),1,0),
    ) %>%
    arrange(Cell.Type,Gene.Symbol)
  

  lung_prot_seq <- {
    mat_lung_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, lung_seq)
    col_fun_3 = colorRamp2(c(-1,0,1), c("#2166ac","white","#b2182b"))
    
  lung_prot_seq <- Heatmap(mat_lung_seq,
                           col = col_fun_3, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                           cluster_rows = F,cluster_columns = T,
                           row_dend_width = unit(0.3, "cm"),
                           row_dend_gp = gpar(lwd = 0.2),
                           column_dend_height  = unit(0.3, "cm"),
                           column_dend_gp = gpar(lwd = 0.2),
                           width = unit(2.5, "cm"),
                           height = unit(2, "cm"),
                           row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                           column_names_rot=45,
                           column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                           #row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",4),rep("All-cause mortality",1)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality")),
                           # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                           #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                           column_title_gp = gpar(fontsize = 5),
                           column_title_side = "bottom",
                           # cell_fun = function(j, i, x, y, width, height, fill) {
                           #   
                           #   if((mat_prot_p)[i, j] < 0.001) {
                           #     grid.text("***", x, y,gp = gpar(fontsize = 4))
                           #   } else if((mat_prot_p)[i, j] < 0.01) {
                           #     grid.text("**", x, y,gp = gpar(fontsize = 4))
                           #   } else if((mat_prot_p)[i, j] < 0.05) {
                           #     grid.text("*", x, y,gp = gpar(fontsize = 4))
                           #   }
                           #   
                           #   
                           # },
                           # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                           #                                  show_legend = FALSE,
                           #                                  width= unit(1, "mm")
                           # ),
                           heatmap_legend_param = list(
                             at = c(-1, 0, 1),
                             #labels = c("0", "10%", "20%", "30%", "40%"),
                             title =  "Mean expression",
                             legend_height = unit(1, "cm"),
                             grid_width = unit(0.1, "mm"),
                             legend_width =unit(0.1, "mm"),
                             title_position = "lefttop-rot",
                             labels_gp = gpar(col = "black", fontsize = 4),
                             title_gp=gpar(col = "black", fontsize = 4)
                             
                           ),
                           top_annotation = columnAnnotation("Lung cancer"=lung_seq %>% filter(Cell.Type=="Type II pneumocyte"& Gene.Symbol!="CCL18") %>% dplyr::select(`Lung cancer`) %>% as.matrix(),
                                                             # `Axon development`=lung_seq %>% filter(Cell.Type=="stem cell") %>% dplyr::select(`Axon development`)%>% as.matrix(),
                                                             # "Lung aging"=rep(1,10),
                                                             col=list("Lung cancer"=colorRamp2(c(0,1), c("white","darkgreen"))
                                                                      # `Axon development`=colorRamp2(c(0,1), c("white","black")),
                                                                      # "Lung aging"=colorRamp2(c(0,1), c("white","darkgreen"))
                                                             ),
                                                             gp = gpar(col = "black",lwd=0.2),
                                                             border=F,
                                                             simple_anno_size = unit(1.5, "mm"),
                                                             show_legend = FALSE,
                                                             annotation_name_gp= gpar(fontsize = 5)
                           ),
                           
                           # right_annotation=rowAnnotation(
                           #   " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                           #                     gp = gpar(fill =col_panel_prot, 
                           #                               col = col_panel_prot,lwd=0.2),bar_width = 0.8,
                           #                     axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                           #   ),
                           #   
                           #   # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                           #   #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                           #   # 
                           #   width = unit(4, "cm"), annotation_name_gp= gpar(fontsize = 6)
                           #   # annotation_width=unit(c(1, 2), c("null", "null")),
                           #   # gap = unit(3, "mm")
                           #   
                           # ),
                           # bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                           #                                     col= list(" "=c(
                           #                                       "BCAN"="#1b9e77",    "C1QL2"="#d95f02",   "CA14"="#7570b3",    "CDH15"="#e7298a",   "CNTN1"="#66a61e",   
                           #                                       "CNTN2"="#e6ab02","COL28A1"="#a6761d", "FGFBP2"="#666666",  "FLT3"="#a6cee3",    "GFAP"="#1f78b4",    
                           #                                       "IGSF21"="#b2df8a",  "KLK6"="#33a02c", "MEPE"="#fb9a99",    "MOG"="#e31a1c",     "NEFL"="#fdbf6f",    
                           #                                       "OMG"="#ff7f00",     "PTPRR"="#cab2d6",   "RBFOX3"="#6a3d9a","SPOCK1"="#ffff99",  "SYT1"="#b15928")),
                           #                                     simple_anno_size = unit(2, "mm"),
                           #                                     show_legend = FALSE), 
                           
                           row_title =NULL
  )
  }
  lung_prot_seq

  #Changes with aging and lung CA in featured proteins of lung aging and cancer 
  library(Matrix)
  library(Seurat)
  lungca_seq <- readRDS(file = fig2_path("FIG2_LUNG_CANCER_RDS"))
    attributes(lungca_seq)
    # https://doi.org/10.1016/j.ccell.2022.10.008
    
    lung_d <-   {
      
      lungca_rna <- lungca_seq@assays[["RNA"]]@data
      lungca_summary <- lungca_seq@meta.data
      
      unique(lungca_summary$cell_type)
      
      lungca_gene_list <- lungca_seq@assays[["RNA"]]@meta.features
      
      lungca_gene_list <- lungca_gene_list %>% filter(feature_name %in% c(top_proteins_list$Lung))
      
      rename_info <- as.list(setNames(row.names(lungca_gene_list),lungca_gene_list$feature_name))
      
      lungca_summary_filtered <- lungca_summary %>%  
        filter(grepl(paste(c("multi-ciliated epithelial cell",
                             "epithelial cell of lung",
                             "club cell",
                             "type I pneumocyte",
                             "type II pneumocyte",
                             "macrophage",
                             "leukocyte"
        ), collapse = "|"), cell_type)) %>% filter(cell_type!="alveolar macrophage")
      
      lungca_summary_filtered$disease %>% unique()
      
      lungca_rna_filtered <- lungca_rna[row.names(lungca_gene_list), row.names(lungca_summary_filtered)] %>% t() %>% as.data.frame()
      
      lungca_filtered_all <- cbind(lungca_summary_filtered,lungca_rna_filtered) %>% 
        rename_with(~ as.character(lungca_gene_list$feature_name), as.character(row.names(lungca_gene_list))) %>% 
        group_by(donor_id,cell_type) %>% 
        mutate(across(as.character(lungca_gene_list$feature_name), ~ log(. + 1))) %>% 
        mutate(across(as.character(lungca_gene_list$feature_name), ~ (.x), .names = "{.col}_sum"))  %>% ungroup() %>% 
        group_by(donor_id) %>% dplyr::select(donor_id,
                                             cell_type,
                                             disease,
                                             ends_with("_sum"),
                                             age,
                                             sex,ever_smoker
                                             
        ) %>% group_by(donor_id,cell_type) %>% 
        mutate(disease_ca=case_when(disease=="normal"~"0",
                                    disease=="chronic obstructive pulmonary disease"~"COPD",
                                    T~"1"
        ),
        disease_copd=ifelse(disease_ca=="COPD","1","0")
        )
      
      # gene_list <- c("SFTPD","ALPP","SFTPA2","SCGB3A2","SFTPA1","SCGB1A1")
      gene_list <- top_proteins_list$Lung
      cell_list <- c("multi-ciliated epithelial cell",
                     "epithelial cell of lung",
                     "club cell",
                     "type I pneumocyte",
                     "type II pneumocyte",
                     "macrophage")
      
      
      run_by_cell_ca <- function(data, cell_list, genelist) {
        data <- data %>% filter(disease_ca!="COPD")
        results <- data.frame()
        for (cell in cell_list) {
          #data <- brain_filtered_all
          if (cell=="all") {subset_data <- data 
          
          }else{subset_data <- data %>% filter(cell_type==cell)}
          
          for (gene in genelist) {
            x <- subset_data$disease_ca
            y <- subset_data[[paste0(gene, "_sum")]]
            
            # model <- glm(y ~ x, family = binomial(link = "logit"))
            model <- lm(y ~ x+age+sex+ever_smoker,data = subset_data) %>% summary()
            
            coef <- model$coefficients["x1","Estimate"]
            se <- model$coefficients["x1","Std. Error"]
            p_value <- model$coefficients["x1","Pr(>|t|)"]
            
            result_row <- data.frame(
              cell=cell,
              gene=gene,
              cell_gene = paste0(cell,"_",gene),
              coef = coef,
              se = se,
              p_value = p_value
            )
            results <- rbind(results, result_row)
          }
        }
        
        results$q <- p.adjust(results$p_value, method = "fdr")
        return(results)
      }
      
      run_by_cell_copd <- function(data, cell_list, genelist) {
        data <- data %>% filter(disease_ca!="1")
        results <- data.frame()
        for (cell in cell_list) {
          #data <- brain_filtered_all
          if (cell=="all") {subset_data <- data 
          
          }else{subset_data <- data %>% filter(cell_type==cell)}
          
          for (gene in genelist) {
            x <- subset_data$disease_copd
            y <- subset_data[[paste0(gene, "_sum")]]
            
            # model <- glm(y ~ x, family = binomial(link = "logit"))
            model <- lm(y ~ x+age+sex+ever_smoker,data = subset_data) %>% summary()
            
            coef <- model$coefficients["x1","Estimate"]
            se <- model$coefficients["x1","Std. Error"]
            p_value <- model$coefficients["x1","Pr(>|t|)"]
            
            result_row <- data.frame(
              cell=cell,
              gene=gene,
              cell_gene = paste0(cell,"_",gene),
              coef = coef,
              se = se,
              p_value = p_value
            )
            results <- rbind(results, result_row)
          }
        }
        
        results$q <- p.adjust(results$p_value, method = "fdr")
        return(results)
      }
      
      
      lungca_coef_list <- run_by_cell_ca(lungca_filtered_all,cell_list, gene_list)
      copd_coef_list <- run_by_cell_copd(lungca_filtered_all,cell_list, gene_list)
      
      result_list_lung_lungca<-cox(data=hallmark_base_prot %>% filter(ca_base!=1),
                                   exposures=as.character(top_proteins_list$Lung),
                                   outcomes=c("Lung Cancer_reg_death"),
                                   covar=c("age","sex","ethnic","tdi","center","smok","pa"))
      
      lungca_change <- result_list_lung_lungca$result_list %>% 
        dplyr::select(Exposure,Outcome,Beta,FDR_P_Value) %>% 
        mutate(Outcome=sub("_reg_death$", "", Outcome)) %>% 
        rbind(.,linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                       a=as.character(top_proteins_list$Lung),
                       b="age",covar = c("age","sex","ethnic","edu","tdi","center","pa","smok")) %>% 
                mutate(q=p.adjust(p, method = "fdr")) %>% dplyr::select(var_a,var_b,beta,q) %>% 
                dplyr::rename("Exposure"=var_a,"Outcome"=var_b,"Beta"=beta,"FDR_P_Value"=q)
        ) %>% 
        rbind(.,lungca_coef_list %>%dplyr::select(cell,gene,coef,q) %>% 
                dplyr::rename("Exposure"=gene,"Outcome"=cell,"Beta"=coef,"FDR_P_Value"=q)
        ) %>% mutate(Outcome=factor(Outcome,levels=c("age","Lung Cancer",
                                                     "epithelial cell of lung","multi-ciliated epithelial cell",
                                                     "type I pneumocyte","type II pneumocyte","club cell","macrophage"
        )))
      
      
      mat_lungca_change <- xtabs(Beta ~ Exposure+Outcome, lungca_change) %>% t() 
      mat_lungca_change_p <- xtabs(FDR_P_Value ~Exposure+Outcome, lungca_change) %>% t() 
      row.names(mat_lungca_change) <- firstup(row.names(mat_lungca_change))
      
      col_fun_4 = colorRamp2(c(-0.2,0,0.2), c("#2166ac","white","#b2182b"))
      lungca <- Heatmap(mat_lungca_change,
                        col = col_fun_4, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                        cluster_rows = F,cluster_columns = F,
                        row_dend_width = unit(0.3, "cm"),
                        row_dend_gp = gpar(lwd = 0.2),
                        column_dend_height  = unit(0.3, "cm"),
                        column_dend_gp = gpar(lwd = 0.2),
                        width = unit(2.5, "cm"),
                        height = unit(1.5, "cm"),
                        row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                        column_names_rot=45,
                        column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                        row_split=factor(c("Plasma","Plasma",rep("Brain",6)),levels = c("Plasma", "Brain")),
                        row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                        column_split =factor(c(0,1,1,1,0,1,1,1,1,1,0)),
                        column_title = NULL,
                        column_title_gp = gpar(fontsize = 6),
                        column_title_side = "top",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          
                          if((mat_lungca_change_p)[i, j] < 0.001) {
                            grid.text("***", x, y,gp = gpar(fontsize = 4))
                          } else if((mat_lungca_change_p)[i, j] < 0.01) {
                            grid.text("**", x, y,gp = gpar(fontsize = 4))
                          } else if((mat_lungca_change_p)[i, j] < 0.05) {
                            grid.text("*", x, y,gp = gpar(fontsize = 4))
                          }
                          
                          
                        },
                        # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                        #                                  show_legend = FALSE,
                        #                                  width= unit(1, "mm")
                        # ),
                        heatmap_legend_param = list(
                          at = c(-0.2, 0, 0.2),
                          #labels = c("0", "10%", "20%", "30%", "40%"),
                          title =  "beta",
                          legend_height = unit(1, "cm"),
                          grid_width = unit(0.1, "mm"),
                          #legend_width =unit(1, "mm"),
                          #title_position = "lefttop-rot",
                          labels_gp = gpar(col = "black", fontsize = 5),
                          title_gp=gpar(col = "black", fontsize = 5)
                          
                        ),
                        row_title ="Lung cancer"
      )
      
      ###Plasma protein and COPD/HC or aging
      
      result_list_lung_copd<-cox(data=hallmark_base_prot %>% filter(ca_base!=1),
                                 exposures=as.character(top_proteins_list$Lung),
                                 outcomes=c("COPD"), #"acd","ad","vd","pd","ms","scz","dep","anx","sd","death"
                                 covar=c("age","sex","ethnic","tdi","center","smok","pa"))
      
      copd_change <- result_list_lung_copd$result_list %>% 
        dplyr::select(Exposure,Outcome,Beta,FDR_P_Value) %>% 
        rbind(.,linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                       a=as.character(top_proteins_list$Lung),
                       b="age",covar = c("age","sex","ethnic","edu","tdi","center","pa","smok")) %>% 
                mutate(q=p.adjust(p, method = "fdr")) %>% dplyr::select(var_a,var_b,beta,q) %>% 
                dplyr::rename("Exposure"=var_a,"Outcome"=var_b,"Beta"=beta,"FDR_P_Value"=q)
        ) %>% 
        rbind(.,copd_coef_list %>%dplyr::select(cell,gene,coef,q) %>% 
                dplyr::rename("Exposure"=gene,"Outcome"=cell,"Beta"=coef,"FDR_P_Value"=q)
        ) %>% mutate(Outcome=factor(Outcome,levels=c("age","COPD",
                                                     "epithelial cell of lung","multi-ciliated epithelial cell",
                                                     "type I pneumocyte","type II pneumocyte","club cell","macrophage"
        )))
      
      
      mat_copd_change <- xtabs(Beta ~ Exposure+Outcome, copd_change) %>% t() 
      mat_copd_change_p <- xtabs(FDR_P_Value ~Exposure+Outcome, copd_change) %>% t() 
      
      mat_copd_change <- mat_copd_change[-1,]
      mat_copd_change_p <- mat_copd_change_p[-1,]
      row.names(mat_copd_change) <- firstup(row.names(mat_copd_change))
      
      col_fun_4 = colorRamp2(c(-0.2,0,0.2), c("#2166ac","white","#b2182b"))
      copd <- Heatmap(mat_copd_change,show_heatmap_legend = F,
                      col = col_fun_4, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                      cluster_rows = F,cluster_columns = F,
                      row_dend_width = unit(0.3, "cm"),
                      row_dend_gp = gpar(lwd = 0.2),
                      column_dend_height  = unit(0.3, "cm"),
                      column_dend_gp = gpar(lwd = 0.2),
                      width = unit(2.5, "cm"),
                      height = unit(1.5, "cm"),
                      row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                      column_names_rot=45,
                      column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                      row_split=factor(c("Plasma",rep("Brain",6)),levels = c("Plasma", "Brain")),
                      row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                      column_split =factor(c(0,1,1,1,0,1,1,1,1,1,0)),
                      column_title = "COPD",
                      column_title_gp = gpar(fontsize = 6),
                      column_title_side = "top",
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        
                        if((mat_copd_change_p)[i, j] < 0.001) {
                          grid.text("***", x, y,gp = gpar(fontsize = 4))
                        } else if((mat_copd_change_p)[i, j] < 0.01) {
                          grid.text("**", x, y,gp = gpar(fontsize = 4))
                        } else if((mat_copd_change_p)[i, j] < 0.05) {
                          grid.text("*", x, y,gp = gpar(fontsize = 4))
                        }
                        
                        
                      },
                      # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                      #                                  show_legend = FALSE,
                      #                                  width= unit(1, "mm")
                      # ),
                      heatmap_legend_param = list(
                        at = c(-0.2, 0, 0.2),
                        #labels = c("0", "10%", "20%", "30%", "40%"),
                        title =  "beta",
                        legend_height = unit(1, "cm"),
                        grid_width = unit(0.1, "mm"),
                        #legend_width =unit(1, "mm"),
                        #title_position = "lefttop-rot",
                        labels_gp = gpar(col = "black", fontsize = 5),
                        title_gp=gpar(col = "black", fontsize = 5)
                        
                      ),
                      row_title ="COPD"
      )
      
      
                  plot_grid(grid.grabExpr(draw(lung_prot_seq,heatmap_legend_side = "left")),
                         grid.grabExpr(draw(lungca %v% copd,ht_gap = unit(1, "mm"),auto_adjust = T))
      )
      
    }
    
    rm(lungca_seq)
    
    f2_lung <- plot_grid(lung_b,ncol=1) %>% 
      plot_grid(.,plot_grid(plot_grid(lung_c,lung_a,nrow=1,rel_widths = c(1.2,1),labels=c("b","c"),label_size = 8),
                            lung_d,ncol=1,rel_heights = c(1,1.5),labels = c("","d"),label_size = 8),
                ncol=2,nrow=1,rel_widths = c(1,1.5),labels = c("a",""),label_size = 8)
    
    #Liver cancer----------
    #Individual disease and liver cancer
    result_hallmark_disease_liver <- cox(
      data=hallmark_base %>% filter(ca_base==0) %>% left_join(hallmark_all_final %>% dplyr::select(eid,date_cols_hesin %>% setdiff("Brain Cancer_date"))) %>% 
        mutate(`Liver Cancer_reg_death_date`=pmax(`Liver Cancer_date`,`Liver Cancer_reg_death_date`))
        ,
      exposures=c(
        paste0(hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% 
                 filter(if_any(c(GI, TA, EA, DNS, MD, CS, SCE, AIC), ~ . == 1)) %>% distinct(Disease) %>% pull(),"_base22")
      ),
      outcomes=c(
        "Liver Cancer_reg_death"
      ), 
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi")
    )
    
    aging_liver <- result_hallmark_disease_liver$result_list %>% filter(FDR_P_Value<0.05)
     
    liver_a <- ggplot(aging_liver %>% mutate(Exposure=gsub("_base22","",Exposure)) %>% 
                       dplyr::arrange(desc(-FDR_P_Value)) %>%
                       dplyr::slice(1:10),
                     aes(x = HR, y=(reorder(Exposure, HR)))) +
      geom_vline(xintercept = 1, color = "grey",linewidth=0.2) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
      geom_point(size = 1,shape=15,position=position_dodge(0.7))+
      #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
      scale_y_discrete()+
      labs(x = "HR (95% CI)", y = NULL, title = NULL) +
      theme_classic(base_size = 5,base_line_size = 0.2)+
      theme(            legend.position = "top",
                        legend.title =element_blank(),
                        legend.key.width = unit(0.1, "cm"),
                        axis.text.x = element_text(color="black",size = 5),
                        axis.text.y = element_text(color="black",size = 5))
    
    
    select <- dplyr::select
    #Individual age-realted proteins with liver cancer
    prot_in_liver <- top_proteins_list[c("Artery", "Organismal", "Immune", "Liver", "Stomach", "Kidney")] %>% unlist() %>%  as.character() %>% unique()
    
    hallmark_base_prot <- hallmark_base %>%
      left_join(hallmark_all_final %>% dplyr::select(eid,date_cols_hesin %>% setdiff("Brain Cancer_date"))) %>% 
      mutate(`Liver Cancer_reg_death_date`=pmax(`Liver Cancer_date`,`Liver Cancer_reg_death_date`)) %>% 
      left_join(hallmark_all_final %>% select(c("eid",c(organ_pac,"met_res","phenoage_res","kdm_res")))) %>%  
      left_join(olink_full %>% dplyr::select(all_of((prot_in_liver)), "eid"), by = "eid") %>%
      mutate(across(all_of((prot_in_liver)), ~ scale(as.numeric(.))))
    
    prot_liver_result <- cox(data=hallmark_base_prot %>% filter(!is.na(Organismal_res)) %>% filter(ca_base!=1) ,
                            exposures=prot_in_liver,
                            outcomes=c("Liver Cancer_reg_death"),
                            covar=c("age","sex","ethnic","tdi","center","pa","smok","bmi"))
    
    
    add_age_results <-linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                             b=prot_in_liver,
                             a="age",covar = c("sex","ethnic","edu","tdi","center","pa","smok","bmi"))
    
    
    prot_liver <- prot_liver_result$result_list %>% 
      mutate(Cat=case_when(Exposure %in% top_proteins_list$Artery~"Artery",
                           Exposure %in% top_proteins_list$Liver~"Liver",
                           Exposure %in% top_proteins_list$Immune~"Immune",
                           Exposure %in% top_proteins_list$Stomach~"Stomach",
                           Exposure %in% top_proteins_list$Kidney~"Kidney",
                           Exposure %in% top_proteins_list$Organismal~"Organismal"
      )) %>% 
      left_join(
        add_age_results %>% dplyr::select(var_b, beta),
        by = c("Exposure" = "var_b")
      ) %>%
      mutate(
        aging = ifelse(Beta / beta > 0, "Synergistic", "Antagonistic"),
        aging = factor(aging)
      )
    
    liver_b <- ggplot(prot_liver %>% filter(Cat!="Stomach"&Cat!="Kidney") %>% filter(Exposure!="CDHR2"), aes(x = HR, y = -log10(FDR_P_Value), color = Cat,shape = aging)) +
      geom_point(size = 1) +
      theme_bw() +
      ggtitle("Liver cancer")+
      labs(x = "HR", y = expression(-log[10](italic(q)~value)), color = NULL) +
      theme(legend.position = "right") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black",lwd=0.2) + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",lwd=0.2)+
      geom_text_repel(data = prot_liver %>% filter(Cat!="Stomach"&Cat!="Kidney") %>% filter(Exposure!="CDHR2") %>%  filter(FDR_P_Value<0.05) %>% arrange(FDR_P_Value) %>% head(15), aes(label = Exposure,color=Cat), size =2 ,segment.size = 0.2, max.overlaps = 15, show.legend = FALSE)+
      theme_classic(base_size = 6,base_line_size = 0.2)+
      theme(            legend.position = c(.01, .9), legend.justification = c("left", "top"),
                        plot.title = element_text(size = 7, face = "bold"),
                        legend.title =element_blank(),
                        legend.key.width = unit(0.1, "cm"),
                        axis.text.x = element_text(color="black",size = 5,angle = 0,vjust = 0.5),
                        axis.text.y = element_text(color="black",size = 6),
                        legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                        legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                        legend.key.spacing.y = unit(0, "pt"),
                        legend.box.spacing = margin(0),
                        # legend.margin = margin(t = 1, b =1, unit = "pt"),
                        # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                        strip.background = element_blank(),
                        # panel.spacing = unit(0.01, "cm")
                        
      )+scale_color_manual(values=c("Organismal"="black","Artery"="#A52A2A","Liver"="#D2B48C","Immune"="#feb24c"
                                    # "Kidney"="#808000","Stomach"="#008000"
                                      ))+
      scale_shape_manual(values = c("Synergistic" = 17, "Antagonistic" = 16))  
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    go_enrichment <- enrichGO(gene = c(prot_liver %>% filter(Cat %in% c("Artery","Organismal","Immune","Liver")) %>% filter(FDR_P_Value<0.05) %>%  pull(Exposure),
                                       c("SCARF2","EDA2R","FBLN2","CDCP1")), 
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",   
                              ont = "BP",
                              pvalueCutoff = 0.05
    )
    
    enrichment_bp_liver <- as.data.frame(go_enrichment)
    dotplot(go_enrichment, showCategory =5)+theme_bw()
    barplot(go_enrichment, showCategory =5)+theme_bw()
    
    
    # go_enrichment@result$qvalue <- round(go_enrichment@result$qvalue, 2)
    
    liver_c <- {
      res2 <- go_enrichment@result %>%
        as.data.frame() %>%
        filter(Description %in% c("blood coagulation","regulation of hemostasis","regulation of fibrinolysis",
                                  "regulation of body fluid levels",
                                  "extracellular matrix organization"
                                  )) %>% 
        dplyr::arrange(qvalue) %>%
        dplyr::slice(1:5) %>%
        dplyr::mutate(
          GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        ) %>%
        dplyr::arrange(GeneRatio_num) %>%
        dplyr::mutate(
          Description = paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description))),
          Description = stringr::str_wrap(Description, width = 28),
          Description = factor(Description, levels = Description)
        )
      
      q_breaks <- signif(seq(min(res2$qvalue), max(res2$qvalue), length.out = 3), 2)
      
      ggplot(res2, aes(x = GeneRatio_num, y = Description)) +
        geom_point(aes(size = Count, color = qvalue)) +
        scale_size_continuous(range = c(0.5, 1.5), breaks = c(2, 4, 6),name = "Gene count") +
        scale_color_gradient(
          low = "#023858", high = "#d0d1e6",expression(italic(q)~value),
          breaks = q_breaks,
          labels = format(q_breaks, scientific = TRUE, digits = 1)
        ) +
        labs(x = "Gene ratio", y = NULL) +
        guides(
          color = guide_colorbar(
            direction = "horizontal",
            barwidth = unit(0.75, "cm"),
            barheight = unit(0.06, "cm"),
            ticks = FALSE
          ),
          size = guide_legend(
            direction = "horizontal",
            override.aes = list(color = "black"),
            keywidth = unit(0.1, "cm"),
            keyheight = unit(0.1, "cm")
          )
        ) +
        theme_bw(base_size = 5, base_line_size = 0.1) +
        theme(
          panel.grid = element_blank(),
          legend.position = c(0.72, 0.15),
          legend.direction = "horizontal",
          legend.title = element_text(size = 3.5),
          legend.text = element_text(size = 3.5),
          axis.text.x = element_text(size = 5, color = "black"),
          axis.text.y = element_text(size = 5, color = "black"),
          legend.spacing.x = unit(0.01, "cm"),
          legend.spacing.y = unit(0.1, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(0.08, "cm"),
          legend.key.height = unit(0.08, "cm"),
          legend.background = element_rect(fill = "white", colour = NA)
        )
    }
    
    liver_c
    
    
    #Protein seq in Vasculature and liver
    vas_seq <-  read.csv(file = fig2_path("FIG2_VASCULAR_RNA_CSV"),skip=9) %>%
      filter(Gene.Symbol!="FRZB") %>% 
      filter(Cell.Type %in% c(
        "vascular associated smooth muscle cell",
        "mononuclear cell",
        "endothelial cell of artery",
        "macrophage",
        "lymphocyte",
        "fibroblast"
      )) %>%
      mutate(Cell.Type=case_when(
        # Cell.Type=="aggregated"~"Brain",
        Cell.Type=="vascular associated smooth muscle cell"~"SMC",
        Cell.Type=="endothelial cell of artery"~"Endothelial cell",
        T~ firstup(Cell.Type))) %>% 
      # mutate(Cell.Type=factor(Cell.Type,levels=c("T cell","NK cell","B cell","Myeloid cell","Progenitor cell"))) %>%
      mutate(Tissue.Composition=gsub("%", "", Tissue.Composition)) %>% mutate(Tissue.Composition=as.numeric(Tissue.Composition) / 100) %>% 
      mutate(Tissue.Composition=ifelse(Tissue.Composition==0,0.001,Tissue.Composition)) %>%
      mutate(expression_all=ifelse(!is.na(Tissue.Composition)&!is.na(Expression),Expression*(Tissue.Composition),NA)) %>% 
      group_by(Gene.Symbol) %>% mutate(expression_nom=scale(expression_all)) %>% ungroup() %>% 
      # mutate(`Synaptic signaling`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[1]] %>% unlist(),1,0),
      #        `Axon development`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[2]] %>% unlist(),1,0)) %>%
      mutate(Aging=ifelse(Gene.Symbol %in% (prot_liver %>% filter(Cat=="Artery") %>% dplyr::select(Exposure) %>% pull),"Artery","Organismal")) %>%
      mutate(`Liver cancer`=ifelse(Gene.Symbol %in% (prot_liver %>% filter(P_Value<0.05) %>% dplyr::select(Exposure) %>% pull),"Yes","No")) %>%
      arrange(Cell.Type,Gene.Symbol)

    mat_vas_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, vas_seq)
    # mat_liver_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, liver_seq)
    
    
    
    col_fun_3 = colorRamp2(c(-1,0,1), c("#2166ac","white","#b2182b"))
    ht_opt$TITLE_PADDING = unit(c(2, 2), "points")
    
    vas_seq_heat <-  Heatmap(mat_vas_seq,
                             col = col_fun_3, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                             cluster_rows = F,cluster_columns = T,
                             row_dend_width = unit(0.3, "cm"),
                             row_dend_gp = gpar(lwd = 0.2),
                             column_dend_height  = unit(0.3, "cm"),
                             column_dend_gp = gpar(lwd = 0.2),
                             # width = unit(4.5, "cm"),
                             # height = unit(2, "cm"),
                             width = ncol(mat_vas_seq)*unit(2, "mm"), 
                             height =nrow(mat_vas_seq)*unit(1.8, "mm"),
                             row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                             column_names_rot=45,
                             column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                             column_split=vas_seq %>% filter(Cell.Type=="SMC") %>% mutate(Aging=factor(Aging,levels=c("Artery","Organismal"))) %>% dplyr::select(Aging) %>% as.matrix(),
                             # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                             #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                             column_title_gp = gpar(fontsize = 4.5),
                             column_title_side = "top",
                             # cell_fun = function(j, i, x, y, width, height, fill) {
                             #   
                             #   if((mat_prot_p)[i, j] < 0.001) {
                             #     grid.text("***", x, y,gp = gpar(fontsize = 4))
                             #   } else if((mat_prot_p)[i, j] < 0.01) {
                             #     grid.text("**", x, y,gp = gpar(fontsize = 4))
                             #   } else if((mat_prot_p)[i, j] < 0.05) {
                             #     grid.text("*", x, y,gp = gpar(fontsize = 4))
                             #   }
                             #   
                             #   
                             # },
                             # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                             #                                  show_legend = FALSE,
                             #                                  width= unit(1, "mm")
                             # ),
                             heatmap_legend_param = list(
                               at = c(-1, 0, 1),
                               #labels = c("0", "10%", "20%", "30%", "40%"),
                               title =  "Mean expression",
                               legend_height = unit(1, "cm"),
                               grid_width = unit(0.1, "mm"),
                               legend_width =unit(0.1, "mm"),
                               title_position = "lefttop-rot",
                               labels_gp = gpar(col = "black", fontsize = 5),
                               title_gp=gpar(col = "black", fontsize = 5)
                               
                             ),
                             top_annotation = columnAnnotation(" "=vas_seq %>% filter(Cell.Type=="SMC") %>% dplyr::select(Aging) %>% as.matrix(),
                                                               "  "=vas_seq %>% filter(Cell.Type=="SMC") %>% dplyr::select(`Liver cancer`) %>% as.matrix(),
                                                               # `Axon development`=lung_seq %>% filter(Cell.Type=="stem cell") %>% dplyr::select(`Axon development`)%>% as.matrix(),
                                                               # "Immune aging"=rep(1,20),
                                                               col=list(" " = c("Artery" = "#A52A2A", "Organismal" = "black"),
                                                                        "  " = c("Yes" = "darkgreen", "No" = "white")
                                                                        ),
                                                               gp = gpar(col = "black",lwd=0.2),
                                                               border=F,
                                                               simple_anno_size = unit(1.5, "mm"),
                                                               show_legend = FALSE,
                                                               annotation_name_gp= gpar(fontsize = 5)
                             ),
                             
                             # right_annotation=rowAnnotation(
                             #   " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                             #                     gp = gpar(fill =col_panel_prot, 
                             #                               col = col_panel_prot,lwd=0.2),bar_width = 0.8,
                             #                     axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                             #   ),
                             #   
                             #   # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                             #   #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                             #   # 
                             #   width = unit(4, "cm"), annotation_name_gp= gpar(fontsize = 6)
                             #   # annotation_width=unit(c(1, 2), c("null", "null")),
                             #   # gap = unit(3, "mm")
                             #   
                             # ),
                             # bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                             #                                     col= list(" "=c(
                             #                                       "BCAN"="#1b9e77",    "C1QL2"="#d95f02",   "CA14"="#7570b3",    "CDH15"="#e7298a",   "CNTN1"="#66a61e",   
                             #                                       "CNTN2"="#e6ab02","COL28A1"="#a6761d", "FGFBP2"="#666666",  "FLT3"="#a6cee3",    "GFAP"="#1f78b4",    
                             #                                       "IGSF21"="#b2df8a",  "KLK6"="#33a02c", "MEPE"="#fb9a99",    "MOG"="#e31a1c",     "NEFL"="#fdbf6f",    
                             #                                       "OMG"="#ff7f00",     "PTPRR"="#cab2d6",   "RBFOX3"="#6a3d9a","SPOCK1"="#ffff99",  "SYT1"="#b15928")),
                             #                                     simple_anno_size = unit(2, "mm"),
                             #                                     show_legend = FALSE), 
                             row_title_gp = gpar(fontsize = 4.5),
                             row_title ="scRNA\nexpression",
                             # column_title =NULL
    )
    
    #Proteins interaction based on STING database
    liver_d <- {
      library(STRINGdb)
      library(igraph)
    #Conduct interaction analyses on: https://string-db.org/
    raw_data <- fread(input = fig2_path("FIG2_VASCULAR_INTERACTIONS_TSV"))
    # Check the structure of the raw_data data frame
    str(raw_data)
    
    unique(raw_data$node1)
    
    # Create an igraph object from the data frame
    g <- graph_from_data_frame(d = raw_data[, c("node1", "node2")], directed = FALSE)
    E(g)$width <- raw_data$combined_score *1.1
    
    
    node_color_mapping <- list(
      "ELN" = "#A52A2A",
      "LTBP2" = "#A52A2A",
      "CRLF1"= "#A52A2A",
      "COMP" = "#A52A2A",
      "TNFRSF11B" = "#A52A2A",
      "FRZB" = "#A52A2A",
      "ITGBL1" = "#A52A2A",
      "EFEMP1" = "#A52A2A",
      
      "FBLN2" = "black",
      "GDF15" = "black",
      "CCDC80"= "black"
    )
    
    # Create a vector of colors for nodes based on the mapping
    node_colors <- sapply(V(g)$name, function(node) {
      if (node %in% names(node_color_mapping)) {
        return(node_color_mapping[[node]])
      } else {
        return("gray")  # Assign gray color if node is not found in the mapping
      }
    })
    
    
    vertex_label_colors <- rep("black", length(V(g)))  # Set all labels to dark red
    
    # Plot the network
    
    par(mar = c(0, 0, 0, 0) + 0.1)
    vas_interact <- plot(g,
                         vertex.label = V(g)$name,  # Remove automatic labeling
                         vertex.color = node_colors,
                         vertex.label.color = vertex_label_colors,
                         vertex.label.cex = 0.5, 
                         vertex.label.dist = 1.5,
                         vertex.size = 8,
                         edge.width = E(g)$width,
                         vertex.label.family = "sans"
    )
    
    c(raw_data$node1,raw_data$node2) %>% unique()
    
    
    library(ggplotify)
    par(mar = c(0, 0, 0, 0) + 0.1)
    plot_grid(grid.grabExpr(draw(vas_seq_heat)),
                         
                         base2grob(~plot(g,
                                         vertex.label = V(g)$name,  # Remove automatic labeling
                                         vertex.color = node_colors,
                                         vertex.label.color = vertex_label_colors,
                                         vertex.label.cex = 0.4, 
                                         vertex.label.dist = 2,
                                         vertex.size = 3,
                                         edge.width = E(g)$width,
                                         vertex.label.family = "sans")),
                         
                         labels = c("","I"),label_size = 6
    )
    }
    
    f2_liver <- plot_grid(liver_b,ncol=1) %>% 
      plot_grid(.,plot_grid(plot_grid(liver_c,liver_a,nrow=1,rel_widths = c(1.2,1),labels=c("f","g"),label_size = 8),
                            liver_d,ncol=1,rel_heights = c(1,1.5),labels = c("","h"),label_size = 8),
                ncol=2,nrow=1,rel_widths = c(1,1.5),labels = c("e",""),label_size = 8)
    
    #NHL cancer----------
    #Individual disease and NHL
    result_hallmark_disease_NHL <- cox(
      data=hallmark_base %>% filter(ca_base==0),
      exposures=c(
        paste0(hesin_all_clean %>% filter(Category!="Cancers") %>% distinct(Disease,.keep_all = T) %>% 
                 filter(if_any(c(GI,CS), ~ . == 1)) %>% distinct(Disease) %>% pull(),"_base22")
      ),
      outcomes=c("NHL_reg_death"), 
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok","bmi")
    )
    
    aging_NHL<- result_hallmark_disease_NHL$result_list %>% filter(FDR_P_Value<0.05)
    
    
    nhl_a <- ggplot(aging_NHL %>% mutate(Exposure=gsub("_base22","",Exposure)) %>% 
                       dplyr::arrange(desc(-FDR_P_Value)) %>%
                       dplyr::slice(1:10),
                     aes(x = HR, y=(reorder(Exposure, HR)))) +
      geom_vline(xintercept = 1, color = "grey",linewidth=0.2) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
      geom_point(size = 1,shape=15,position=position_dodge(0.7))+
      #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
      scale_y_discrete()+
      labs(x = "HR (95% CI)", y = NULL, title = NULL) +
      theme_classic(base_size = 5,base_line_size = 0.2)+
      theme(            legend.position = "top",
                        legend.title =element_blank(),
                        legend.key.width = unit(0.1, "cm"),
                        axis.text.x = element_text(color="black",size = 5),
                        axis.text.y = element_text(color="black",size = 5))
    
    prot_in_NHL<- top_proteins_list[c("Muscle", "Organismal", "Kidney","Immune","Heart")] %>% unlist() %>%  as.character() %>% unique() 
    
    
    hallmark_base_prot <- hallmark_base %>%
      left_join(hallmark_all_final %>% select(c("eid",c(organ_pac,"met_res","phenoage_res","kdm_res")))) %>%  
      left_join(olink_full %>% dplyr::select(all_of((prot_in_NHL)), "eid"), by = "eid") %>%
      mutate(across(all_of((prot_in_NHL)), ~ scale(as.numeric(.))))
      
    
    prot_nhl_result <- cox(data=hallmark_base_prot %>% filter(!is.na(Organismal_res)) %>% filter(ca_base!=1) ,
                             exposures=prot_in_NHL,
                             outcomes=c("NHL_reg_death"),
                             covar=c("age","sex","ethnic","tdi","center","pa","smok","bmi"))
    
    
    add_age_results <-linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                             b=prot_in_NHL,
                             a="age",covar = c("sex","ethnic","edu","tdi","center","pa","smok","bmi"))
    
    prot_nhl<- prot_nhl_result$result_list %>% 
      mutate(Cat=case_when(Exposure %in% top_proteins_list$Muscle~"Muscle",
                           Exposure %in% top_proteins_list$Kidney~"Kidney",
                           Exposure %in% top_proteins_list$Immune~"Immune",
                           Exposure %in% top_proteins_list$Heart~"Heart",
                           Exposure %in% top_proteins_list$Organismal~"Organismal"
      )) %>% 
      left_join(
        add_age_results %>% dplyr::select(var_b, beta),
        by = c("Exposure" = "var_b")
      ) %>%
      mutate(
        aging = ifelse(Beta / beta > 0, "Synergistic", "Antagonistic"),
        aging = factor(aging)
      )

    nhl_b <- ggplot(prot_nhl, aes(x = HR, y = -log10(FDR_P_Value), color = Cat,shape = aging)) +
      geom_point(size = 1) +
      theme_bw() +
      ggtitle("NHL")+
      labs(x = "HR", y = expression(-log[10](italic(q)~value)), color = NULL) +
      theme(legend.position = "right") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black",lwd=0.2) + 
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",lwd=0.2)+
      geom_text_repel(data = prot_nhl  %>%  filter(FDR_P_Value<0.05) %>% arrange(FDR_P_Value) %>% head(15), aes(label = Exposure,color=Cat), size =2 ,segment.size = 0.2, max.overlaps = 15, show.legend = FALSE)+
      theme_classic(base_size = 6,base_line_size = 0.2)+
      theme(            legend.position = c(.01, .9), legend.justification = c("left", "top"),
                        plot.title = element_text(size = 7, face = "bold"),
                        legend.title =element_blank(),
                        legend.key.width = unit(0.1, "cm"),
                        axis.text.x = element_text(color="black",size = 5,angle = 0,vjust = 0.5),
                        axis.text.y = element_text(color="black",size = 6),
                        legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                        legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                        legend.key.spacing.y = unit(0, "pt"),
                        legend.box.spacing = margin(0),
                        # legend.margin = margin(t = 1, b =1, unit = "pt"),
                        # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                        strip.background = element_blank(),
                        # panel.spacing = unit(0.01, "cm")
                        
      )+scale_color_manual(values=c("Organismal"="black","Muscle"="#0000FF","Heart"="#FF7F50","Immune"="#feb24c",
                                      "Kidney"="#808000"
      ))+
      scale_shape_manual(values = c("Synergistic" = 17, "Antagonistic" = 16)) +
      scale_x_continuous(limits = c(0.75, NA))
    
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    go_enrichment <- enrichGO(gene = prot_nhl %>% filter(Cat=="Immune"&FDR_P_Value<0.05) %>% select(Exposure) %>% pull, 
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",   
                              ont = "BP",
                              pvalueCutoff = 0.05
    )
    dotplot(go_enrichment, showCategory =5)+theme_bw()
    barplot(go_enrichment, showCategory =5)+theme_bw()
    
    
    # go_enrichment@result$qvalue <- round(go_enrichment@result$qvalue, 2)
    
    nhl_c <- {
      res2 <- go_enrichment@result %>%
        as.data.frame() %>%
        mutate(
          Description = dplyr::recode(
            Description,
            "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" = "Adaptive immune receptor recombination"
          )
        ) %>% 
        dplyr::arrange(qvalue) %>%
        dplyr::slice(1:5) %>%
        dplyr::mutate(
          GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        ) %>%
        dplyr::arrange(GeneRatio_num) %>%
        dplyr::mutate(
          Description = paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description))),
          Description = stringr::str_wrap(Description, width = 28),
          Description = factor(Description, levels = Description)
        )
      
      q_breaks <- signif(seq(min(res2$qvalue), max(res2$qvalue), length.out = 3), 2)
      
      ggplot(res2, aes(x = GeneRatio_num, y = Description)) +
        geom_point(aes(size = Count, color = qvalue)) +
        scale_size_continuous(range = c(0.5, 1.5), breaks = c(2, 4, 6),name = "Gene count") +
        scale_color_gradient(
          low = "#023858", high = "#d0d1e6",expression(italic(q)~value),
          breaks = q_breaks,
          labels = format(q_breaks, scientific = TRUE, digits = 1)
        ) +
        labs(x = "Gene ratio", y = NULL) +
        guides(
          color = guide_colorbar(
            direction = "horizontal",
            barwidth = unit(0.75, "cm"),
            barheight = unit(0.06, "cm"),
            ticks = FALSE
          ),
          size = guide_legend(
            direction = "horizontal",
            override.aes = list(color = "black"),
            keywidth = unit(0.1, "cm"),
            keyheight = unit(0.1, "cm")
          )
        ) +
        theme_bw(base_size = 5, base_line_size = 0.1) +
        theme(
          panel.grid = element_blank(),
          legend.position = c(0.72, 0.15),
          legend.direction = "horizontal",
          legend.title = element_text(size = 3.5),
          legend.text = element_text(size = 3.5),
          axis.text.x = element_text(size = 5, color = "black"),
          axis.text.y = element_text(size = 5, color = "black"),
          legend.spacing.x = unit(0.01, "cm"),
          legend.spacing.y = unit(0.1, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(0.08, "cm"),
          legend.key.height = unit(0.08, "cm"),
          legend.background = element_rect(fill = "white", colour = NA)
        )
    }
    
    nhl_c
    

    #sc-RNAseq of NHL (annotated with immune aging and NHL)
    top_proteins_list$Immune 
    
    #Spleen RNA seq
    immune_seq <-  read.csv(file = fig2_path("FIG2_IMMUNE_RNA_CSV"),skip=9) %>%
      filter(Cell.Type %in% c(
        "T cell",
        "natural killer cell",
        "B cell",
        "myeloid cell",
        "progenitor cell"
      )) %>%
      mutate(Cell.Type=case_when(
        # Cell.Type=="aggregated"~"Brain",
        Cell.Type=="natural killer cell"~"NK cell",
        Cell.Type=="professional antigen presenting cell"~"Professional APCs",
        T~ firstup(Cell.Type))) %>% 
      mutate(Cell.Type=factor(Cell.Type,levels=c("T cell","NK cell","B cell","Myeloid cell","Progenitor cell"))) %>%
      mutate(Tissue.Composition=gsub("%", "", Tissue.Composition)) %>% mutate(Tissue.Composition=as.numeric(Tissue.Composition) / 100) %>% 
      mutate(Tissue.Composition=ifelse(Tissue.Composition==0,0.001,Tissue.Composition)) %>%
      mutate(expression_all=ifelse(!is.na(Tissue.Composition)&!is.na(Expression),Expression*(Tissue.Composition),NA)) %>% 
      group_by(Gene.Symbol) %>% mutate(expression_nom=scale(expression_all)) %>% ungroup() %>% 
      # mutate(`Synaptic signaling`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[1]] %>% unlist(),1,0),
      #        `Axon development`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[2]] %>% unlist(),1,0)) %>%
      mutate(NHL=ifelse(Gene.Symbol %in% (prot_nhl %>% filter(Cat=="Immune" &FDR_P_Value<0.05) %>% dplyr::select(Exposure) %>% pull),1,2)) %>%
      arrange(Cell.Type,Gene.Symbol)
    
    unique(immune_seq$Gene.Symbol)
    
    mat_immune_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, immune_seq)
    col_fun_3 = colorRamp2(c(-1,0,1), c("#2166ac","white","#b2182b"))
    
    immune_seq_ht <- Heatmap(mat_immune_seq,
                             col = col_fun_3, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                             cluster_rows = F,cluster_columns = T,
                             row_dend_width = unit(0.3, "cm"),
                             row_dend_gp = gpar(lwd = 0.2),
                             column_dend_height  = unit(0.3, "cm"),
                             column_dend_gp = gpar(lwd = 0.2),
                             # width = unit(4.5, "cm"),
                             # height = unit(2, "cm"),
                             width = ncol(mat_immune_seq)*unit(2, "mm"), 
                             height =nrow(mat_immune_seq)*unit(2, "mm"),
                             row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                             column_names_rot=45,
                             column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                             column_split=immune_seq %>% filter(Cell.Type=="B cell") %>% mutate(NHL=factor(NHL,levels=c(1,2))) %>% dplyr::select(NHL) %>% as.matrix(),
                             # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                             #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                             column_title_gp = gpar(fontsize = 6),
                             column_title_side = "bottom",
                             # cell_fun = function(j, i, x, y, width, height, fill) {
                             #   
                             #   if((mat_prot_p)[i, j] < 0.001) {
                             #     grid.text("***", x, y,gp = gpar(fontsize = 4))
                             #   } else if((mat_prot_p)[i, j] < 0.01) {
                             #     grid.text("**", x, y,gp = gpar(fontsize = 4))
                             #   } else if((mat_prot_p)[i, j] < 0.05) {
                             #     grid.text("*", x, y,gp = gpar(fontsize = 4))
                             #   }
                             #   
                             #   
                             # },
                             # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                             #                                  show_legend = FALSE,
                             #                                  width= unit(1, "mm")
                             # ),
                             heatmap_legend_param = list(
                               at = c(-1, 0, 1),
                               #labels = c("0", "10%", "20%", "30%", "40%"),
                               title =  "Mean expression",
                               legend_height = unit(1, "cm"),
                               grid_width = unit(0.1, "mm"),
                               legend_width =unit(0.1, "mm"),
                               title_position = "lefttop-rot",
                               labels_gp = gpar(col = "black", fontsize = 5),
                               title_gp=gpar(col = "black", fontsize = 5)
                               
                             ),
                             top_annotation = columnAnnotation("NHL"=immune_seq %>% filter(Cell.Type=="B cell") %>% dplyr::select(`NHL`) %>% as.matrix(),
                                                               # `Axon development`=lung_seq %>% filter(Cell.Type=="stem cell") %>% dplyr::select(`Axon development`)%>% as.matrix(),
                                                               # "Immune aging"=rep(1,20),
                                                               col=list("NHL"=colorRamp2(c(2,1), c("white","darkgreen"))
                                                                        # `Axon development`=colorRamp2(c(0,1), c("white","black")),
                                                                        # "Immune aging"=colorRamp2(c(0,1), c("white","darkgreen"))
                                                               ),
                                                               gp = gpar(col = "black",lwd=0.2),
                                                               border=F,
                                                               simple_anno_size = unit(1.5, "mm"),
                                                               show_legend = FALSE,
                                                               annotation_name_gp= gpar(fontsize = 5)
                             ),
                             row_title_gp = gpar(fontsize = 4.5),
                             row_title ="Spleen\nscRNA",column_title =NULL
    )
    immune_seq_ht
    
    
    #Changes with aging and NHL in featured proteins of immune aging
    library(Seurat)
    nhl_seq <- readRDS(file = fig2_path("FIG2_NHL_RDS"))
    attributes(nhl_seq)
    # https://doi.org/10.1016/j.xcrm.2023.101158
    
    nhl_d <-  {
    nhl_rna <- nhl_seq@assays[["RNA"]]@data
    nhl_summary <- nhl_seq@meta.data
    
    
    unique(nhl_summary$cell_type)
    nhl_gene_list <- nhl_seq@assays[["RNA"]]@meta.features
    
    # lungca_gene_list <- lungca_gene_list %>% filter(feature_name %in% c("SFTPD","ALPP","SFTPA2","SCGB3A2","SFTPA1","SCGB1A1"))
    nhl_gene_list <- nhl_gene_list %>% filter(feature_name %in% c(top_proteins_list$Immune))
    
    rename_info <- as.list(setNames(row.names(nhl_gene_list),nhl_gene_list$feature_name))
    
    
    nhl_summary_filtered <- nhl_summary 
    # %>%  
    #   filter(grepl(paste(c(
    #                        
    #                        "macrophage",
    #                        "naive B cell",
    #                        "memory B cell",
    #                        "precursor B cell",
    #                        "plasma cell",
    #                        "hematopoietic stem cell",
    #                        "natural killer cell",
    #                        "lymphocyte"
    #                        ), collapse = "|"), cell_type))
    
    
    
    nhl_rna_filtered <- nhl_rna[row.names(nhl_gene_list), row.names(nhl_summary_filtered)] %>% as.matrix() %>%  t() %>% as.data.frame()
    
    nhl_filtered_all <- cbind(nhl_summary_filtered,nhl_rna_filtered) %>% 
      rename_with(~ as.character(nhl_gene_list$feature_name), as.character(row.names(nhl_gene_list))) %>% 
      group_by(donor_id,cluster) %>% 
      mutate(across(as.character(nhl_gene_list$feature_name), ~ log(. + 1))) %>% 
      mutate(across(as.character(nhl_gene_list$feature_name), ~ (.x), .names = "{.col}_sum"))  %>% ungroup() %>% 
      group_by(donor_id) %>% dplyr::select(donor_id,
                                           cluster,
                                           cell_type,
                                           disease,
                                           ends_with("_sum"),
                                           `Age (years)`,
                                           sex,
                                           self_reported_ethnicity,
                                           analysis_group,
                                           therapy) %>%
      ungroup() %>% 
      mutate(disease_nhl=case_when(disease=="normal"~"0",T~"1"),
             age=as.numeric(`Age (years)`),
             age=ifelse(is.na(age),61.36,age)
      ) %>% 
      mutate(
        cluster=as.character(cluster),
        cluster2=case_when(     str_detect(cluster,"CD8 T")~"CD8 T cell",
                                str_detect(cluster,"CD4 T")~"CD4 T cell",
                                str_detect(cluster,"NK")~"NK cell",
                                str_detect(cluster, "C1[1-5]") ~ "B cell",
                                str_detect(cluster, "C1[6-9]|c20") ~ "Myeloid cell",
                                str_detect(cluster, "C2[3-5]") ~ "Stem and progenitor cell",
                                T~NA
        )) %>% filter(!is.na(cluster2))
    
    describe(nhl_filtered_all$cluster)
    
    gene_list <- top_proteins_list$Immune
    cell_list <-   unique(nhl_filtered_all$cluster2) %>% as.character()
    
    
    run_by_cell_ca <- function(data, cell_list, genelist) {
      data <- data
      results <- data.frame()
      for (cell in cell_list) {
        #data <- brain_filtered_all
        if (cell=="all") {subset_data <- data 
        
        }else{subset_data <- data %>% filter(cluster2==cell)}
        
        for (gene in genelist) {
          
          
          x <- subset_data$disease_nhl
          y <- subset_data[[paste0(gene, "_sum")]]
          
          # model <- glm(y ~ x, family = binomial(link = "logit"))
          model <- lm(y ~ x+age+sex+self_reported_ethnicity+analysis_group+therapy,data = subset_data) %>% summary()
          
          coef <- model$coefficients["x1","Estimate"]
          se <- model$coefficients["x1","Std. Error"]
          p_value <- model$coefficients["x1","Pr(>|t|)"]
          
          result_row <- data.frame(
            cell=cell,
            gene=gene,
            cell_gene = paste0(cell,"_",gene),
            coef = coef,
            se = se,
            p_value = p_value
          )
          results <- rbind(results, result_row)
        }
      }
      
      results$q <- p.adjust(results$p_value, method = "fdr")
      return(results)
    }
    
    ###Plasma protein and NHL/HC or aging
    nhl_coef_list <- run_by_cell_ca(nhl_filtered_all,cell_list, gene_list) %>% dplyr::select(cell,gene,coef,q)
    
    result_list_immune_nhl<-cox(data=hallmark_base_prot %>% filter(ca_base!=1),
                                exposures=as.character(top_proteins_list$Immune),
                                outcomes=c("NHL_reg_death"),
                                covar=c("age","sex","ethnic","edu","tdi","center","smok","pa"))
    
    nhl_change <- prot_nhl %>% 
      filter(Cat=="Immune") %>% 
      dplyr::mutate(Outcome = gsub("_reg_death", "", Outcome)) %>% 
      dplyr::select(Exposure,Outcome,Beta,FDR_P_Value) %>% 
      rbind(.,linear(hallmark_base_prot %>% mutate(age=scale(age) %>% as.numeric()),
                     b=as.character(top_proteins_list$Immune),
                     a="age",covar = c("sex","ethnic","edu","tdi","center","pa","smok",as.character(top_proteins_list$Immune))) %>% 
              mutate(q=p.adjust(p, method = "fdr")) %>% dplyr::select(var_a,var_b,beta,q) %>% 
              rename("Exposure"=var_b,"Outcome"=var_a,"Beta"=beta,"FDR_P_Value"=q)
      ) %>% 
      rbind(.,nhl_coef_list %>%dplyr::select(cell,gene,coef,q) %>% 
              rename("Exposure"=gene,"Outcome"=cell,"Beta"=coef,"FDR_P_Value"=q)
      )  %>% 
      mutate(Outcome=ifelse(Outcome=="age","Immune aging",Outcome),Outcome=firstup(Outcome)) %>% 
      arrange(Exposure) %>% mutate(NHL=ifelse(Outcome=="NHL"&FDR_P_Value<0.05,1,2) %>% factor(.,levels=c(2,1))) %>% 
      mutate(Outcome=factor(Outcome,levels=c(  "Immune aging",
                                               "NHL","CD8 T cell","CD4 T cell","NK cell","B cell","Myeloid cell" ,"Stem and progenitor cell"))) %>% 
      mutate(group="NHL vs normal")
    
    
    # immune_seq %>% select(Cell.Type,Gene.Symbol,)
    
    mat_nhl <- xtabs(Beta ~ Outcome+Exposure, nhl_change) 
    mat_nhl_q <- xtabs(FDR_P_Value ~ Outcome+Exposure, nhl_change) 
    
    col_fun_4 = colorRamp2(c(-0.25,0,0.25), c("#2166ac","white","#b2182b"))
    nhl_seq_ht <- Heatmap(mat_nhl,
                          col = col_fun_4, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                          cluster_rows = F,cluster_columns = F,
                          row_dend_width = unit(0.3, "cm"),
                          row_dend_gp = gpar(lwd = 0.2),
                          column_dend_height  = unit(0.3, "cm"),
                          column_dend_gp = gpar(lwd = 0.2),
                          # width = unit(4, "cm"),
                          # height = unit(2, "cm"),
                          width = ncol(mat_nhl)*unit(2, "mm"), 
                          height =nrow(mat_nhl)*unit(2, "mm"),
                          row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                          column_names_rot=45,
                          column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                          row_split=factor(c("Plasma\nprotein","Plasma\nprotein",rep("Bone marrow\nscRNA",6)),levels = c("Plasma\nprotein", "Bone marrow\nscRNA")),
                          row_title_gp = gpar(fontsize = 4.5),
                          column_split = nhl_change %>% filter(Outcome=="NHL") %>% dplyr::select(NHL),
                          column_title = NULL,
                          column_title_gp = gpar(fontsize = 6),
                          column_title_side = "bottom",
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            
                            if((mat_nhl_q)[i, j] < 0.001) {
                              grid.text("***", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_nhl_q)[i, j] < 0.01) {
                              grid.text("**", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_nhl_q)[i, j] < 0.05) {
                              grid.text("*", x, y,gp = gpar(fontsize = 4))
                            }
                            
                            
                          },
                          heatmap_legend_param = list(
                            at = c(-0.2, 0, 0.2),
                            #labels = c("0", "10%", "20%", "30%", "40%"),
                            title =  "beta",
                            legend_height = unit(1, "cm"),
                            grid_width = unit(0.1, "mm"),
                            #legend_width =unit(1, "mm"),
                            #title_position = "lefttop-rot",
                            labels_gp = gpar(col = "black", fontsize = 5),
                            title_gp=gpar(col = "black", fontsize = 5)
                            
                          ),
                          # row_title =T
    )

    grid.grabExpr(draw(immune_seq_ht %v% nhl_seq_ht,ht_gap = unit(0, "mm"),auto_adjust = T)) %>% plot_grid()
    
    }
    
    rm(nhl_seq)
    
    f2_nhl <- plot_grid(nhl_b,ncol=1) %>% 
      plot_grid(.,plot_grid(plot_grid(nhl_c,nhl_a,nrow=1,rel_widths = c(1.2,1),labels=c("k","l"),label_size = 8),
                            nhl_d,ncol=1,rel_heights = c(1,1.5),labels = c("","m"),label_size = 8),
                ncol=2,nrow=1,rel_widths = c(1,1.5),labels = c("j",""),label_size = 8)
    
    
    #Figure 3 output----
    cairo_pdf(fig2_path("FIG2_OUTPUT_PDF"), 
              width = 7.4, height = 8.2)
    
    plot_grid(f2_lung,f2_liver,f2_nhl,ncol=1)+
      theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
    
    dev.off()  
    
    ### Saving working image--------
    save.image(fig2_path("FIG2_OUTPUT_WORKSPACE"))
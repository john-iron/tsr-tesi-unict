#!/usr/bin/env Rscript
#library(ggplot2)
library(PharmacoGx)
library(ragg)
source("Functions/getScenarioInput_tumor_vs_normal.R")
source("Functions/getSensitivitySignature.R")
source("Functions/getSignificanceSameDirection.R")
source("Functions/get_Alternative_CMAP_score.R")
source("Functions/cmap_score_new.R")
source("Functions/get_Alternative_Connectivity_Score_random_distribution.R")

  ### Load data: ###
  scenarioInput <- getScenarioInput_tumor_vs_normal()
  TCGA_project_to_PRISM_tissue_final<-readRDS("Output/1.process_PRISM_data/TCGA_project_to_PRISM_tissue_enriched.Rds")
  processing_tumor<-TCGA_project_to_PRISM_tissue_final$TCGA_tumor
  
  scenarioInput <- subset(scenarioInput, scenarioInput$log2_FC_input==processing_tumor)
  log2_FC_input <- readRDS("Output/3.Tumor_DEG/combined_tumor_vs_normal_DEG.Rds")
  log2_FC_input <- subset(log2_FC_input, Tumor == processing_tumor)
  log2_FC_input$t.value <- log2_FC_input$log2_FC_estimate / log2_FC_input$log2_FC_estimate_SE

  # PRISM:
  drug_ranef_and_sd <- readRDS(file = "Output/1.process_PRISM_data/drug_ranef_and_sd.Rds")

  # LINCS:
  gene_meta <- readRDS("External_input/gene_meta.Rds")
  DEG_LINCS <- readRDS("Output/2.LINCS_preprocessing/PRISM_combined_DEG_drugs_combined_incl.corrected.Rds")
  
  DEG_LINCS <- merge(x = DEG_LINCS, y = gene_meta, by.x = "gene", by.y = "Gene name") 
  DEG_LINCS$gene <- DEG_LINCS$`Gene stable ID`
  DEG_LINCS <- DEG_LINCS[DEG_LINCS$gene %in% unique(log2_FC_input$Gene),]
  DEG_LINCS <- DEG_LINCS[DEG_LINCS$drug %in% drug_ranef_and_sd$drug_name,]
  length(unique(DEG_LINCS$drug)) # N = 507
  LINCS_genes <- unique(DEG_LINCS$gene) # 952 genes left
  
  drug_ranef_and_sd <- drug_ranef_and_sd[drug_ranef_and_sd$drug_name %in% unique(DEG_LINCS$drug),]
  log2_FC_input <- log2_FC_input[log2_FC_input$Gene %in% LINCS_genes,]
  #Creo la cartella di Output
 
  abs_path <- getwd()
  scenario_path <-'Output/4.AUC_versus_connectivity_score/A549_LUNG'
  newpath <- file.path(abs_path, scenario_path)
  dir.create(newpath)
  ### Start loop: ###
  #Gli scenari generati saranno sempre 16 per ogni Linea Cellulare
  for(i in 1:16){
    print(paste0("Scenario INPUT = ", i))
    drugs <- as.data.frame(table(drug_ranef_and_sd$drug_name), stringsAsFactors = F)
    colnames(drugs) <- c("drug_name","number_of_cellLines_tested")
    drugs$gene_selection <- scenarioInput[i,]$gene_selection
    drugs$DEG_condition <- scenarioInput[i,]$DEG_condition
    drugs$log2_FC_input <- scenarioInput[i,]$log2_FC_input
  
    log2_FC_input <- log2_FC_input[order(abs(log2_FC_input$t.value), decreasing = T),]
  
  if (scenarioInput[i,]$gene_selection == "50_most_significant") {
    
    tumorSignature <- log2_FC_input[1:50,]
    
  } else if (scenarioInput[i,]$gene_selection == "100_most_significant"){
    
    tumorSignature <- log2_FC_input[1:100,]
    
  } else if (scenarioInput[i,]$gene_selection == "150_most_significant"){
    
    tumorSignature <- log2_FC_input[1:150,]
    
  } else if (scenarioInput[i,]$gene_selection == "Bin_Chen_original"){
    
    tumorSignature <- 
      subset(log2_FC_input, abs(log2_FC_estimate) > 1.5 & adj.P.Val < 0.001, select = c("Gene","t.value"))
    
  }

  if (nrow(tumorSignature) >= 10){
    
    alt_cs_distribution <- 
      get_Alternative_Connectivity_Score_random_distribution(
        n_genes_up = sum(tumorSignature$t.value > 0), 
        n_genes_down = sum(tumorSignature$t.value < 0), 
        n_genes = length(unique(DEG_LINCS$gene)), # N = 952
        n_permutations = 10^5
      )
    
    rownames(tumorSignature) <- tumorSignature$Gene
    tumorSignature <- tumorSignature[,"t.value", drop=F]
    colnames(tumorSignature) <- "estimate"
    
    for (d in 1:nrow(drugs)){
      #507 farmaci
      #print(d)
      print(paste0("FARMACI  = ",drugs[d,"drug_name"] ))
      drugSignature <- subset(DEG_LINCS, drug == drugs[d,"drug_name"])
      
      if (nrow(drugSignature) > 0){
        
        drugSignature <- drugSignature[!duplicated(drugSignature$gene),]
        drugSignature$t.value <- drugSignature[,paste0("T.value_",scenarioInput[i,]$DEG_condition)]
        rownames(drugSignature) <- drugSignature$gene
        
        ### Prepare signatures for connectivity scoring: ###
        drugSignature <- drugSignature[,"t.value",drop=F]
        colnames(drugSignature) <- "estimate"
        
        cs <- PharmacoGx::connectivityScore(
          x = drugSignature,
          y = tumorSignature,
          method = "fgsea",
          nperm = 10^4
        )
        
        drugs[d,c("connectivity_score_with_LINCS_profile","connectivity_score_with_LINCS_profile.pvalue")] <- cs
        
        alt_cs <- 
          get_Alternative_CMAP_score(
            tumorSignature, 
            drugSignature, 
            alt_cs_distribution
          )
        
        drugs[d,c("alt_connectivity_score_with_LINCS_profile","alt_connectivity_score_with_LINCS_profile.pvalue")] <- alt_cs
        
      }
    }
    
    saveRDS(drugs, file = paste0(scenario_path,"/","scenario",i,".Rds"))
    
  }
  
}

net_drug_effect_accross_tumor_type <- readRDS("Output/1.process_PRISM_data/drug_ranef_and_sd.Rds")
net_drug_effect_by_tumor_type <- readRDS("Output/1.process_PRISM_data/net_drug_effect_by_tumor_type.Rds")
files <- list.files(scenario_path)
scenarios <- as.integer(substr(files, 9, nchar(files)-4)) # Missing PAAD with Bin Chen method because of no significant genes
scenarios <- scenarios[order(scenarios)]


for (i in scenarios){
  tmp <- readRDS(paste0(scenario_path,"/","scenario",i,".Rds"))
  scenarioInput <- getScenarioInput_tumor_vs_normal()
  scenarioInput <- subset(scenarioInput, scenarioInput$log2_FC_input==processing_tumor)
  scenarioInput$scenario_index <- i
  
  tmp$connectivity_score_with_LINCS_profile.fdr <- 
    p.adjust(tmp$connectivity_score_with_LINCS_profile.pvalue, method="BH")
  
  scenarioInput$fraction_positive_connectivity_score_below.P_5pct <-
    nrow(subset(tmp, connectivity_score_with_LINCS_profile.pvalue < 0.05 & connectivity_score_with_LINCS_profile > 0)) /
    nrow(tmp)
  
  scenarioInput$fraction_negative_connectivity_score_below.P_5pct <-
    nrow(subset(tmp, connectivity_score_with_LINCS_profile.pvalue < 0.05 & connectivity_score_with_LINCS_profile < 0)) /
    nrow(tmp)
  
  tmp <- merge(x = tmp, y = subset(net_drug_effect_by_tumor_type, TCGA_project == scenarioInput$log2_FC_input), by = "drug_name")
  
  #tmp <- merge(x = tmp, y = net_drug_effect_by_tumor_type)
  
  #predicted_log_AUC <- subset(net_drug_effect_by_tumor_type, TCGA_project == paste0("Everything_except_",scenarioInput$log2_FC_input))
  
  predicted_log_AUC<-net_drug_effect_by_tumor_type
  
  predicted_log_AUC <- predicted_log_AUC[,c("drug_name","net_drug_effect")]
  colnames(predicted_log_AUC)[2] <- "predicted_net_drug_effect"
  tmp <- merge(x = tmp, y = predicted_log_AUC, by = "drug_name")
  tmp$net_drug_effect_centered <- tmp$net_drug_effect - mean(tmp$net_drug_effect)
  tmp$predicted_net_drug_effect_centered <- tmp$predicted_net_drug_effect - mean(tmp$predicted_net_drug_effect)
  
  linear_model_null_model <- lm(net_drug_effect_centered ~ 1, data = tmp)
  linear_model_connectivity_scores <- lm(net_drug_effect_centered ~ 1 + alt_connectivity_score_with_LINCS_profile, data = tmp)
  linear_model_AUC_other_cellLines <- lm(net_drug_effect_centered ~ 1 + predicted_net_drug_effect_centered, data = tmp)
  linear_model_combined <- lm(net_drug_effect_centered ~ 1 + alt_connectivity_score_with_LINCS_profile + predicted_net_drug_effect_centered, data = tmp)
  
  coef_cs_model <- as.data.frame(coef(summary(linear_model_connectivity_scores)))
  coef_cs_model$coef_name <- rownames(coef_cs_model)
  
  if (coef_cs_model[coef_cs_model$coef_name == "alt_connectivity_score_with_LINCS_profile","Estimate"] < 0){
    
    P.value_connectivity_scores <- 1
    
  } else {
    
    P.value_connectivity_scores <- anova(linear_model_null_model, linear_model_connectivity_scores)$`Pr(>F)`[2]
    
  }
  
  coef_combined_model <- as.data.frame(coef(summary(linear_model_combined)))
  coef_combined_model$coef_name <- rownames(coef_combined_model)
  
  if (coef_combined_model[coef_combined_model$coef_name == "alt_connectivity_score_with_LINCS_profile","Estimate"] < 0){
    
    P.value_connectivity_scores_added_to_AUC_other_celllines <- 1
    
  } else {
    
    P.value_connectivity_scores_added_to_AUC_other_celllines <- anova(linear_model_AUC_other_cellLines, linear_model_combined)$`Pr(>F)`[2]
    
  }
  
  scenarioInput[,c(
    "R2_connectivity_score",
    "R2_AUC_other_celllines",
    "R2_benefit_adding_connectivity_score",
    "P.value_connectivity_scores", 
    "P.value_AUC_other_celllines", 
    "P.value_connectivity_scores_added_to_AUC_other_celllines"
  )] <-
    c(
      signif(100 * summary(linear_model_connectivity_scores)$r.squared, digits = 3),
      signif(100 * summary(linear_model_AUC_other_cellLines)$r.squared, digits = 3),
      signif(100 * (summary(linear_model_combined)$r.squared - summary(linear_model_AUC_other_cellLines)$r.squared), digits = 3),
      signif(P.value_connectivity_scores, digits = 3),
      signif(anova(linear_model_null_model, linear_model_AUC_other_cellLines)$`Pr(>F)`[2], digits = 3),
      signif(P.value_connectivity_scores_added_to_AUC_other_celllines, digits = 3)
    )
  
  spearman_result <- cor.test(x = tmp$alt_connectivity_score_with_LINCS_profile, y = tmp$net_drug_effect, method = "spearman")
  scenarioInput$Rho <- spearman_result$estimate
  scenarioInput$Rho.p.value <- spearman_result$p.value
  if(i==1){
  scenarioOutput<-scenarioInput
  }
  else{
    scenarioOutput<-rbind(scenarioOutput,scenarioInput)
  }
  
}

scenarioOutput$P_significant <- scenarioOutput$Rho.p.value < 0.05
scenarioOutput$Rho_positive <- scenarioOutput$Rho > 0

Figure_4_plotInput <- subset(scenarioOutput, DEG_condition %in% c("6h_24h_combined","6h_24h_combined_corrected"))

Figure_4_plotInput$`Number of genes in tumor signature` <- 
  factor(
    Figure_4_plotInput$gene_selection, 
    levels = c("50_most_significant","100_most_significant","150_most_significant","Bin_Chen_original"),
    labels = c("50 most statistically significant DE tumor genes","100 most statistically significant DE tumor genes","150 most statistically significant DE tumor genes",">1.5 logFC, adj. P < 0.001"),
    ordered = T
  )
Figure_4_plotInput$`Drug signatures used` <-
  factor(
    Figure_4_plotInput$DEG_condition,
    levels = c("6h_24h_combined","6h_24h_combined_corrected"),
    labels = c("Uncorrected","Corrected")
  )


### Figure 4: ###
Figure_4 <- ggplot(data = Figure_4_plotInput, aes(x = `Drug signatures used`, y = Rho)) + 
  geom_hline(yintercept = 0, color = "red") +
  geom_boxplot() + facet_grid(~ `Number of genes in tumor signature`) + labs(
    y = "Spearman correlation between connectivity score and mnAUC of cell lines belonging to the tissue"
  ) + scale_y_continuous(limits = c(-0.25,0.33))

ragg::agg_png("Publication/Figures/Figure_4.png", width = 12.80, height = 7.20, units = "in", res = 300)
print(Figure_4)
dev.off()



# P-values:
wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "50_most_significant" & `Drug signatures used` == "Uncorrected")$Rho)
wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "50_most_significant" & `Drug signatures used` == "Corrected")$Rho)

wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "100_most_significant" & `Drug signatures used` == "Uncorrected")$Rho)
wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "100_most_significant" & `Drug signatures used` == "Corrected")$Rho)

wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "150_most_significant" & `Drug signatures used` == "Uncorrected")$Rho)
wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "150_most_significant" & `Drug signatures used` == "Corrected")$Rho)

wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "Bin_Chen_original" & `Drug signatures used` == "Uncorrected")$Rho)
wilcox.test(x = subset(Figure_4_plotInput, gene_selection == "Bin_Chen_original" & `Drug signatures used` == "Corrected")$Rho)

# Tables:
Figure_4_plotInput$Label <- paste0(round(Figure_4_plotInput$Rho, digits = 2)," (",Figure_4_plotInput$Rho.p.value,")")
Table_1 <- subset(Figure_4_plotInput, `Drug signatures used` == "Uncorrected", select = c("log2_FC_input","gene_selection","Label"))
write.csv(Table_1, file = "Publication/Tables/Table_1.csv", row.names = F)

Table_2 <- subset(Figure_4_plotInput, `Drug signatures used` == "Corrected", select = c("log2_FC_input","gene_selection","Label"))
write.csv(Table_2, file = "Publication/Tables/Table_2.csv", row.names = F)

Table_3 <- subset(Figure_4_plotInput, `Drug signatures used` == "Uncorrected" & gene_selection == "150_most_significant")
Table_3 <- Table_3[,c("log2_FC_input","R2_connectivity_score","P.value_connectivity_scores","R2_AUC_other_celllines","P.value_AUC_other_celllines","R2_benefit_adding_connectivity_score","P.value_connectivity_scores_added_to_AUC_other_celllines")]
write.csv(Table_3, file = "Publication/Tables/Table_3.csv", row.names = F)
median(Table_3$R2_connectivity_score)
median(Table_3$R2_AUC_other_celllines)
median(Table_3$R2_benefit_adding_connectivity_score)

Table_S4 <- subset(scenarioOutput, DEG_condition == "6h" & gene_selection == "150_most_significant")
Table_S4 <- Table_S4[,c("log2_FC_input","R2_connectivity_score","P.value_connectivity_scores","R2_benefit_adding_connectivity_score","P.value_connectivity_scores_added_to_AUC_other_celllines")]
write.csv(Table_S4, file = "Publication/Tables/Table_S4.csv", row.names = F)
median(Table_S4$R2_connectivity_score)
median(Table_S4$R2_benefit_adding_connectivity_score)

Table_S5 <- subset(scenarioOutput, DEG_condition == "24h" & gene_selection == "150_most_significant")
Table_S5 <- Table_S5[,c("log2_FC_input","R2_connectivity_score","P.value_connectivity_scores","R2_benefit_adding_connectivity_score","P.value_connectivity_scores_added_to_AUC_other_celllines")]
write.csv(Table_S5, file = "Publication/Tables/Table_S5.csv", row.names = F)
median(Table_S5$R2_connectivity_score)
median(Table_S5$R2_benefit_adding_connectivity_score)




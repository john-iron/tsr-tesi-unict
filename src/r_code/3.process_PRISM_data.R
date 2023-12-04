#!/usr/bin/env Rscript
#library(ggplot2)
library(lme4)
library(egg)
library(sets)
source(snakemake@input[['fun_loadprism']])
source(snakemake@input[['fun_add_ass']])

#importo il dataset modificato
input <- 
  read.csv(
    snakemake@input[['input_file']],
    stringsAsFactors = F
  )

cell.line.info <- read.csv(
    snakemake@input[['cell_lineinfo']],
    stringsAsFactors = F
)



### Calculate average AUC over all cell lines: ###
allData_Prism_data_and_model <- load_Prism_data_and_model()

saveRDS(allData_Prism_data_and_model$drug_ranef, file = snakemake@output[['drug_ranef_and_sd']])
saveRDS(allData_Prism_data_and_model$number_of_cellLines_tested_per_drug, file = snakemake@output[['number_of_cellLines_tested_per_drug']])

### Calculate average AUC over by tumor type: ###
colnames(input)[colnames(input) == "name"] <- "drug_name"
colnames(input)[colnames(input) == "depmap_id"] <- "cellLine_id"


### Calculate average AUC over by tumor type: ###
summary(input)
sum(is.finite(input$ic50)) / nrow(input)
min(input$ec50)
max(input$ec50)


input<-subset(input, !is.na(cellLine_id))
input<-subset(input, !is.na(drug_name))
input<-merge(x=input, y=cell.line.info, by.x="cellLine_id", by.y="depmap_id")


#importo le signature dei tumori
TCGA_project_to_PRISM_tissue <- read.csv(snakemake@input[['TCGA_project_to_PRISM_tissue']], row.names=NULL)

# -> Tumore... che tumore è?
# Identifico dai tessuti

prim=input[1,]$primary_tissue
seco=input[1,]$secondary_tissue
tert=input[1,]$tertiary_tissue
for (i in 1:(nrow(TCGA_project_to_PRISM_tissue))){
if(toString(TCGA_project_to_PRISM_tissue[i,]$PRISM_primary_tissue)==toString(prim))
  TCGA_project_to_PRISM_tissue<-TCGA_project_to_PRISM_tissue[i,]
else if(toString(TCGA_project_to_PRISM_tissue[i,]$PRISM_secondary_tissue)==toString(seco))
  TCGA_project_to_PRISM_tissue<-TCGA_project_to_PRISM_tissue[i,]
else if(toString(TCGA_project_to_PRISM_tissue[i,]$PRISM_tertiary_tissue)==toString(tert))
  TCGA_project_to_PRISM_tissue<-TCGA_project_to_PRISM_tissue[i,]
}

original_nrow <- nrow(TCGA_project_to_PRISM_tissue)

for (i in 1:(nrow(TCGA_project_to_PRISM_tissue))){
  
  print(i)
  
 
  tmp <- TCGA_project_to_PRISM_tissue[i,]
  label <- tmp$TCGA_tumor

  
  if (tmp$PRISM_primary_tissue != "-"){
    
    PRISM_tissue <- paste0("PRISM_primary_tissue_", tmp$PRISM_primary_tissue)
    
    if (i > original_nrow){
      
      input_subset <- subset(input, primary_tissue != tmp$PRISM_primary_tissue)
      
    } else {
      
      input_subset <- subset(input, primary_tissue == tmp$PRISM_primary_tissue)
      
    }
    
  } else if (tmp$PRISM_secondary_tissue != "-"){
    
    PRISM_tissue <- paste0("PRISM_secondary_tissue_", tmp$PRISM_secondary_tissue)
    
    if (i > original_nrow){
      
      input_subset <- subset(input, secondary_tissue != tmp$PRISM_secondary_tissue)
      
    } else {
      
      input_subset <- subset(input, secondary_tissue == tmp$PRISM_secondary_tissue)
      
    }
    
  } else if (tmp$PRISM_tertiary_tissue != "-"){
    
    PRISM_tissue <- paste0("PRISM_tertiary_tissue_", tmp$PRISM_tertiary_tissue)
    
    if (i > original_nrow){
      
      input_subset <- subset(input, tertiary_tissue != tmp$PRISM_tertiary_tissue)
      
    } else {
      
      input_subset <- subset(input, tertiary_tissue == tmp$PRISM_tertiary_tissue)
      
    }
    
  } else {
    
    next
    
  }
  
  print('')
  print('')
  TCGA_project_to_PRISM_tissue[i,"TCGA_tumor"] <- label
  TCGA_project_to_PRISM_tissue[i,"number_of_experiments"] <- nrow(input_subset)
  TCGA_project_to_PRISM_tissue[i,"number_of_different_cellLines"] <- length(unique(input_subset$cellLine_id))
  TCGA_project_to_PRISM_tissue[i,"number_of_different_drugs"] <- length(unique(input_subset$drug_name))
  
  if (length(unique(input_subset$cellLine_id)) > 1){
    
    model_output <- 
      lmer(log(auc) ~ (1 | drug_name) + (1 |cellLine_id), data = input_subset)
    
  } else {
    
    model_output <- 
      lmer(log(auc) ~ (1 | drug_name), data = input_subset)
    
  }
  
  output <- ranef(model_output)$`drug_name`
  colnames(output)[1] <- "random_effect"
  output$drug_name <- rownames(output)
  output$TCGA_project <- label
  output$net_drug_effect <- output$random_effect + coef(summary(model_output))[1,1]
  output <- output[,c("TCGA_project","drug_name","random_effect","net_drug_effect")]
  
  if (i == 1){
    
    combined_output <- output
    
  } else {
    
    combined_output <- 
      rbind(combined_output, output)
    
  }
  
  # Open a connection to a file for writing
  file_conn <- file(snakemake@input[['tumor_type']], open = "w")
  
  # Write the value to the file
  cat(label, file = file_conn)
  
  # Close the file connection
  close(file_conn)
  
}

#end for

TCGA_project_to_PRISM_tissue <- TCGA_project_to_PRISM_tissue[!is.na(TCGA_project_to_PRISM_tissue$number_of_experiments),]
saveRDS(combined_output, snakemake@output[['net_drug_effect_by_tumor_type']])
saveRDS(TCGA_project_to_PRISM_tissue, snakemake@output[['TCGA_project_to_PRISM_tissue_enriched']])

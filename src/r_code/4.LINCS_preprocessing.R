#!/usr/bin/env Rscript
library("readr")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cmapR")
BiocManager::install("PharmacoGx")
library("cmapR")
library("limma")
library("PharmacoGx")
library("statmod")
library("parallel")
library("lme4")
library("metafor")
#library("ggplot2")

#Qui vengono istanziate delle variabili con il nome uguale al nome file, escluso il prefisso del file
source(snakemake@input[['load_LINCS_meta']])

#Lincs 
#https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit?pli=1#heading=h.5lvc853bqeqc
#Level 3 - INF_mlr12k - Gene expression (GEX, Level 2) levels that have been normalized 
#to invariant gene set curves, and quantile normalized across each plate, 
#and inferred values based on those normalized values

# 65 GB file, download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 seperately and specify path:
gctx_path <- snakemake@input[["gctx_path"]]

drug_ranef_and_sd <- readRDS(snakemake@input[['drug_ranef_and_sd']])
PRISM_drugs <- unique(drug_ranef_and_sd$drug_name)
# Using 6h experiments: 48,075 obs, 908 rna plates and 67 cell lines, 572 unique pert_iname
# Using 24h experiments: 41,952 obs, 1,022 rna plates and 28 cell lines, 514 unique pert_iname

for (pert_time_subset in c("6","24")){
  loadLINCSmeta()

  gene_info <- subset(gene_info, pr_is_lm == 1)
  gene_info$pr_gene_id <- as.character(gene_info$pr_gene_id)
  #978 geni
  landmark_genes <- gene_info$pr_gene_id
  
  tumor_cell_lines <- 
    subset(cell_info, sample_type=="tumor")$cell_id
  
  test <- subset(inst_info, pert_type == "trt_cp" & pert_dose_unit != "-666") # N = 672,106
  test <- test[test$pert_iname %in% PRISM_drugs,] # N = 105,641
  #View(table(factor(test$pert_dose))) # 66,184 of 105,641 (63%) are tested at 10 microM
  test_10um <- subset(test, pert_dose == 10)
  
  #View(table(factor(test_10um$pert_time)))
  test_10um_subset <- subset(test_10um, pert_time == as.integer(pert_time_subset))
  
  experiments_per_drug_freq <- as.data.frame(table(test_10um_subset$pert_iname))
  cell_lines_per_drug_freq <- as.data.frame(table(unique(test_10um_subset[,c("cell_id","pert_iname")])$pert_iname),stringsAsFactors = FALSE)
  
  LINCS_experiments_and_cell_lines_per_drug <- list(experiments_per_drug_freq = experiments_per_drug_freq, cell_lines_per_drug_freq = cell_lines_per_drug_freq)
  saveRDS(LINCS_experiments_and_cell_lines_per_drug, file = paste0(snakemake@output[['lincs_experiments']],pert_time_subset,"h.Rds"))
  
  sum(cell_lines_per_drug_freq$Freq >= 5) / 669 # 85% for 6h
  cell_lines_per_drug_freq <- cell_lines_per_drug_freq[cell_lines_per_drug_freq$Freq >= 5,]
  summary(cell_lines_per_drug_freq$Freq)
  
  compounds <- cell_lines_per_drug_freq$Var1 
  inst_info <- inst_info[inst_info$rna_plate %in% unique(test_10um_subset$rna_plate),]
  pert_iname_freq <- as.data.frame(table(inst_info$pert_iname))
  pert_iname_freq <- pert_iname_freq[order(pert_iname_freq$Freq, decreasing = T),]
  nrow(subset(pert_iname_freq, Freq >= 50))
  
  inst_info <- 
    inst_info[
      (inst_info$pert_type == "ctl_vehicle" & inst_info$pert_time == pert_time_subset)
      | 
        (
          inst_info$pert_iname %in% compounds &
            inst_info$pert_type == "trt_cp" &
            inst_info$pert_dose == 10 &
            inst_info$pert_time == pert_time_subset
        ),
    ] 

  table(as.character(subset(inst_info, pert_type == "ctl_vehicle")$pert_iname))
  drug_names <- unique(inst_info$pert_iname)
  drug_names <- c("DMSO",drug_names[drug_names != "DMSO"]) 
  inst_info$pert_iname <- factor(inst_info$pert_iname, levels = drug_names)
  
  # 1 run through loop costs 3 min per gene on Macbook Pro 2014 16 GB memory; 
  # processo i 978 landmark genes
  for (i in 1:978){
    
    print(i)
    input <- inst_info
    
    data_subset <-
      parse_gctx(
        fname = gctx_path,
        rid = gene_info$pr_gene_id[i],
        cid = input$inst_id 
      )@mat  


    input$expr <- data_subset[1,]

    model_output <- 
      lmer(expr ~ pert_iname + (1|cell_id) + (1|rna_plate), data = input)
    
    drug_coefs <- as.data.frame(coef(summary(model_output)))
    drug_coefs$drug <- rownames(drug_coefs)
    drug_coefs <- drug_coefs[2:nrow(drug_coefs),]
    drug_coefs$drug <- substr(drug_coefs$drug, 11, nchar(drug_coefs$drug))
    drug_coefs$gene <- gene_info$pr_gene_symbol[i]
    
    drug_coefs <- drug_coefs[,c("drug","gene","Estimate","Std. Error","t value")]
    rownames(drug_coefs) <- NULL

    if (i == 1){
      
      combined_DEG_drugs <- 
        drug_coefs
      
    } else {
      
      combined_DEG_drugs <- 
        rbind(combined_DEG_drugs, drug_coefs)
      
    }
  }
  
  saveRDS(combined_DEG_drugs, file = paste0(snakemake@output[['combined_deg_drug']],pert_time_subset,"h.Rds"))
  print('arrivato2')
}

### Combine 6h and 24h using meta-analysis: ###

PRISM_combined_DEG_drugs_6h <- readRDS(paste0(snakemake@output[['combined_deg_drug']],"6h.Rds"))
PRISM_combined_DEG_drugs_24h <- readRDS(paste0(snakemake@output[['combined_deg_drug']],"24h.Rds"))

shared_drugs <- intersect(unique(PRISM_combined_DEG_drugs_6h$drug), unique(PRISM_combined_DEG_drugs_24h$drug))
PRISM_combined_DEG_drugs_6h <- PRISM_combined_DEG_drugs_6h[PRISM_combined_DEG_drugs_6h$drug %in% shared_drugs,]
PRISM_combined_DEG_drugs_24h <- PRISM_combined_DEG_drugs_24h[PRISM_combined_DEG_drugs_24h$drug %in% shared_drugs,]
PRISM_combined_DEG_drugs_combined <- merge(x = PRISM_combined_DEG_drugs_6h, y = PRISM_combined_DEG_drugs_24h, by = c("gene","drug"))
colnames(PRISM_combined_DEG_drugs_combined)[3:5] <- c("Estimate_6h","Std.error_6h","T.value_6h")
colnames(PRISM_combined_DEG_drugs_combined)[6:8] <- c("Estimate_24h","Std.error_24h","T.value_24h")

for (i in 1:nrow(PRISM_combined_DEG_drugs_combined)){
  
  if (i %% 1000 == 0) print(i)
  
  tmp <- PRISM_combined_DEG_drugs_combined[i,]
  meta_output <- rma(yi = c(tmp$Estimate_6h, tmp$Estimate_24h), sei = c(tmp$Std.error_6h, tmp$Std.error_24h))
  
  PRISM_combined_DEG_drugs_combined[i,c("Estimate_6h_24h_combined", "Std.error_6h_24h_combined","T.value_6h_24h_combined","P.value_6h_24h_combined","Heterogeneity.p.value")] <-
    c(meta_output$b[,1], meta_output$se, meta_output$zval, meta_output$pval, meta_output$QEp)
  
}
saveRDS(PRISM_combined_DEG_drugs_combined, file = snakemake@output[['PRISM_combined_DEG_drugs']])

### Filtering out down-stream effects of reduced cell viability: ###
PRISM_combined_DEG_drugs_combined <- readRDS(snakemake@output[['PRISM_combined_DEG_drugs']])
drug_ranef_and_sd <- readRDS(snakemake@output[['drug_ranef_and_sd']])
PRISM_combined_DEG_drugs_combined <- merge(
  x = PRISM_combined_DEG_drugs_combined, y = drug_ranef_and_sd, by.x = "drug", by.y = "drug_name"
)
PRISM_combined_DEG_drugs_combined <- PRISM_combined_DEG_drugs_combined[order(PRISM_combined_DEG_drugs_combined$gene,PRISM_combined_DEG_drugs_combined$drug),]

input_columns <- c("Estimate_6h", "Estimate_24h","Estimate_6h_24h_combined")
for (i in 1:length(input_columns)){
  
  print(paste0("Busy with ",input_columns[i]))
  
  PRISM_combined_DEG_drugs_combined$DE_expr <- PRISM_combined_DEG_drugs_combined[,input_columns[i]]
  
  
  gene_association_with_toxicity <- data.frame(
    Input = input_columns[i],
    Gene = unique(PRISM_combined_DEG_drugs_combined$gene),
    Rho = NA,
    Rho.P.value = NA
  )
  
  for (j in 1:nrow(gene_association_with_toxicity)){
    
    indices <- which(PRISM_combined_DEG_drugs_combined$gene == gene_association_with_toxicity$Gene[j])
    model_input <- PRISM_combined_DEG_drugs_combined[indices,]
    
    if (input_columns[i] != "Estimate_6h_24h_combined"){
      
      model_output <- lm(DE_expr ~ drug_ranef, data = model_input)
      PRISM_combined_DEG_drugs_combined[indices,paste0(input_columns[i],"_corrected")] <- resid(model_output)
      
    }
    
    spearman_output <- cor.test(x = model_input$drug_ranef, y = model_input$DE_expr, method = "spearman")
    gene_association_with_toxicity[j,c("Rho","Rho.P.value")] <- c(spearman_output$estimate, spearman_output$p.value)
    
  }
  
  gene_association_with_toxicity$Rho.FDR.value <-
    p.adjust(gene_association_with_toxicity$Rho.P.value, method = "BH")
  
  if (i == 1){
    
    gene_association_with_toxicity_combined <- gene_association_with_toxicity
    
  } else {
    
    gene_association_with_toxicity_combined <- 
      rbind(gene_association_with_toxicity_combined, gene_association_with_toxicity)
    
  }
  
  if (input_columns[i] != "Estimate_6h_24h_combined"){
    
    suffix <- substr(input_columns[i], 10,nchar(input_columns[i]))
    
    PRISM_combined_DEG_drugs_combined[,paste0("T.value_",suffix,"_corrected")] <-
      PRISM_combined_DEG_drugs_combined[,paste0(input_columns[i],"_corrected")] /
      PRISM_combined_DEG_drugs_combined[,paste0("Std.error_",suffix)]
    
  }
  
}

saveRDS(gene_association_with_toxicity_combined,snakemake@output[['gene_association_with_toxicity']])


### Calculate meta-analysis combined estimate of corrected 6h and 24h estimates: ###

for (i in 1:nrow(PRISM_combined_DEG_drugs_combined)){
  
  if (i %% 1000 == 0) print(i)
  
  tmp <- PRISM_combined_DEG_drugs_combined[i,]
  meta_output <- rma(yi = c(tmp$Estimate_6h_corrected, tmp$Estimate_24h_corrected), sei = c(tmp$Std.error_6h, tmp$Std.error_24h))
  
  PRISM_combined_DEG_drugs_combined[i,c("Estimate_6h_24h_combined_corrected", "Std.error_6h_24h_combined_corrected","T.value_6h_24h_combined_corrected","P.value_6h_24h_combined_corrected","Heterogeneity.p.value_corrected")] <-
    c(meta_output$b[,1], meta_output$se, meta_output$zval, meta_output$pval, meta_output$QEp)
  
}
# Example for drug == "YM-155" & gene == "GADD45A":
# Standard error decreases from 0.3 uncorrected to 0.1 for combined corrected data, because of no more heterogenity in the different estimates from 6h and 24h
saveRDS(PRISM_combined_DEG_drugs_combined, file = snakemake@output[['PRISM_combined_DEG_drugs_corrected']])
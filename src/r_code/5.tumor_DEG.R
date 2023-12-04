#!/usr/bin/env Rscript
library(lme4)
library(edgeR)
library(metafor)
library(ggplot2)
library(egg)
library(readr)
library(matrixStats)
library(ragg)

# Load files:
project_info <- as.data.frame(read_delim(snakemake@input[['project_info']],"\t", escape_double = FALSE, trim_ws = TRUE))
gene_expr <- as.data.frame(read_delim(snakemake@input[['gene_expr']],"\t", escape_double = FALSE, trim_ws = TRUE))
gene_info <- readRDS(snakemake@input[['drug_ranef_and_sd']])


# Project info 
project_info <- subset(project_info, program_code.project == "TCGA")
project_info$sample_short <- substr(project_info$sample, 1, 15)
project_info$tissue_type <- "Other"
project_info$tissue_type[endsWith(project_info$sample_short, "01")] <- "Tumor"
project_info$tissue_type[endsWith(project_info$sample_short, "11")] <- "Normal"
project_info <- subset(project_info, tissue_type != "Other")
project_info <- project_info[project_info$sample %in% colnames(gene_expr),]

Normal_tissue_type_frequency <- as.data.frame(table(subset(project_info, tissue_type == "Normal")[,c("disease_code.project")]), stringsAsFactors = F)
Tumor_tissue_type_frequency <- as.data.frame(table(subset(project_info, tissue_type == "Tumor")[,c("disease_code.project")]), stringsAsFactors = F)

Normal_tissue_type_frequency <- subset(Normal_tissue_type_frequency, Freq >= 2)
Tumor_tissue_type_frequency <- Tumor_tissue_type_frequency[Tumor_tissue_type_frequency$Var1 %in% Normal_tissue_type_frequency$Var1,]
combined_Normal_Tumor_frequency <- merge(x = Tumor_tissue_type_frequency, y = Normal_tissue_type_frequency, by = "Var1")
colnames(combined_Normal_Tumor_frequency) <- c("Tumor_type","N_tumor_samples","N_normal_samples")
combined_Normal_Tumor_frequency <- combined_Normal_Tumor_frequency[order(combined_Normal_Tumor_frequency$Tumor_type),]
write_csv(combined_Normal_Tumor_frequency, file = "Publication/Tables/Table_S1.csv")

split_genes <- unlist(strsplit(gene_expr$Ensembl_ID,".", fixed=TRUE))
split_genes <- split_genes[seq(from=1, to=length(split_genes), by=2)]
rownames(gene_expr) <- split_genes
protein_coding_genes <- unique(subset(gene_info, protein_coding==T)$`Gene stable ID`)
gene_expr <- gene_expr[rownames(gene_expr) %in% protein_coding_genes,]
gene_expr <- gene_expr[,colnames(gene_expr) %in% project_info$sample]


gene_stats <- 
  data.frame(
    Gene = rownames(gene_expr),
    Mean = rowMeans(gene_expr),
    SDs = rowSds(as.matrix(gene_expr))
  )



# From Mean > 6 there's a downward trend visible
filtered_gene_stats <- subset(gene_stats, Mean > 6)
gene_expr <- gene_expr[rownames(gene_expr) %in% filtered_gene_stats$Gene,] # 13,386 genes left

for (i in 1:ncol(gene_expr)) gene_expr[,i] <- (2 ^ gene_expr[,i]) - 1
project_info <- project_info[order(project_info$sample),]
gene_expr <- gene_expr[,order(colnames(gene_expr))]

#tumorType <- snakemake@input[['tumor_type']]
tumorType = read_file(snakemake@input[['tumor_type_r']], header = FALSE)


  print(paste0("Processing ", tumorType, "; file = ", i))
  
  project_info_subset <- subset(project_info, disease_code.project == tumorType)
  project_info_subset$tissue_type <-
    factor(project_info_subset$tissue_type, levels = c("Normal","Tumor"))
    
  gene_expr_subset <- gene_expr[,colnames(gene_expr) %in% project_info_subset$sample]
  
  design <- model.matrix(~ tissue_type, data = project_info_subset)
  dge  <- DGEList(gene_expr_subset, remove.zeros = TRUE)
  dge <- calcNormFactors(dge, method = 'upperquartile') 
  v <- voom(dge, design, plot = T) # Looks good

  fitvoom <- lmFit(v, design)
  fitvoom <- eBayes(fitvoom)
  SE <- sqrt(fitvoom$s2.post) * fitvoom$stdev.unscaled
  
  toptable_result <- topTable(fitvoom,coef=2, number = 10^6)
  toptable_result$Gene <- rownames(toptable_result)
  colnames(toptable_result)[1] <- "log2_FC_estimate"
  
  standard_errors <- as.data.frame(SE[order(rownames(SE)),])
  standard_errors$Gene <- rownames(standard_errors)
  standard_errors <- standard_errors[,c("Gene","tissue_typeTumor")]
  colnames(standard_errors)[2] <- "log2_FC_estimate_SE"
  
  toptable_result <- merge(x = toptable_result, y = standard_errors, by = "Gene")
  toptable_result$Tumor <- tumorType
  toptable_result <- toptable_result[,c("Tumor","Gene","log2_FC_estimate","log2_FC_estimate_SE","t","P.Value","adj.P.Val")]
  
  if (i == 1){
    
    combined_tumor_vs_normal_DEG <- toptable_result
    
  } else {
    
    combined_tumor_vs_normal_DEG <- rbind(combined_tumor_vs_normal_DEG, toptable_result)
    
  }



saveRDS(combined_tumor_vs_normal_DEG, file = snakemake@output[["tumor_vs_normal"]])

#### Combine all tumor type DEG into 1 meta-analysis estimate each: ###
tumor_types <- unique(combined_tumor_vs_normal_DEG$Tumor)

average_log2_fold_change_vs_normal <- data.frame(
  Gene = unique(combined_tumor_vs_normal_DEG$Gene),
  Mean_log2_fold_change = NA,
  Mean_log2_fold_change_SE = NA,
  Mean_log2_fold_change_Tau = NA,
  Mean_log2_fold_change_p.value = NA
)

for (i in 1:nrow(average_log2_fold_change_vs_normal)){
  
  input <- subset(combined_tumor_vs_normal_DEG, Gene == average_log2_fold_change_vs_normal[i,"Gene"])
  rma_output <- rma(yi = log2_FC_estimate, sei = log2_FC_estimate_SE, data = input)
  
  average_log2_fold_change_vs_normal[i, 2:5] <-
    c(rma_output$b[,1], rma_output$se, sqrt(rma_output$tau2) ,rma_output$pval)
  
}


average_log2_fold_change_vs_normal$Mean_log2_fold_change_div_Tau <-
  average_log2_fold_change_vs_normal$Mean_log2_fold_change / 
  average_log2_fold_change_vs_normal$Mean_log2_fold_change_Tau

ggplot(data = average_log2_fold_change_vs_normal, aes(x = Mean_log2_fold_change_div_Tau)) +
  geom_density()

saveRDS(average_log2_fold_change_vs_normal, file = snakemake@out[["average_log2fold"]])

getSensitivitySignature <- function(drug, fit_method, sensivity.measure, genexp_normalization_method){
  
  data_subset <- 
    subset(mean_log_auc, drug_name == drug)
  
  logCPMs_subset <-
    logCPMs[, colnames(logCPMs) %in% data_subset$cellLine_id]

  if (genexp_normalization_method != "none"){
    
    cellLine_info_subset <-
      cellLine_info[cellLine_info$DepMap_ID %in% colnames(logCPMs_subset),]
    
    if(sum(cellLine_info_subset$DepMap_ID == colnames(logCPMs_subset)) != nrow(cellLine_info_subset)){
      
      stop("cellLine_info_subset IDs not matching genexp IDs")
      
    }
    
    for (i in 1:ncol(logCPMs_subset)){
      
      logCPMs_subset[,i] <-
        log(
          exp(logCPMs_subset[,i]) / 
          cellLine_info_subset[i,paste0(genexp_normalization_method,"_norm.factor")]
        )
      
    }
    
  }
  
  
  data_subset <- data_subset[order(data_subset$cellLine_id),]
  logCPMs_subset <- logCPMs_subset[,order(colnames(logCPMs_subset))]
  
  output <-
    data.frame(
      Drug = drug,
      Gene = rownames(logCPMs_subset),
      estimate = NA,
      p.value = NA
    )
  
  for (i in 1:nrow(logCPMs_subset)){
    
    data_subset$logCPM <- logCPMs_subset[i,]

    library(ggplot2)
    ggplot(data = data_subset, aes(x = logCPM, y = exp(uncorrected_log_auc))) + geom_point() + geom_smooth(method = "lm")
    
    if (fit_method == "spearman"){
      
      cor_result <- 
        cor.test(x = data_subset$logCPM, y = data_subset[,sensivity.measure], method = "spearman")
      
      output[i,c("estimate","p.value")] <-
        c(cor_result$estimate, cor_result$p.value)
      
    } else if (fit_method == "lm_with_weights"){
      
      data_subset$ss_corrected_sigma_cellLine <-
        data_subset$sigma_cellLine / sqrt(data_subset$sample_size)
      
      cor_result <- 
        coef(summary(lm(
          formula = as.formula(paste0(sensivity.measure, " ~ logCPM")), 
          data = data_subset,
          weights = (1 / data_subset$ss_corrected_sigma_cellLine^2)
        )))
      
      output[i,c("estimate","p.value")] <-
        c(cor_result[2,3], cor_result[2,4])
      
    } else if (fit_method == "lm_without_weights"){
      
      cor_result <- 
        coef(summary(lm(
          formula = as.formula(paste0(sensivity.measure, " ~ logCPM")), 
          data = data_subset
        )))
      
      output[i,c("estimate","p.value")] <-
        c(cor_result[2,3], cor_result[2,4])
      
    }
    
  }
  
  output$fdr.value <-
    p.adjust(output$p.value, method = "BH")
  
  return(output)
}

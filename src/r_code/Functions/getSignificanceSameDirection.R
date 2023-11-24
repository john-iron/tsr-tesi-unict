getSignificanceSameDirection <- function(drugSensitivitySignature, drugSignature){
  
  output <- data.frame(
    total_number_of_genes_overlap = NA,
    total_number_of_genes_same_direction = NA,
    total_same_direction.p.value = NA,
    upregulated_number_of_genes_overlap = NA,
    upregulated_number_of_genes_same_direction = NA,
    upregulated_same_direction.p.value = NA,
    downregulated_number_of_genes_overlap = NA,
    downregulated_number_of_genes_same_direction = NA,
    downregulated_same_direction.p.value = NA
  )
  
  input <- merge(
    x = subset(drugSensitivitySignature, fdr.value < 0.1),
    y = subset(drugSignature, fdr < 0.1),
    by.x = "Gene", by.y = "gene"
  )
  
  if (nrow(input) > 10){
    
    input$same_direction <-
      (input$estimate.x > 0 & input$estimate.y > 0) |
      (input$estimate.x < 0 & input$estimate.y < 0)
    
    upregulated <- subset(input, estimate.y > 0)
    downregulated <- subset(input, estimate.y < 0)
    
    output <- data.frame(
      total_number_of_genes_overlap = nrow(input),
      total_number_of_genes_same_direction = sum(input$same_direction),
      total_same_direction.p.value = NA,
      upregulated_number_of_genes_overlap = nrow(upregulated),
      upregulated_number_of_genes_same_direction = sum(upregulated$same_direction),
      upregulated_same_direction.p.value = NA,
      downregulated_number_of_genes_overlap = nrow(downregulated),
      downregulated_number_of_genes_same_direction = sum(downregulated$same_direction),
      downregulated_same_direction.p.value = NA
    )
    
    output$total_same_direction.p.value <-
      prop.test(
        p = 0.5, 
        alternative = "greater", 
        x = output$total_number_of_genes_same_direction, 
        n = output$total_number_of_genes_overlap
      )$p.value  
    
    
    if (nrow(upregulated) >= 10){
      
      output$upregulated_same_direction.p.value <-
        prop.test(
          p = 0.5, 
          alternative = "greater", 
          x = output$upregulated_number_of_genes_same_direction, 
          n = output$upregulated_number_of_genes_overlap
        )$p.value  
      
    }

    
    if (nrow(downregulated) >= 10){
      
      output$downregulated_same_direction.p.value <-
        prop.test(
          p = 0.5, 
          alternative = "greater", 
          x = output$downregulated_number_of_genes_same_direction, 
          n = output$downregulated_number_of_genes_overlap
        )$p.value  
      
    }

  }
  
  return(output)

} 

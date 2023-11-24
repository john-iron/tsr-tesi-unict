add_association_between_residuals_and_cellLine_ranef <- function(input){
  
  input$drug_ranef$Resid_vs_cellLine_ranef_rho <- NA
  input$drug_ranef$Resid_vs_cellLine_ranef_Spearman_rho.p.value <- NA
  
  for (i in 1:nrow(input$drug_ranef)){
    
    resid_and_fitted_subset <- 
      subset(input$resid_and_fitted, drug_name == input$drug_ranef[i,"drug_name"])
    
    resid_and_fitted_subset <-
      merge(x = resid_and_fitted_subset, y = input$cellLine_ranef, by = "cellLine_id")
    
    ggplot(data = resid_and_fitted_subset, aes(x = cellLine_ranef, y = model_resid)) + geom_point() + geom_smooth() + geom_abline(slope=-1, linetype = "dashed")
    
    cor_result <- 
      cor.test(x = resid_and_fitted_subset$cellLine_ranef, y = resid_and_fitted_subset$model_resid, method = "spearman")
    
    input$drug_ranef[i,c("Resid_vs_cellLine_ranef_rho","Resid_vs_cellLine_ranef_Spearman_rho.p.value")] <-
      c(cor_result$estimate,cor_result$p.value)
    
  }
  
  plot <- 
    ggplot(data = input$drug_ranef, aes(x = Resid_vs_cellLine_ranef_rho, y = sd_resid_drug)) + geom_point() + geom_smooth()
  
  print(plot)
  
  
  return(input)
}
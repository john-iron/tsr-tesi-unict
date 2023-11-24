load_Prism_data_and_model <- function(
  drugs_subset = NULL, 
  cellLines_subset = NULL, 
  estimated_cellLine_ranef = NULL,
  sd_resid_drug = NULL
)
{
  
  input <- 
    read.csv(
      "External_input/DepMap_Public_21Q2/secondary-screen-dose-response-curve-parameters_reduced.csv",
      stringsAsFactors = F
    )
  
  #rinomino le colonne
  colnames(input)[colnames(input) == "name"] <- "drug_name"
  colnames(input)[colnames(input) == "depmap_id"] <- "cellLine_id"
  
  #escludo i valori vuoti
  input <- subset(input, !is.na(cellLine_id))
  input <- subset(input, !is.na(drug_name))

  cell.line.info <- read.csv(
    "External_input/DepMap_Public_21Q2/secondary-screen-cell-line-info.csv",
    stringsAsFactors = F
  )
  
  #merge dei valori tra i due dataset
  input <- merge(x = input, y = cell.line.info, by.x = "cellLine_id", by.y = "depmap_id")
  
  table(input$primary_tissue)


  #salto if
  if (!is.null(drugs_subset)){
    
    input <- input[input$drug_name %in% drugs_subset,]
    print('ciao, sono in drug subset')
  }
  
  #salto if
  if (!is.null(cellLines_subset)){
    
    input <- input[input$cellLine_id %in% cellLines_subset,]
    print('ciao sono in celline subset')
  }
  
  #salto if
  if (!is.null(estimated_cellLine_ranef)){
    
    input <- 
      merge(
        x = input, 
        y = estimated_cellLine_ranef[,c("cellLine_id","cellLine_ranef")], 
        by = "cellLine_id"
      )
    
  }
  
  #salto if
  if (!is.null(sd_resid_drug)){
    
    input <- 
      merge(
        x = input, 
        y = sd_resid_drug[,c("drug_name","sd_resid_drug")], 
        by = "drug_name"
      )
    
  }
  
  #preparo i dati per il plot
  number_of_cellLines_tested_per_drug <- as.data.frame(table(unique(input[,c("drug_name","cellLine_id")])$drug_name), stringsAsFactors = F)
  number_of_drugs_tested_per_cellLine <- as.data.frame(table(unique(input[,c("drug_name","cellLine_id")])$cellLine_id), stringsAsFactors = F)
  
  plot1 <- ggplot(data = number_of_cellLines_tested_per_drug, aes(x = Freq)) + geom_histogram() + labs(
    x = "Number of cell lines tested per drug", y = "Count"
  ) 
  
  plot2 <- ggplot(data = number_of_drugs_tested_per_cellLine, aes(x = Freq)) + geom_histogram() + labs(
    x = "Number of drugs tested per cell line", y = "Count"
  ) 
  #entro dentro il blocco
  if (is.null(sd_resid_drug)){
    #entro dentro il blocco
    if (is.null(estimated_cellLine_ranef)){
      #creo il modello lineare
      #model_output <- 
      #  lmer(log(auc) ~ (1 | drug_name) + (1 | cellLine_id), data = input)
      model_output <- 
        lmer(log(auc) ~ (1 | drug_name), data = input)
      
      
    } else {
      
      model_output <- 
        lmer(log(auc) ~ cellLine_ranef + (1 + cellLine_ranef | drug_name), data = input)
      
    }
    
  } else {
    
    if (is.null(estimated_cellLine_ranef)){
      
      model_output <- 
        lmer(log(auc) ~ (1 | drug_name) + (1 | cellLine_id), weights = 1 / (sd_resid_drug ^ 2), data = input)
      
    } else {
      
      model_output <- 
        lmer(log(auc) ~ cellLine_ranef + (1 + cellLine_ranef | drug_name), weights = 1 / (sd_resid_drug ^ 2), data = input)
      
    }
    
  }
  
  

  
  print(summary(model_output))
  # The average drug in the average cell line reduces AUC by exp(-0.078613) = -8%. 
  
  resid_and_fitted <- model_output@frame
  resid_and_fitted$cellLine_id <- input$cellLine_id
  resid_and_fitted$model_resid <- resid(model_output)
  resid_and_fitted$fitted <- fitted(model_output)
  
  drug_ranef <- ranef(model_output)$drug_name
  colnames(drug_ranef)[1] <- "drug_ranef"
  drug_ranef$drug_name <- rownames(drug_ranef)
  sd_residuals_drug <- aggregate(model_resid ~ drug_name, FUN = sd, data = resid_and_fitted)
  colnames(sd_residuals_drug)[2] <- "sd_resid_drug"
  drug_ranef <- merge(x = drug_ranef, y = sd_residuals_drug, by = "drug_name")
  drug_ranef$net_drug_effect <- drug_ranef$drug_ranef + coef(summary(model_output))[1,1]
  
  plot3 <- ggplot(data = drug_ranef, aes(x = net_drug_effect, y = sd_resid_drug)) + 
    geom_point() + geom_smooth() + labs(
      x = "Net drug effect", y = expression(sigma[resid]*' across cell lines, same drug')
    ) + scale_y_log10()
  
  output <- list()
  output$input <- input
  output$resid_and_fitted <- resid_and_fitted
  output$model_output <- model_output
  output$drug_ranef <- drug_ranef
  output$number_of_cellLines_tested_per_drug <- number_of_cellLines_tested_per_drug
  output$number_of_drugs_tested_per_cellLine <- number_of_drugs_tested_per_cellLine
  


  
  return(output)
  
}
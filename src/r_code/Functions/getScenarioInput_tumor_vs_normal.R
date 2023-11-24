getScenarioInput_tumor_vs_normal <- function(){
  
  gene_selections <- c("50_most_significant","100_most_significant","150_most_significant","Bin_Chen_original")
  DEG_conditions <- c("6h","24h","6h_24h_combined","6h_24h_combined_corrected")
  log2_FC_inputs <- c(
    "BLCA", "BRCA", "CHOL", "COAD", "ESCA", "GBM", "KICH", "KIRC", 
    "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "SARC", "STAD", 
    "THCA", "UCEC"
  )
  
  output <- data.frame(
    gene_selection = rep(NA, times=length(gene_selections)*length(DEG_conditions)),
    DEG_condition = NA,
    log2_FC_input = NA
  )
  
  counter <- 1
  for (i in 1:length(gene_selections)){
    
    for (j in 1:length(DEG_conditions)){
      
      for (k in 1:length(log2_FC_inputs)){
        
        output[counter,] <-
          c(
            gene_selections[i],
            DEG_conditions[j],
            log2_FC_inputs[k]
          )
        
        counter <-
          counter + 1
        
      }

    }
  }
  
  return(output)
  
}
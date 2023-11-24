get_Alternative_Connectivity_Score_random_distribution <- function(n_genes_up, n_genes_down, n_genes, n_permutations){
  
  output <- numeric(n_permutations)
  
  for (i in 1:length(output)){
    
    #if (i %% 10000 == 0) print(i)
    
    drug_signature <- data.frame(
      ids = 1:n_genes,
      rank = sample(1:n_genes, replace=F)
    )
    
    DEG_genes <- sample(1:(n_genes_up+n_genes_down), replace=F)
    sig_up <- data.frame(GeneID = DEG_genes[1:n_genes_up])
    sig_down <- data.frame(GeneID = DEG_genes[(n_genes_up+1):length(DEG_genes)])

    output[i] <- cmap_score_new(
      sig_up[,"GeneID",drop=F],
      sig_down[,"GeneID",drop=F],
      drug_signature[,c("ids","rank")]
    )
    
  }
  
  return(output)
  
}
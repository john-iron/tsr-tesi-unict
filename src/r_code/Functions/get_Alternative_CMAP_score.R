get_Alternative_CMAP_score <- function(drugSensitivitySignature, drugSignature, alt_cs_distribution){
  
  sig_up <- subset(drugSensitivitySignature, estimate > 0)
  sig_up$GeneID <- as.character(as.integer(substr(rownames(sig_up), 5, nchar(sig_up))))
  sig_up <- sig_up[order(sig_up$estimate, decreasing = T),]
  
  sig_down <- subset(drugSensitivitySignature, estimate < 0)
  sig_down$GeneID <- as.character(as.integer(substr(rownames(sig_down), 5, nchar(sig_down))))
  sig_down <- sig_down[order(sig_down$estimate, decreasing = F),]
  
  drug_signature <- drugSignature
  drug_signature$ids <- as.character(as.integer(substr(rownames(drug_signature), 5, nchar(drug_signature))))
  drug_signature$rank <- rank(-1 * drug_signature$estimate, ties.method="random")
  
  cs <-
    cmap_score_new(
      sig_up[,"GeneID",drop=F],
      sig_down[,"GeneID",drop=F],
      drug_signature[,c("ids","rank")]
    )
  
  output <- c(
    cs, 
    sum(alt_cs_distribution < cs) / length(alt_cs_distribution)
  )

  names(output) <- c("connectivity_score","single_sided_p.value")
  return(output)
  
}









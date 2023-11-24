loadLINCSmeta <- function(){
  
  dir <- "External_input/LINCS_meta"
  
  for (f in list.files(dir)){
    
    var_name <- 
      strsplit(strsplit(f,"GSE92742_Broad_LINCS_")[[1]][2],".txt")[[1]][1]
    
    print(paste0("Loading ", var_name,"..."))
    
    df <- 
      read_delim(
        paste0(dir,"/",f), 
        "\t", 
        escape_double = FALSE, 
        trim_ws = TRUE
      )
    
    assign(var_name,df, envir=.GlobalEnv)
    
  }

}
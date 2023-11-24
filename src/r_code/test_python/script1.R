#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default sample computing
  args[2] = ""
}
print('INIZIO SCRIPT1')
Sys.sleep(5)
linea_cellulare <- c(args[2],args[1])
campioni <-args[3]
print(paste(linea_cellulare, collapse="_"))
print(campioni)
print('FINE SCRIPT1')

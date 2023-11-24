#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default sample computing
  args[2] = ""
}
print('INIZIO Main')
linea_cellulare <- c(args[2],args[1])
campioni <-args[3]
print(paste(linea_cellulare, collapse="_"))

newpath1 <- file.path(args[4], 'script1.R')
newpath2 <- file.path(args[4], 'script2.R')
newpath3 <- file.path(args[4], 'script3.R')
newpath4 <- file.path(args[4], 'script4.R')
newpath5 <- file.path(args[4], 'script5.R')

source(newpath1)
source(newpath2)
source(newpath3)
source(newpath4)
source(newpath5)
print('FINE Main')
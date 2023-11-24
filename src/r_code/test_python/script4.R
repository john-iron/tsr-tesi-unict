#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default sample computing
  args[2] = ""
}
print('INIZIO SCRIPT4')
Sys.sleep(5)

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])
print('Fine SCRIPT4')

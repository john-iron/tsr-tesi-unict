#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(sva)
library(drc)

sec_screencellline_info = "External_input/DepMap_Public_21Q2/secondary-screen-dose-response-curve-parameters.csv"
secondary_drc_gen = "External_input/DepMap_Public_21Q2/secondary_drc_table.csv"

#Seleziono la row_name corrispondente tramite sec_drc[1,8]
#Esempio per A549_LUNG row_name == ACH-000681

sec_drc = data.table::fread(secondary_drc_gen)

sec_screencellline_info_reduced = data.table::fread(sec_screencellline_info) %>%
  dplyr::filter(row_name==toString(sec_drc[1,8]))


sec_drc_distinct = distinct(sec_drc)

DATA2 = dplyr::left_join(x=sec_drc_distinct,y=sec_screencellline_info_reduced,by="broad_id")

df <- subset(DATA2, select = -c(upper_limit, lower_limit,slope,r2,auc.y,ec50,ic50,row_name.y))
colnames(df)[colnames(df) == "EC50"] <- "ec50"
colnames(df)[colnames(df) == "LowerLimit"] <- "lower_limit"
colnames(df)[colnames(df) == "auc.x"] <- "auc"
colnames(df)[colnames(df) == "R2"] <- "r2"
colnames(df)[colnames(df) == "row_name.x"] <- "row_name"
colnames(df)[colnames(df) == "ccle_name.x"] <- "ccle_name"


write.csv(df, "External_input/DepMap_Public_21Q2/secondary-screen-dose-response-curve-parameters_reduced.csv", row.names=FALSE)

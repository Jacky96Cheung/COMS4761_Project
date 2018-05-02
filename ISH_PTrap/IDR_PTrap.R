# IDR for Control and IP PhosphoTrap
library(dplyr)
library(idr)
# get to the correct directory 
setwd("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/results")

sample_ids <- dir(file.path("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/results"))
reps <- sample_ids[grepl("IP_Sucrose", sample_ids)]

for(rep in reps){
  filepath <- paste("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/results/" , rep, "/abundance.tsv", sep = "")
  temp <- read.table(filepath, sep = '\t', header = TRUE)
  t[rep] <- temp$tpm
}


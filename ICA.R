# Run Independent Component Analysis on PhosphoTrap data -- IP and CR

library(dplyr)
library(fastICA)

# Get to the correct directory 
setwd("/Users/elizabethboylesobolik/COMS4671_Project/ISH_Ptrap")

# load the data
ish_pt <- read.csv("output_idr.csv", header = TRUE, sep = ",")
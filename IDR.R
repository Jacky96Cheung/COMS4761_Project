# Run IDR on Unionized Allen Brain Atlas and Homogenate RNA Seq gene expression data
library(dplyr)
library(idr)

# get to the correct directory 
setwd("/Users/elizabethboylesobolik/COMS4671_Project/ISH_Ptrap")

# load the data 
ish_pt <- read.csv("output_idr.csv", header = TRUE, sep = ",")
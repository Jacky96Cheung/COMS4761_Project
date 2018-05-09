# Exploratory plots for relationship between ISH and PhosphoTrap mRNA seq data
# These are specifically for transcript per million and expression level metrics. 
library(dplyr)

# relevant file names and paths
p <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap"
dir <- "/exp_plots/" 
et <- "energy_tpm.jpg"
elt <- "energy_logtpm.jpg"
lelt <- "logenergy_logtpm.jpg"

setwd(p)

# read in the input file
ish_pt_raw <- read.csv("output.csv", header = TRUE, sep = ",")

# filter the data so that no values are undefined
ish_pt_raw <- filter(ish_pt_raw, Avg_TPM < 500, log(Avg_TPM) > 0)
ish_pt <- filter(ish_pt_raw, log(Avg_TPM) > 0 , log(Energy) > 0 )

# make plots!
filename <- paste(p, dir, et, sep = "")
jpeg(filename)
scatter.smooth(x = ish_pt_raw$Energy, y = ish_pt_raw$Avg_TPM, main = "Energy ~ Avg_TPM")
dev.off()

filename <- paste(p, dir, elt, sep = "")
jpeg(filename)
scatter.smooth(x = ish_pt_raw$Energy, y = log(ish_pt_raw$Avg_TPM), main = "Energy ~ Avg_TPM")
dev.off()

filename <- paste(p, dir, lelt, sep = "")
jpeg(filename)
scatter.smooth(x = log(ish_pt$Energy), y = log(ish_pt$Avg_TPM), main = "Energy ~ Avg_TPM")
dev.off()



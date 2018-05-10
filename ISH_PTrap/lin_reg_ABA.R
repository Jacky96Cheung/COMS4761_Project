# Plot simple scatter plots for linear regression analysis between
# ISH and mRNA Seq. metrics of gene expression 

library(dplyr)
d <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap"

setwd(d)

ish_pt_raw <- read.csv("output.csv", header = TRUE, sep = ",")

#Filter the results so that log plot does not produce undefined values
ish_pt <- filter(ish_pt_raw, Avg_TPM > 0 , Energy > 0 )

folder <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap/exploratory_scatter/"
jpeg(filename = paste(folder, "energy_count.jpg", sep = ""))
scatter.smooth(x = ish_pt$Energy, y = log(ish_pt$Avg_Count), main = "Energy ~ Avg_Count")
dev.off()

jpeg(filename = paste(folder, "density_count.jpg", sep = ""))
scatter.smooth(x = ish_pt$Density, y = log(ish_pt$Avg_Count), main = "Density ~ Avg_Count")
dev.off()

jpeg(filename = paste(folder, "energy_tpm.jpg", sep = ""))
scatter.smooth(x = ish_pt$Avg_TPM , y = ish_pt$Energy, main = "Energy ~ Avg_TPM")
dev.off()

jpeg(filename = paste(folder, "density_tpm.jpg", sep = ""))
scatter.smooth(x = log(ish_pt$Density), y = log(ish_pt$Avg_TPM), main = "Density ~ Avg_TPM")
dev.off()

jpeg(filename = paste(folder, "intensity_density.jpg", sep = ""))
scatter.smooth(x = log(ish_pt$Intensity), y = log(ish_pt$Density), main = "Intensity ~ Density")
dev.off()

jpeg(filename = paste(folder, "intensity_energy.jpg", sep = ""))
scatter.smooth(x = ish_pt$Intensity, y = ish_pt$Energy, main = "Intensity ~ Energy")

jpeg(filename = paste(folder, "density_energy.jpg", sep = ""))
scatter.smooth(x = ish_pt$Density, y = ish_pt$Energy, main = "Density ~ Energy")
dev.off()

jpeg(filename = paste(folder, "count_tpm.jpg", sep = ""))
scatter.smooth(x = log(ish_pt$Avg_Count), y = log(ish_pt$Avg_TPM), main = "Count ~ TPM")
dev.off()
# Plot simple scatter plots for linear regression analysis between
# ISH and mRNA Seq. metrics of gene expression 

library(dplyr)
d <- "/Users/elizabethboylesobolik/Desktop/"

setwd(d)
ish_pt_raw <- read.csv("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap/output.csv", header = TRUE, sep = ",")

ish_pt <- filter(ish_pt_raw, Avg_TPM > 0 , Energy > 0 )

jpeg("energy_count.jpg")
scatter.smooth(x = ish_pt$Energy, y = log(ish_pt$Avg_Count), main = "Energy ~ Avg_Count")
dev.off()

jpeg("density_count.jpg")
scatter.smooth(x = ish_pt$Density, y = log(ish_pt$Avg_Count), main = "Density ~ Avg_Count")
dev.off()

jpeg("energy_tpm.jpg")
scatter.smooth(x = ish_pt$Avg_TPM , y = ish_pt$Energy, main = "Energy ~ Avg_TPM")
dev.off()

jpeg("density_tpm.jpg")
scatter.smooth(x = log(ish_pt$Density), y = log(ish_pt$Avg_TPM), main = "Density ~ Avg_TPM")
dev.off()

jpeg("intensity_density.jpg")
scatter.smooth(x = log(ish_pt$Intensity), y = log(ish_pt$Density), main = "Intensity ~ Density")
dev.off()

jpeg("intensity_energy.jpg")
scatter.smooth(x = ish_pt$Intensity, y = ish_pt$Energy, main = "Intensity ~ Energy")

jpeg("density_energy.jpg")
scatter.smooth(x = ish_pt$Density, y = ish_pt$Energy, main = "Density ~ Energy")
dev.off()

jpeg("count_tpm.jpg")
scatter.smooth(x = log(ish_pt$Avg_Count), y = log(ish_pt$Avg_TPM), main = "Count ~ TPM")
dev.off()
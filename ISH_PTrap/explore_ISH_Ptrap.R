# Exploratory plots for relationship between ISH and PhosphoTrap mRNA seq data
setwd("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap")

ish_pt_raw <- read.csv("output.csv", header = TRUE, sep = ",")
ish_pt_raw <- filter(ish_pt_raw, Avg_TPM < 500, log(Avg_TPM) > 0)
ish_pt <- filter(ish_pt_raw, log(Avg_TPM) > 0 , log(Energy) > 0 )

jpeg("exp/energy_tpm.jpg")
scatter.smooth(x = ish_pt_raw$Energy, y = ish_pt_raw$Avg_TPM, main = "Energy ~ Avg_TPM")
dev.off()

jpeg("exp/energy_logtpm.jpg")
scatter.smooth(x = ish_pt_raw$Energy, y = log(ish_pt_raw$Avg_TPM), main = "Energy ~ Avg_TPM")
dev.off()

jpeg("exp/logenergy_logtpm.jpg")
scatter.smooth(x = log(ish_pt$Energy), y = log(ish_pt$Avg_TPM), main = "Energy ~ Avg_TPM")
dev.off()



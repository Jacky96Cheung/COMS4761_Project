# Run IDR on Unionized Allen Brain Atlas and Homogenate RNA Seq gene expression data
# Only looking at TPM, Counts, and Expression Level
library(dplyr)
library(idr)

# get to the correct directory 
d <- "/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap"
folder <- "/IDR_test/" # for saving the plots later
setwd(d)

# load the data 
e <- exp(1)
ish_pt_raw <- read.csv("all_data.csv", header = TRUE, sep = ",")

rank_tpm <- ish_pt_raw$TPM_Rank
rank_count <- ish_pt_raw$Count_Rank
rank_exp <- ish_pt_raw$Exp_Rank

# Specify percentage thresholds 
t <- seq(0.01, 0.99, by = 1/100)

tpm_count <- get.correspondence(rank_tpm, rank_count, t)
tpm_exp <- get.correspondence(rank_tpm, rank_exp, t)
count_exp <- get.correspondence(rank_count, rank_exp, t)

ish_pt_raw <- filter(ish_pt_raw, Exp_Level > 0, Avg_TPM > 1)
jpeg(filename = paste(d, folder,"tpm_exp.jpg", sep = ""))
scatter.smooth(x = log(ish_pt_raw$Exp_Level), y = log(ish_pt_raw$Avg_TPM), main = "Exp_Level ~ Avg_TPM")
dev.off()

jpeg(filename = paste(d, folder,"count_exp.jpg", sep = ""))
scatter.smooth(x = log(ish_pt_raw$Exp_Level), y = log(ish_pt_raw$Avg_Count), main = "Exp_Level ~ Avg_Count")
dev.off()

# plot correspondence curve on the scale of percentage for tpm and count
jpeg(filename = paste(d, folder,"tpm_count.jpg", sep = ""))
plot(tpm_count$psi$t, tpm_count$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_count$psi$t)),
     ylim=c(0, max(tpm_count$psi$value)), cex.lab=2)
lines(tpm_count$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

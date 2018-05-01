# Run IDR on Unionized Allen Brain Atlas and Homogenate RNA Seq gene expression data

library(dplyr)
library(idr)
# get to the correct directory 
setwd("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap")

# load the data 
ish_pt_raw <- read.csv("all_data.csv", header = TRUE, sep = ",")

# a positive control matrix -- rank of TPM and Counts (both RNASeq)
rank_tpm_count <- cbind(ish_pt_raw$TPM_Rank, ish_pt_raw$Count_Rank)
mu <- mean(ish_pt_raw$TPM_Rank)
sigma <- var(ish_pt_raw$TPM_Rank)
rho <- 1
p <- 1
eps <- 0.001 # default according to idr package documentation
max.ite <- 30 # default according to idr package documentation

rank_tpm <- ish_pt_raw$TPM_Rank
rank_count <- ish_pt_raw$Count_Rank
rank_exp <- ish_pt_raw$Exp_Rank

t <- seq(0.01, 0.99, by = 1/28)

tpm_count <- get.correspondence(rank_tpm, rank_count, t)
tpm_exp <- get.correspondence(rank_tpm, rank_exp, t)
count_exp <- get.correspondence(rank_count, rank_exp, t)

# plot correspondence curve on the scale of percentage for tpm and count
jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/tpm_count.jpg")
plot(tpm_count$psi$t, tpm_count$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_count$psi$t)),
     ylim=c(0, max(tpm_count$psi$value)), cex.lab=2)
lines(tpm_count$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

# plot correspondence curve on the scale of percentage for tpm and exp
jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/tpm_exp.jpg")
plot(tpm_exp$psi$t, tpm_exp$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_exp$psi$t)),
     ylim=c(0, max(tpm_exp$psi$value)), cex.lab=2)
lines(tpm_exp$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

# plot correspondence curve on the scale of percentage for count and exp
jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/count_exp.jpg")
plot(count_exp$psi$t, count_exp$psi$value, xlab="t", ylab="psi", xlim=c(0, max(count_exp$psi$t)),
     ylim=c(0, max(count_exp$psi$value)), cex.lab=2)
lines(count_exp$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()
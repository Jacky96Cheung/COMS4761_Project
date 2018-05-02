# Run IDR on Unionized Allen Brain Atlas and Homogenate RNA Seq gene expression data

library(dplyr)
library(idr)
# get to the correct directory 
setwd("/Users/elizabethboylesobolik/Desktop/COMS4761_Project/ISH_PTrap")

# load the data 
e <- exp(1)
ish_pt_raw <- read.csv("all_data.csv", header = TRUE, sep = ",")

#######################################################
#################### ESTIMATE IDR  ####################
#######################################################
# total_count <- sum(ish_pt_raw$Avg_Count)
# total_exp <- sum(ish_pt_raw$Exp_Level)
# total_tpm <- sum(ish_pt_raw$Avg_TPM)
# 
# norm_tpm <- ish_pt_raw$Avg_TPM/total_tpm
# norm_exp <- ish_pt_raw$Exp_Level/total_exp
# norm_count <- ish_pt_raw$Avg_Count/total_count
# 
# tpm <- log(1+ e^(ish_pt_raw$Avg_TPM))
# exp <- log(1 + e^(norm_exp))  
# count <- log(1+ e^(ish_pt_raw$Avg_Count))
# 
# tpm_count <- data.frame("tpm" = tpm, "count" = count)
# mu <- 10
# sigma <- 0.1
# rho <- 0.8
# p <- 0.7
# eps <- 0.001 # default according to idr package documentation
# maxit <- 30 # default according to idr package documentation
# idr.out <- est.IDR(x = tpm_count, mu, sigma, rho, p, eps, max.ite = 10)
# names(idr.out)
############################################################
#################### GET CORRESPONDENCE ####################
############################################################
rank_tpm <- ish_pt_raw$TPM_Rank
rank_count <- ish_pt_raw$Count_Rank
rank_exp <- ish_pt_raw$Exp_Rank

t <- seq(0.01, 0.99, by = 1/28)

tpm_count <- get.correspondence(rank_tpm, rank_count, t)
tpm_exp <- get.correspondence(rank_tpm, rank_exp, t)
count_exp <- get.correspondence(rank_count, rank_exp, t)
ish_pt_raw <- filter(ish_pt_raw, Exp_Level > 0, Avg_TPM > 1)
jpeg("/Users/elizabethboylesobolik/Desktop/tpm_exp.jpg")
scatter.smooth(x = log(ish_pt_raw$Exp_Level), y = log(ish_pt_raw$Avg_TPM), main = "Exp_Level ~ Avg_TPM")
dev.off()

jpeg("/Users/elizabethboylesobolik/Desktop/count_exp.jpg")
scatter.smooth(x = log(ish_pt_raw$Exp_Level), y = log(ish_pt_raw$Avg_Count), main = "Exp_Level ~ Avg_Count")
dev.off()

# plot correspondence curve on the scale of percentage for tpm and count
# jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/tpm_count.jpg")
# plot(tpm_count$psi$t, tpm_count$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_count$psi$t)),
#      ylim=c(0, max(tpm_count$psi$value)), cex.lab=2)
# lines(tpm_count$psi$smoothed.line, lwd=4)
# abline(coef=c(0,1), lty=3)
# dev.off()

# plot correspondence curve on the scale of percentage for tpm and exp
# jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/tpm_exp.jpg")
# plot(tpm_exp$psi$t, tpm_exp$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_exp$psi$t)),
#      ylim=c(0, max(tpm_exp$psi$value)), cex.lab=2)
# lines(tpm_exp$psi$smoothed.line, lwd=4)
# abline(coef=c(0,1), lty=3)
# dev.off()

# plot correspondence curve on the scale of percentage for count and exp
# jpeg("/Users/elizabethboylesobolik/Desktop/IDR_test/count_exp.jpg")
# plot(count_exp$psi$t, count_exp$psi$value, xlab="t", ylab="psi", xlim=c(0, max(count_exp$psi$t)),
#      ylim=c(0, max(count_exp$psi$value)), cex.lab=2)
# lines(count_exp$psi$smoothed.line, lwd=4)
# abline(coef=c(0,1), lty=3)
# dev.off()
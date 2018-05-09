# IDR for Control and IP PhosphoTrap

library(dplyr)
library(idr)

# get to the correct directory
d <- "/Users/elizabethboylesobolik/Desktop"
setwd(d)

cr_raw <- read.csv("tpma.csv", header = TRUE, sep = ",")

# load the rank data from the raw input file
rank_a <- cr_raw$rank_A
rank_b <- cr_raw$rank_B
rank_c <- cr_raw$rank_C

# specify the different percentage thresholds
t <- seq(0.01, 0.99, by = 1/100)

# Get Correspondence!
tpm_ab <- get.correspondence(rank_a, rank_b, t)
tpm_ac <- get.correspondence(rank_a, rank_c, t)
tpm_bc <- get.correspondence(rank_b, rank_c, t)

# Plot the percentage threshold vs. psi value for each condition
ab <- "/tpm_ab.jpg"
ac <- "/tpm_ab.jpg"
bc <- "/tpm_ab.jpg"

filename <- paste(d, ab, sep = "")
jpeg(filename)
plot(tpm_ab$psi$t, tpm_ab$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_ab$psi$t)),
     ylim=c(0, max(tpm_ab$psi$value)), cex.lab=2)
lines(tpm_ab$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

filename <- paste(d, ac, sep = "")
jpeg(filename)
plot(tpm_ac$psi$t, tpm_ac$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_ac$psi$t)),
     ylim=c(0, max(tpm_ac$psi$value)), cex.lab=2)
lines(tpm_ac$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

filename <- paste(d, bc, sep = "")
jpeg(filename)
plot(tpm_bc$psi$t, tpm_bc$psi$value, xlab="t", ylab="psi", xlim=c(0, max(tpm_bc$psi$t)),
     ylim=c(0, max(tpm_bc$psi$value)), cex.lab=2)
lines(tpm_bc$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
dev.off()

########################################################################
# Chapter 7
(n_wd <- getwd())               # Confirmation of working directory
source('myfunc/myfunc.R')       # Loading self-made functions
library(cmdstanr)            # load and attach add-on package cmdstanr
library(posterior)           # load and attach add-on package posterior

# Table 6.1 English Grades of Experimental and Control Groups
# Experimental group
x1 <- c(49,66,69,55,54,72,51,76,40,62,66,51,59,68,66,57,53,66,58,57)
# Control group
x2 <- c(41,55,21,49,53,50,52,67,54,69,57,48,31,52,56,50,46,38,62,59)
print(x1); print(x2)


# Execution of MCMC
out0601 <- G2Ind(x1, x2, EQU = 0, prior = F)
#save(out0601,  file="./script/obje/out0601")
#load(file="./script/obje/out0601"); # Load the results of a pre-run of MCMC

gqcal(out0601$ext[,1:6])

# Extracting values as MCMC sample names are long.
mu1    <- out0601$mu1;   mu2    <- out0601$mu2;
sigma1 <- out0601$sigma1; sigma2 <- out0601$sigma2;
xaste1 <- out0601$xaste1; xaste2 <- out0601$xaste2;
dmu    <- mu1 - mu2;  
dsigma <- sigma1 - sigma2; 

# 7.1.1 Within-group Standard Deviation
sigma_within <- sqrt((sigma1^2 + sigma2^2) / 2);  # Equation (7.1)

# Table 7.1 Numerical Summary of the Posterior Distribution of the Standardized Mean Difference
delta  <- dmu / sigma_within;        # Equation (7.2), standardized mean difference
delta1 <- dmu / sigma1;              # Equation (7.3)
delta2 <- dmu / sigma2;              # Equation (7.3)
gqcal(cbind(sigma_within, delta, delta1, delta2))


# Figure 7.3 Posterior Distribution of Within-Group Standard Deviation - Checking the Positional Relationship of the Posterior Distribution of Delta
par(mfrow = c(1, 2))
# Left figure
hist(sigma_within, breaks = 100, col = gray(1.0), cex.axis = 2.0, cex.lab = 2.0,
     main = "", xlab = "Within-Group Standard Deviation", ylab = "") 
# Right figure
plot(density(delta), xlim = c(-1, 3), ylim = c(0, 1.3), cex.axis = 2.0, cex.lab = 2.0,
     main = "", xlab = "", ylab = "", lwd = 2.0) 
par(new = T)
plot(density(delta1), xlim = c(-1, 3), ylim = c(0, 1.3), cex.axis = 2.0, cex.lab = 2.0,
     main = "", xlab = "", ylab = "", lwd = 2.0) 
par(new = T)
plot(density(delta2), xlim = c(-1, 3), ylim = c(0, 1.3), cex.axis = 2.0, cex.lab = 2.0,
     main = "", xlab = "Standardized Mean Difference", ylab = "", lwd = 2.0) 
text(0.2, 1.0, "delta2", cex = 2.0)
text(1.25, 1.1, "delta", cex = 2.0)
text(1.9, 0.3, "delta1", cex = 2.0)
par(mfrow = c(1, 1))
#dev.copy2pdf(file="../English_tex/fig/chap07/fig07_03.pdf")


# 7.1.8 Consideration using PHC
# Table 7.2 PHC Table of Standardized Mean Difference
PHC01(seq(0.0, 0.5, 0.1), delta, 0, cc = "gtc", byoga = "no"); # Equation (7.2)
PHC01(seq(0.0, 0.5, 0.1), delta1, 0, cc = "gtc", byoga = "no"); # Equation (7.3)
PHC01(seq(0.0, 0.5, 0.1), delta2, 0, cc = "gtc", byoga = "no"); # Equation (7.3)
# ROPE
PHC01(seq(0.0, 0.5, 0.1), delta, 0, cc = "rope", byoga = "no"); #
PHC01(seq(0.0, 0.5, 0.1), delta1, 0, cc = "rope", byoga = "no"); #
PHC01(seq(0.0, 0.5, 0.1), delta2, 0, cc = "rope", byoga = "no"); #

# Figure 7.4 PHC Curves of Standardized Mean Difference
par(mfrow = c(2, 1))
PHC01(seq(0.0, 2.5, 0.1), delta, 0, cc = "gtc", byoga = "yes", lwd = 2.0); par(new = T)
PHC01(seq(0.0, 2.5, 0.1), delta1, 0, cc = "gtc", byoga = "yes", lwd = 2.0, lty = 3); par(new = T)
PHC01(seq(0.0, 2.5, 0.1), delta2, 0, cc = "gtc", byoga = "yes", lwd = 2.0, lty = 3);
text(1.2, 0.5, "delta1", cex = 2.0)
text(0.6, 0.5, "delta2", cex = 2.0)
# ROPE
PHC01(seq(0.0, 2.5, 0.1), delta, 0, cc = "rope", byoga = "yes", lwd = 2.0); par(new = T)
PHC01(seq(0.0, 2.5, 0.1), delta1, 0, cc = "rope", byoga = "yes", lwd = 2.0, lty = 3); par(new = T)
PHC01(seq(0.0, 2.5, 0.1), delta2, 0, cc = "rope", byoga = "yes", lwd = 2.0, lty = 3);
text(1.2, 0.5, "delta1", cex = 2.0)
text(0.6, 0.5, "delta2", cex = 2.0)
par(mfrow = c(1, 1))
#dev.copy2pdf(file = "../English_tex/fig/chap07/fig07_04.pdf")


# Table 7.3: Estimation Results of Posterior Distribution of Measure of Nonoverlap
U1   <- pnorm(mu1, mu2, sigma_within); #(7.6) equation, measure of nonoverlap, denominator is within-group standard deviation
U2   <- pnorm(mu2, mu1, sigma_within); #(7.7) equation, measure of nonoverlap, denominator is within-group standard deviation
U1a  <- pnorm(mu1, mu2, sigma2);       #(7.9) equation, measure of nonoverlap, denominator is standard deviation of group 2
U2a  <- pnorm(mu2, mu1, sigma1);       #(7.10) equation, measure of nonoverlap, denominator is standard deviation of group 1
gqcal(cbind(U1, U2, U1a, U2a))


# Figure 7.7: Comparison of Positional Relationships of Posterior Distribution of Measure of Nonoverlap
plot(density(U1), xlim=c(0, 1), ylim=c(0, 5), cex.axis=2.0, cex.lab=2.0, main="", xlab="", ylab="", lwd=2.0); par(new=T)
plot(density(U1a), xlim=c(0, 1), ylim=c(0, 5), cex.axis=2.0, cex.lab=2.0, main="", xlab="", ylab="", lwd=2.0); par(new=T)
plot(density(U2), xlim=c(0, 1), ylim=c(0, 5), cex.axis=2.0, cex.lab=2.0, main="", xlab="", ylab="", lwd=2.0); par(new=T)
plot(density(U2a), xlim=c(0, 1), ylim=c(0, 5), cex.axis=2.0, cex.lab=2.0, main="", xlab="", ylab="", lwd=2.0) 
text(0.92, 3.5, "U1", cex=2.0)
text(0.68, 3.5, "U1*", cex=2.0)
text(0.3, 3.5, "U2", cex=2.0)
text(0.03, 3.5, "U2*", cex=2.0)
#dev.copy2pdf(file="../English_tex/fig/chap07/fig07_07.pdf")


# Table 7.4: PHC Table of Measure of Nonoverlap / PHC Table of ROPE, Effectively 0.5
PHC01(seq(0.5, 0.75, 0.05), U1, 0, cc="gtc", byoga="no"); #
PHC01(seq(0.5, 0.75, 0.05), U1a, 0, cc="gtc", byoga="no"); #
PHC01(seq(0.25, 0.5, 0.05), U2, 0, cc="ltc", byoga="no"); #
PHC01(seq(0.25, 0.5, 0.05), U2a, 0, cc="ltc", byoga="no"); #
# ROPE
PHC01(seq(0.0, 0.25, 0.05), U1, 0.5, cc="rope", byoga="no"); #
PHC01(seq(0.0, 0.25, 0.05), U1a, 0.5, cc="rope", byoga="no"); #
PHC01(seq(0.0, 0.25, 0.05), U2, 0.5, cc="rope", byoga="no"); #
PHC01(seq(0.0, 0.25, 0.05), U2a, 0.5, cc="rope", byoga="no"); #


# Figure 7.8: PHC Curves of Measure of Nonoverlap / PHC Curves of ROPE, Effectively 0.5
par(mfrow=c(3, 1))
# Upper row
PHC01(seq(0.5, 1.0, 0.01), U1, 0, cc="gtc", byoga="yes", lwd=2.0); par(new=T)
PHC01(seq(0.5, 1.0, 0.01), U1a, 0, cc="gtc", byoga="yes", lwd=2.0);
text(0.85, 0.6, "U1", cex=2.0)
text(0.70, 0.6, "U1*", cex=2.0)
# Middle row
PHC01(seq(0.0, 0.5, 0.01), U2, 0, cc="ltc", byoga="yes", lwd=2.0); par(new=T)
PHC01(seq(0.0, 0.5, 0.01), U2a, 0, cc="ltc", byoga="yes", lwd=2.0);
text(0.3, 0.6, "U2", cex=2.0)
text(0.1, 0.6, "U2*", cex=2.0)
# Lower row
PHC01(seq(0.0, 0.5, 0.01), U1, 0.5, cc="rope", byoga="yes", lwd=2.0); par(new=T)
PHC01(seq(0.0, 0.5, 0.01), U2, 0.5, cc="rope", byoga="yes", lwd=2.0); par(new=T)
PHC01(seq(0.0, 0.5, 0.01), U1a, 0.5, cc="rope", byoga="yes", lwd=2.0); par(new=T)
PHC01(seq(0.0, 0.5, 0.01), U2a, 0.5, cc="rope", byoga="yes", lwd=2.0);
text(0.33, 0.6, "U1 U2", cex=2.0)
text(0.20, 0.6, "U1*", cex=2.0)
text(0.45, 0.6, "U2*", cex=2.0)
par(mfrow=c(1, 1))
#dev.copy2pdf(file="../English_tex/fig/chap07/fig07_08.pdf")


# Figure 7.9: Visualization of Probability Beyond Threshold
nob <- 1000
plot(xaste1[1:nob], xaste2[1:nob], xlim=c(0, 100), ylim=c(0, 100), cex.axis=2.0,
     cex.lab=1.3, main="", xlab="Experimental Group", ylab="Control Group", lwd=2.0) 
abline(0, 1, lwd=2, lty=1)
abline(-3, 1, lwd=2, lty=2)
abline(-5, 1, lwd=2, lty=3)
abline(-10, 1, lwd=2, lty=4)
arrows(25, 25, 35, 10, lwd=4)
text(39, 6, "72%", cex=4)
#dev.copy2pdf(file="../English_tex/fig/chap07/fig07_09.pdf")

# ikijyori2( ): Function Returning Curves and Tables for Probability Beyond Threshold
# seq01: Vector (or scalar) of reference points c
# out  : MCMC object containing random vectors of mu1, mu2, sigma1, sigma2
# If byoga="yes", it outputs the curve of probability beyond threshold; otherwise, it outputs the table of probability beyond threshold


# Table 7.5: Numerical Summary of Posterior Distribution of Probability Beyond Threshold
ikijyori2(c(0, 3, 5, 10), out0601, byoga="no")


# Figure 7.10: EAP of Probability Beyond Threshold and Upper & Lower Limits of 95% Confidence Interval
ikijyori2(seq(-10, 30, 1), out0601, byoga="yes")
#dev.copy2pdf(file="../English_tex/fig/chap07/fig07_10.pdf")


# Appendix: Probability that the Judgment 'The Probability of Scoring More Than 3 Points is Greater Than 0.7' is Correct
ikijyo3 <- pnorm((mu1 - mu2 - 3) / sqrt(sigma1^2 + sigma2^2), 0, 1)
PHC01(0.7, ikijyo3, 0, cc="gtc", byoga="no");     #


########################################################################
# Chapter 7 Practical Assignment

# Table 6.8: Measurement Results of Perception Time (seconds) under Auditory Conditions
# Data input for auditory condition
y1 <- c(32.30, 34.24, 28.10, 33.40, 37.71, 31.62, 31.37, 35.85, 32.33, 34.04,
        34.96, 31.43, 35.28, 30.19, 35.09, 33.38, 31.49, 28.44, 32.12, 31.81)

# Data input for control condition
y2 <- c(31.43, 31.09, 33.38, 30.49, 29.62, 35.40, 32.58, 28.96, 29.43, 28.52,
        25.39, 32.68, 30.51, 30.15, 32.33, 30.43, 32.50, 32.07, 32.35, 31.57)

out0602 <- G2Ind(y1, y2, EQU=0, prior=F)
#save(out0602,  file="./script/obje/out0602")
#load(file="./script/obje/out0602"); # Load the results of a pre-run of MCMC
gqcal(out0602$ext[,1:6])

# Other script examples are omitted

########################################################################







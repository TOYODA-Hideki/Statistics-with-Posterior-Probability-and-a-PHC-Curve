########################################################################
# Chapter 10
(n_wd <- getwd())              # Check working directory
source('myfunc/myfunc.R')      # Load custom functions
library(cmdstanr)           # Load cmdstanr package
library(posterior)          # Load posterior package

# Table 10.2 Time and Number of Deviations for the 15th Trial of the Mirror Drawing Experiment by Group
(mirror_drawing01 <- read.csv("dat/mirror_tracing01.csv", header = TRUE))
(x1 <- mirror_drawing01[mirror_drawing01$group == 1, 1])
(x2 <- mirror_drawing01[mirror_drawing01$group == 2, 1])
(x3 <- mirror_drawing01[mirror_drawing01$group == 3, 1])
n1 <- length(x1); n2 <- length(x2); n3 <- length(x3);

# Figure 10.2 Box plot showing the distribution differences among the three groups
y <- c(x1, x2, x3)
group <- c(rep("1", n1), rep("2", n2), rep("3", n3))
boxplot(y ~ group, cex.axis = 2.0, cex.lab = 2.0, ylab = "", xlab = "Group")
#dev.copy2eps(file = "z0503.eps", family = "Japan1")

# Table 10.3 Mean, Variance, and Standard Deviation
round(apply(cbind(y, x1, x2, x3), 2, mean), 2)
round(apply(cbind(y, x1, x2, x3), 2, var), 2)
round(sqrt(apply(cbind(y, x1, x2, x3), 2, var)), 2)

# Inference for One-Factor Experiments
A <- c(rep(1, n1), rep(2, n2), rep(3, n3))  # Specifying levels
out1001 <- E1Ind(y, A, prior = T, mL = 0, mH = 100, sL = 0, sH = 50)
#save(out1001, file = "./script/obje/out1001")
#load(file = "./script/obje/out1001"); # Load pre-run MCMC results

# Extracting Random Variables
muA <- out1001$muA; sigmaA <- out1001$sigmaA; sigmaE <- out1001$sigmaE;
mu <- out1001$mu; aj <- out1001$aj; U2 <- out1001$U2; eta2 <- out1001$eta2

# Table 10.4 Estimation Results of Parameters
gqcal(cbind(muA, sigmaE), 1)

# Table 10.5 Estimation Results of Total Mean and Level Effects
gqcal(cbind(mu, aj), 1)

# Figure 10.3 Boxplot of the Posterior Distribution of Level Effects
size<-nrow(mu);             # Number of random numbers
yy <- aj[, 1]; for (i in 2:3) {yy <- c(yy, aj[, i])};
boxplot(yy~rep(c("Control","Non-Dominant Hand","Dominant Hand"),each=size),
         ylab = "", xlab = "", cex.axis = 2.0, cex.lab = 2.0)
abline(h = 0.0, lwd = 1.0, lty = 3)
#dev.copy2eps(file = "z0501.eps", family = "Japan1")

# 10.3.2 Decomposition of Variance in Data
(sy2 <- var(mirror_drawing01$time));                         # Total variance
(tm <- mean(mirror_drawing01$time));                         # Total mean

# Factor Variance / Between-Group Variance               Equation (10.16)
(sa2 <- ((mean(x1) - tm)^2 + (mean(x2) - tm)^2 + (mean(x3) - tm)^2) / 3); 

# Error Variance / Within-Group Variance                 Equations (10.17) & (10.18)
(se2 <- (var(x1) + var(x2) + var(x3)) / 3);

# Check that the sum of Factor Variance and Error Variance equals Total Variance  Equations (10.19) & (10.20)
sa2 + se2;  

# Error Variance / Within-Group Variance can also be calculated as follows
(se2 <- (sum((x1 - mean(x1))^2) + sum((x2 - mean(x2))^2) + sum((x3 - mean(x3))^2)) / 90)

# 10.4 Consideration of Factors and Levels
eta2 <- sigmaA^2 / (sigmaA^2 + sigmaE^2);       # Variance Explained / Coefficient of Determination

# Table 10.6 Numerical Summary and PHC Table for Between-Group Standard Deviation
gqcal(sigmaA)
PHC01(seq(1.0, 4.0, 0.5), sigmaA, 0, cc = "gtc", byoga = "no");
PHC01(seq(1.0, 4.0, 0.5), sigmaA, 0, cc = "rope", byoga = "no"); #rope

# Table 10.7 Numerical Summary and PHC Table for Variance Explained
gqcal(eta2)
PHC01(seq(0, 0.06, 0.01), eta2, 0, cc = "gtc", byoga = "no");
PHC01(seq(0, 0.06, 0.01), eta2, 0, cc = "rope", byoga = "no"); #rope

# Figure 10.4 PHC Curves for Between-Group Standard Deviation and Variance Explained
par(mfrow = c(2, 1))
PHC01(seq(1.0, 8.0, 0.05), sigmaA, 0, cc = "gtc", byoga = "yes", xlab = "Between-Group Standard Deviation");
PHC01(seq(0, 0.25, 0.01), eta2, 0, cc = "gtc", byoga = "yes", xlab = "Proportion of Variance Explained");
par(mfrow = c(1, 1))
#dev.copy2eps(file = "z1004as.eps", family = "Japan1")


# PHC02(c=0, ext, cc="gtc", digits=3) Matrix of probabilities that there is a difference greater than c between multiple levels
# c   Scalar, reference point
# ext A matrix with random variables in rows and parameters in columns
# cc="gtc" means row parameter - column parameter > c,
#  "rope" means abs(row parameter - column parameter) < c,
#  otherwise, it means row parameter - column parameter < c
# Table 10.8 Probability that level j in the row is greater than level j′ in the column by more than c
PHC02(c = 0.0, muA)
PHC02(c = 1.0, muA)
PHC02(c = 2.0, muA)

# Table 10.9 Analysis that provides insights into the difference between Group 1 and Group 2 (Hypothesis A)
E1between(out1001, 2, 1, 2, cr1 = c(0, 0.5, 1.0))

# Table 10.10 Analysis that provides insights into the difference between Group 2 and Group 3 (Hypothesis B)
E1between(out1001, 2, 2, 3, cr1 = c(0, 0.5, 1.0))

# conjunctive(): Probability that two research hypotheses are simultaneously true under the reference vectors c1, c2
# x: A matrix[, a] composed of random variables in columns a (not 0/1 but random numbers)
# IJ: A 2*2 matrix, evaluating the probability that IJ[1, 1] - IJ[1, 2] and IJ[2, 1] - IJ[2, 2] are both true
# c1, c2: Vectors of each reference point
# Table 10.11 Probability that the conjunctive proposition related to Hypothesis C is true
conjunctive(muA, rbind(c(1, 2), c(2, 3)), c1 = seq(0, 2.5, 0.5), c2 = seq(0, 2.5, 0.5), 3)

### One-Way ANOVA F-test
print(summary(aov(y ~ A)))

########################################################################
# Chapter 10 10.6 Practical Assignment

# Inputting "Perception Time" Data for Table 10.12
# Control
x1 <- c(31.43, 31.09, 33.38, 30.49, 29.62, 35.40, 32.58, 28.96, 29.43, 28.52,
        25.39, 32.68, 30.51, 30.15, 32.33, 30.43, 32.50, 32.07, 32.35, 31.57)
# Auditory
x2 <- c(32.30, 34.24, 28.10, 33.40, 37.71, 31.62, 31.37, 35.85, 32.33, 34.04,
        34.96, 31.43, 35.28, 30.19, 35.09, 33.38, 31.49, 28.44, 32.12, 31.81)
# Reading Aloud
x3 <- c(31.62, 37.04, 33.76, 30.01, 34.18, 33.08, 28.77, 33.90, 28.06, 37.54,
        33.89, 32.23, 35.95, 36.68, 33.57, 30.87, 32.20, 29.98, 33.08, 35.12)
n1 <- length(x1); n2 <- length(x2); n3 <- length(x3);
y <- c(x1, x2, x3)
A <- c(rep(1, n1), rep(2, n2), rep(3, n3))  # Specifying levels
out1002 <- E1Ind(y, A, prior = T, mL = 0, mH = 100, sL = 0, sH = 50)
#save(out1002, file = "./script/obje/out1002")
#load(file = "./script/obje/out1002"); # Load pre-run MCMC results

# The scripts for the following examples are omitted.

########################################################################
# Self-Study Chapter 10 10.7 Practical Assignment
# "Mouse Weight Gain Data (g)" in Table 10.13
# 1) Report the order of conditions in which weight gain is most likely, based on point estimates.
# 2) What percentage of the difference in weight gain is explained by the difference in breeding conditions?
# 3) What is the probability that there is a difference of 0g, 2g, 4g or more in the population mean?
# 4) Describe in detail the difference between the largest and smallest groups. Calculate the threshold rate for 0g, 2g, 4g.
# 5) What is the probability that the conjunctive proposition that there is a difference in the population mean according to the order of point estimates is true?

x1 <- c(05.02, 06.67, 08.17, 02.79, 08.13, 06.34, 06.32, 03.97)
x2 <- c(09.89, 09.58, 11.20, 09.05, 12.33, 09.39, 10.88, 09.37, 17.40)
x3 <- c(10.20, 07.29, 07.57, 03.42, 05.82, 10.92, 05.21, 13.47, 08.64, 06.05)
y <- c(x1, x2, x3)
A <- c(rep(1, 8), rep(2, 9), rep(3, 10))
out1003 <- E1Ind(y, A, prior = T, mL = 0, mH = 50, sL = 0, sH = 50)
#save(out1003, file = "./script/obje/out1003")
#load(file = "./script/obje/out1003"); # Load pre-run MCMC results

# 1) 2)
gqcal(out1003$ext[,1:17])

# 3)
# Table 10.8 Probability that level j in the row is greater than level j′ in the column by more than c
PHC02(c = 0.0, out1003$muA)
PHC02(c = 2.0, out1003$muA)
PHC02(c = 4.0, out1003$muA)

# 4)
E1between(out1003, 2, 2, 1, cr1 = c(0.0, 2.0, 4.0))

# 5)
conjunctive(out1003$muA, rbind(c(2, 3), c(3, 1)),
            c1 = seq(0, 0.5, 0.1), c2 = seq(0, 0.5, 0.1), 2)

########################################################################


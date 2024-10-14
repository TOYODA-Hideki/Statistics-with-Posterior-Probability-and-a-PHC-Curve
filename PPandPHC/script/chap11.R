########################################################################
# Chapter 11
(n_wd <- getwd())              # Check working directory
source('myfunc/myfunc.R')      # Load custom functions
library(cmdstanr)           # Load cmdstanr package
library(posterior)          # Load posterior package

# Table 11.1 Data of Subject S1 in the "Psychological Perception Experiment" with two factors (unit: seconds)
y1 <- c(
  28.10, 31.67, 33.60, 25.37, 32.03, 32.21, 29.62, 29.69, 29.65, 28.86,
  31.35, 30.10, 30.64, 32.65, 34.80, 32.24, 31.27, 30.31, 30.49, 38.30,
  33.32, 31.82, 31.94, 34.10, 31.34, 29.52, 34.38, 30.54, 32.96, 30.75,
  33.95, 30.16, 29.60, 30.10, 27.39, 28.50, 26.07, 28.08, 30.59, 30.18,
  34.48, 28.44, 28.95, 30.33, 28.61, 28.68, 28.34, 28.00, 29.74, 29.81,
  28.17, 31.10, 29.84, 30.06, 32.11, 33.85, 36.45, 30.64, 36.35, 29.72,
  29.58, 34.12, 27.92, 26.39, 27.98, 32.27, 25.27, 31.28, 31.19, 28.81,
  31.92, 30.81, 31.02, 33.46, 36.87, 31.53, 28.56, 30.16, 32.52, 31.14)
A <- rep(1:4, each = 10, times = 2) # Specifying levels of Factor A
B <- rep(1:2, each = 40)           # Specifying levels of Factor B

# Two-Factor Experiment Inference Analysis of Subject S1
out1101 <- E2Ind(y1, A, B, prior = T, mL = 0, mH = 100, sL = 0, sH = 50)
#save(out1101, file = "./script/obje/out1101")
#load(file = "./script/obje/out1101"); # Load pre-run MCMC results

############### 11.2 Inference Example 1 (Analysis of Interaction) #######

# Table 11.2 Mean and Standard Deviation of Subject S1's Data
# Figure 11.1 Mean Value Plot of Subject S1's Data
# mean_sd_plot() is a function to plot cell means, standard deviations, and means for a two-factor experiment
# y: Vector of measurement values
# A: Specification of levels (more than 2 levels)
# B: Specification of levels (2 levels)
mean_sd_plot(y1, A, B, ylab = "Perceived 30 Seconds", xlab = "Control・Auditory・Reading Aloud・Exercise")

# Table 11.3 Estimated Results of Parameters (Subject S1)
gqcal(out1101$mu,        2);    # Equation (11.1)
gqcal(out1101$muA,       2);    # Equation (11.2)
gqcal(out1101$muB[,1],   2);    # Equation (11.3)
gqcal(out1101$muAB[,1:4],2);    # Equation (11.4)

# Table 11.4 Estimated Results of Generated Quantities Related to the Magnitude of Effects
gqcal(out1101$sigmaA, 3);  # Equation (11.15)
gqcal(out1101$sigmaB, 3);  # Equation (11.15)
gqcal(out1101$sigmaAB, 3); # Equation (11.15)
gqcal(out1101$sigmaE, 2);  # Equation (11.15)
gqcal(out1101$eta2A,  3);  # Explained Variance of the Main Effect of Factor A
gqcal(out1101$eta2B,  3);  # Explained Variance of the Main Effect of Factor B
gqcal(out1101$eta2AB, 3);#Explained Variance of the Interaction of Factors A&B
gqcal(out1101$eta2T,  3);  # Explained Variance of All Factors

# Table 11.5 Estimated Results of Cell Means (Order of output differs from the textbook) Equations (11.8) and (11.9)
gqcal(out1101$cellmean,2)

# PHC02(c=0, ext, cc="gtc", digits=3) Matrix of probabilities that there is a difference greater than c between multiple levels
# c   Scalar, reference point
# ext A matrix with random variables in rows and parameters in columns
# cc="gtc" means row parameter - column parameter > c,
#  "rope" means abs(row parameter - column parameter) < c,
#  and otherwise, it means row parameter - column parameter < c
# Table 11.6 Checking Necessary Conditions for a Substantial Difference  PHC(μjk ? μj′k′ > 0)
PHC02(c=0.0,out1101$cellmean,cc="gtc",3)

# E2between() Comparison of two cells (threshold rate reference points specified as a vector)
# x: An object of class 'E2Ind'
# I=1: Integer, cell with the higher mean value
# J=2: Integer, cell with the lower mean value 
# cr1: Vector specifying reference points for threshold rate
# digits=3 : Rounding of decimals
# Table 11.7 Estimated Results of the Difference between 'Auditory Condition' and 'Exercise Condition' under 'Normal Condition'
colnames(out1101$cellmean)
E2between(out1101,I=2,J=4,cr1=c(0,1,2),digits=3)


############ 11.4  ###############################

# Table 11.8 Data of Subject S2 in the "Psychological Perception Experiment" with two factors (unit: seconds)
y2 <- c(
  33.00, 31.92, 27.85, 30.99, 31.61, 30.50, 31.06, 29.48, 33.08, 32.55,
  34.89, 32.33, 31.98, 34.22, 35.48, 27.54, 35.53, 28.57, 29.57, 33.94,
  36.30, 30.51, 33.70, 28.77, 32.60, 35.11, 29.35, 31.99, 33.22, 33.04,
  27.34, 33.40, 24.63, 27.24, 28.27, 31.85, 25.05, 29.04, 31.29, 30.20,
  29.38, 29.80, 31.28, 31.54, 31.83, 30.50, 27.22, 30.06, 29.61, 35.33,
  35.69, 35.02, 31.94, 33.98, 30.38, 26.21, 34.91, 33.36, 31.37, 33.49,
  30.00, 31.23, 30.09, 34.73, 35.43, 33.44, 34.65, 33.35, 30.41, 33.41,
  26.35, 28.61, 29.52, 27.63, 27.05, 33.32, 29.85, 34.44, 28.52, 32.40)

# Two-Factor Experiment Inference Analysis of Subject S2
out1102 <- E2Ind(y2, A, B, prior = T, mL = 0, mH = 100, sL = 0, sH = 50)
#save(out1102, file = "./script/obje/out1102")
#load(file = "./script/obje/out1102"); # Load pre-run MCMC results

# Mean and Standard Deviation of Subject S2's Data
# Figure 11.2 Mean Value Plot of Subject S2's Data
mean_sd_plot(y2, A, B, ylab = "Perceived 30 Seconds", xlab = "Control・Auditory・Reading Aloud・Exercise")

# Table 11.9 Estimated Results of Parameters (Subject S2)
gqcal(out1102$mu,        2);    # Equation (11.1)
gqcal(out1102$muA,       2);    # Equation (11.2)
gqcal(out1102$muB[,1],   2);    # Equation (11.3)
gqcal(out1102$muAB[,1:4],2);    # Equation (11.4)

# Table 11.10 Estimated Results of Generated Quantities Related to the Magnitude of Effects (Subject S2)
gqcal(out1102$eta2A,  3);  # Explained Variance of the Main Effect of Factor A
gqcal(out1102$eta2B,  3);  # Explained Variance of the Main Effect of Factor B
gqcal(out1102$eta2AB, 3);#Explained Variance of the Interaction of Factors A&B
gqcal(out1102$eta2T,  3);  # Explained Variance of All Factors

# Table 11.11 PHC(aj - aj′ > c)  Identifying differences between levels of Factor A
PHC02(c = 0.0, out1102$muA, cc = "gtc", 3)
PHC02(c = 1.0, out1102$muA, cc = "gtc", 3)

# Table 11.12 PHC Table for Differences between Control Condition and Other Conditions
PHC01(seq(0, 1.2, 0.2),out1102$muA[,2],out1102$muA[,1],cc="gtc",byoga = "no")
PHC01(seq(0, 1.2, 0.2),out1102$muA[,3],out1102$muA[,1],cc="gtc",byoga = "no")
PHC01(seq(0, 1.2, 0.2),out1102$muA[,1],out1102$muA[,4],cc="gtc",byoga = "no")

############### 11.5 Inference Example 3 (Analysis of Main Effects of Two Levels) ####################

# Table 11.13 Data of Subject S3 in the "Psychological Perception Experiment" with two factors (unit: seconds)
y3 <- c(
  35.24, 31.29, 31.00, 32.23, 30.13, 26.27, 30.71, 29.88, 31.11, 26.99,
  28.73, 31.34, 32.09, 33.68, 31.60, 32.17, 29.54, 28.75, 29.57, 26.70,
  29.59, 30.51, 32.05, 28.99, 28.10, 31.29, 28.48, 30.80, 33.88, 33.52,
  34.34, 29.29, 30.00, 33.31, 27.98, 29.23, 32.34, 30.54, 31.82, 31.85,
  27.88, 30.26, 29.96, 24.01, 29.48, 34.10, 30.69, 27.32, 26.67, 33.97,
  27.20, 30.14, 31.62, 29.46, 29.19, 29.18, 31.26, 28.61, 26.66, 25.21,
  28.52, 29.13, 26.07, 30.27, 27.26, 26.49, 30.62, 29.11, 29.28, 30.80,
  32.03, 32.48, 28.64, 24.34, 29.19, 25.09, 35.23, 24.04, 30.39, 31.51)

# Two-Factor Experiment Inference Analysis of Subject S3
out1103 <- E2Ind(y3, A, B, prior = T, mL = 0, mH = 100, sL = 0, sH = 50)
#save(out1103, file = "./script/obje/out1103")
#load(file = "./script/obje/out1103"); # Load pre-run MCMC results

# Mean and Standard Deviation of Subject S3's Data
# Figure 11.3 Mean Value Plot of Subject S3's Data
mean_sd_plot(y3, A, B, ylab = "Perceived 30 Seconds", xlab = "Control・Auditory・Reading Aloud・Exercise")

# Table 11.14 Estimated Results of Parameters (Subject S3)
gqcal(out1103$mu,        2);    # Equation (11.1)
gqcal(out1103$muA,       2);    # Equation (11.2)
gqcal(out1103$muB[, 1],  2);    # Equation (11.3)
gqcal(out1103$muAB[,1:4],2);    # Equation (11.4)

# Table 11.15 Estimated Results of Generated Quantities Related to the Magnitude of Effects (Subject S3) (Numerical Summary of Explained Variance)
gqcal(out1103$eta2B,  3);

# Table 11.15 Estimated Results of Generated Quantities Related to the Magnitude of Effects (Subject S3) (PHC Table for Factor B)
PHC01(seq(0, 1.2, 0.2),out1103$muB[,1],out1103$muB[,2],cc="gtc", byoga = "no")


########################################################################

# Self-Study
# Table 5.16 Pitcher E's Pitch Type Speed With and Without Runners
A <- c(rep(1, 49), rep(2, 49))  # Factor A (presence of runners)
B <- c(rep(1, 10), rep(2, 8), rep(3, 7), rep(4, 9), rep(5, 8), rep(6, 7),
       rep(1, 10), rep(2, 8), rep(3, 7), rep(4, 9), rep(5, 8), rep(6, 7))  # Factor B (pitch type)
y <- c(140, 146, 149, 136, 147, 147, 143, 143, 143, 141,
       139, 136, 136, 140, 135, 132, 140, 134,
       123, 127, 131, 130, 138, 128, 129,
       115, 120, 118, 118, 121, 124, 129, 119, 128,
       128, 124, 123, 121, 122, 126, 131, 122,
       121, 121, 120, 116, 117, 113, 118,
       143, 141, 142, 145, 149, 145, 143, 141, 142, 155,
       138, 134, 142, 136, 135, 136, 131, 133,
       131, 128, 128, 128, 127, 130, 130,
       117, 125, 132, 122, 119, 122, 129, 117, 127,
       117, 120, 124, 122, 122, 122, 118, 122,
       119, 125, 122, 116, 119, 113, 122)
a <- length(unique(A))  # Number of levels in Factor A
b <- length(unique(B))  # Number of levels in Factor B

out1104 <- E2Ind(y, A, B, prior = T, mL = 0, mH = 1000, sL = 0, sH = 500)
#save(out1104, file = "./script/obje/out1104")
#load(file = "./script/obje/out1104"); # Load pre-run MCMC results

# 1) Determine the cell means and standard deviations for the two-factor experiment, plot the means, and interpret them.
#    Also, provide an estimate of the explanatory power of the factors.
#    Note that A (a=2) and B (b=6) are specified in the opposite way compared to the perceived time.
# 2) Interpret the point estimates for μ, a_{1}, b_{1}, ..... b_{6}, ab_{15}, ab_{25}
#    (the 5th is a slider).

mean_sd_plot(y, B, A, ylab="Pitch Speed",
             xlab="Fastball, Cutball, Forkball, Changeup, Slider, Curveball")

gqcal(out1104$ext[,1:41],3)

# 3) Interpret η_{A}^{2}, η_{B}^{2}, η_{AB}^{2}, and determine which of the three analysis patterns this data falls into.
#    The output is included in gqcal(out1104$ext[,1:41],3).

# 4) Calculate the probability of a difference between levels. Use 5 km/h and 10 km/h as reference points.
# PHC(aj - aj′ > c)  Identifying differences between levels of Factor B
PHC02(c = 5.0, out1104$muB, cc = "gtc", 3)
PHC02(c = 10.0, out1104$muB, cc = "gtc", 3)

# 5) Interpret the difference between a changeup and a curveball when there are no runners, in terms of mean difference,
# delta, non-overlap rate, and threshold rate (0 km/h, 1 km/h, 2 km/h, 3 km/h).
#    (Comparison of two cells of particular interest)
colnames(out1104$cellmean);              #Compare columns 8 and 12
E2between(out1104,digits=3,I=8,J=12,cr1=c(0,1,2,3))

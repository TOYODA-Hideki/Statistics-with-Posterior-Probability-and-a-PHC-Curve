####### Chapter 14 Simple Linear Regression Analysis
# Confirmation of working directory
(n_wd <- getwd())                # Checking the working directory
source('myfunc/myfunc.R')        # Loading custom functions
library(cmdstanr)                # Load and attach add-on package cmdstanr
library(posterior)               # Load and attach add-on package posterior
prob <- c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points


########## Functions to be used in Chapter 14 - Calculate initial values for simple linear regression using the least squares method
initi01 <- function(y, x) {      # y: dependent variable, x: predictor variable
  lmm <- lm(y ~ x)             # Store the result of the least squares method
  a <- as.numeric(coef(lmm)[1])   # Extract the intercept, the first element
  b <- as.numeric(coef(lmm)[2])   # Extract the regression coefficient, the second element
  sigma <- as.numeric(summary(lmm)$sigma)  # Extract σe from Equation (1.5)
  list(a=a, b=b, sigma=sigma) #Combine them into a list as the return value
}

########## 14.1 Galton (1886) Parent-Child Height Data

########  Table 14.2 - Partial "Height Data" (in cm)
parent_child_height <- read.csv("dat/parent_child_height.csv", header = TRUE)
head(parent_child_height, n = 14)  # View the first 14 rows
tail(parent_child_height, n = 14)  # View the last 14 rows

####### Table 1.3 - Numerical Summary of "Height Data"
hsd <- function(x) {sqrt(mean((x - mean(x))^2))}  # Function to calculate standard deviation using sample variance
round(apply(parent_child_height, 2, mean), 1)    # Mean
round(apply(parent_child_height, 2, hsd), 2)    # Standard deviation
round(apply(parent_child_height, 2, quantile), 1)  # Quantiles

####### Figure 14.1 - Scatterplot of "Height Data" + Random Noise
set.seed(1234);                               # Set the random seed
R_parent_child_height <- parent_child_height + rnorm(928 * 2, 0, 0.8)  # Add random noise
plot(R_parent_child_height, type = "p", xlim = c(155, 190), ylim = c(155, 190), pch = 1,
   cex.axis = 1.5, cex.lab = 1.3)
#dev.copy2pdf(file = "../English_tex/fig/chap14/fig14_01.pdf")

####### Figure 14.3 - Scatterplot with Regression Line (Overlay on Figure 1.1)
abline(0, 1, lty = 2, lwd = 2.0)               # Intercept 0, Slope 1 Line
abline(60.856, 0.646, lty = 1, lwd = 2.0)      # Regression line based on point estimates
#dev.copy2pdf(file = "../English_tex/fig/chap14/fig14_03.pdf")

############################### Stan Code for Simple Linear Regression Model
reg_01 <- '
// Simple Linear Regression Model
data {
  int<lower=0> n; // 01
  vector[n] y; // 02 Definition of the number of data points on page 2, line 5
  vector[n] X; // 03 Definition of the dependent variable
}
parameters { // 05
  real a; // 06 Definition of the intercept (Equation 14.3)
  real b; // 07 Definition of the regression coefficient (Equation 14.3)
  real<lower=0> sigma; // 08 Definition of the error standard deviation (Equation 14.5)
}
transformed parameters { // 09
  vector[n] yhat; // 10 Definition of predicted values
  for (i in 1 : n) { // 11 Loop from 1 to n
    yhat[i] = a + X[i] * b; // 12 Predicted values (Equation 14.2)
  } 
}
model { // 13
  for (i in 1 : n) { // 14 Loop from 1 to n
    y[i] ~ normal(yhat[i], sigma); // 15 Normal distribution (Equation 14.6)
  } 
}
'

par <- c("a", "b", "sigma")  # Parameters
dataSet_tankaiki01 <- list(n = nrow(parent_child_height), y = parent_child_height[, 2], X = parent_child_height[, 1])

########### Execution with cmdstanr
modfile_tankaiki01 <- write_stan_file(reg_01)   # Temporary file for writing
mod_tankaiki01 <- cmdstan_model(modfile_tankaiki01)  # Compile
csrfit_tankaiki01<- mod_tankaiki01$sample(data= dataSet_tankaiki01, chains = 5,
  init=function(){initi01(parent_child_height[, 2],parent_child_height[, 1]) },
  iter_sampling=20000,iter_warmup=1000,parallel_chains=5,seed=1234)  # MCMC
ext_tankaiki01<-as_draws_df(csrfit_tankaiki01$draws(par))
(colnames(ext_tankaiki01) <- gsub("\\[|\\]|,", "", colnames(ext_tankaiki01)))

# Save the 'fit_tankaiki01' object to a file
# save(ext_tankaiki01, file = "./script/obje/out1401")

# Load a previously saved 'fit_tankaiki01' object
# load(file = "./script/obje/out1401")

# Table 14.4 - Summary of Posterior Distributions for Parameters, R-squared, and Correlation Coefficients (First Half)
gqcal(ext_tankaiki01$a)     # Calculate summary statistics for the intercept
gqcal(ext_tankaiki01$b)     # Calculate summary statistics for the regression coefficient
gqcal(ext_tankaiki01$sigma) # Calculate summary statistics for sigma (error standard deviation)

# Plot histograms for the posterior distribution of the intercept (not in the textbook)
hist(ext_tankaiki01$a, breaks = 100)

# Plot histograms for the posterior distribution of the regression coefficient (not in the textbook)
hist(ext_tankaiki01$b, breaks = 100)

# Execute simple linear regression analysis using a function
fit_tankaiki02 <- Reg(parent_child_height[, 2], parent_child_height[, 1], 
                      cha = 5, ite = 20000)

# Save the 'fit_tankaiki02' object to a file
# save(fit_tankaiki02, file = "./script/obje/out1402")

# Load a previously saved 'fit_tankaiki02' object
# load(file = "./script/obje/out1402")

# Table 14.4 - Summary of Posterior Distributions for Parameters, R-squared, and Correlation Coefficients (All)
gqcal(fit_tankaiki02$ext[,1:5])
out2 <- print(fit_tankaiki02)

# 14.4.2 "Probability that the Research Hypothesis is True" (PHC)
bp <- as.vector(fit_tankaiki02$b)
r2p <- as.vector(fit_tankaiki02$r2)
round(mean(bp > 0.6), 3)  # PHC ("Regression coefficient is greater than 0.6") Equation (1.21) on page 9, line 10
round(mean(r2p > 0.2), 3)  # PHC ("R-squared is greater than 0.2") Equation (1.21) on page 9, line 11

####### Table 14.5 PHC table of regression coefficients and coefficient of determination
PHC01(seq(0.54, 0.60, 0.01), bp, cc="gtc", byoga="no")
PHC01(seq(0.14, 0.20, 0.01), r2p, cc="gtc", byoga="no")

####### Figure 14.4 PHC curve of regression coefficients and coefficient of determination
par(mfrow=c(2, 1));     # Figure 1.4
PHC01(seq(0.5, 0.75, 0.01), bp, cc="gtc", byoga="yes", xlab="Regression coefficient")
PHC01(seq(0.1, 0.30, 0.01), r2p, cc="gtc", byoga="yes", xlab="Coefficient of determination")
par(mfrow=c(1, 1))
#dev.copy2pdf(file="../English_tex/fig/chap14/fig14_04.pdf")

####### 14.4.4 Regression effect
####### Table 14.6 Summary of the posterior distribution of predicted values
out3 <- print(fit_tankaiki02, degits=1, Xnew=seq(155, 190, 2))

######Figure 14.5 Regression line, confidence interval, and prediction interval
plot(R_parent_child_height, type="p",xlim=c(155, 190), ylim=c(155, 190), pch=1,
     cex.axis=1.5, cex.lab=1.3)
abline(60.856, 0.646, lty=1, lwd=2.0)
lines(seq(155, 190, 2), out3$yhatc[,"0.025"], lty=2, lwd=2.0)
lines(seq(155, 190, 2), out3$yhatc[,"0.975"], lty=2, lwd=2.0)
lines(seq(155, 190, 2), out3$yastc[,"0.025"], lty=3, lwd=2.0)
lines(seq(155, 190, 2), out3$yastc[,"0.975"], lty=3, lwd=2.0)
#dev.copy2pdf(file="../English_tex/fig/chap14/fig14_05.pdf");       

####### Table 14.7 Summary of the posterior predictive distribution of the reference variable
Xnew = seq(155, 190, 5)
out4 <- print(fit_tankaiki02, degits=1, Xnew=Xnew)

####### Figure 14.6 Residual plot (scatter plot of predictor variables and residuals)
plot(parent_child_height[,1], out2$resi, type="p", xlim=c(155, 190), pch=1,
     cex.axis=1.5, cex.lab=1.3, xlab='Parents', ylab='Residuals')
abline(h=0, lwd=2);                                          # End of Figure 14.6
#dev.copy2pdf(file="../English_tex/fig/chap14/fig14_06.pdf")

############## Homework
# Input 'Pasta' data
x1 <- c(130, 268, 104, 185, 128, 147, 162, 68, 142, 175); # Actual measurement
x2 <- c(110, 232, 176, 207, 122, 202, 191, 124, 193, 250); # Estimated
x <- cbind(x1, x2)


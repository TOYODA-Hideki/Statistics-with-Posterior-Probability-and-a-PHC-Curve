########################################################################
###The working directory should be 'PPandPHC'
#chapter16
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)           #Load and attach add-on package cmdstanr
library(posterior)          #Load and attach add-on package posterior
library(psych)
library(rgl)
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) #Definition of probability points
par<-c("bs","bs12","r2","r","indi","tota");#Since the parameters are common, only this chapter is declared first.

####################### Functions used in Chapter 16 start here
### Initial values for multiple and simple regression
initi <- function(y, x1, x2) {
   lm2 <- lm(y ~ x1 + x2)
   lm1 <- lm(x2 ~ x1)
   ay <- as.numeric(coef(lm2)[1])
   a2 <- as.numeric(coef(lm1)[1])
   b <- as.vector(coef(lm2)[c(2, 3)])
   b12 <- as.numeric(coef(lm1)[2])
   sigma_y <- as.numeric(summary(lm2)$sigma)
   sigma_2 <- as.numeric(summary(lm1)$sigma)
   list(ay=ay, a2=a2, b=b, b12=b12, sigma_y=sigma_y, sigma_2=sigma_2)
}

#### Print function
gqcalc3 <- function(ext) {
   print('Partial regression coefficient y2'); print(gqcal(ext$bs2))
   print('Partial regression coefficient 21'); print(gqcal(ext$bs12))
   print('Direct effect y1 '); print(gqcal(ext$bs1))
   print('Indirect effect '); print(gqcal(ext$indi))
   print('Total effect '); print(gqcal(ext$tota))
   print('Coefficient of determination '); print(gqcal(ext$r2))
}
####################### Functions used in Chapter 16 end here

############### Stan code for multiple and simple regression models
mu_si_reg <- '
data { 
  int<lower=0>  n;                     // Number of data 
  vector[n]     y;                     // Reference variable   
  matrix[n,2]   X;                     // Predictor variable matrix   
}
transformed data {
  real      sdy;                       // Standard deviation of y
  vector[2] sdX;                       // Standard deviations of x1, x2
  sdy    = sqrt((variance(y)*(n-1))/n);      // 01  Standard deviation of y
  sdX[1] = sqrt((variance(X[,1])*(n-1))/n);  // 02  Standard deviation of x1
  sdX[2] = sqrt((variance(X[,2])*(n-1))/n);  // 03  Standard deviation of x2
}
parameters {
  real        ay;                      // Intercept (multiple regression)
  real        a2;                      // Intercept (simple regression)
  vector[2]   b;                       // Partial regression coefficient vector  
  real        b12;                     // Partial regression coefficient
  real<lower=0> sigma_y;               // Error standard deviation (multiple regression)
  real<lower=0> sigma_2;               // Error standard deviation (simple regression)
}
transformed parameters {
  vector[n] yhat;                      // yhat
  vector[n] x2hat;                     // x2hat
  yhat = ay + X*b;                     // 04  Multiple regression equation (equation(16.4))
  x2hat= a2 + X[,1]*b12;               // 05  Simple regression equation (equation(16.5))
}
model {
  y    ~ normal(yhat,  sigma_y);       // 06  Multiple regression model (equation(15.3))
  X[,2]~ normal(x2hat, sigma_2);       // 07  Simple regression model (equation(14.6))
}
generated quantities{
  vector[2]         bs;                
  real              bs12;              
  real<lower=0> vyhat;                 
  real<lower=0,upper=1>    r2;         
  real<lower=0,upper=1>    r;          
  real                     indi;       
  real                     tota;       
  vyhat =(variance(yhat)*(n-1))/n;     // 08  Variance of yhat
  r2 = vyhat/(vyhat+sigma_y^2);        // 09  Coefficient of determination (equation(14.15), (15.8))
  r =sqrt(r2);                         // 10  Multiple correlation coefficient (equation(15.9))
  bs[1] = (b[1]*sdX[1])/sdy;           // 11  Standardized partial regression coefficient 1 (equation(15.7))
  bs[2] = (b[2]*sdX[2])/sdy;           // 12  Standardized partial regression coefficient 2 (equation(15.7))
  bs12 = (b12*sdX[1])/sdX[2];          // 13  Standardized regression coefficient (equation(15.7))
  indi = bs12*bs[2];                   // 14  Indirect effect (equation(16.6))
  tota = indi+bs[1];                   // 15  Total effect (equation(16.6))
}
';

###### 16.1 Difficulties in interpreting partial regression coefficients when there are many predictor variables
###### 16.2 Direct, Indirect, and Total Effects

###### 15.1 Analysis of Rectal Cancer Data (Continuation from Chapter 15, hence classification (3) is at the beginning)
(can <- read.csv("dat/rectal_cancer.csv", header=T))
dataSet03 <- list(n=nrow(can), y=can[,3], X=can[,c(2,1)]);      # Input
ini <- function() { initi(y=can[,3], x1=can[,2], x2=can[,1]) }; #Initial values

########### Run with cmdstanr
modfile03 <- write_stan_file(mu_si_reg)             # Write to a temporary file
mod03 <- cmdstan_model(modfile03)                   # Compile
csrfit03 <- mod03$sample(data = dataSet03, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext03<-as_draws_df(csrfit03$draws(par))
(colnames(ext03) <- gsub("\\[|\\]|,", "", colnames(ext03)))
#save(ext03, file = "./script/obje/out1603")
#load(file="./script/obje/out1603");#Load object previously created with stan()

gqcalc3(ext = ext03);                       # Estimated values



####### 16.3.1 (1) When a suppressor variable exists
d01Biology <- read.csv("dat/ch02/d01Biology_raw.csv", header = T)
####### Figure 16.4 Multivariate scatter plot (1)
pairs.panels(d01Biology, smooth = F, density = F, 
             ellipses = F, pch = 16, rug = F, cex.axis = 1.0)
dataSet01 <- list(n = nrow(d01Biology), y = d01Biology[, 3], 
             X = d01Biology[, c(1, 2)]) # Input
ini <- function() { initi(y = d01Biology[, 3], x1 = d01Biology[, 1], 
             x2 = d01Biology[, 2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit01 <- mod03$sample(data = dataSet01, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext01<-as_draws_df(csrfit01$draws(par))
(colnames(ext01) <- gsub("\\[|\\]|,", "", colnames(ext01)))
#save(ext01, file = "./script/obje/out1601")
#load(file="./script/obje/out1601");#Load object previously created with stan()

######## Figure 16.6 Path diagram (1)
gqcalc3(ext = ext01)
######## Figure 16.5 3D scatter plot regression plane (1)
plot3d(d01Biology[, 1], d01Biology[, 2], d01Biology[, 3], size = 3, type = "p"
 , xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = -0.410, b = 0.749, c = -1, d = 0, alpha = 0.1)



####### 16.3.2 (2) When the direct effect from x1 to y is positive, and the correlation & indirect effect is negative
d02prejudice <- read.csv("dat/ch02/d02prejudice_raw.csv", header=T)
######## Figure 16.7 Multivariate scatter plot (2)
pairs.panels(d02prejudice,smooth=F,density=F,ellipses=F,
             pch=16,rug=F,cex.axis=1.0)
dataSet02 <- list(n=nrow(d02prejudice), y=d02prejudice[,3], 
             X=d02prejudice[,c(1,2)]) # Input
ini <- function(){ initi(y=d02prejudice[,3], x1=d02prejudice[,1], 
             x2=d02prejudice[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit02 <- mod03$sample(data = dataSet02, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext02<-as_draws_df(csrfit02$draws(par))
(colnames(ext02) <- gsub("\\[|\\]|,", "", colnames(ext02)))
#save(ext02, file = "./script/obje/out1602")
#load(file="./script/obje/out1602");#Load object previously created with stan()

######## Figure 16.9 Path diagram (2)
gqcalc3(ext = ext02)
######## Figure 16.8 3D scatter plot regression plane (2)
plot3d(d02prejudice[,1], d02prejudice[,2], d02prejudice[,3], size=3, type ="p",
       xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = 0.341, b = 0.875, c = -1, d = 0, alpha = 0.1)



####### 16.3.3 (4) When the indirect effect from x1 to y is positive, and the correlation & direct effect is negative
d04Impression <- read.csv("dat/ch02/d04Impression_raw.csv", header=T)
######## Figure 16.11 3D scatter plot regression plane (4)
pairs.panels(d04Impression, smooth=F, density=F, ellipses=F, 
             pch=16, rug=F, cex.axis=1.0)
dataSet04 <- list(n=nrow(d04Impression), y=d04Impression[,3], 
             X=d04Impression[,c(1,2)])                           # Input
ini <- function(){ initi(y=d04Impression[,3], x1=d04Impression[,1], 
             x2=d04Impression[,2]) }                    # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit04 <- mod03$sample(data = dataSet04, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext04<-as_draws_df(csrfit04$draws(par))
(colnames(ext04) <- gsub("\\[|\\]|,", "", colnames(ext04)))
#save(ext04, file = "./script/obje/out1604")
#load(file ="./script/obje/out1604");#Load object previously created with stan

######## Figure 16.12 Path diagram (4)
gqcalc3(ext = ext04)
######## Figure 16.11 3D scatter plot regression plane (4)
plot3d(d04Impression[,1],d04Impression[,2],d04Impression[,3], size=3, type="p",
       xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = -0.720, b = 0.952, c = -1, d = 0, alpha = 0.1)



####### 16.3.4 (5) When the indirect effect from x1 to y is negative, and the correlation & direct effect is positive
d05Yield <- read.csv("dat/ch02/d05Yield_raw.csv", header=T)
######## Figure 16.13 Multivariate scatter plot (5)
pairs.panels(d05Yield,smooth=F,density=F,ellipses=F,pch=16,rug=F, cex.axis=1.0)
dataSet05 <- list(n=nrow(d05Yield), y=d05Yield[,3], X=d05Yield[,c(1,2)]) #Input
ini <- function(){ initi(y=d05Yield[,3], x1=d05Yield[,1], 
        x2=d05Yield[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit05 <- mod03$sample(data = dataSet05, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext05<-as_draws_df(csrfit05$draws(par))
(colnames(ext05) <- gsub("\\[|\\]|,", "", colnames(ext05)))
#save(ext05, file = "./script/obje/out1605")
#load(file="./script/obje/out1605");#Load object previously created with stan

######## Figure 16.15 Path diagram (5)
gqcalc3(ext = ext05)
######## Figure 16.14 3D scatter plot regression plane (5)
plot3d(d05Yield[,1], d05Yield[,2], d05Yield[,3], size = 3, type = "p",
       xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = 0.460, b = 0.540, c = -1, d = 0, alpha = 0.1)



####### 16.3.5 (6) When there is a positive additive effect
d06Grade8 <- read.csv("dat/ch02/d06Grade8_raw.csv", header=T)
######## Figure 16.16 Multivariate scatter plot (6)
pairs.panels(d06Grade8,smooth=F,density=F,ellipses=F,pch=16,rug=F,cex.axis=1.0)
dataSet06 <- list(n=nrow(d06Grade8),y=d06Grade8[,3],X=d06Grade8[,c(1,2)])#Input
ini <- function(){ initi(y=d06Grade8[,3], x1=d06Grade8[,1], 
       x2=d06Grade8[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit06 <- mod03$sample(data = dataSet06, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext06<-as_draws_df(csrfit06$draws(par))
(colnames(ext06) <- gsub("\\[|\\]|,", "", colnames(ext06)))
#save(ext06, file = "./script/obje/out1606")
#load(file ="./script/obje/out1606");#Load object previously created with stan

######## Figure 16.18 Path diagram (6)
gqcalc3(ext = ext06)
######## Figure 16.17 3D scatter plot regression plane (6)
plot3d(d06Grade8[,1], d06Grade8[,2], d06Grade8[,3], size = 3, type = "p",
       xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = 0.468, b = 0.579, c = -1, d = 0, alpha = 0.1)



####### 16.3.6 (7) When there is a negative additive effect
d07Loneliness <- read.csv("dat/ch02/d07Loneliness_raw.csv", header=T)
######## Figure 16.19 Multivariate scatter plot (7)
pairs.panels(d07Loneliness, smooth=F, density=F, ellipses=F, 
             pch=16, rug=F, cex.axis=1.0)
dataSet07 <- list(n=nrow(d07Loneliness), y=d07Loneliness[,3], 
             X=d07Loneliness[,c(1,2)]) # Input
ini <- function(){ initi(y=d07Loneliness[,3], x1=d07Loneliness[,1], 
             x2=d07Loneliness[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit07 <- mod03$sample(data = dataSet07, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext07<-as_draws_df(csrfit07$draws(par))
(colnames(ext07) <- gsub("\\[|\\]|,", "", colnames(ext07)))
#save(ext07, file = "./script/obje/out1607")
#load(file="./script/obje/out1607");#Load object previously created with stan

######## Figure 16.21 Path diagram (7)
gqcalc3(ext = ext07)
######## Figure 16.20 3D scatter plot regression plane (7)
plot3d(d07Loneliness[,1], d07Loneliness[,2], d07Loneliness[,3], size=3, 
       type="p", xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = -0.657, b = 0.433, c = -1, d = 0, alpha = 0.1)



####### 16.3.7 (8) When multicollinearity occurs
d08Height <- read.csv("dat/ch02/d08Height_raw.csv", header=T)
######## Figure 16.22 Multivariate scatter plot (8)
pairs.panels(d08Height,smooth=F,density=F,ellipses=F,pch=16,rug=F,cex.axis=1.0)
dataSet08 <-list(n=nrow(d08Height),y=d08Height[,3],X=d08Height[,c(1,2)]) #Input
ini <- function(){ initi(y=d08Height[,3], x1=d08Height[,1], 
           x2=d08Height[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit08 <- mod03$sample(data = dataSet08, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext08<-as_draws_df(csrfit08$draws(par))
(colnames(ext08) <- gsub("\\[|\\]|,", "", colnames(ext08)))
#save(ext08, file ="./script/obje/out1608")
#load(file="./script/obje/out1608");#Load object previously created with stan

######## Figure 16.24 Path diagram (8)
gqcalc3(ext = ext08)
######## Figure 16.23 3D scatter plot regression plane (8)
plot3d(d08Height[,1], d08Height[,2], d08Height[,3], size = 3, type = "p",
       xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = 4.859, b = -4.641, c = -1, d = 0, alpha = 0.1)



##### 16.3.8 (9) When the correlation between predictor variables is close to 0
d09WritingTest <- read.csv("dat/ch02/d09WritingTest_raw.csv", header=T)
######## Figure 16.25 Multivariate scatter plot (9)
pairs.panels(d09WritingTest, smooth=F, density=F, ellipses=F, 
             pch=16, rug=F, cex.axis=1.0)
dataSet09 <- list(n=nrow(d09WritingTest), y=d09WritingTest[,3], 
             X=d09WritingTest[,c(1,2)]) # Input
ini <- function(){ initi(y=d09WritingTest[,3], x1=d09WritingTest[,1], 
             x2=d09WritingTest[,2]) } # Initial values

########### Run with cmdstanr #### Use 'mod03' for compilation
csrfit09 <- mod03$sample(data = dataSet09, chains = 5, iter_sampling = 20000,
    init = ini, iter_warmup = 1000, parallel_chains = 5, seed = 1234)  # MCMC
ext09<-as_draws_df(csrfit09$draws(par))
(colnames(ext09) <- gsub("\\[|\\]|,", "", colnames(ext09)))
#save(ext09, file = "./script/obje/out1609")
#load(file ="./script/obje/out1609");#Load object previously created with stan

######## Figure 16.27 Path diagram (9)
gqcalc3(ext = ext09)
######## Figure 16.26 3D scatter plot regression plane (9)
plot3d(d09WritingTest[,1], d09WritingTest[,2], d09WritingTest[,3], size = 3, 
       type = "p", xlab = "z1", ylab = "z2", zlab = "z3") 
planes3d(a = 0.434, b = 0.626, c = -1, d = 0, alpha = 0.1)

########## Chapter 19: Inference with Independent Single Factors by Various Distributions
(n_wd<-getwd())                # Confirmation of working directory
source('myfunc/myfunc.R')      # Loading self-made functions
library(cmdstanr)              # load and attach add-on package cmdstanr
library(posterior)             # load and attach add-on package posterior
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points

################### 19.1 Log-normal Distribution

##### Figure 19.1 Probability Density Function of the Log-normal Distribution
par(mfrow=c(1,2))
for (i in c(0.5,1.0,1.8)){
  d1<-function(x){dlnorm(x,0,i)}
  curve(d1,0,3,ylim=c(0,1.2),xlab="",ylab="",cex.axis=2.0)
  par(new=T)
}
text(1.3,0.8,"σ=0.5",cex=2)
text(0.8,0.5,"σ=1.0",cex=2)
text(0.4,1.1,"σ=1.8",cex=2)
text(2.3,0.8,"μ=0",cex=2.5);                      # End of first plot
par(new=F)
for (i in c(0.0,0.2,0.4,0.6)){
  d1<-function(x){dlnorm(x,i,1)}
  curve(d1,0,4,ylim=c(0,0.7),xlab="",ylab="",cex.axis=2.0)
  par(new=T)
}
text(1.0,0.6,"μ=0.0",cex=2)
text(1.2,0.5,"μ=0.2",cex=2)
text(1.4,0.4,"μ=0.4",cex=2)
text(0.8,0.3,"μ=0.6",cex=2)
text(3.0,0.50,"σ=1.0",cex=2.5)
par(new=F)
par(mfrow=c(1,1));                                    # End of Figure 19.1
#dev.copy2pdf(file="../English_tex/fig/chap19/z1901.pdf")

################### 19.1.1 Numerical Example

### Figure 19.2 Distribution of Annual Income in a Population
hist(rlnorm(10000,15.3,0.5)/10000,xlim=c(0,2000),breaks=30,
     xlab="Ten thousand yen",ylab="",cex.axis=2.0,cex.lab=2.0)

exp(15.3+0.5*0.5^2);                                  # Mean (Formula 19.3)
exp(15.3);                                            # Median (Formula 19.5)
exp(15.3-0.5^2);                                      # Mode (Formula 19.6)
qlnorm(0.95, 15.3,0.5);                               # Percentile points
qlnorm(0.99, 15.3,0.5)
qlnorm(0.999,15.3,0.5)


################### 19.1.2 Estimation of the Log-normal Distribution
################### Table 19.1 Annual Income of 30 People in Occupation I (in Ten Thousand Yen)
x1<-c(367, 1611, 263, 405, 754, 415, 837, 341, 396, 342,
      282,  296, 886, 412, 572, 471, 781, 531, 757, 723,
      870,  933, 392, 612, 394, 343, 372, 747, 280, 941)


##################### Stan code for Estimation of the Log-normal Distribution
lognormal01<-'
data {
  int<lower=0> n; // Number of data
  array[n] real<lower=0> x; // Data
}
parameters {
  real mu; // Parameter μ (Formula 19.1)
  real<lower=0> sigma; // Parameter σ (Formula 19.1)
}
model {
  for (i in 1 : n) {
    x[i] ~ lognormal(mu, sigma);
  } // Log-normal distribution (Formulas 19.2 and 19.8)
}
generated quantities {
  real heikin; // Mean value
  real tyuouti; // Median
  real saihin; // Mode
  real sigma0; // Standard deviation
  heikin = exp(mu + 0.5 * sigma ^ 2); // (Formula 19.3)
  tyuouti = exp(mu); // (Formula 19.5)
  saihin = exp(mu - sigma ^ 2); // (Formula 19.6)
  sigma0 = sqrt(exp(2 * mu + sigma ^ 2) * (exp(sigma ^ 2) - 1)); // (Formula 19.4)
}
'

par<-c("mu","sigma","heikin","tyuouti","saihin","sigma0")    # Parameters
dataSet_1 <-list(n=length(x1),x=x1)                         # Input


########### Execution with cmdstanr
modfile_1 <- write_stan_file(lognormal01)          # Write out a temporary file
mod_1 <- cmdstan_model(modfile_1)                  # Compile
csrfit_1 <- mod_1$sample(data = dataSet_1,chains = 5,iter_sampling = 20000,
                     iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
ext_1<-as_draws_df(csrfit_1$draws(par))
(colnames(ext_1) <- gsub("\\[|\\]|,", "", colnames(ext_1)))

#save(ext_1,  file="./script/obje/out1901")
#load(file="./script/obje/out1901");#Load an object created earlier with stan()

###### Table 19.2 Summary of the Posterior Distribution of Parameters and Generated Quantities for Occupation I
gqcal(ext_1$mu)
gqcal(ext_1$sigma)
gqcal(ext_1$heikin ,0)
gqcal(ext_1$tyuouti,0)
gqcal(ext_1$saihin ,0)
gqcal(ext_1$sigma0 ,0)


################## 19.1.3 Comparison of Two Log-normal Distributions
################## Table 19.3 Annual Income of 30 People in Occupation II (in Ten Thousand Yen)
x2<-c(697, 681, 330, 176,  307, 1668, 641, 1150, 301, 365,
      1033, 373, 836, 367,  265,  217, 941,  278, 504, 299,
      865, 504, 790, 273, 2292,  385, 277,  622, 475,1205)

################### Stan code for Estimation of Two Log-normal Distributions
lognormal02<-'
data {
  int<lower=0> n1; // Data count I
  int<lower=0> n2; // Data count II
  array[n1] real<lower=0> x1; // Data I
  array[n2] real<lower=0> x2; // Data II
}
parameters {
  array[2] real mu; // Parameters μ1, μ2 (Formula 19.8)
  array[2] real<lower=0> sigma; // Parameters σ1, σ2 (Formula 19.8)
}
transformed parameters {

}
model {
  for (i in 1 : n1) {
    x1[i] ~ lognormal(mu[1], sigma[1]);
  } // Log-normal distribution
  for (i in 1 : n2) {
    x2[i] ~ lognormal(mu[2], sigma[2]);
  } // (Formula 19.8)
}
generated quantities {
  array[2] real heikin; // Two mean values
  array[2] real tyuouti; // Two median values
  array[2] real saihin; // Two mode values
  array[2] real sigma0; // Two standard deviations
  heikin[1] = exp(mu[1] + 0.5 * sigma[1] ^ 2); // (Formula 19.3) I
  tyuouti[1] = exp(mu[1]); // (Formula 19.5) I
  saihin[1] = exp(mu[1] - sigma[1] ^ 2); // (Formula 19.6) I
  sigma0[1] = sqrt(exp(2 * mu[1] + sigma[1] ^ 2) * (exp(sigma[1] ^ 2) - 1)); // (Formula 19.4) I
  heikin[2] = exp(mu[2] + 0.5 * sigma[2] ^ 2); // (Formula 19.3) II
  tyuouti[2] = exp(mu[2]); // (Formula 19.5) II
  saihin[2] = exp(mu[2] - sigma[2] ^ 2); // (Formula 19.6) II
  sigma0[2] = sqrt(exp(2 * mu[2] + sigma[2] ^ 2) * (exp(sigma[2] ^ 2) - 1)); // (Formula 19.4) II
}
';

par<-c("mu","sigma","heikin","tyuouti","saihin","sigma0")  # Parameters
dataSet_2 <-list(n1=length(x1),n2=length(x2),x1=x1,x2=x2) # Input

########### Execution with cmdstanr
modfile_2 <- write_stan_file(lognormal02)       # Write out a temporary file
mod_2 <- cmdstan_model(modfile_2)               # Compile
csrfit_2 <- mod_2$sample(data = dataSet_2,chains = 5,iter_sampling = 20000,
                     iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
ext_2<-as_draws_df(csrfit_2$draws(par))
(colnames(ext_2) <- gsub("\\[|\\]|,", "", colnames(ext_2)))

#save(ext_2, file="./script/obje/out1902")
#load(file="./script/obje/out1902");#Load an object created earlier with stan()

###### Table 19.4 Summary of the Posterior Distribution of Parameters and Generated Quantities for Occupation II
gqcal(ext_2$mu2)
gqcal(ext_2$sigma2)
gqcal(ext_2$heikin2 ,0)
gqcal(ext_2$tyuouti2,0)
gqcal(ext_2$saihin2 ,0)
gqcal(ext_2$sigma02 ,0)

# Table 19.5 Summary of the Posterior Distribution of the Difference in Location between Two Log-normal Distributions: Differences in Mean, Median, and Mode
gqcal(ext_2$heikin2-ext_2$heikin1)
gqcal(ext_2$tyuouti2-ext_2$tyuouti1)
gqcal(ext_2$saihin2-ext_2$saihin1)

###### Table 19.6 Probability of Difference in Mean, Median, and Mode Values
sequ66<-c(0,10,20,30)
PHC01(sequ66,ext_2$heikin2 ,ext_2$heikin1 ,cc="gtc",byoga="no")
PHC01(sequ66,ext_2$tyuouti2,ext_2$tyuouti1,cc="gtc",byoga="no")
PHC01(sequ66,ext_2$saihin1 ,ext_2$saihin2 ,cc="gtc",byoga="no")

###### Table 19.7 Probability of Difference in Mean, Median, and Mode Values being within ROPE (Range of Practical Equivalence)
sequ67<-c(0,50,100,120)
PHC01(sequ67,ext_2$heikin2 ,ext_2$heikin1 ,cc="rope",byoga="no")
PHC01(sequ67,ext_2$tyuouti2,ext_2$tyuouti1,cc="rope",byoga="no")
PHC01(sequ67,ext_2$saihin1 ,ext_2$saihin2 ,cc="rope",byoga="no")

######################## 19.2 One-factor Experiment with Log-normal Distribution
##### Table 19.9 Time and Deviation Number for the 15th Trial of the Mirror Tracing Task by Group
(mirror_tracing01<-read.csv("dat/mirror_tracing01.csv", header = TRUE))

##### Figure 19.3 Time by Group in the Mirror Tracing Task (seconds)
boxplot(mirror_tracing01$time~mirror_tracing01$group,             # Figure 19.3
        xlab="group", ylab="",cex.axis=2.0,cex.lab=2.0)
#dev.copy2pdf(file="../English_tex/fig/chap19/z1903.pdf")
tapply(mirror_tracing01$time,mirror_tracing01$group,mean); # Mean time by group p.82

############## Stan code for One-factor Experiment with Log-normal Distribution
one_factor_lognormal<-'
data {
  int<lower=0> n; // Total data count
  int<lower=0> J; // Number of groups
  vector[n] y; // Dependent variable
  array[n] int<lower=0> k; // Categorical variable
}
parameters {
  vector[J] mu; // mu for each group     J items
  vector<lower=0>[J] sigma; // sigma for each group  J items
}
model {
  for (i in 1 : n) {
    y[i] ~ lognormal(mu[k[i]], sigma[k[i]]);
  } // (Formula 19.14)
}
generated quantities {
  vector[J] mu0; // Mean for each group     J items
  vector[J] medi0; // Median for each group   J items
  vector[J] mode0; // Mode for each group   J items
  vector<lower=0>[J] sigma0; // Standard deviation for each group J items
  for (j in 1 : J) {
    mu0[j] = exp(mu[j] + (0.5 * sigma[j] ^ 2)); // (Formula 19.3)
    medi0[j] = exp(mu[j]); // (Formula 19.5)
    mode0[j] = exp(mu[j] - sigma[j] ^ 2); // (Formula 19.6)
    sigma0[j] = sqrt(exp(2 * mu[j] + sigma[j] ^ 2) * (exp(sigma[j] ^ 2) - 1)); // (Formula 19.4)
  }
}
';
par<-c("mu","sigma","mu0","medi0","mode0","sigma0")               # Parameters
dataSetlogn <-list(n=length(mirror_tracing01$time), 
                   J=max(mirror_tracing01$group),y=mirror_tracing01$time, k=mirror_tracing01$group ) # Input

########### Execution with cmdstanr
modfilelogn <- write_stan_file(one_factor_lognormal)#Write out a temporary file
modlogn <- cmdstan_model(modfilelogn)                   # Compile
csrfitlogn <- modlogn$sample(data=dataSetlogn,chains=5,iter_sampling = 20000,
                   iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
extlogn<-as_draws_df(csrfitlogn$draws(par))
(colnames(extlogn) <- gsub("\\[|\\]|,", "", colnames(extlogn)))

#save(extlogn, file="./script/obje/out1903")
#load(file="./script/obje/out1903");#Load an object created earlier with stan()

##### Table 19.10 Posterior Distribution of Parameters of the One-factor Log-normal Distribution Model
gqcal(extlogn[,1:18])

##### Figure 19.4 Distribution of Time by Group in the Mirror Tracing Task (Log-normal Distribution)
lognor01<-function(x){dlnorm(x,3.153, 0.685)}
lognor02<-function(x){dlnorm(x,2.981, 0.449)}
lognor03<-function(x){dlnorm(x,2.822, 0.488)}
curve(lognor01,0,80,ylim=c(0,0.055),xlab="",ylab="",lty=1,lwd=2)
par(new=T)
curve(lognor02,0,80,ylim=c(0,0.055),xlab="",ylab="",lty=2,lwd=2)
par(new=T)
curve(lognor03,0,80,ylim=c(0,0.055),xlab="time",
      ylab="Probability Density",lty=3,lwd=2)
text(60,0.01,"Group 1",cex=1.6)
text(30,0.04,"Group 2",cex=1.6)
text( 5,0.05,"Group 3",cex=1.6)                       # End of Figure 19.4

# Table 19.11 Probability that the Mean of Level j is Greater than that of Level j′
PHC02(0,extlogn[,7:9])

# Probability that the Rest Group Takes Longer than the Non-dominant Hand Group, and the Non-dominant Hand Group Takes Longer than the Dominant Hand Group
round(mean((extlogn$mu01>extlogn$mu02)&
           (extlogn$mu02>extlogn$mu03)),3)   # (Formulas 19.17 and 19.18)

######################## 19.3 One-factor Experiment with Poisson Distribution
##### Table 19.9 Number of Deviations and Time for the 15th Trial of the Mirror Tracing Task by Group
(mirror_tracing01<-read.csv("dat/mirror_tracing01.csv", header = TRUE))

##### Figure 19.5 Bar Graph of the Number of People per Deviation Number in Each Group for the Mirror Tracing Task
barplot(table(mirror_tracing01$deviation, mirror_tracing01$group),
        xlab="Number of Deviations per Group", ylab="Number of People", beside=T, legend.text=0:5 )

# Mean number of deviations by group p.84,bl.8
tapply(mirror_tracing01$deviation, mirror_tracing01$group, mean)

################ Stan Code for One-factor Experiment with Poisson Distribution
one_factor_Poisson<-'
data {
  int<lower=0> n; // Total data count
  int<lower=0> J; // Number of groups
  array[n] int<lower=0> y; // Dependent variable
  array[n] int<lower=0> k; // Categorical variable
}
parameters {
  vector<lower=0>[J] lambda; // lambda for each group J items
}
model {
  for (i in 1 : n) {
    y[i] ~ poisson(lambda[k[i]]);
  } // (Formula 19.19)
}
';

par<-c("lambda")                                                   # Parameters
dataSetpo <-list(n=length(mirror_tracing01$deviation),
                 J=max(mirror_tracing01$group), y=mirror_tracing01$deviation, k=mirror_tracing01$group ) # Input

########### Execution with cmdstanr
modfilepo <- write_stan_file(one_factor_Poisson)   # Write out a temporary file
modpo <- cmdstan_model(modfilepo)                  # Compile
csrfitpo <- modpo$sample(data = dataSetpo, chains = 5, iter_sampling = 20000,
               iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extpo<-as_draws_df(csrfitpo$draws(par))
(colnames(extpo) <- gsub("\\[|\\]|,", "", colnames(extpo)))

#save(extpo, file="./script/obje/out1904")
#load(file="./script/obje/out1904");#Load an object created earlier with stan()


#### Table 19.12 Posterior Distribution of Parameters of the One-factor Poisson Distribution Model
gqcal(extpo[,1:3])

##### Figure 19.6 Distribution of Deviation Numbers by Group in the Mirror Tracing Experiment (Poisson Distribution)
x<-0:5
plot(x, dpois(x,1.100), type="o", ylim=c(0,0.7), xlab="", ylab="", lty=1, lwd=2)
par(new=T)
plot(x, dpois(x,0.666), type="o", ylim=c(0,0.7), xlab="", ylab="", lty=2, lwd=2)
par(new=T)
plot(x, dpois(x,0.401), type="o", ylim=c(0,0.7),
     xlab="Number of Deviations", ylab="Probability", lty=3, lwd=2)
legend(3, 0.6, paste("Group", 1:3), lty=1:3)          # End of Figure 19.6

## Table 19.13 Probability that the Mean of Level j is Greater than that of Level j'
PHC02(0, extpo[,1:3])

## Probability that the Rest Group Takes Longer than the Non-dominant Hand Group, and the Non-dominant Hand Group Takes Longer than the Dominant Hand Group
round(mean((extpo$lambda1>extpo$lambda2) &
           (extpo$lambda2>extpo$lambda3)),3)    # (Formula 19.20)


################### 19.4 One-factor Experiment with Binomial Distribution
#### Table 19.14 Number of Basic and Expressed Words in Training for Conversation Facilitation by Group
(conversation_facilitation<-read.csv("dat/conversation_facilitation.csv", header = TRUE))

#### Figure 19.7 Expression Rate by Group in Conversation Facilitation Training
boxplot(conversation_facilitation$expressed_words/conversation_facilitation$basic_words~conversation_facilitation$group,  
        xlab="group", ylab="", cex.axis=2.0, cex.lab=2.0)
#dev.copy2pdf(file="../English_tex/fig/chap19/z1907.pdf")

##### Mean expression rate by group
round(tapply(conversation_facilitation$expressed_words/conversation_facilitation$basic_words, conversation_facilitation$group, mean), 3)

########################### Stan Code for One-factor Experiment with Binomial Distribution
one_factor_binom<-'
data {
  int<lower=0> n; // Total data count
  int<lower=0> J; // Number of groups
  array[n] int<lower=0> y; // Dependent variable
  array[n] int<lower=0> m; // Number of trials
  array[n] int<lower=0> k; // Categorical variable
}
parameters {
  vector<lower=0, upper=1>[J] p; // Probability for each group J items
}
model {
  for (i in 1 : n) {
    y[i] ~ binomial(m[i], p[k[i]]);
  } // (Formula 19.21)
}
';

par<-c("p")                                                       # Parameters
dataSetbi <-list(n=length(conversation_facilitation$expressed_words), J=max(conversation_facilitation$group),
                 y=conversation_facilitation$expressed_words, m=conversation_facilitation$basic_words, k=conversation_facilitation$group ) # Input

########### Execution with cmdstanr
modfilebi <- write_stan_file(one_factor_binom)   # Write out a temporary file
modbi <- cmdstan_model(modfilebi)                # Compile
csrfitbi <- modbi$sample(data = dataSetbi, chains = 5, iter_sampling = 20000,
                 iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extbi<-as_draws_df(csrfitbi$draws(par))
(colnames(extbi) <- gsub("\\[|\\]|,", "", colnames(extbi)))

#save(extbi, file="./script/obje/out1905")
#load(file="./script/obje/out1905");#Load an object created earlier with stan()

###### Table 19.15 Posterior Distribution of Parameters of the One-factor Binomial Distribution Model
gqcal(extbi[,1:3])

##### Figure 19.8 Binomial Distribution for Each Group Based on EAP Estimates
x<-0:10
plot(x, dbinom(x, 10, 0.338), type="o", ylim=c(0,0.3), xlab="", ylab="", lty=1, lwd=2)
par(new=T)
plot(x, dbinom(x, 10, 0.510), type="o", ylim=c(0,0.3), xlab="", ylab="", lty=2, lwd=2)
par(new=T)
plot(x, dbinom(x, 10, 0.613), type="o", ylim=c(0,0.3),
     xlab="Number of Expressions", ylab="Probability", lty=3, lwd=2)
legend(8, 0.27, paste("Group", 1:3), lty=1:3);          # End of Figure 19.8


### Table 19.16 Probability that the Rate of Level j is Greater than that of Level j' 
PHC02(0,extbi[,1:3])

# Probability that the Conjunction Proposition is True (Corresponding to Hypothesis C, not in textbook)
round(mean((extbi$p1<extbi$p2)&(extbi$p2<extbi$p3)), 3)


############ 19.5 One-factor Experiment with Normal Distribution (Variance Model)
### Table 19.17 Results of 10 Trials for 6 Athletes in the 100m Track Event
(M100m_run<-read.csv("dat/M100m_run.csv", header = TRUE))

#### Figure 19.9 Box Plot of the Times for 10 Trials of 6 Athletes in the 100m Track Event
boxplot(M100m_run$time~M100m_run$athlete, xlab="Athlete Number", ylab="time"); # Figure 19.9
#dev.copy2pdf(file="../English_tex/fig/chap19/z1909.pdf")

############## Stan Code for One-factor Experiment with Normal Distribution (Variance Model)
one_factor_normal<-'
data {
  int<lower=0> n; // Total data count
  int<lower=0> J; // Number of groups
  vector[n] y; // Dependent variable
  array[n] int<lower=0> k; // Categorical variable
}
parameters {
  vector[J] mu; // Mean of each group J items
  real<lower=0> sigma; // Error SD
  real<lower=0> s_mu; // Factor SD
  real mu0; // Overall mean
}
model {
  y ~ normal(mu[k], sigma); // (Formula 19.24)
  mu ~ normal(mu0, s_mu); // (Formula 19.25)
}
generated quantities {
  real<lower=0, upper=1> eta2; // Explained variance
  eta2 = (s_mu ^ 2) / ((s_mu ^ 2) + (sigma ^ 2)); // (Formula 19.27)
}
';

par<-c("mu","sigma","s_mu","mu0","eta2")     # Parameters
dataSetnorm <-list(n=length(M100m_run$time), 
      J=max(M100m_run$athlete), y=M100m_run$time, k=M100m_run$athlete ) # Input

########### Execution with cmdstanr
modfilenorm <- write_stan_file(one_factor_normal) # Write out a temporary file
modnorm <- cmdstan_model(modfilenorm)             # Compile
csrfitnorm <- modnorm$sample(data=dataSetnorm, chains=5, iter_sampling = 20000,
                   iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extnorm<-as_draws_df(csrfitnorm$draws(par))
(colnames(extnorm) <- gsub("\\[|\\]|,", "", colnames(extnorm)))

#save(extnorm, file="./script/obje/out1906")
#load(file="./script/obje/out1906");#Load an object created earlier with stan()

##### Table 19.18 Summary of Posterior Distribution of Parameters of the Variance Model
gqcal(extnorm[,1:10])

##### Table 19.19 Probability that the Mean of One Athlete is Greater than that of Another Athlete
PHC02(0,extnorm[,1:6])


################# Practical Assignment
###  1. Create the PHC curve and PHC table for "The median income is greater than c yen".

#
#     Correct answer omitted
#

###  2. Create the PHC curve and PHC table for "The modal income of Occupation II is c yen higher than that of Occupation I".

#
#     Correct answer omitted
#

###  3. Create the PHC table for formulas (19.15) and (19.16).

#
#     Correct answer omitted
#

###  4.

#
#     Correct answer omitted
#


###  5.
lambda_sa12 <- extpo$lambda1-extpo$lambda2
lambda_sa23 <- extpo$lambda2-extpo$lambda3
PHC01(seq(0,1,0.05),lambda_sa12,cc="gtc",byoga="yes",
      xlab="Difference between rest group and non-dominant hand group")
PHC01(seq(0,1,0.05),lambda_sa23,cc="gtc",byoga="yes",
      xlab="Difference between non-dominant hand group and dominant hand group")
PHC01(seq(0,0.3,0.05),lambda_sa12,cc="gtc",byoga="no")
PHC01(seq(0,0.2,0.05),lambda_sa23,cc="gtc",byoga="no")


###  6.

#
#     Correct answer omitted
#


###  7.
hiritu_sa21 <- extbi$p2-extbi$p1
hiritu_sa32 <- extbi$p3-extbi$p2
PHC01(seq(0,0.3,0.01),hiritu_sa21,cc="gtc",byoga="yes",
      xlab="Difference between 2 months and 1 month")
PHC01(seq(0,0.2,0.01),hiritu_sa32,cc="gtc",byoga="yes",
      xlab="Difference between 3 months and 2 months")
PHC01(seq(0,0.15,0.03),hiritu_sa21,cc="gtc",byoga="no")
PHC01(seq(0,0.1,0.02),hiritu_sa32,cc="gtc",byoga="no")


###  8.
eta2<-extnorm$eta2
PHC01(seq(0.1,0.8,0.05),eta2,cc="gtc",byoga="yes",
      xlab="Explanation rate")
PHC01(seq(0,3,0.02),extnorm$mu1-extnorm$mu2,cc="gtc",byoga="yes",
      xlab="Difference between athlete 1 and athlete 2")
PHC01(seq(0,0.15,0.03),eta2,cc="gtc",byoga="no")
PHC01(seq(0,1,0.2),extnorm$mu1-extnorm$mu2,cc="gtc",byoga="no")


############################### Mirror Tracing Experiment Data
# Time
rest_time <- c(40.03, 35.04, 32.38, 23.25, 24.45,
               30.66, 32.26, 21.18, 25.53, 37.66,
               22.51, 33.46, 14.87, 23.52, 26.63,
               35.27, 54.38, 17.59, 39.75, 29.63,
               36.72, 26.77, 24.67, 61.77, 31.28,
               27.83, 27.78, 33.15)
non_dominant_hand_time <- c(21.39, 13.89, 24.76, 24.20, 17.77,
                            24.90, 24.15, 13.56, 21.86, 17.58,
                            34.03, 20.41, 19.09, 27.90, 29.91,
                            20.88, 26.44, 21.29, 10.71, 21.92,
                            17.09, 31.98, 31.18, 18.37, 16.53,
                            27.40, 26.03)
dominant_hand_time <- c(08.63, 12.40, 15.11, 23.68, 16.28,
                        29.56, 25.40, 12.53, 19.00, 21.59,
                        21.93, 23.03, 17.09, 11.37, 23.84,
                        22.03, 18.98, 12.40, 25.77, 12.62,
                        13.78, 10.43, 24.49, 34.28, 13.74,
                        18.22, 18.56, 21.22)

# Number of Deviations
rest_freq <- c(1, 1, 0, 2, 1, 0, 2, 2, 2, 2, 0, 4, 2, 0, 1, 2, 2, 1, 6, 5, 1, 2, 3, 1, 0, 1, 2, 2)
non_dominant_hand_freq <- c(2, 0, 0, 0, 1, 2, 0, 3, 0, 1, 2, 1, 0, 0, 1, 0, 1, 2, 0, 0, 2, 1, 1, 0, 2, 4, 1)
dominant_hand_freq <- c(0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)


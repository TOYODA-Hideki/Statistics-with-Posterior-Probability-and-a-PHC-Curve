########################################################################
# Chapter 12
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
   library(cmdstanr)           #load and attach add-on package cmdstanr
   library(posterior)          #load and attach add-on package posterior

######### 12.3 Estimation of Ratios (One Binomial Distribution) #########

# Governor's approval rate
  x1 <- 305                               # Number of positive responses
  n1 <- 500                               # Number of data points

# Estimation of ratio (One binomial distribution)
out1201<-Bi01(x1,n1)
#save(out1201,  file="./script/obje/out1201")
#load(file="./script/obje/out1201"); # Load the results of a pre-run of MCMC

# Table 12.1 Numerical summary of the posterior and predictive distribution of ratio
# Table 12.2 Numerical summary of the posterior distribution of the odds of the governor's approval rate
gqcal(out1201$ext[,1:3])

# Figure 12.1 PHC curves for population ratio and odds
par(mfrow=c(2,1))
PHC01(seq(0.5,0.7,0.01),out1201$theta,0,cc="gtc",byoga="yes",xlab="Population ratio");
PHC01(seq(1.2,1.9,0.01),out1201$Odds ,0,cc="gtc",byoga="yes",xlab="Odds");
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap12/z1201.pdf")

# Table 12.3 PHC table for population ratio and odds
PHC01(seq(0.55,0.61,0.01),out1201$theta,0,cc="gtc",byoga="no");
PHC01(seq(1.28,1.40,0.02),out1201$Odds ,0,cc="gtc",byoga="no");

########## 12.4 Inference on a 2x2 Contingency Table (Product of Two Binomial Distributions) #################

# Brand recognition problem 1 
# Table 12.4 Results of brand recognition survey (number of people)
  x2 <- c(85,31)                          # Number of positive responses
  n2 <- c(123,121)                        # Number of data points

# Inference on a 2x2 contingency table (Product of two binomial distributions)
out1202<-Bi02(x2,n2)
#save(out1202,  file="./script/obje/out1202")
#load(file="./script/obje/out1202") # Load the results of a pre-run of MCMC

# Table 12.5 "Brand Recognition Issue 1": Summary of posterior distribution of population size and generated quantity
gqcal(out1202$ext[,1:9])

# PHC curve, not in the textbook
par(mfrow=c(3,1))
PHC01(seq(0.3,0.55,0.02),out1202$p_sa,0,cc="gtc",byoga="yes",xlab="Difference in population ratio")
PHC01(seq(1.8,3.8 ,0.02),out1202$p_hi,0,cc="gtc",byoga="yes",xlab="Ratio of population ratios")
PHC01(seq(3,12.0,0.02),out1202$Odds_hi,0,cc="gtc",byoga="yes",xlab="Odds ratio")
par(mfrow=c(1,1))

# Table 12.6 PHC Table for difference in ratios, ratio of ratios, odds ratio
PHC01(seq(0.28,0.40,0.02),out1202$p_sa    ,0,cc="gtc",byoga="no");
PHC01(seq(1.8,2.4,0.1),   out1202$p_hi    ,0,cc="gtc",byoga="no");
PHC01(seq(3.0,4.5,0.25),  out1202$Odds_hi ,0,cc="gtc",byoga="no");


########12.5 Prediction of g×2 cross-table (product of g binomial distributions)#####################

# Table 12.7 Whether there was a New Year's gift or not (Number of people)
x3 <- c(42,31,29,20)                     # Number of positive responses
n3 <- c(51,49,50,48)                     # Number of data points
g <- length(x3)                          # Number of groups

# Prediction of g×2 cross-table (product of g binomial distributions)
out1203<-Bi03(x3,n3)
#save(out1203,  file="./script/obje/out1203")
#load(file="./script/obje/out1203"); # Load the results of a pre-run of MCMC

# Table 12.8 Numerical summary of the posterior distribution of parameters for the "New Year's gift problem"
colnames(out1203$p) <-paste("p",1:4,sep="")
gqcal(out1203$p)

# Table 12.9 Probability that category in row i has a larger ratio than category in column j
PHC02(0,out1203$p,cc="gtc",3)

# Creation of generated quantities
p_sa13  <-out1203$p[,1]-out1203$p[,3]
p_hi13  <-out1203$p[,1]/out1203$p[,3]
Odds_hi13<-(out1203$p[,1]*(1-out1203$p[,3]))/(out1203$p[,3]*(1-out1203$p[,1]))

# Table 12.10 Comparison between high school students and university students (later)
gqcal(p_sa13  );  # Difference in ratios
gqcal(p_hi13  );  # Ratio of ratios
gqcal(Odds_hi13); # Odds ratio

# Figure 12.2 PHC curves for difference in population ratios, ratio of population ratios, Odds ratio
par(mfrow=c(3,1))
PHC01(seq(0.05,0.45,0.02),p_sa13   ,0,cc="gtc",byoga="yes",xlab="Difference in population ratios");
PHC01(seq(1.00,2.00,0.05),p_hi13   ,0,cc="gtc",byoga="yes",xlab="Ratio of population ratios");
PHC01(seq(1.00,8.00,0.25),Odds_hi13,0,cc="gtc",byoga="yes",xlab="Odds ratio");
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap12/z1202.pdf")

########################################################################
# Self-study Chapter 12 12.9 Practical Problems Pfizer's COVID Vaccine

# Table 12.12 Overview of the clinical trial
x4 <- c(8,162)                          # Number of cases
n4 <- c(2214,2222)                      # Number of subjects

# Inference about two binomial distributions
out1204<-Bi02(x4,n4)
#save(out1204,  file="./script/obje/out1204")
#load(file="./script/obje/out1204"); # Load the results of a pre-run of MCMC

# 1) Calculate and discuss the posterior distribution of the population ratio for the group that received the vaccine and the group that did not, and calculate PHC.
gqcal(out1204$ext[,1:9])
PHC01(seq(0.001,0.01,0.0005),out1204$p[,1],0,cc="ltc",byoga="yes",
      xlab="Experimental group population ratio");
PHC01(seq(0.001,0.01,0.001),out1204$p[,1],0,cc="ltc",byoga="no",
      xlab="Experimental group population ratio");
PHC01(seq(0.05,0.1,0.001),out1204$p[,2],0,cc="ltc",byoga="yes",
      xlab="Control group population ratio");
PHC01(seq(0.05,0.1,0.01),out1204$p[,2],0,cc="ltc",byoga="no",
      xlab="Control group population ratio");

# 2) How much does the risk reduce by taking the vaccine? Calculate and discuss the posterior distribution of the risk ratio and PHC.
PHC01(seq(0.01,0.1,0.001),out1204$p_hi,0,cc="ltc",byoga="yes",xlab="Risk ratio");
PHC01(seq(0.01,0.1,0.01),out1204$p_hi,0,cc="ltc",byoga="no",xlab="Risk ratio");
PHC01(seq(1,41,2),1/out1204$p_hi,0,cc="gtc",byoga="yes",xlab="Inverse of risk ratio");
PHC01(seq(2,20,2),1/out1204$p_hi,0,cc="gtc",byoga="no",xlab="Inverse of risk ratio");

# 3) As an indicator of the vaccine's effectiveness, the efficacy rate (=1-risk ratio) is often used, which is a medical indicator not explained in introductory statistics textbooks.
#    The efficacy rate indicates "what percentage of cases that would have occurred without vaccination are prevented by vaccination."
#    Construct the efficacy rate as a generated quantity, calculate and discuss its posterior distribution and PHC.
efficacy_rate <- 1-out1204$p_hi 
gqcal(efficacy_rate)
PHC01(seq(0.9,1.0,0.05),efficacy_rate,0,cc="gtc",byoga="yes",xlab="Efficacy rate");
PHC01(seq(0.9,1.0,0.01),efficacy_rate,0,cc="gtc",byoga="no",xlab="Efficacy rate");



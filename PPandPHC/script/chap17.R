############ Chapter 17: Logistic Regression/Meta-Analysis
(n_wd<-getwd())                # Confirmation of working directory
source('myfunc/myfunc.R')      # Loading self-made functions
library(cmdstanr)              # load and attach add-on package cmdstanr
library(posterior)             # load and attach add-on package posterior
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points

############# 17.1 Logistic Regression (Bernoulli Distribution)
### Table 4.1 Whether a Manager or Not and Age
(promotion<-read.csv("dat/promotion.csv", header = TRUE))

### Table 17.2 Number of Managers by Age
table(promotion);                                       # Table 17.2

### Figure 17.1 Scatter Plot and Regression Line
temp<-promotion[,2]+runif(200,-0.5,0.5); # Add randomness to make the figure clearer
plot(temp,promotion[,1],ylim=c(-0.2,1.2),xlim=c(30,46),
  xlab="Age",ylab="Manager");
  abline(a= -3.06335, b=0.09377, lwd=1.8);
  text(40,0.6,"Regression Line",cex=1.5);               # End of Figure 17.1

### 17.1.2 Logistic Transformation
# Odds conversion - Transformation from all real numbers to probability (Equation 17.4)
odds<-function(prob){prob/(1-prob)}
plot(odds,0,1,xlab="Probability",main="Odds Conversion")

# Logit transformation (logarithmic odds conversion) (Equation 17.5)
log_odds<-function(prob){log(prob/(1-prob))}
plot(log_odds,0,1,xlab="Probability",main="Logit Conversion")

# Logistic transformation (inverse of logit transformation) (Equation 17.7)
logistic<-function(yhat){1/(1+exp((-1)*yhat))}

### Figure 17.2 Inverse Logit Function logit-1( )
plot(logistic,-3,3,xlab="yhat",main="Logistic Transformation")

############### Logistic Regression (Bernoulli Distribution) in Stan code
logistic_reg_01<-'
data {
  int<lower=0> n; // Number of data 
  array[n] int<lower=0, upper=1> y; // Criterion variable (0-1 data)   
  vector[n] X; // Predictor variable vector   
}
parameters {
  real a; // Intercept
  real b; // Partial regression coefficient 
}
transformed parameters {
  vector<lower=0, upper=1>[n] prob;
  for (i in 1 : n) {
    prob[i] = inv_logit(a + X[i] * b);
  } // Response probability (Equation 17.10)
}
model {
  for (i in 1 : n) {
    y[i] ~ bernoulli(prob[i]);
  } // Bernoulli distribution model (Equation 17.8)
}
';

par<-c("a","b")                        #Parameters
dataSet01 <-list(n=nrow(promotion),y=promotion[,1], X=promotion[,2] ) # Input

########### Execution by cmdstanr
modfile01 <- write_stan_file(logistic_reg_01)       # Write temporary file
mod01 <- cmdstan_model(modfile01)                   # Compile
csrfit01 <- mod01$sample(data = dataSet01,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
ext<-as_draws_df(csrfit01$draws(par))
(colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))

#save(ext,  file="./script/obje/out1701")
#load(file="./script/obje/out1701"); # Load pre-run stan( ) object


##### Histogram of the posterior distribution (not in the textbook)
hist(ext$a,breaks=100)
hist(ext$b,breaks=100)

############### Generated Quantities
p0_5<- (-1)*ext$a/ext$b;         # Age at which there's a 50% chance of becoming a manager
oddshi<-exp(ext$b);              # Odds ratio when age changes by one year (Equation 17.13)
p35<- 1/(1+exp((-1)*(ext$a+ext$b*35))); # Probability of promotion at age 35 (Equation 17.14)

### Table 17.3 Summary of Posterior Distributions of Parameters and Generated Quantities
gqcal(ext$a);
gqcal(ext$b);
gqcal(p0_5);  
gqcal(oddshi);
gqcal(p35)

### Figure 17.3 Change in Promotion Probability with Age
plot(temp,promotion[,1],ylim=c(-0.2,1.2),xlim=c(30,46),
    xlab="Age",ylab="Manager",main="Logistic Regression")
  par(new=T)
  curve(1/(1+exp((-1)*((-18.679)+0.492*x))),30,46,lwd=1.8,
    ylim=c(-0.2,1.2),xlim=c(30,46),xlab="",ylab="");      # End of Figure 17.3


##############  17.2  Logistic Regression (Binomial Distribution)

### Table 17.4 Aphasia Conversation Facilitation Data
(aphasia_training<-read.csv("dat/AphasiaTraining.csv", header = TRUE));# Table 17.4

### Figure 17.4 Line Graph of Training Count and Correct Answer Rate
(ab<-coef(glm(aphasia_training[,2]/aphasia_training[,1]~aphasia_training[,3])))
plot(aphasia_training[,3],aphasia_training[,2]/aphasia_training[,1],
      ylim=c(0,1),type="o",xlab="Training Count",ylab="Accuracy Rate")
      abline(a=ab[1],b=ab[2]);                             # End of Figure 17.4


############ Figure 17.5   Viewing Each Probability Function of the Binomial Distribution Individually
size<-10; 
for (prob in (1:9)/10 ){
  x<-c(0:size)
  plot(x, dbinom(x, size, prob), type="o", 
        main=paste("Binomial Distribution size:",size,", prob:",prob) )
  locator(); # Right-click to "Stop"
}

############ Figure 17.5 Probability Function of the Binomial Distribution
size<-10; x<-c(0:size)
plot(x, dbinom(x, size, 0.1), type="o",ylim=c(0,0.4),xlab="",ylab="") 
text(1,0.3,paste("p=",prob,sep=""))
for (prob in c(0.3,0.5,0.7,0.9)){
  par(new=T)
  plot(x, dbinom(x, size, prob), type="o",ylim=c(0,0.4),xlab="",ylab="") 
  text(prob*10,0.3,paste("p=",prob,sep=""))
}

########### Case of Large Number of Trials (Becomes Indistinguishable from Normal Distribution) - Not in textbook
size<-100; 
for (prob in (1:9)/10 ){
  x<-c(0:size)
  plot(x, dbinom(x, size, prob), type="o", 
        main=paste("Binomial Distribution size:",size,", prob:",prob) )
  locator(); # Right-click to "Stop"
}

################ Logistic Regression (Binomial Distribution) Stan Code
binom_reg<-'
data {
  int<lower=0> n; // Number of data 
  array[n] int<lower=0> t; // Number of trials 
  array[n] int<lower=0> y; // Number of successes 
  vector[n] X; // Predictor variable vector   
}
parameters {
  real a; // Intercept
  real b; // Partial regression coefficient
}
transformed parameters {
  vector<lower=0, upper=1>[n] prob;
  for (i in 1 : n) {
    prob[i] = inv_logit(a + X[i] * b);
  } // Response probability (Equation 17.17)
}
model {
  for (i in 1 : n) {
    y[i] ~ binomial(t[i], prob[i]);
  } // Binomial distribution model (Equations 17.16 and 17.17)
}
';

par<-c("a","b")                                       # Parameters
dataSet02 <-list(n=nrow(aphasia_training),t=aphasia_training[,1],y=aphasia_training[,2],
                 X=aphasia_training[,3] )             # Input

########### Execution by cmdstanr
modfile02 <- write_stan_file(binom_reg)              # Write temporary file
mod02 <- cmdstan_model(modfile02)                    # Compile
csrfit02 <- mod02$sample(data = dataSet02,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
ext<-as_draws_df(csrfit02$draws(par))
(colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))

#save(ext,  file="./script/obje/out1702")
#load(file="./script/obje/out1702"); # Load pre-run stan( ) object


############### Posterior Distributions and Generated Quantities
p0_5<- (-1)*ext$a/ext$b;      # Number of trainings to express half the words 
oddshi<-exp(ext$b);           # Odds ratio for one more training session 
p35<- 1/(1+exp((-1)*(ext$a+ext$b*35)));    # Substitute 35 in (Equation 17.14) 
p35aste<-rbinom(length(p35),15,p35);       # (Equation 17.19)

##### Table 17.5 Summary of Posterior Distributions of Parameters and Generated Quantities
gqcal(ext$a  ,3);             # Intercept
gqcal(ext$b  ,3);             # Slope
gqcal(p0_5   ,1);             # Number of trainings to express half the words
gqcal(oddshi ,2);             # Odds ratio for one more training session
gqcal(p35    ,3);             # Correct answer rate after 5 more trainings
gqcal(p35aste,1);             # Posterior distribution of the number of words expressed in the 35th session with a base word count of 15

##### Figure 17.6 Logistic Regression Analysis (Binomial Distribution)
plot(aphasia_training[,3],aphasia_training[,2]/aphasia_training[,1],         
    ylim=c(0,1),xlim=c(0,36),type="o",xlab="Training Count",ylab="Correct Answer Rate")
  par(new=T)
  curve(1/(1+exp((-1)*((-1.803)+0.099*x))),0,36,lwd=1.8,
    ylim=c(0,1),xlim=c(0,36),xlab="",ylab="");
    lines(c(35,35),c(0.747,0.905),lwd=1.8)
    text(33,0.95,"95% Credible Interval");                          # End of Figure 17.6

##### Figure 17.7 Posterior Predictive Distribution of 'Number of Expressions' when x = 35, n = 15
hist(p35aste)

##### 17.3 Meta-analysis

######17.3.1 Prospective Studies
###### Table 17.10 Review of Prospective Studies on the Effectiveness of BCG
(prospective<-read.csv("dat/pre-tuberculosisBCG.csv", header = TRUE))

round(risk_actual_prospective<-prospective[1,3]/(prospective[1,3]+prospective[1,4]),4);   # Risk (Equation 17.20)
round(risk_control_prospective<-prospective[1,5]/(prospective[1,5]+prospective[1,6]),4);   # (Equation 17.21)
round(risk_actual_prospective/risk_control_prospective,3);                           # Risk ratio
round((prospective[1,3]*prospective[1,6])/(prospective[1,4]*prospective[1,5]),3); # Odds ratio (Equation 17.22)

######17.3.2 Retrospective Studies
##### Table 17.12 Review of Retrospective Studies on the Effectiveness of BCG
(retrospective<-read.csv("dat/post-tuberculosisBCG.csv", header = TRUE))   # Table 17.12
round(risk_actual_retro<-retrospective[1,3]/(retrospective[1,3]+retrospective[1,4]),4);   # Risk (Equation 17.23)
round(risk_control_retro<-retrospective[1,5]/(retrospective[1,5]+retrospective[1,6]),4);   # (Equation 17.24)
round((retrospective[1,3]*retrospective[1,6])/(retrospective[1,4]*retrospective[1,5]),3); # Odds ratio (Equation 17.25)

########################### Meta-analysis (Odds Ratio) Stan Code
meta_ana<-'
data {
  int<lower=0> m; // Number of studies 
  array[m] int<lower=0> u1; // Success count in the numerator of odds ratio 
  array[m] int<lower=0> n1; // Trial count in the numerator of odds ratio 
  array[m] int<lower=0> u0; // Success count in the denominator of odds ratio 
  array[m] int<lower=0> n0; // Trial count in the denominator of odds ratio 
}
parameters {
  array[m] real a; // Intercept
  array[m] real b; // Coefficient
  real mu_a; // Mean of intercept
  real mu_b; // Mean of coefficient
  real sigma_a; // SD of intercept
  real sigma_b; // SD of coefficient
}
model {
  for (i in 1 : m) {
    a[i] ~ normal(mu_a, sigma_a); // Prior distribution of intercept (Equation 17.30)
    b[i] ~ normal(mu_b, sigma_b); // Prior distribution of coefficient (Equation 17.31)
    u1[i] ~ binomial(n1[i], inv_logit(a[i] + b[i])); // (Equations 17.26, 17.28, 17.34)
    u0[i] ~ binomial(n0[i], inv_logit(a[i])); // (Equations 17.26, 17.29, 17.35)
  }
}
generated quantities {
  array[m] real Odds_r;
  real Odds_r_mu;
  Odds_r = exp(b); // Odds ratio for each study (Equation 17.13)
  Odds_r_mu = exp(mu_b); // Mean of odds ratio
}
';

##########17.3.4 Meta-analysis of Prospective Studies

par<-c("Odds_r","Odds_r_mu","mu_b","sigma_b")           # Parameters
dataSetProspective <-list(m=nrow(prospective),
      u1=prospective$u11,n1=prospective$u11+prospective$u12,
      u0=prospective$u21,n0=prospective$u21+prospective$u22)          # Input

########### Execution by cmdstanr
modfileProspective <- write_stan_file(meta_ana)        # Write temporary file  
modProspective <- cmdstan_model(modfileProspective)                   # Compile
csrfitProspective <- modProspective$sample(data = dataSetProspective,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
extProspective<-as_draws_df(csrfitProspective$draws(par))
(colnames(extProspective) <- gsub("\\[|\\]|,", "", colnames(extProspective)))

#save(extProspective,  file="./script/obje/out1703")
#load(file="./script/obje/out1703"); # Load pre-run stan( ) object

##### Table 17.10 Odds/EAP
gqcal(extProspective[,1:13])

##### Histogram of Posterior Distribution of the Mean of Odds Ratio (Not in textbook)
hist(extProspective$Odds_r_mu,breaks=200,xlim=c(0,1.0),main="")

##### Histogram of Posterior Predictive Distribution of Odds Ratio (Equation 17.33) (Not in textbook)
Odds_r_aste<-exp(rnorm(length(extProspective$mu_b),extProspective$mu_b,extProspective$sigma_b))
hist(Odds_r_aste,breaks=1000,xlim=c(0,4.0),main="")

#### Table 17.11 Summary of Posterior Distribution and Posterior Predictive Distribution of Average Odds in Prospective Studies
gqcal(exp(extProspective$mu_b),3);  # Table 17.11
gqcal(Odds_r_aste,3)

#### PHC
mean(1.0<Odds_r_aste)
mean(Odds_r_aste<0.8)


########### 17.3.5 Meta-analysis of Retrospective Studies

par<-c("Odds_r","Odds_r_mu","mu_b","sigma_b")                # Parameters  
dataSetRetro <-list(m=nrow(retrospective),
      u1=retrospective$u11,n1=retrospective$u11+retrospective$u21,
      u0=retrospective$u12,n0=retrospective$u12+retrospective$u22)                      # Input

########### Execution by cmdstanr
                                             # Compilation is "modProspective"
csrfitRetro <- modProspective$sample(data = dataSetRetro,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  # MCMC
extRetro<-as_draws_df(csrfitRetro$draws(par))
(colnames(extRetro) <- gsub("\\[|\\]|,", "", colnames(extRetro)))

#save(extRetro,  file="./script/obje/out1704")
#load(file="./script/obje/out1704"); # Load pre-run stan( ) object

#### Table 17.12 Odds/EAP
gqcal(extRetro[,1:10])
Odds_r_aste<-exp(rnorm(length(extRetro$mu_b),extRetro$mu_b,extRetro$sigma_b))

#### Table 17.13 Summary of Posterior Distribution and Posterior Predictive Distribution of Average Odds in Retrospective Studies
gqcal(exp(extRetro$mu_b),3)
gqcal(Odds_r_aste,3)

#####PHC
mean(1.0<Odds_r_aste)
mean(Odds_r_aste<0.8)

#### Figure 17.8 Posterior Distribution of Odds Ratio for Each Retrospective Study
aa<-extRetro[[1]]; for (i in 2:10){aa<-c(aa,extRetro[[i]])}
boxplot(aa~rep(1:10,each=100000),outline=F,horizontal=T,cex.axis=1.5);
abline(v=1.0,lwd=2.0)


################# Practical Assignment
###  1.
#load(file="./script/obje/out1701"); # Load pre-run stan( ) object
p0_5<- (-1)*ext$a/ext$b;       # Age at which there's a 50% chance of becoming a manager
oddshi<-exp(ext$b);              # Odds ratio when age changes by one year
p35<- 1/(1+exp((-1)*(ext$a+ext$b*35))); # Probability of promotion at age 35
gqcal(p0_5)
gqcal(oddshi)

#
#  Parts using function PHC01 are omitted
#

###  2.
#load(file="./script/obje/out1702"); # Load pre-run stan( ) object
p35<- 1/(1+exp((-1)*(ext$a+ext$b*35))); # Correct answer rate after 5 more trainings
gqcal(p35)

#
#  Parts using function PHC01 are omitted
#

###  3.
#load(file="./script/obje/out1704"); # Load pre-run stan( ) object
gqcal(extRetro$Odds_r_mu) 
PHC01(seq(0.25,0.70,0.01),extRetro$Odds_r_mu,cc="gtc",byoga="yes",
      xlab="Average of Odds")
PHC01(seq(0.25,0.70,0.05),extRetro$Odds_r_mu,cc="gtc",byoga="no")
Odds_r_aste<-exp(rnorm(length(extRetro$mu_b),extRetro$mu_b,extRetro$sigma_b))
gqcal(Odds_r_aste,3)
PHC01(seq(0.05,2.1,0.01),extRetro$Odds_r_mu,cc="gtc",byoga="yes",
      xlab="Posterior Predictive Distribution of Odds")
PHC01(seq(0.05,0.3,0.05),extRetro$Odds_r_mu,cc="gtc",byoga="no")



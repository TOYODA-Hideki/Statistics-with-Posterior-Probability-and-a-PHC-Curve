# Chapter 20: Analysis of Covariance/Trend Score
# (n_wd<-getwd())                # Confirmation of working directory
source('myfunc/myfunc.R')      # Loading self-made functions
library(cmdstanr)              # load and attach add-on package cmdstanr
library(posterior)             # load and attach add-on package posterior
library(pROC)
library(Matching)
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points


####################### Functions used in Chapter 20 start here
# Function to draw scatter plots of bivariate data for two groups, changing pch
plot2g<-function(y,x,a,ba=5,xlab="pretest"){
  plot(x,y,pch=ifelse(a-1,1,16),xlab=xlab,ylab="posttest",
       cex.lab=1.5,cex.axis=1.5,xlim=c(min(x)-ba,max(x)+ba),
       ylim=c(min(y)-ba,max(y)+ba))
}
# Function to draw scatter plots for two groups and regression lines with the same slope for both
plot2gequ<-function(y,x,a,ext,ba=5,xlab="pre test"){
  plot(x,y,pch=ifelse(a-1,1,16),xlab=xlab,ylab="post test",
     cex.lab=1.5,cex.axis=1.5,xlim=c(min(x)-ba,max(x)+ba),
     ylim=c(min(y)-ba,max(y)+ba))
  abline(mean(ext$a1),mean(ext$b))     
  abline(mean(ext$a2),mean(ext$b))     
}
# Sample variance
bunsan<-function(x){mean((x-mean(x))^2)}

# Function to draw scatter plot with two regression lines of different slopes
plot2gdef<-function(y,x,a,ext,ba=3){
  plot(x,y,pch=ifelse(a-1,1,16),xlab="pre test",ylab="post test",
     cex.lab=1.5,cex.axis=1.5,xlim=c(min(x)-ba,max(x)+ba),
     ylim=c(min(y)-ba,max(y)+ba))
  abline(mean(ext$a1),mean(ext$b1))     
  abline(mean(ext$a2),mean(ext$b2))     
}
# Function to draw EAP and 95% credible interval of the difference in predicted values of two regression lines with different slopes
plot2gdefdef<-function(x,ext){
  leng=30
  xlabe<-seq(min(x),max(x),length=leng)
  yhat_sa<-matrix(0,nrow(ext),leng)
  yhat_sa_HL<-matrix(0,leng,2)
  for (i in (1:leng)){
    yhat_sa[,i]<-ext$a1-ext$a2+(ext$b1-ext$b2)*xlabe[i]}
  for (i in (1:leng)){
    yhat_sa_HL[i,]<-quantile(yhat_sa[,i],prob=c(0.025,0.975))}
  plot(xlabe,colMeans(yhat_sa),ylim=c(min(yhat_sa_HL),max(yhat_sa_HL)),
    type="l",ylab="",xlab="",lwd=2.0)
  abline(a=0,b=0,lwd=2.0)
  par(new=T)
  plot(xlabe,yhat_sa_HL[,1],ylim=c(min(yhat_sa_HL),max(yhat_sa_HL)),
    type="l",ylab="",xlab="",lty=2,lwd=2.0)
  par(new=T)
  plot(xlabe,yhat_sa_HL[,2],ylim=c(min(yhat_sa_HL),max(yhat_sa_HL)),
    type="l",ylab="postの予測値の差",xlab="pre値",lty=2,lwd=2.0)
}
####################### Functions used in Chapter 7 end here

##### 20.2 Analysis of Covariance with Common Slope

##### 20.2.1 Analysis of "Data1"

##### Table 20.1 Average Values of "Data1"
learning01<-read.csv("dat/learning01.csv", header = TRUE)
round(tapply(learning01$post, learning01$subject, mean), 1)
round(tapply(learning01$pre, learning01$subject, mean), 1)

##### Figure 20.1 Scatter Plot of "Data1"
plot2g(y=learning01$post, x=learning01$pre, a=learning01$subject, ba=3)
#dev.copy2eps(file="z2001.eps", family='Japan1')

############################# Analysis of Covariance (Model with Common Slope) in Stan Code
ANCOVAequ<-'
data {
  int<lower=0> n; // Number of data points 
  vector[n] y; // Dependent variable   
  vector[n] x; // Independent variable   
  array[n] int<lower=0, upper=2> k; // Classification variable 
}
parameters {
  vector[2] a; // Intercept
  real b; // Regression coefficient 
  real<lower=0> sigma; // Error standard deviation
}
transformed parameters {
  vector[n] yhat; // Predicted values yhat  
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b;
  } // Equations (20.4), (20.5)
}
model {
  y ~ normal(yhat, sigma); // Equations (20.1), (20.2), (20.3)
}
';

par<-c("a","b","sigma");                                                  # Parameters
dataSetequ01<-list(n=nrow(learning01), y=learning01[,2], x=learning01[,1], k=learning01[,3]); # Input

########### Execution with cmdstanr
modfileequ01 <- write_stan_file(ANCOVAequ)         # Write to temporary file
modequ01 <- cmdstan_model(modfileequ01)            # Compile
csrfitequ01<-modequ01$sample(data=dataSetequ01, chains=5, iter_sampling=20000,
                    iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extequ01<-as_draws_df(csrfitequ01$draws(par))
(colnames(extequ01) <- gsub("\\[|\\]|,", "", colnames(extequ01)))

#save(extequ01,  file="./script/obje/out2001")
#load(file="./script/obje/out2001");#Load the object previously madewith stan()

##### Table 20.2 Summary of Posterior Distribution of Parameters for "Data1" with Common Slope Model
gqcal(extequ01$a1,2)
gqcal(extequ01$a2,2)
gqcal(extequ01$b, 2)
gqcal(extequ01$sigma,2)

##### Figure 20.2 Regression Lines of "Data1" (Common Slope)
plot2gequ(y=learning01$post, x=learning01$pre,
          a=learning01$subject, ext=extequ01, ba=3);
#dev.copy2eps(file="z2002.eps", family='Japan1')

a_sa01<-extequ01$a1-extequ01$a2; # Posterior distribution of the difference in intercepts
gqcal(a_sa01,1);                                      

##### Figure 20.3 Posterior Distribution of the Difference in Intercepts for "Data1"
hist(a_sa01, breaks=100, ylab="", xlab="Difference in intercepts", main="", cex.lab=2, cex.axis=1.5)
#dev.copy2eps(file="z2003.eps", family='Japan1')

####### 20.2.2 Analysis of "Data2"
##### Table 20.3 Average Values of "Data2"
learning02<-read.csv("dat/learning02.csv", header = TRUE)
round(tapply(learning02$post, learning02$subject, mean), 1)
round(tapply(learning02$pre, learning02$subject, mean), 1)

##### Figure 20.4 Scatter Plot of "Data2"
plot2g(y=learning02$post, x=learning02$pre, a=learning02$subject, ba=3)
#dev.copy2eps(file="z2004.eps", family='Japan1')

par<-c("a","b","sigma");                                                  # Parameters
dataSetequ02<-list(n=nrow(learning02), y=learning02[,2], x=learning02[,1], k=learning02[,3]); # Input

########### Execution with cmdstanr
### Using the compilation result of 'modequ01' (Model with Common Slope)
csrfitequ02<-modequ01$sample(data=dataSetequ02, chains=5, iter_sampling=20000,
             iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extequ02<-as_draws_df(csrfitequ02$draws(par))
(colnames(extequ02) <- gsub("\\[|\\]|,", "", colnames(extequ02)))

#save(extequ02, file="./script/obje/out2002")
#load(file="./script/obje/out2002");#Load the object previously made withstan()


##### Table 20.4 Summary of Posterior Distribution of Parameters for "Data2" with Common Slope Model
gqcal(extequ02$a1,2)
gqcal(extequ02$a2,2)
gqcal(extequ02$b, 2)
gqcal(extequ02$sigma,2)

##### Figure 20.5 Regression Lines of "Data2" (Common Slope)
plot2gequ(y=learning02$post, x=learning02$pre, a=learning02$subject, ext=extequ02, ba=3)
#dev.copy2eps(file="z2005.eps", family='Japan1')

a_sa02<-extequ02$a1-extequ02$a2;  # Posterior distribution of the difference in intercepts
gqcal(a_sa02,1);                                         

##### Figure 20.6 Posterior Distribution of the Difference in Intercepts for "Data2"
hist(a_sa02, breaks=100, ylab="", xlab="Difference in intercepts",   
     main="", cex.lab=2, cex.axis=1.5)
#dev.copy2eps(file="z2006.eps", family='Japan1')

mean(abs(a_sa02)<1.0);                                # Probability within ROPE
mean(abs(a_sa02)<1.5);                                

########## 20.3 Analysis of Covariance with Different Slopes

###### 20.3.1 Analysis of "Data3"
###### Table 20.5 Average Values of "Data3"
learning03<-read.csv("dat/learning03.csv", header = TRUE)
round(tapply(learning03$post, learning01$subject, mean), 1)
round(tapply(learning03$pre, learning01$subject, mean), 1)

###### Figure 20.7 Scatter Plot of "Data3"
plot2g(y=learning03$post, x=learning03$pre, a=learning02$subject, ba=3); # Figure 20.7
#dev.copy2eps(file="z2007.eps", family='Japan1')

############################# Analysis of Covariance (Model with Different Slopes) in Stan Code
ANCOVAdef<-'
// Analysis of Covariance (Model with Different Slopes)
data {
  int<lower=0> n; // Number of data points 
  vector[n] y; // Dependent variable   
  vector[n] x; // Independent variable   
  array[n] int<lower=0, upper=2> k; // Classification variable 
}
parameters {
  vector[2] a; // Intercept
  vector[2] b; // Regression coefficients 
  real<lower=0> sigma; // Error standard deviation
}
transformed parameters {
  vector[n] yhat; // Predicted values yhat  
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b[k[i]];
  } // Equations (20.7), (20.8)
}
model {
  y ~ normal(yhat, sigma); // Equation (20.9)
}
';

par<-c("a","b","sigma")                       
dataSetdef03 <-list(n=nrow(learning03), y=learning03[,2], x=learning03[,1], k=learning03[,3])

########### Execution with cmdstanr
modfiledef03 <- write_stan_file(ANCOVAdef)          # Write to temporary file
moddef03 <- cmdstan_model(modfiledef03)             # Compile
csrfitdef03<-moddef03$sample(data=dataSetdef03, chains=5, iter_sampling=20000,
           iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extdef03<-as_draws_df(csrfitdef03$draws(par))
(colnames(extdef03) <- gsub("\\[|\\]|,", "", colnames(extdef03)))

#save(extdef03, file="./script/obje/out2003")
#load(file="./script/obje/out2003");#Load the object previously made withstan()


###### Table 20.6 Summary of Posterior Distribution of Parameters for "Data3" with Different Slope Model
gqcal(extdef03$a1,2)
gqcal(extdef03$a2,2)
gqcal(extdef03$b1,2)
gqcal(extdef03$b2,2)
gqcal(extdef03$sigma,2)

###### Figure 20.8 Regression Lines of "Data3" (Different Slopes)
plot2gdef(y=learning03$post, x=learning03$pre,       
          a=learning03$subject, ext=extdef03, ba=3)
#dev.copy2eps(file="z2008.eps", family='Japan1')

###### Figure 20.9 Point Estimate and 95% Prediction Interval of the Difference in Predicted Values for "Data3"
plot2gdefdef(x=learning03[,1], extdef03);                  # Figure 20.9
#dev.copy2eps(file="z2009.eps", family='Japan1')

######## 20.3.2 Reanalysis of "Data1" and "Data2"

# Reanalysis of "Data1" with Different Slope Model
########### Execution with cmdstanr
### Using the compilation result of 'moddef03' and data 'dataSetequ01'
csrfitdef01<-moddef03$sample(data=dataSetequ01, chains=5, iter_sampling=20000,
                    iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extdef01<-as_draws_df(csrfitdef01$draws(par))
(colnames(extdef01) <- gsub("\\[|\\]|,", "", colnames(extdef01)))

#save(extdef01, file="./script/obje/out2004")
#load(file="./script/obje/out2004");#Load the object previously made withstan()

#### Summary of Posterior Distribution of Parameters for "Data1" with Different Slope Model (not in the textbook)
gqcal(extdef01$a1,2)
gqcal(extdef01$a2,2)
gqcal(extdef01$b1,2)
gqcal(extdef01$b2,2)
gqcal(extdef01$sigma,2)

###### Figure 20.10 Regression Lines of "Data1" (Different Slopes)
plot2gdef(y=learning01$post, x=learning01$pre, a=learning01$subject, ext=extdef01)
#dev.copy2eps(file="z2010.eps", family='Japan1')

###### Figure 20.11 Point Estimate and 95% Prediction Interval of the Difference in Predicted Values for "Data1"
plot2gdefdef(x=learning01[,1], extdef01); abline(h=4.8)
#dev.copy2eps(file="z2011.eps", family='Japan1')

# Reanalysis of "Data2" with Different Slope Model
########### Execution with cmdstanr
### Using the compilation result of 'moddef03' and data 'dataSetequ02'
csrfitdef02<-moddef03$sample(data=dataSetequ02, chains=5, iter_sampling=20000,
                   iter_warmup = 1000, parallel_chains = 5, seed=1234)  # MCMC
extdef02<-as_draws_df(csrfitdef02$draws(par))
(colnames(extdef02) <- gsub("\\[|\\]|,", "", colnames(extdef02)))

#save(extdef02, file="./script/obje/out2005")
#load(file="./script/obje/out2005");#Load the object previously made withstan()

###### Figure 20.12 Regression Lines of "Data2" (Different Slopes)
plot2gdef(y=learning02$post, x=learning02$pre, a=learning02$subject, ext=extdef02)
#dev.copy2eps(file="z2012.eps", family='Japan1')

###### Figure 20.13 Point Estimate and 95% Prediction Interval of the Difference in Predicted Values for "Data2"
plot2gdefdef(x=learning02[,1], extdef02)
#dev.copy2eps(file="z2013.eps", family='Japan1')

# Self-study 1. From Figure 20.12, it is visually apparent that the two regression lines for "Data2" in the model with different slopes are almost identical. As a precaution, show the summary of the posterior distribution corresponding to Table 20.6.

######## 20.3.3 Posterior Distribution of the Difference in Slopes
###### Table 20.7 Summary of the Posterior Distribution of the Difference in Slopes
b_sa01<-extdef01$b1-extdef01$b2
gqcal(b_sa01)
b_sa02<-extdef02$b1-extdef02$b2
gqcal(b_sa02)
b_sa03<-extdef03$b1-extdef03$b2
gqcal(b_sa03)

# Self-study 2. In the summary of the posterior distribution of the difference in slopes in Table 20.7 of 20.3.3, the summary of the posterior distribution of the difference in slopes is shown. To visually confirm the relationship of these distributions with 0, draw histograms for all three.

###### 20.4 Propensity Score
###### 20.4.2 Data Description
###### Table 20.8 Data from an Experiment Testing the Effectiveness of a New Teaching Method Using Existing Classes at a Cram School
learning04<-read.csv("dat/propen.csv", header = TRUE)
head(learning04)
tail(learning04)

###### Frequency of each variable 
table(learning04$group)
table(learning04$high_school)
table(learning04$income)
table(learning04$mother_education)
table(learning04$father_education)

###### Table 20.9 Average Values of the Experimental and Control Groups
round(tapply(learning04$score, learning04$group, mean), 1)
round(tapply(learning04$high_school, learning04$group, mean), 1)
round(tapply(learning04$income, learning04$group, mean), 1)
round(tapply(learning04$mother_education, learning04$group, mean), 1)
round(tapply(learning04$father_education, learning04$group, mean), 1)

################################# Analysis with Propensity Scores in Stan Code
propensity01<-'
data {
  int<lower=0> n; // Number of data points 
  int<lower=0> p; // Number of predictor variables 
  vector[n] y; // Dependent variable   
  matrix[n, p] x; // Predictor variables   
  array[n] int<lower=1, upper=2> k; // Classification variable 1 for experimental group, 2 for control group
}
transformed data {
  array[n] int<lower=0, upper=1> kk; // Classification variable 
  for (i in 1 : n) {
    kk[i] = k[i] * (-1) + 2;
  }
}
parameters {
  vector[2] a; // Intercept
  real b; // Regression coefficient 
  real alpha; // Logistic regression intercept
  vector[p] beta; // Logistic regression coefficients
  real<lower=0> sigma; // Error standard deviation
}
transformed parameters {
  vector[n] yhat; // yhat  
  vector<lower=0, upper=1>[n] propensity; // Propensity score
  for (i in 1 : n) {
    propensity[i] = inv_logit(alpha + x[i,  : ] * beta);
    yhat[i] = a[k[i]] + propensity[i] * b;
  }
}
model {
  for (i in 1 : n) {
    kk[i] ~ bernoulli(propensity[i]); // Bernoulli distribution model
    y[i] ~ normal(yhat[i], sigma); // Normal distribution model
  }
}
generated quantities {
  
}
';

par<-c("a","b","sigma","alpha","beta","propensity")    # Parameters
dataSetpropen01 <-list(n=nrow(learning04), p=4, y=learning04[,2],
                       x=learning04[,c(3:6)], k=learning04[,1]) # Input

########### Execution with cmdstanr
modfilepropen01 <- write_stan_file(propensity01)    # Write to temporary file
modpropen01 <- cmdstan_model(modfilepropen01)       # Compile
csrfitpropen01 <- modpropen01$sample(data = dataSetpropen01, chains = 5,
  iter_sampling=20000, iter_warmup=1000, parallel_chains= 5, seed=1234)  # MCMC
extpropen01<-as_draws_df(csrfitpropen01$draws(par))
(colnames(extpropen01) <- gsub("\\[|\\]|,", "", colnames(extpropen01)))

#save(extpropen01, file="./script/obje/out2006")
#load(file="./script/obje/out2006");#Load the object previously made withstan()


###### Table 20.10 Summary of Posterior Distribution of Parameters for Analysis of Covariance with Propensity Score as Covariate
gqcal(extpropen01[,1:13],2)


###### Figure 20.14 Histogram of Propensity Scores by Group
propen01<-apply(extpropen01[,10:159], 2, mean);  # EAP of propensity scores
classes<-seq(0, 1, 0.1)
barplot(table(learning04$group, cut(propen01, classes, include.lowest = T, right = F)),
        xlab="propensity score", axisnames = F, cex.lab=2.0, cex.axis=1.5)
#dev.copy2eps(file="z2014.eps", family='Japan1')

###### Figure 20.15 ROC Curve
(ROC1<-roc(learning04$group ~ propen01, ci=T))
plot(ROC1, cex.lab=2.0, cex.axis=1.5)
#dev.copy2eps(file="z2015.eps", family='Japan1')

###### Figure 20.16 Analysis of Covariance with Propensity Score as a Covariate
plot2gequ(y=learning04$score, x=propen01, a=learning04$group,
          ext=extpropen01, ba=0, xlab="propensity score")
#dev.copy2eps(file="z2016.eps", family='Japan1')

###### Figure 20.17 Posterior Distribution of the Difference in Intercepts
a_sa04<-extpropen01$a1-extpropen01$a2
gqcal(a_sa04, 1);
hist(a_sa04, br=100, cex.lab=2.0, cex.axis=1.5, xlab="a1 - a2", ylab="", main="")
#dev.copy2eps(file="z2017.eps", family='Japan1')

####### Probability of 'Difference in Intercepts being within ROPE (1.0)' 
mean(abs(a_sa04)<1)

# The following practice will not be performed
set.seed(1234)
match<-Match(Tr= learning04$group*-1+2, X=propen01, M=1, caliper=0.2, replace=F)
match

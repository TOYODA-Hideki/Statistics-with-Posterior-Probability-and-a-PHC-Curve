###The working directory should be 'PPandPHC'
# Chapter 18
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)              #Load and attach add-on package cmdstanr
library(posterior)             #Load and attach add-on package posterior
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points


#################18.1 Poisson Distribution

###Figure 18.1 Probability Function of the Poisson Distribution
par(mfrow=c(1,2))
  u<-c(0:5);p061<-dpois(u, 0.61);names(p061)<-u
  barplot(p061); text(4,0.3,"λ=0.61",cex=2.0)
  x<-c(0:10)
  plot(x, dpois(x, 2.0), type="o",xlab="",ylab="",xlim=c(0,10),ylim=c(0,0.3))
  for (lambda in seq(3.0,6.0,1.0)){
    par(new=T)
    plot(x,dpois(x,lambda),type="o",xlab="",ylab="",xlim=c(0,10),ylim=c(0,0.3))
  }
  text(c(2,3,4,5,6),c(0.29,0.25,0.22,0.19,0.18),
    paste("λ=",2:6,sep=""),cex=1.5)
par(mfrow=c(1,1));                                    #End of Figure 18.1
#dev.copy2pdf(file="fig1801.pdf",family="Japan1")

####Approximates normal distribution as lambda grows (not in textbooks)
size<-100; x<-c(0:size);
plot(x, dpois(x, 50), type="o",main=paste("Poisson Distribution　lambda:",50))

###18.1.1 Numerical Example (Number of Soldiers Killed by Horse Kicks)
###Table 18.1: Number of Soldiers Killed by Horse Kicks in the Prussian Army. 
Horse<-c(dpois(0, 0.61),dpois(1, 0.61),dpois(2, 0.61),
         dpois(3, 0.61),dpois(4, 0.61),dpois(5, 0.61));
print(0:5);print(c(109, 65, 22, 3, 1, 0));round(Horse,4);round(Horse*200,1)

###18.1.2 Numerical Example (Probability of Winning Prizes)
ice<-c(dpois(0, 0.25),dpois(1, 0.25),dpois(2, 0.25),           #Equation (18.7)
   dpois(3, 0.25),dpois(4, 0.25));round(ice,3)

#################18.2 Estimation of the Poisson Distribution
####18.2.1 Number of Bomb Hits and Number of Sections
bomb<-c(rep(0,229),rep(1,211),rep(2,93),rep(3,35),rep(4,7),rep(5,1))

###Table 18.2 First half of number of bomb hits and posterior predictive checks in London
table(bomb)

barplot(table(bomb),xlab="Number of Hits",ylab="Number of Sections"); #Histogram (not in textbook)
round(mean(bomb),3);                              #mean　
round(mean((bomb-mean(bomb))^2),3);               #variance 　


################# stan code for estimation of Poisson distribution
Poisson<-'
data {
  int<lower=0> n;          //Number of data
  array[n] int<lower=0> u; //count data 
}
parameters {
  real lambda;             //parameter
}
model {
  for (i in 1 : n) {
    u[i] ~ poisson(lambda);
  } //Poisson distribution (18.5) Equation
}
generated quantities {
  int<lower=0> uaste;           //predictive distribution  
  uaste = poisson_rng(lambda); //Poisson random number (18.8)Equation
}
';

par<-c("lambda","uaste")                  #parameter
dataSetbomb <-list(n=length(bomb),u=bomb) #input

########### Execution by cmdstanr
modfilebomb <- write_stan_file(Poisson)                     #Temporary file to be exported
modbomb <- cmdstan_model(modfilebomb)                       #compile
csrfitbomb <- modbomb$sample(data = dataSetbomb,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)       #MCMC
extbomb<-as_draws_df(csrfitbomb$draws(par))
(colnames(extbomb) <- gsub("\\[|\\]|,", "", colnames(extbomb)))

#save(extbomb,  file="./script/obje/out1801")
#load(file="./script/obje/out1801");       #Load object created with stan( ) in advance

########## Table 18.3: Summary of the Posterior Distribution of Parameters and Generated Quantities
sdbomb<-sqrt(extbomb$lambda);                             # Standard Deviation
gqcal(extbomb$lambda,3);                         # Mean Posterior Distribution 
gqcal(sdbomb   ,3);                                # sd Posterior Distribution 

###Table 18.2 Second half of number of bomb hits and posterior predictive checks in London
round(table(extbomb$uaste)/length(extbomb$uaste),3);        #probability
round(table(extbomb$uaste)/length(extbomb$uaste)*576,1);    #posterior predictive value

#############18.2.2　Number of Absentees in Class
#####Table 18.4: Number of Absentees in "Psychology Statistics" Classes
(absent<-c(0,1,2,0,0,1,2,0,2,0,1,0,0,1))

round(mean(absent),3);                
round(mean((absent-mean(absent))^2),3); 

par<-c("lambda","uaste")                        #parameter
dataSetabsent <-list(n=length(absent),u=absent) #input

###########Execution by cmdstanr
#    Use compiled "modbomb".
csrfitabsent  <- modbomb$sample(data = dataSetabsent,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)               #MCMC
extabsent<-as_draws_df(csrfitabsent$draws(par))
(colnames(extabsent) <- gsub("\\[|\\]|,", "", colnames(extabsent)))

#save(extabsent, file="./script/obje/out1802")
#load(file="./script/obje/out1802"); #Load object created with stan( ) in advance

######Table 18.5: Summary of the Posterior Distribution of Parameters and Generated Quantities
sdabsent<-sqrt(extabsent$lambda)
gqcal(extabsent$lambda,3)
gqcal(sdabsent        ,3)

#####Table 18.6: Posterior Predictive Checking for Absentees
round(table(extabsent$uaste)/length(extabsent$uaste),3);        #Table 18.6
round(table(extabsent$uaste)/length(extabsent$uaste)*14,1)

#Table 18.7 Summary of the Posterior Distribution of Parameters and Generated Quantities
#           (Comparison of the distribution of posterior predictions)
preprob2bomb<-dpois(2,extbomb$lambda)
preprob2absent<-dpois(2,extabsent$lambda)
prevalue2bomb<-preprob2bomb*576
pre2absent<-preprob2absent*14
gqcal(preprob2bomb,3)
gqcal(prevalue2bomb,1)
gqcal(preprob2absent,3)
gqcal(pre2absent,2)


############# 18.3   Comparison of Two Poisson Distributions

####### Table 18.8: Distribution of Deviations in Mirror Drawing Task
u1<-c(1,0,2,3,2,3,3,1,4,3,1,0,0,1,2,1,3,2,1,4,1,1,1,1,3,2,3,1,1,5)
u2<-c(1,3,2,7,4,4,2,1,0,2,1,2,1,0,1,5,5,3,4,2,1,3,1,1,1,2,3,1,1,3,0)
mean(u1);mean(u2)
table(u1);table(u2)

################# stan code for two Poisson distribution guesses
Poisson2<-'
data {
  int<lower=0> n1;            //Number of data 
  int<lower=0> n2;            //Number of data 
  array[n1] int<lower=0> u1; //count data 
  array[n2] int<lower=0> u2; //count data 
}
parameters {
  real lambda1; //parameter1
  real lambda2; //parameter2
}
model {
  for (i in 1 : n1) {
    u1[i] ~ poisson(lambda1);
  } //Equation (18.10)
  for (i in 1 : n2) {
    u2[i] ~ poisson(lambda2);
  } //Equation (18.10)
}
';
par<-c("lambda1","lambda2")                                    #parameter
dataSetmirror2 <-list(n1=length(u1),n2=length(u2),u1=u1,u2=u2) #input

###########Execution by cmdstanr
modfilemirror2 <- write_stan_file(Poisson2)                #Temporary file to be exported
modmirror2 <- cmdstan_model(modfilemirror2)                #compile
csrfitmirror2 <- modmirror2$sample(data = dataSetmirror2,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)      #MCMC
extmirror2<-as_draws_df(csrfitmirror2$draws(par))
(colnames(extmirror2) <- gsub("\\[|\\]|,", "", colnames(extmirror2)))

#save(extmirror2, file="./script/obje/out1803")
#load(file="./script/obje/out1803"); #Load object created with stan( ) in advance

#######Table 18.9: Summary of the Posterior Distribution of Parameters 
#                  of the Two Group  Poisson Model
gqcal(extmirror2$lambda1,3)
gqcal(extmirror2$lambda2,3)

##### Table 18.10: PHC within ROPE, PHC not within ROPE
sequ<-seq(0,1.2,0.2);a<-extmirror2$lambda2;b<-extmirror2$lambda1
  PHC01(sequ,a,b,cc="rope",byoga="no")
1-PHC01(sequ,a,b,cc="rope",byoga="no")

#### Table 18.11: PHC indicating that λ2 is larger. 
PHC01(sequ,a,b,cc="gtc", byoga="no")


####################18.4 Poisson Regression

###### Table 18.12: Number of Deviations in the Mirror Drawing Task
(X<-data.frame(
  NoTrials=c(1:15),
  Nodeviations=c(6,4,3,2,2,2,1,1,1,0,0,2,1,0,0)))

##### Figure 18.3: Number of Deviations over 15 Trials
plot(X[,1],X[,2],type="o",xlab="Number of Trials",ylab="Number of Deviations")
  abline(a=coef(glm(X[,2]~X[,1]))[[1]],coef(glm(X[,2]~X[,1]))[[2]])
#dev.copy2pdf(file="fig1803.pdf")

##### 18.4.1 Exponential Transformation   
#####        Conversion from whole real numbers to positive real numbers  
##### Figure 18.4: Exponential Transformation exp( )
plot(exp,-2,2,xlab="a+bx",ylab="λ",main="", lwd=2.0,cex.lab=1.5)

####################Poisson regression stan code
Poisson_reg<-'
data {
  int<lower=0> n;         //Number of data 
  array[n] int<lower=0> y; //criterion variable (count data)   
  vector[n] X;             //Predictor variable vector  
}
parameters {
  real a; //intercept
  real b; //partial regression coefficient
}
transformed parameters {
  vector<lower=0>[n] lambda; //lambda 
  for (i in 1 : n) {
    lambda[i] = exp(a + X[i] * b);
  } //Equation (18.14), (18.19)
}
model {
  for (i in 1 : n) {
    y[i] ~ poisson(lambda[i]);
  } //Poisson distribution model Equation (18.13),(18.15)
}
generated quantities {
  array[n] int<lower=0> yaste; //predictive distribution 
  for (i in 1 : n) {
    yaste[i] = poisson_rng(lambda[i]);
  } //Poisson distribution model (18.20) Equation
}
';

par<-c("a","b","lambda","yaste")                         #parameter
dataSetregression <-list(n=nrow(X),y=X[,2], X=X[,1] )    #input

###########Execution by cmdstanr
modfileregression <- write_stan_file(Poisson_reg)            #Temporary file to be exported
modregression <- cmdstan_model(modfileregression)            #compile
csrfitregression <- modregression$sample(data = dataSetregression, chains=5,iter_sampling=20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)        #MCMC
extregression<-as_draws_df(csrfitregression$draws(par))
(colnames(extregression) <- gsub("\\[|\\]|,", "", colnames(extregression)))

#save(extregression,  file="./script/obje/out1804")
#load(file="./script/obje/out1804");       #Load object created with stan( ) in advance

#### generated quantities ：
#### Rate of change of parameter when predictor variable is changed by one unit
expb<-exp(extregression$b);    #Equation (18.21)

###### Table 18.13: Summary of the Posterior Distribution of Parameters 
######              and Generated Quantities for the Poisson Regression Model
gqcal(extregression$a,3)
gqcal(extregression$b,3)
gqcal(extregression$lambda7,3)
gqcal(extregression$yaste7,3)
gqcal(expb,3)

###### Figure 18.5: Poisson Regression of the Number of Deviations
lambda_i<-apply(extregression[,3:17],2,mean)
plot(X[,1],X[,2],type="o",xlab="Number of Trials",ylab="Number of Deviations",xlim=c(1,15),ylim=c(0,6))
   par(new=T);                                        
plot(X[,1],lambda_i,type="l",xlab="",ylab="",lwd=2,xlim=c(1,15),ylim=c(0,6))
#dev.copy2pdf(file="fig1805.pdf")

############# 18.5 Log-linear Model
############# 18.5.1 Structuring Frequencies
############# Table 18.14: Survey Results Investigating the Advertising Effectiveness of Ad B
(u<-matrix(c(29,21,50,145),2,2,byrow=T))

##################Analysis of contingency table 01, decomposition of frequencies (2*2) stan code
Poisson_Cross01<-'
data {
  array[2, 2] int<lower=0> u; //criterion variable (count data)  
}
parameters {
  real mu; //total mean 
  real A1; //main effect of factor A
  real B1; //main effect of factor B
  real AB11; //interaction
}
transformed parameters {
  matrix<lower=0>[2, 2] lambda;            //lambda 
  lambda[1, 1] = exp(mu + A1 + B1 + AB11); //Equation (18.23), (18.24)
  lambda[2, 1] = exp(mu - A1 + B1 - AB11);
  lambda[1, 2] = exp(mu + A1 - B1 - AB11);
  lambda[2, 2] = exp(mu - A1 - B1 + AB11);
}
model {
  u[1, 1] ~ poisson(lambda[1, 1]); //Equation (18.25)
  u[2, 1] ~ poisson(lambda[2, 1]);
  u[1, 2] ~ poisson(lambda[1, 2]);
  u[2, 2] ~ poisson(lambda[2, 2]);
}
';

par<-c("mu","A1","B1","AB11","lambda")    #parameter
dataSetcontable01<-list(u=u )             #input

###########Execution by cmdstanr
modfilecontable01 <- write_stan_file(Poisson_Cross01)  #Temporary file to be exported
modcontable01 <- cmdstan_model(modfilecontable01)        #compile
csrfitcontable01 <- modcontable01$sample(data = dataSetcontable01,chains = 5,
  iter_sampling=20000, iter_warmup=1000, parallel_chains=5, seed=1234)
extcontable01<-as_draws_df(csrfitcontable01$draws(par))
(colnames(extcontable01) <- gsub("\\[|\\]|,", "", colnames(extcontable01)))

#save(extcontable01,  file="./script/obje/out1805")
#load(file="./script/obje/out1805"); #Load object created with stan( ) in advance


##### Table 18.15: Summary of the Posterior Distribution of Parameters 
#####              and Generated Quantities for Advertisement Effect Analysis
gqcal(extcontable01$mu)
gqcal(extcontable01$A1)
gqcal(extcontable01$B1)
gqcal(extcontable01$AB11)


############# 18.5.2 Structuring Ratios

############# Table 18.16: Number of Victims among Crew and Passengers of the Titanic
(u<-matrix(c(4,118, 13,154, 106,422, 3,670),4,2,byrow=T));    #Victims
(n<-matrix(c(145,180, 106,179, 196,510, 23,862),4,2,byrow=T));#the number of people

#################Analysis of contingency table 02, decomposition of ratios (a*b) stan code
Poisson_Cross02<-'
functions {
  vector zero_sum_vector(int a, vector m1A) {
    //main effect zero sum
    vector[a] muA;
    for (i in 1 : (a - 1)) {
      muA[i] = m1A[i];
    }                   //Substitute as is up to a-1
    muA[a] = -sum(m1A); //Substitute -1 times the sum for the a-th
    return muA;
  }
  matrix zero_sum_matrix(int a, int b, matrix m1AB) {
    //Interaction Zero sum
    vector[a - 1] m1a; //Vector to insert the sum of rows
    vector[b - 1] m1b; //Vector to insert column sums
    matrix[a, b] muAB; //Returning matrix
    for (i in 1 : (a - 1)) {
      m1a[i] = 0.0;
    } //Zero Clear
    for (j in 1 : (b - 1)) {
      m1b[j] = 0.0;
    } //Zero Clear
    for (i in 1 : (a - 1)) {
      for (j in 1 : (b - 1)) {
        muAB[i, j] = m1AB[i, j];
      }
    } //Assign as is
    for (i in 1 : (a - 1)) {
      for (j in 1 : (b - 1)) {
        m1a[i] = m1a[i] + m1AB[i, j];
      }
    } //make the sum of a column
    for (j in 1 : (b - 1)) {
      for (i in 1 : (a - 1)) {
        m1b[j] = m1b[j] + m1AB[i, j];
      }
    } //make the sum of a row
    for (i in 1 : (a - 1)) {
      muAB[i, b] = (-1) * m1a[i];
    } //Substitute -1 times the sum of the rows in column b
    for (j in 1 : (b - 1)) {
      muAB[a, j] = (-1) * m1b[j];
    } //Substitute -1 times the sum of the columns in row a.
    muAB[a, b] = sum(m1a); //Total assignment to row a, column b elements, note: not negative
    return muAB;
  }
}
data {
  int<lower=0> a; //A level-number
  int<lower=0> b; //B level-number
  array[a, b] int<lower=0> u; //Observation frequency  
  array[a, b] int<lower=0> n; //number of trials  
}
parameters {
  real mu; //全平均
  vector[a - 1] m1A; //A average one less
  vector[b - 1] m1B; //B average one less
  matrix[a - 1, b - 1] m1AB; //One less interaction
}
transformed parameters {
  matrix<lower=0>[a, b] lambda; //lambda 
  vector[a] A; //A main effect
  vector[b] B; //B main effect
  matrix[a, b] AB; //AB interaction
  matrix[a, b] p; //Ratio
  A = zero_sum_vector(a, m1A); //constraints Equation (18.28)
  B = zero_sum_vector(b, m1B); //constraints Equation (18.29)
  AB = zero_sum_matrix(a, b, m1AB); //constraints Equation (18.30),(18.31)
  for (i in 1 : a) {
    for (j in 1 : b) {
      p[i, j] = exp(mu + A[i] + B[j] + AB[i, j]); //Equation (18.27), experimental design
      lambda[i, j] = n[i, j] * p[i, j]; //Equation(18.26), Structure of lambda
    }
  }
}
model {
  for (i in 1 : a) {
    for (j in 1 : b) {
      u[i, j] ~ poisson(lambda[i, j]); //Equation (18.26), Poisson model
    }
  }
}
';

par<-c("mu","A","B","AB","p")                             #parameter
dataSetcontable02 <-list(a=nrow(u) ,b=ncol(u), u=u ,n=n )   #input

###########Execution by cmdstanr
modfilecontable02 <- write_stan_file(Poisson_Cross02)  #Temporary file to be exported
modcontable02 <- cmdstan_model(modfilecontable02)        #compile
csrfitcontable02 <- modcontable02$sample(data = dataSetcontable02,chains = 5,
  iter_sampling=20000, iter_warmup=1000, parallel_chains=5, seed=1234)
extcontable02<-as_draws_df(csrfitcontable02$draws(par))
(colnames(extcontable02) <- gsub("\\[|\\]|,", "", colnames(extcontable02)))

#save(extcontable02,  file="./script/obje/out1806")
#load(file="./script/obje/out1806"); #Load object created with stan( ) in advance

##Table 18.17: Summary of Posterior Distribution of Parameters 
##             and Generated Quantities for Titanic Incident Victims Analysis
gqcal(extcontable02[1:23])

###Table 18.18: Probability of Levels & Interaction Effects Being Greater (or Smaller) Than Zero
#Table1 18.18 left
  round(mean(extcontable02$A1>0),3)
  round(mean(extcontable02$A2>0),3)
  round(mean(extcontable02$A3>0),3)
  round(mean(extcontable02$A4>0),3)
#Table1 18.18 center
  round(mean(extcontable02$B1>0),3)
#Table 18.18 right
  round(mean(extcontable02$AB11>0),3)
  round(mean(extcontable02$AB21>0),3)
  round(mean(extcontable02$AB31>0),3)
  round(mean(extcontable02$AB41>0),3)


###Table 18.19: PHC of Row Proportions Being Greater Than Column Proportions
PHC02(0,extcontable02[,c(16,20,17,21,18,22,19,23)])


################# Practical Assignment
p_11<-extcontable02$p11
p_31<-extcontable02$p31
gqcal(p_11)
gqcal(p_31)
PHC01(seq(0,1,0.01),p_31-p_11,cc="gtc",byoga="yes",xlab="Mortality Risk Difference")
PHC01(seq(0.38,0.48,0.01),p_31-p_11,cc="gtc",byoga="no")

PHC01(seq(1,30,1),p_31/p_11,cc="gtc",byoga="yes",xlab="Mortality Risk Ratios")
PHC01(seq(5,15,1),p_31/p_11,cc="gtc",byoga="no")


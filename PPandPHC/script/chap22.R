###The working directory should be 'PPandPHC'
# Chapter 22
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)              #Load and attach add-on package cmdstanr
library(posterior)             #Load and attach add-on package posterior
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) # Definition of probability points

###### 22. 1 Models with Distributions Assumed for Intercepts and Regression Coefficients
###### 22. 1. 1 Experiment on Perceived Length

###### Table22.1 Experiment on Perceived Length
(Epl11<-read.csv("dat/Epl11.csv", header = TRUE))

###### Figure22.1 Experiment on Perceived Length (Data for Two Students)
X<-Epl11[((Epl11$Student==17)|(Epl11$Student==1)), ]
plot(X$Actual.Measurement,X$Visual.Estimation,pch=X$Student,
  xlab="Actual Measurement",ylab="Visual Estimation")
  legend(50,170,c("Student1 ","Student17 "),pch=c(1,17),cex=1.5,text.width=20)
#dev.copy2pdf(file="fig2201.pdf")

###### Figure22.2 Experiment on Perceived Length (Data for 25 Students)
plot(Epl11$Actual.Measurement,Epl11$Visual.Estimation,pch=Epl11$Student,
     xlab="Actual Measurement",ylab="Visual Estimation")   # Figure22.2
#dev.copy2pdf(file="fig2202.pdf")

############# stan code for model assuming distribution 
############# for intercept and regression coefficients
HLM01<-'
data {
  int<lower=0> n;//Number of All Data 
  int<lower=0> J;//Number of Groups
  vector[n] y;//Criterion Variable   
  vector[n] x;//Predictor Variable    
  array[n] int<lower=0> k; //Classification Variables 
}
parameters {
  vector[J] a;//Intercept
  vector[J] b;//Regression Coefficient
  real<lower=0> sigma;//Error standard deviation
  real a0;//Intercept mean
  real b0;//Regression Coefficient mean 
  real<lower=0> s_a;//Intercept SD
  real<lower=0> s_b;//Regression Coefficient SD
}
transformed parameters {
  vector[n] yhat;//Predicted value
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b[k[i]];
  }//Equation(22.1)
}
model {
  y ~ normal(yhat, sigma);//Equation(22.2)
  a ~ normal(a0, s_a);//Equation(22.3)
  b ~ normal(b0, s_b);//Equation(22.4)
}
';

par<-c("sigma","a0","b0","s_a","s_b","a","b")                      #parameter
dataSet81 <-list(n=length(Epl11$Visual.Estimation), J=max(Epl11$Student), 
 y=Epl11$Visual.Estimation, x=Epl11$Actual.Measurement, k=Epl11$Student) #input

########### Execution by cmdstanr
modfile81 <- write_stan_file(HLM01)            #Temporary file to be exported
mod81 <- cmdstan_model(modfile81)              #compile
csrfit81 <- mod81$sample(data = dataSet81,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext81<-as_draws_df(csrfit81$draws(par))
(colnames(ext81) <- gsub("\\[|\\]|,", "", colnames(ext81)))

#save(ext81,  file="./script/obje/out2201")
#load(file="./script/obje/out2201");#Load object created with stan()in advance

######Table22.2 Summary of Posterior Distribution of Parameters for Models with
######     Distributions for Intercepts and Regression Coefficients
gqcal(ext81$sigma)
gqcal(ext81$a0   )
gqcal(ext81$b0   )
gqcal(ext81$s_a  )
gqcal(ext81$s_b  )

##### EAP estimates of 25 intercept and regression coefficients 
##### (for plotting purposes and not in the textbook)
av<-apply(ext81[, 6:30],2,mean)
bv<-apply(ext81[,31:55],2,mean)
round(rbind(av,bv),2)

###### Figure 22.3 Shows the 25 regression lines for each student
plot(Epl11$Actual.Measurement,Epl11$Visual.Estimation,pch=Epl11$Student,
    xlab="Actual Measurement",ylab="Visual Estimation",type="n")
    for (i in 1:25) {abline(av[i],bv[i],lwd=1.8)}
#dev.copy2pdf(file="fig2203.pdf")

##### Figure 22.4 displays the posterior distribution of the standard deviation of the regression coefficients
hist(ext81$s_b,breaks=100,main="",ylab="",xlab="Standard deviation of regression coefficients")
#dev.copy2pdf(file="fig2204.pdf")

##### Figure 22.5 Scatter Plot of Intercepts and Regression Coefficients
plot(av,bv,pch=16,xlab="Intercept",ylab="Regression Coefficient")
#dev.copy2pdf(file="fig2205.pdf")

cor(av,bv);     #Correlation coefficient between intercept and slope


###### 22. 2 Model assuming distribution on the intercept
###### 22. 2. 1 High School-Specific Junior High Grades and High School Deviation Values

##### Table 22.3 High School-Specific Junior High Grades and High School Deviation Values
(Chper<-read.csv("dat/Chperformance.csv", header = TRUE))

##### Figure 22.6 Deviation Values in Junior and Senior High School for Schools 2 and 9
X<-Chper[((Chper$School==2)|(Chper$School==9)), ]
plot(X$Junior.High.School,X$High.School,pch=X$School,
   xlab="Junior High School 3rd Year Deviation Value",
   ylab="High School First-Year Deviation Value", xlim=c(25,75),ylim=c(25,75))
legend(30,70,c("High School 2","High School 9"),pch=c(2,9))
#dev.copy2pdf(file="fig2206.pdf")

###Figure 22.7:Deviation Values in Junior and Senior High School for 10 Schools
plot(Chper$Junior.High.School,Chper$High.School,pch=Chper$School,
   xlab="Junior High 3rd Year Deviation Value",
   ylab="High School First-Year Deviation Value", xlim=c(25,75),ylim=c(25,75))
#dev.copy2pdf(file="fig2207.pdf")

###### First, a model assuming a distribution for the intercept and regression coefficients
par<-c("sigma","a0","b0","s_a","s_b","a","b")                      #parameter
dataSet821 <-list(n=length(Chper$High.School),J=max(Chper$School),
      y=Chper$High.School,x=Chper$Junior.High.School, k=Chper$School) #input

########### Execution by cmdstanr
#Use 'mod81' for compilation results
csrfit821 <- mod81$sample(data = dataSet821,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext821<-as_draws_df(csrfit821$draws(par))
(colnames(ext821) <- gsub("\\[|\\]|,", "", colnames(ext821)))

#save(ext821,  file="./script/obje/out2202")
#load(file="./script/obje/out2202");#Load object created with stan()in advance

###### Table9.4 Summary of Posterior Distribution of Parameters 
###### for Models with Distributions for Intercepts and Regression Coefficients
######        (Deviation Value Data)
gqcal(ext821$sigma)
gqcal(ext821$a0   )
gqcal(ext821$b0   )
gqcal(ext821$s_a  )
gqcal(ext821$s_b  )

##### EAP estimates of intercept and regression coefficients 
##### (for the sake of illustration and not in the textbook)
av<-apply(ext821[, 6:15],2,mean)
bv<-apply(ext821[,16:25],2,mean)
round(rbind(av,bv),2)

##### Figure 22.8: Ten Single Regression Lines by High School
plot(Chper$Junior.High.School,Chper$High.School,pch=Chper$School,
     xlab="Junior High Third-Year Deviation Value",
     ylab="High School First-Year Deviation Value", 
     xlim=c(25,75),ylim=c(25,75),type="n")
     for (i in 1:10) {abline(av[i],bv[i],lwd=1.8)}
#dev.copy2pdf(file="fig2208.pdf")

##### Figure 22.9 Posterior distribution of the standard deviation of the regression coefficient
hist(ext821$s_b,breaks=100,main="",ylab="",
     xlab="Standard deviation of regression coefficients")
#dev.copy2pdf(file="fig2209.pdf")

##############stan code for a model assuming a distribution on the intercept
HLM02<-'
data {
  int<lower=0> n; //Number of All Data 
  int<lower=0> J; //Number of Groups
  vector[n] y; //Criterion Variable   
  vector[n] x; //Predictor Variable   
  array[n] int<lower=0> k; //Classification Variables 
}
parameters {
  vector[J] a; //Intercept
  real b; //Regression Coefficient 
  real<lower=0> sigma; //Error standard deviation
  real a0; //Intercept mean 
  real<lower=0> s_a; //Intercept SD
}
transformed parameters {
  vector[n] yhat; //Predicted value
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b;
  } //Equation(22.5),(22.8)
}
model {
  y ~ normal(yhat, sigma); //Equation(22.6)
  a ~ normal(a0, s_a); //Equation(22.7)
}
';

par<-c("sigma","a0","b","s_a","a")                               #parameter
dataSet822 <-list(n=length(Chper$High.School),J=max(Chper$School),
   y=Chper$High.School,x=Chper$Junior.High.School, k=Chper$School) #input

########### Execution by cmdstanr
modfile822 <- write_stan_file(HLM02)           #Temporary file to be exported
mod822 <- cmdstan_model(modfile822)            #compile
csrfit822 <- mod822$sample(data = dataSet822,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext822<-as_draws_df(csrfit822$draws(par))
(colnames(ext822) <- gsub("\\[|\\]|,", "", colnames(ext822)))

#save(ext822,  file="./script/obje/out2203")
#load(file="./script/obje/out2203");#Load object created with stan() in advance

###### Table 22.5 Summary of the posterior distribution of the population of
######  the model with distribution on the intercept
gqcal(ext822$sigma )
gqcal(ext822$a0    )
gqcal(ext822$b     )
gqcal(ext822$s_a   )
print('-------------------------------------')
gqcal(ext822$a1 )
gqcal(ext822$a2 )
gqcal(ext822$a3 )
gqcal(ext822$a4 )
gqcal(ext822$a5 )
gqcal(ext822$a6 )
gqcal(ext822$a7 )
gqcal(ext822$a8 )
gqcal(ext822$a9 )
gqcal(ext822$a10)

###### Table 22.6 Probability that the intercept of high school i is greater than that of high school j
PHC02(0,ext822[,5:14],digits=2)

###### Figure 22.10 Probability that the difference in intercepts between high school 8 and 7 is greater than c
PHC01(seq(2,8,0.2),ext822$a8,ext822$a7,
   xlab="Difference between intercepts for high schools 8 and 7")
#dev.copy2pdf(file="fig2210.pdf")

###### 22.3 Model with Correlated Slope and Intercept
###### 22.3.1 Weight Change Over Four Weeks Since Diet Start

###### Table 22.7 Weight change over four weeks since the start of dieting.
(wlprocess<-read.csv("dat/wlprocess.csv", header = TRUE))

###### Figure 22.11 Four-week weight change (for four individuals)
Y<-wlprocess[wlprocess$female==14,]
plot(Y$Week,Y$weight,pch=Y$female,xlab="Weeks since weight loss started",
  ylab="weight",ylim=c(46,54),type="b")
  for (i in 15:17){
    Y<-wlprocess[wlprocess$female==i,]
    par(new=T)
    plot(Y$Week,Y$weight,pch=Y$female,xlab="",ylab="",
      ylim=c(46,54),type="b")
  }
  legend(2.8,53.5,c("female14","female15","female16","female17"),pch=c(14:17),text.width=1)
#dev.copy2pdf(file="fig2211.pdf")

###### Figure 22.12 Four-week weight change (for all participants)
Y<-wlprocess[wlprocess$female==1,]
plot(Y$Week,Y$weight,xlab="Weeks since weight loss started",
  ylab="weight",ylim=c(46,54),type="l")
  for (i in 2:50){
    Y<-wlprocess[wlprocess$female==i,]
    par(new=T)
    plot(Y$Week,Y$weight,xlab="",ylab="",ylim=c(46,54),type="l")
}

par<-c("sigma","a0","b0","s_a","s_b","a","b")                     #parameter
dataSet91 <-list(n=length(wlprocess$weight), J=max(wlprocess$female), 
    y=wlprocess$weight, x=wlprocess$Week, k=wlprocess$female)            #input
#dev.copy2pdf(file="fig2212.pdf")

########### Execution by cmdstanr 
#Use 'mod81' for compilation results
csrfit91 <- mod81$sample(data = dataSet91,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext91<-as_draws_df(csrfit91$draws(par))
(colnames(ext91) <- gsub("\\[|\\]|,", "", colnames(ext91)))

#save(ext91,  file="./script/obje/out2204")
#load(file="./script/obje/out2204");#Load object created with stan() in advance

###### Table 22.8 Summary of the posterior distribution of the model with distributed slope 
######            and intercept for independent weight change
gqcal(ext91$sigma)
gqcal(ext91$a0   )
gqcal(ext91$b0   )
gqcal(ext91$s_a  )
gqcal(ext91$s_b  )

##### EAP estimates of intercept and regression coefficients
av<-apply(ext91[,c(6:55)],2,mean)
bv<-apply(ext91[,c(56:105)],2,mean)
round(rbind(av,bv),2)

###### Figure 22.13 Scatter plot of intercepts and regression coefficients
plot(av,bv,pch=16,xlab="Intercept",ylab="Regression Coefficient");abline(h=0.0)
#dev.copy2pdf(file="fig2213.pdf")

cor(av,bv);      #Correlation coefficient between intercept and slope


################ stan code for models with correlated intercept and slope
HLM03<-'
data {
  int<lower=0> n; //Number of All Data 
  int<lower=0> J; //Number of Groups
  vector[n] y; //Criterion Variable   
  vector[n] x; //Predictor Variable   
  array[n] int<lower=0> k; //Classification Variables 
}
parameters {
  matrix[J, 2] ab; //Intercept and Regression Coefficient 
  real<lower=0> sigma; //Error standard deviation 
  real a0; //Intercept mean 
  real b0; //Regression Coefficient mean 
  real<lower=0> s_a; //Intercept SD
  real<lower=0> s_b; //Regression Coefficient SD
  real<lower=-1, upper=1> rho; //Correlation between intercept and regression coefficient
}
transformed parameters {
  vector[n] yhat; //Predicted value
  vector[2] mu; //Intercept/Regression Coefficient mean
  cov_matrix[2] Sigma; //Variance-covariance matrix of intercept and regression coefficients
  mu[1] = a0; //Substitute the mean of the intercept into the first element of μ
  mu[2] = b0; //Substitute the mean of the regression coefficients into the second element of μ
  Sigma[1, 1] = pow(s_a, 2); //Variance of intercept 
  Sigma[2, 2] = pow(s_b, 2); //Variance of regression coefficients 
  Sigma[1, 2] = s_a * s_b * rho; //Covariance of intercept and regression coefficients
  Sigma[2, 1] = Sigma[1, 2]; //same as above
  for (i in 1 : n) {
    yhat[i] = ab[k[i], 1] + x[i] * ab[k[i], 2];
  } //Equation(22.10) Predicted Value Structure
}
model {
  //y follows a normal distribution around yhat
  y ~ normal(yhat, sigma); //Equation(22.11)　No for statement in vector notation
  for (j in 1 : J) {
    //Intercept and regression coefficients follow a bivariate normal distribution
    ab[j,  : ] ~ multi_normal(mu, Sigma);
  } //Equation(22.12)
}
';

par<-c("sigma","a0","b0","s_a","s_b","rho","ab")                 #parameter
dataSet92 <-list(n=length(wlprocess$weight), J=max(wlprocess$female), 
    y=wlprocess$weight, x=wlprocess$Week, k=wlprocess$female )          #input

########### Execution by cmdstanr
modfile92 <- write_stan_file(HLM03)         #Temporary file to be exported
mod92 <- cmdstan_model(modfile92)           #compile
csrfit92 <- mod92$sample(data = dataSet92,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext92<-as_draws_df(csrfit92$draws(par))
(colnames(ext92) <- gsub("\\[|\\]|,", "", colnames(ext92)))

#save(ext92,  file="./script/obje/out2205")
#load(file="./script/obje/out2205");#Load object created with stan( )in advance

###### Table 22.9 Summary of the posterior distribution of the model 
######            with correlated slope and intercept
gqcal(ext92$sigma)
gqcal(ext92$a0   )
gqcal(ext92$b0   )
gqcal(ext92$s_a  )
gqcal(ext92$s_b  )
gqcal(ext92$rho  )

##### EAP estimates of intercept and regression coefficients
av<-apply(ext92[, 7: 56],2,mean)
bv<-apply(ext92[,57:106],2,mean)
round(rbind(av,bv),2)

#### Figure 22.14 50 individual regression lines
plot(wlprocess$Week,wlprocess$weight,pch=wlprocess$female,
    xlab="Week",ylab="weight",type="n")
    for (i in 1:50) {abline(av[i],bv[i],lwd=1.8)}
#dev.copy2pdf(file="fig2214.pdf")

##### Figure 22.15 Posterior distribution of the correlation 
#####              between the intercept and regression coefficient
hist(ext92$rho,breaks=100,main="",ylab="",xlab="Correlation between intercept and regression coefficient")
#dev.copy2pdf(file="fig2215.pdf")

cor(av,bv);                  #Correlation coefficient between intercept and slope (not in textbook)


####### 22 4 Model with Level 2 Variables
####### 22. 4. 1 Years of service and annual salary analysis

##############Table 22.10 Company-wise years of service and salary
(yoservice<-read.csv("dat/yoservice.csv", header = TRUE))

############## Table 22.11 Capital scale of companies
(capital  <-read.csv("dat/capital.csv",   header = TRUE))

############### Scatterplot of years of service and annual income for the three firms
############### (not in textbook)
X<-yoservice[((yoservice$company==3)|(yoservice$company==15)|
             (yoservice$company==21)), ]
plot(X$years.of.service,X$salary,pch=X$company,xlab="Years with company",ylab="salary")
    legend(2.5,900,c("company3","company15","company21"),pch=c(3,15,21))

############## Figure 22.16 Scatter plot of years of service and salary
plot(yoservice$years.of.service,yoservice$salary,pch=yoservice$company,
    xlab="Years with company",ylab="Annual Salary");               #Figure 22.16
#dev.copy2pdf(file="fig2216.pdf")

############################ stan code for a model with level 2 variables
HLM04<-'
data {
  int<lower=0> n; //Number of All Data
  int<lower=0> J; //Number of Groups
  vector[n] y; //Criterion Variable   
  vector[n] x; //Predictor Variable   
  vector[J] z; //Level 2 Variables   
  array[n] int<lower=0> k; //Classification Variables 
}
parameters {
  vector[J] a; //Intercept
  vector[J] b; //Regression Coefficient 
  real<lower=0> sigmae; //Error standard deviation
  real a0; //Intercept of intercept 
  real b0; //Intercept of regression coefficient 
  real ca; //intercept coefficient 
  real cb; //Regression coefficient coefficients 
  real<lower=0> sigmaa; //intercept SD
  real<lower=0> sigmab; //regression coefficient SD
}
transformed parameters {
  vector[n] yhat; //Predicted value
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b[k[i]];
  } //Equation (22.14)
}
model {
  y ~ normal(yhat, sigmae); //Equation (22.15)
  a ~ normal(a0 + z * ca, sigmaa); //Equation (22.16)
  b ~ normal(b0 + z * cb, sigmab); //Equation (22.17)
}
';

par<-c("sigmae","a0","b0","ca","cb","sigmaa","sigmab","a","b")#parameter
dataSet93 <-list(n=length(yoservice$salary),J=max(yoservice$company), 
 y=yoservice$salary, x=yoservice$years.of.service, z=capital$classification, k=yoservice$company )#input

########### Execution by cmdstanr
modfile93 <- write_stan_file(HLM04)             #Temporary file to be exported
mod93 <- cmdstan_model(modfile93)               #compile
csrfit93 <- mod93$sample(data = dataSet93,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext93<-as_draws_df(csrfit93$draws(par))
(colnames(ext93) <- gsub("\\[|\\]|,", "", colnames(ext93)))

#save(ext93,  file="./script/obje/out2206")
#load(file="./script/obje/out2206");#Load object created with stan( )in advance

###### Table 22.12 Posterior distribution of parameters for a model with level 2 variables
gqcal(ext93$sigmae)
gqcal(ext93$a0)
gqcal(ext93$b0)
gqcal(ext93$ca)
gqcal(ext93$cb)
gqcal(ext93$sigmaa)
gqcal(ext93$sigmab)

#### EAP estimates of intercept and regression coefficients
av<-apply(ext93[, 8:32],2,mean)
bv<-apply(ext93[,33:57],2,mean)
round(rbind(av,bv),2)

###### Figure 22.17 Regression lines for 25 companies at level 1
plot(yoservice$years.of.service,yoservice$salary,xlab="Years with company",ylab="Annual Salary",type="n")
  for (i in 1:25) {abline(av[i],bv[i],lwd=1.8,lty=(capital[i,2]+3))}
  legend(2.5,950,paste("classification",-2:2,sep=""),lty=1:5,text.width=10)
#dev.copy2pdf(file="fig2217.pdf")
  
######Figure 22.18 Regression line for intercept based on level 2 variable
plot(capital[,2],av,ylab="Level 1 intercept",xlab="Capital Classification",pch=16)
  abline(a=mean(ext93$a0),b=mean(ext93$ca))
#dev.copy2pdf(file="fig2218.pdf")
  
###### Figure 22.19 Regression line for regression coefficient based on level 2 variable
plot(capital[,2],bv,ylab="Level 1 regression coefficients",xlab="Capital Classification",pch=16)
  abline(a=mean(ext93$b0),b=mean(ext93$cb))
#dev.copy2pdf(file="fig2219.pdf")

cor(av,bv);                  #Correlation coefficient between intercept and slope (not in textbook)


###### 22. 5 Model with Common Slope and Level 2 Qualitative Variables

######Table 22.13 Level 2 variable indicating whether special English education is implemented
(Chperformance<-read.csv("dat/Chperformance.csv", header = TRUE))
SpecialEnglish<-c(0,0,1,1,0,0,1,1,1,1)

################stan code for a model with a common slope and level 2 qualitative variables
HLM05<-'
data {
  int<lower=0> n; //Number of All Data 
  int<lower=0> J; //Number of Groups
  vector[n] y; //Criterion Variable   
  vector[n] x; //Predictor Variable   
  vector[J] z; //Level 2 Variables   
  array[n] int<lower=0> k; //Classification Variables 
}
parameters {
  vector[J] a; //Intercept
  real b; //Regression Coefficient 
  real<lower=0> sigmae; //Error standard deviation
  real a0; //mean of intercept 
  real ca; //intercept coefficient 
  real<lower=0> sigmaa; //intercept SD
}
transformed parameters {
  vector[n] yhat; //Predicted value
  for (i in 1 : n) {
    yhat[i] = a[k[i]] + x[i] * b;
  } //Equation (22.18)
}
model {
  y ~ normal(yhat, sigmae); //Equation (22.19)
  a ~ normal(a0 + z * ca, sigmaa); //Equation (22.20),(22.21)
}
';

par<-c("sigmae","a0","ca","sigmaa","b","a")                     #parameter
dataSet94 <-list(n=length(Chperformance$High.School),J=max(Chperformance$School),
   y=Chperformance$High.School,x=Chperformance$Junior.High.School, z=SpecialEnglish, k=Chperformance$School)#input

########### Execution by cmdstanr
modfile94 <- write_stan_file(HLM05)             #Temporary file to be exported
mod94 <- cmdstan_model(modfile94)               #compile
csrfit94 <- mod94$sample(data = dataSet94,chains = 5,iter_sampling = 20000,
    iter_warmup = 1000,parallel_chains = 5,seed=1234)  #MCMC
ext94<-as_draws_df(csrfit94$draws(par))
(colnames(ext94) <- gsub("\\[|\\]|,", "", colnames(ext94)))

#save(ext94,  file="./script/obje/out2207")
#load(file="./script/obje/out2207");#Load object created with stan( )in advance

### Table 22.14 Summary of the posterior distribution of the model 
###             with a common slope and level 2 qualitative variable
gqcal(ext94$sigmae)
gqcal(ext94$a0)
gqcal(ext94$ca)
gqcal(ext94$sigmaa)
gqcal(ext94$b)

#### EAP estimates of intercept and regression coefficients
av<-apply(ext94[,6:15],2,mean)
bv<-mean(ext94$b)
round(rbind(av,bv),2)

###### Figure 22.20 Regression lines differ based on whether special education is conducted
plot(Chperformance$Junior.High.School,Chperformance$High.School,xlab="Middle School 3rd Year Deviation Score",
  ylab="High School 1st Year Deviation Score", xlim=c(25,75),ylim=c(25,75),type="n")
  for (i in 1:10) {
    abline(av[i],bv,lwd=1.9,lty=((SpecialEnglish[i])*(-1)+2) )
  }
  legend(30,70,c("With Special Education","Without Special Education"),lty=c(1,2),text.width=20)
#dev.copy2pdf(file="fig2220.pdf")

mean(ext94$ca>0);   #Probability of effectiveness of special education

##### Figure 22.21 Probability that the average difference in intercepts 
#####             for schools conducting special education is greater than c
PHC01(seq(0,10,0.5),ext94$ca,0,xlab="intercept difference")
#dev.copy2pdf(file="fig2221.pdf")

################################# Practical Assignment

#
#     Correct answer is omitted.
#

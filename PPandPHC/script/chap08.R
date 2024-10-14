########################################################################
###The working directory should be 'PPandPHC'
# Chapter 8
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)              #Load and attach add-on package cmdstanr
library(posterior)             #Load and attach add-on package posterior

#Table 8.1 Change in mental health scores over six months
#Post-assistance:
x1<-c(73,72,56,58,71,42,78,77,75,72,56,71,69,77,84,51,62,88,56,58,84,91,71,82,
81,77,65,78,79,60,66,70,65,57,64,61,56,67,75,64,68,67,80,55,48,85,56,62,65,79)
#Pre-assistance:：
x2<-c(62,54,19,54,47,22,35,77,64,60,27,41,41,44,57,16,42,89,40,67,69,46,74,62,
60,87,32,42,73,25,42,57,31,35,33,38,43,53,55,62,67,56,76,05,31,70,66,65,34,48)
print(x1);print(x2)     #Equation (8.1),Equation (8.2)

#Table 8.2 Numerical summary of 'Mental Health' data
n<-length(x1)                           #Number of data
h1<-mean(x1);round(h1,1)                #mean
h2<-mean(x2);round(h2,1)
va1<-van(x1);round(va1,1)               #variance
va2<-van(x2);round(va2,1)
sd1<-sqrt(va1);round(sd1,2)             #standard deviation
sd2<-sqrt(va2);round(sd2,2)
(me1<-median(x1))                       #median
(me2<-median(x2))
round(quantile(x1,type =2),1)           #percentile
round(quantile(x2,type =2),1)

#Figure 8.1 Boxplot of 'Mental Health' data
x<-c(x1,x2);y<-c(rep("post-assistance group",n),rep("pre-assistance group",n))
boxplot(x~y,cex.axis=2.0,xlab="",ylab="")                 
#dev.copy2pdf(file="fig0801.pdf")

#Figure 8.2 Scatterplot of 'Mental Health' data
plot(x2,x1,cex.axis=2.0,cex.lab=2.0,main="",xlab="pre-assistance",ylab="post-assistance",lwd=2)
abline(0,1,lwd=1.5)         
#dev.copy2pdf(file="fig0802.pdf")

#Table 8.3 Average deviation data of 'Mental Health'
v1<-x1-h1;round(v1,1)     #Equation (8.3)
v2<-x2-h2;round(v2,1)     #Equation (8.4)

#Covariance Equation (8.5),(8.6)
Co<-mean(v1*v2);round(Co,1)

#Figure 8.3 Scatter plot of average deviation data
plot(v2,v1,cex.axis=2.0,cex.lab=2.0,main="",xlab="pre-assistance",ylab="post-assistance",lwd=2)
abline(h=0,lwd=1.5);abline(v=0,lwd=1.5)         #
text(-40,-25,"+",cex=3.5)
text( 35, 20,"+",cex=3.5)
text(-40, 20,"-",cex=3.5)
text( 35,-25,"-",cex=3.5)
#dev.copy2pdf(file="fig0803.pdf")

#Table 8.4 Standardized data of 'Mental Health'
z1<-v1/sd1;round(z1,2)     #(Equation 8.7)
z2<-v2/sd2;round(z2,2)     #Equation (8.8)

#Correlation coefficient Equations (8.9),(8.10)
r<-mean(z1*z2);round(r,2)

#8.2     Bivariate normal distribution
#Difference inference between paired two groups
x<-cbind(x1,x2)
out0801<-G2pair(x,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0801,  file="./script/obje/out0801")
#load(file="./script/obje/out0801"); #Load the results of a pre-run of MCMC

#Table 8.5 Numerical summary of the posterior distribution of parameters 
#          and predictive distribution
gqcal(out0801$ext[,1:7])

#Extraction of random numbers
mu1<-out0801$mu1;mu2<-out0801$mu2;
sigma1<-out0801$sigma1;sigma2<-out0801$sigma2;
xaste1<-out0801$xaste1;xaste2<-out0801$xaste2;
rho<-out0801$rho

# Composition of generated quantities
dmu     <-mu1-mu2;                          # mean difference 
dsigma  <-sigma1-sigma2;                    # Difference in standard deviations
sigma_within <-sqrt((sigma1^2 + sigma2^2)/2)#Equation (7.1), within group standard deviation
delta   <-dmu/sigma_within;                 #Equation (7.2), standardized mean difference
delta1  <-dmu/sigma1;                       #Equation (7.3)
delta2  <-dmu/sigma2;                       #Equation (7.4)

#Equation (7.6), measure of nonoverlap , denominator is within-group standard deviation
U1      <-pnorm(mu1,mu2,sigma_within);
#Equation (7.7), measure of nonoverlap , denominator is within-group standard deviation
U2      <-pnorm(mu2,mu1,sigma_within);
#Equation (7.9), measure of nonoverlap , denominator is second group standard deviation
U1a     <-pnorm(mu1,mu2,sigma2);
#Equation (7.10), measure of nonoverlap , denominator is first group standard deviation
U2a     <-pnorm(mu2,mu1,sigma1);

#Table 8.6 Numerical summary of the posterior distribution of generated quantities
gqcal(cbind(dmu,dsigma,sigma_within,delta,delta1,delta2,U1,U2,U1a,U2a))

#Figure 8.7 PHC curve
par(mfrow=c(4,1))
PHC01(seq(10,25,0.5),     dmu,   0,cc="gtc",byoga="yes",xlab="mean difference dμ");#
PHC01(seq(0.5,1.3,0.05),  delta2,0,cc="gtc",byoga="yes",
      xlab="standardized mean difference δ2");#
PHC01(seq(0.65,0.95,0.02),U1a   ,0,cc="gtc",byoga="yes",xlab="measure of nonoverlap  U1*");#
ikijyori2(seq(0,30,5),    out0801, byoga="yes",xlab="probability beyond threshold πc");#
par(mfrow=c(1,1))
#dev.copy2pdf(file="fig0807.pdf")

#Table 8.7 PHC table
PHC01(seq(5,30,5),     dmu,   0,cc="gtc",byoga="no");  #Difference in means
PHC01(seq(0.6,1.1,0.1),  delta2,0,cc="gtc",byoga="no");#standardized mean difference 
PHC01(seq(0.7,0.95,0.05),U1a   ,0,cc="gtc",byoga="no");#measure of nonoverlap
ikijyori2(seq(5,30,5),    out0801,byoga="no");         #probability beyond threshold

########################################################################
#Chapter 8 Section 8.5 Practical assignment

#Table 8.8 Raw data of estimated and actual lengths (mm)
y1<- c(110,232,176,207,122,202,191,124,193,250); #Estimation group
y2<- c(130,268,104,185,128,147,162, 68,142,175); #Actual measurement group
y<-cbind(y1,y2)

out0802<-G2pair(y,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0802,  file="./script/obje/out0802")
#load(file="./script/obje/out0802"); #Load the results of a pre-run of MCMC

#Subsequent script examples are omitted


########################################################################
#Self-study
#To investigate the effect of a certain diet method, 
#20 women were recruited to participate. The body weight before 
#and after program participation is as follows. 
#Is this diet method effective?
#x1 is before, and x2 is after. Assuming equal variances, 
#respond to the following.
#1) Consider the posterior distribution and PHC for the difference in means.
#2) Consider the posterior distribution and PHC for the standardized mean difference.
#3) Consider the posterior distribution and PHC for measure of nonoverlap.
#4) Interpret the probability beyond threshold.

w1<-c(53.1,51.5,45.5,55.5,49.6,50.1,59.2,54.7,53.0,48.6,
      55.3,52.6,51.7,48.6,56.4,42.9,50.3,42.4,51.2,39.1)
w2<-c(48.3,45.2,46.6,56.6,41.2,44.6,51.9,55.5,45.4,47.6,
      50.6,54.5,49.0,43.9,53.8,40.1,52.8,35.3,55.6,38.0)

w<-cbind(w1,w2)
out0803<-G2pair(w,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0803,  file="./script/obje/out0803")
#load(file="./script/obje/out0803"); #Load the results of a pre-run of MCMC

#Table 8.5 Numerical summary of the posterior distribution of parameters and predictive distribution
print(out0803$fit)

#Extraction of random numbers
mu1<-out0803$mu1;mu2<-out0803$mu2;
sigma1<-out0803$sigma1;sigma2<-out0803$sigma2;
xaste1<-out0803$xaste1;xaste2<-out0803$xaste2;
rho<-out0803$rho

#Posterior distribution of generated quantities
dmu     <-mu1-mu2;                          #Difference in means
sigma_within <-sqrt((sigma1^2 + sigma2^2)/2)#Equation (7.1), within group standard deviation
delta   <-dmu/sigma_within;                 #Equation (7.2), standardized mean difference 
U1      <-pnorm(mu1,mu2,sigma_within);  #Equation (7.6), measure of nonoverlap , denominator is within group standard deviation 
# 1)
gqcal(dmu   )
PHC01(seq(0.5,5.0,0.5),dmu, 0,cc="gtc",byoga="yes"); #mean difference
PHC01(seq(0.5,5.0,0.5),dmu, 0,cc="gtc",byoga="no");  #mean difference

# 2)
gqcal(delta )
PHC01(seq(0.0,1.0,0.1),delta,0,cc="gtc",byoga="yes");#standardized mean difference 
PHC01(seq(0.0,1.0,0.2),delta,0,cc="gtc",byoga="no"); #standardized mean difference 

# 3)
gqcal(U1    )
PHC01(seq(0.4,0.9,0.02),U1,0,cc="gtc",byoga="yes");#measure of nonoverlap
PHC01(seq(0.4,0.9,0.05),U1,0,cc="gtc",byoga="no" );#measure of nonoverlap

# 4)
ikijyori2(seq(0,3,0.2),    out0803,byoga="yes");    #probability beyond threshold curve
ikijyori2(seq(0,3,0.5),    out0803,byoga="no")      #probability beyond threshold table

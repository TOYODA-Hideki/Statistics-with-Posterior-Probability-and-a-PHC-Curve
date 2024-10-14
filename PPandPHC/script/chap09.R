########################################################################
###The working directory should be 'PPandPHC'
# Chapter 9
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

#Difference inference between paired two groups
x<-cbind(x1,x2)
out0801<-G2pair(x,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0801,  file="./script/obje/out0801")
#load(file="./script/obje/out0801"); #Load the results of a pre-run of MCMC

#Table 8.5 Numerical summary of the posterior distribution of parameters
# and predictive distribution
#(Review of previous chapter)
gqcal(out0801$ext[,1:7])

#Extraction of random numbers
mu1<-out0801$mu1;mu2<-out0801$mu2;                #平均値
sigma1<-out0801$sigma1;sigma2<-out0801$sigma2;    #標準偏差
xaste1<-out0801$xaste1;xaste2<-out0801$xaste2;    #予測分布
rho<-out0801$rho;                                 #相関係数
dmu     <-mu1-mu2;                                #平均値の差
dsigma  <-sigma1-sigma2;                          #標準偏差の差

#The posterior predictive distribution of difference scores     (9.1)式
ds<-xaste1-xaste2
gqcal(ds)

#The posterior distribution of the mean of difference scores   (9.2)式
gqcal(dmu);

#The generated quantities of the standard deviation of difference scores (9.3)式
sigmaD  <-sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2);

#Figure 9.1 The posterior predictive distribution of difference scores 
#      and the posterior distribution of the standard deviation of difference scores
par(mfrow=c(1,2))
#Left figure
hist(ds ,breaks=100,col=gray(1.0),cex.axis=1.2,cex.lab=1.2,
      main="",xlab="The posterior predictive distribution of difference scores",ylab="") 
      arrows(18.6,5000,18.6,0,lwd=2)
      text(22.0,5150,"Average: 18.6 points",cex=1.5)
      arrows(18.6,1500,18.6-15.6,1500,lwd=4,lty=1)
      text(60,1550,"Standard Deviation (SD) ",cex=1)
      text(60,1350,"of Difference Scores: 15.7",cex=1)
      arrows(0,1000,18.6,1000,lwd=4,lty=1)
      text(-25,1000,"Effect of Assistance",cex=1)
#Right figure
hist(sigmaD, breaks=100,col=gray(1.0),cex.axis=1.1,cex.lab=1.1,
      main="",xlab="The Posterior Distribution of the Standard Deviation of Difference Scores",ylab="")
par(mfrow=c(1,1))
#dev.copy2pdf(file="fig0901.pdf")

#Table 9.1 Numerical summary and PHC table of the posterior distribution 
#          of the standard deviation of the difference score
gqcal(sigmaD,1)
PHC01(seq(15,21,1),sigmaD,cc="ltc",byoga="no")

#Figure 9.2 PHC curve of standard deviation of difference scores
PHC01(seq(10,24,0.2),sigmaD,cc="ltc",byoga="yes",xlab="the standard deviation of the difference score")
#dev.copy2pdf(file="fig0902.pdf")

#Table 9.2 the Standardized Mean of the Difference Score
deltaD  <-dmu/sigmaD;                 #(9.9)式、(9.10)式
gqcal(deltaD)
PHC01(seq(0.6,1.2,0.1),deltaD,cc="gtc",byoga="no")

#Figure 9.3 Posterior Distribution and PHC Curve of the Standardized Mean of the Difference Score
par(mfrow=c(1,2))
#Left figure
hist(deltaD, breaks=100,col=gray(1.0),cex.axis=2.0,cex.lab=2.0,main="",
      xlab="Posterior Distribution of the Standardized Mean of the Difference Score",ylab="")
#Right figure
PHC01(seq(0.5,2.0,0.05),deltaD,cc="gtc",byoga="yes",
      xlab="Standardized Mean of the Difference Score")
par(mfrow=c(1,1))
#dev.copy2pdf(file="fig0903.pdf")

#Figure 9.4: Visualization of the Probability Beyond Threshold for Difference Score
nob<-1000
plot(xaste1[1:nob],xaste2[1:nob],xlim=c(0,100),ylim=c(0,100),cex.axis=2.0,
      cex.lab=1.3,main="",xlab="Post-assistance",ylab="Pre-assistance",lwd=2.0) 
      abline( 0 ,1,lwd=2,lty=1)
      abline(-3 ,1,lwd=2,lty=2)
      abline(-5 ,1,lwd=2,lty=3)
      abline(-10,1,lwd=2,lty=4)
      arrows(25,25,35,10,lwd=4)
      text(39,6,"88.5％",cex=2.5)
#dev.copy2pdf(file="fig0904.pdf")

#9.3.4 "probability beyond threshold" appearing in the relationship 
#      between probability beyond threshold and probability beyond threshold of difference scores.
ikijyori2(c(0,3,5,10),out0801,byoga="no")

#Table 9.3: Estimation Results of the Posterior Distribution of the Difference Score’s Probability Beyond Threshold 
Dikijyori2(c(0,3,5,10),out0801,byoga="no")

#Figure 9.5: EAP and the upper and lower limits of the 95% Credible Interval of the Probability Beyond Threshold
Dikijyori2(seq(-10,30,2),out0801,byoga="yes")
#dev.copy2pdf(file="fig0905.pdf")

########################################################################
#Chapter 9 9.7 Practical Problems

#Table 8.8 Raw data of estimated and actual lengths (mm)
y1<- c(110,232,176,207,122,202,191,124,193,250); #Estimation group
y2<- c(130,268,104,185,128,147,162, 68,142,175); #Actual measurement group

y<-cbind(y1,y2)
out0802<-G2pair(y,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0802,  file="./script/obje/out0802")
#load(file="./script/obje/out0802"); #Load the results of a pre-run of MCMC

#Subsequent script examples are omitted.


########################################################################
#Appendix: The probability that the judgment "the likelihood of scoring 10 points higher 
#          after assistance being greater than 60%" is correct.
Dikijyo10<-pnorm((mu1-mu2-10)/sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2),0,1)
PHC01(0.6,Dikijyo10,0,cc="gtc",byoga="no")



########################################################################
# Self-Study
# Investigating the Effectiveness of a Diet Method
# To evaluate the effectiveness of a particular diet method, 
# 20 women participated in the program.
# The weights of the participants before and after the program are recorded. 
# The weights before are denoted as x1, and after as x2.
# Assuming equal variances, answer the following:
#1) Compare and discuss the standardized mean difference, 
#   the posterior distribution of the difference score's mean, and PHC.
#2) Compare and discuss probability beyond threshold and  
#   probability beyond threshold of difference scores.
#3) Conduct a test for the mean difference of the two related groups.

w1<-c(53.1,51.5,45.5,55.5,49.6,50.1,59.2,54.7,53.0,48.6,
      55.3,52.6,51.7,48.6,56.4,42.9,50.3,42.4,51.2,39.1)
w2<-c(48.3,45.2,46.6,56.6,41.2,44.6,51.9,55.5,45.4,47.6,
      50.6,54.5,49.0,43.9,53.8,40.1,52.8,35.3,55.6,38.0)

w<-cbind(w1,w2)
out0803<-G2pair(w,EQU=0,prior=T,mL=0, mH=250, sL=0, sH=125)
#save(out0803,  file="./script/obje/out0803")
#load(file="./script/obje/out0803"); #Load the results of a pre-run of MCMC

gqcal(out0803$ext[,1:7])

# Extracting random numbers and creating generation quantities
dmu     <-out0803$mu1-out0803$mu2;             # mean difference
sigma1<-out0803$sigma1;sigma2<-out0803$sigma2; # Difference in standard deviations
rho<-out0803$rho;                              # correlation coefficient 
sigma_within <-sqrt((sigma1^2 + sigma2^2)/2)   # Equation (7.1), within group standard deviation
sigmaD  <-sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2);# Equation (9.3),Standard deviation of difference scores  
delta   <-dmu/sigma_within;                    #Equation (7.2), standardized mean difference
deltaD  <-dmu/sigmaD;                          #Equation (9.9),(9.10)

#1) Posterior distribution of standardized mean differences, mean of difference scores and PHC
gqcal(delta)
PHC01(seq(0.0,1.2,0.1),delta, cc="gtc",byoga="no")
gqcal(deltaD)
PHC01(seq(0.0,1.2,0.1),deltaD,cc="gtc",byoga="no")

#2) probability beyond threshold and difference scores probability beyond threshold 
ikijyori2(seq(0,3,0.5),out0803,byoga="no",dedits=2)
Dikijyori2(seq(0,3,0.5),out0803,byoga="no",dedits=2)

#3) Test for difference of means of two corresponding groups
t.test(w1,w2,paired=T)  #t-test


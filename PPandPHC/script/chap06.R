########################################################################
###The working directory should be 'PPandPHC'
#Chapter 6
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
   library(cmdstanr)           #load and attach add-on package cmdstanr
   library(posterior)          #load and attach add-on package posterior

#Table 6.1 English performance of experimental group and control group
#Experimental group
x1<-c(49,66,69,55,54,72,51,76,40,62,66,51,59,68,66,57,53,66,58,57)
#Control group
x2<-c(41,55,21,49,53,50,52,67,54,69,57,48,31,52,56,50,46,38,62,59)
print(x1);print(x2)

#Table 6.2 Summary statistic for 'English performance data'
(n1<-length(x1));(n2<-length(x2))         #The number of data
mean(x1);mean(x2)                         #Mean
van<-function(x){mean((x-mean(x))^2)}     #Function to calculate variance
van(x1);van(x2)                           #Variance
sqrt(van(x1));sqrt(van(x2))               #Standard deviation
sort(x1);sort(x2)                         #Sort
median(x1);median(x2)                     #Median
quantile(x1,type =2);quantile(x2,type =2) #Percentile

#Table 6.3 Sorted measurements
sort(x1);sort(x2)

#Fig. 6.1 Box plot of 'English performance data'
x<-c(x1,x2);y<-c(rep("Experimental",n1),rep("Control",n2))
boxplot(x~y,cex.axis=2.0,xlab="",ylab="", border=2, col="lightblue",lwd=1.8)
#dev.copy2pdf(file="../English_tex/fig/chap06/fig6_1.pdf")

#G2Ind, Inferring the difference between two independent groups
# x1: vector. Group 1 data (Group with large mean value)
# x2: vector. Group 2 data (Group with small mean value)
# EQU: logical. If 1, two variances are treated as equal.
#               If 0, two variances are treated separately
# prior: logical. If it is TRUE, specify the range of the prior distribution.
#                 If it is FALSE, stan default is used.
# mL, mH, sL, sH: Prior distribution parameters
#    Upper and lower bounds of the mean (sd) value of a uniform distribution.
out0601<-G2Ind(x1,x2,EQU=0,prior=F)

#save(out0601,  file="./script/obje/out0601")
#load(file="./script/obje/out0601"); #Load the results of a pre-run of MCMC

#Table 6.4 Summary statistic of the posterior distribution of the parameters
gqcal(out0601$ext[,1:6])

# Take out the MCMC sample.
mu1   <-out0601$mu1;      mu2<-out0601$mu2;
sigma1<-out0601$sigma1;sigma2<-out0601$sigma2;
xaste1<-out0601$xaste1;xaste2<-out0601$xaste2;

#Fig. 6.2 Comparison of posterior distributions of experimental and control groups
par(mfrow=c(1,2))
#figure on the left 
plot(density(mu1),xlim=c(40,70),ylim=c(0,0.2),cex.axis=2.0,
      cex.lab=2.0,main="",xlab="",ylab="",lwd=2.0)
      text(55.0,0.17,"experimental",cex=2.0)
par(new=T)
plot(density(mu2),xlim=c(40,70),ylim=c(0,0.2),cex.axis=2.0,
      cex.lab=2.0,main="",xlab="mean",ylab="",lwd=2.0)
      text(44.5,0.15,"control",cex=2.0)
#figure on the right
plot(density(sigma1),xlim=c(5,25),ylim=c(0,0.28)
,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="",lwd=2.0)  
text( 12.5,0.25,"experimental",cex=2.0)
par(new=T)
plot(density(sigma2),,xlim=c(5,25),ylim=c(0,0.28)
,cex.axis=2.0,cex.lab=2.0,main="",xlab="standard deviation",ylab="",lwd=2.0)  
text(16.0,0.17,"control",cex=2.0)
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap06/fig6_2.pdf")


#Table 6.5 Summary statistic of the posterior distribution of the generated quantity
dmu     <-mu1-mu2;       #Difference between means
dsigma  <-sigma1-sigma2; #Difference between standard deviations
gqcal(cbind(dmu,dsigma),2)

#Fig. 6.3 Graphical summary of the posterior distribution of generated quantity
par(mfrow=c(1,2))
hist(dmu ,breaks=100,border=2, col="lightblue",lwd=1.8
,cex.axis=2.0,cex.lab=2.0,main="",xlab="Difference between means",ylab="") 
hist(dsigma ,breaks=100,border=2, col="lightblue",lwd=1.8
,cex.axis=2.0,cex.lab=2.0,main="",xlab="Difference between standard deviations",ylab="") 
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap06/fig6_3.pdf")

#Table 6.6 PHC table for the hypothesis that there is a difference
PHC01(c(0,3,5,10,15,20),dmu,cc="gtc",byoga="no")
PHC01(c(-20,-15,-10,-5,-3,0),dsigma,cc="ltc",byoga="no")

#Fig. 6.4  PHC curve for the hypothesis that there is a difference
par(mfrow=c(2,1))
PHC01(seq(0,20,1.0),dmu,cc="gtc",byoga="yes",xlab="Difference between means")
PHC01(seq(-20,5,1.0),dsigma,cc="ltc",byoga="yes",xlab="Difference between standard deviations")
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap06/fig6_4.pdf")

#Table 6.7 PHC table for ROPE's hypothesis that there is virtually no difference
PHC01(c(0,3,5,10,15,20),dmu,   cc="rope",byoga="no")
PHC01(c(0,3,5,10,15,20),dsigma,cc="rope",byoga="no")

#Fig. 6.5 PHC curve for ROPE's hypothesis that there is virtually no difference
par(mfrow=c(2,1))
PHC01(seq(0,20,1.0),dmu,   cc="rope",byoga="yes",xlab="Difference between means")
PHC01(seq(0,20,1.0),dsigma,cc="rope",byoga="yes",xlab="Difference between standard deviations")
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap06/fig6_5.pdf")

#traditional t-test
t.test(x1,x2,var.equal = F)

########################################################################
#Chapter 6 6.5 Practical assignment

#Table 6.8 Measurement results of perceived time under listening conditions(sec.)
#listening condition
y1<-c(32.30,34.24,28.10,33.40,37.71,31.62,31.37,35.85,32.33,34.04,
     34.96,31.43,35.28,30.19,35.09,33.38,31.49,28.44,32.12,31.81)

#control condition
y2<-c(31.43,31.09,33.38,30.49,29.62,35.40,32.58,28.96,29.43,28.52,
     25.39,32.68,30.51,30.15,32.33,30.43,32.50,32.07,32.35,31.57)

out0602<-G2Ind(y1,y2,EQU=0,prior=F)
#save(out0602,  file="./script/obje/out0602")
#load(file="./script/obje/out0602"); #Load the results of a pre-run of MCMC
gqcal(out0602$ext[,1:6])

#The script after this is omitted.


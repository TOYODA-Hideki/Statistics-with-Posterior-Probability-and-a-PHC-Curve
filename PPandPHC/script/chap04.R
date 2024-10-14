########################################################################
###The working directory should be 'PPandPHC'
#Chapter 4
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)              #Load and attach add-on package cmdstanr
library(posterior)             #Load and attach add-on package posterior

#Table 1.1 Measurement Results (Days)
#Antigen test data
x<-c(
11.5, 14.0, 15.0, 10.0, 14.0, 14.5, 12.5, 12.5, 12.5, 12.0, 
13.0, 12.0, 12.5, 13.5, 15.0, 13.5, 12.5, 12.0, 12.0, 17.0)

#Inference on the normal distribution
out0301 <-G1mean(x,prior=F)
#save(out0301,  file="script/obje/out0301")
#load(file="script/obje/out0301"); #Load the results of a pre-run of MCMC

#Extraction of random numbers
mu<-out0301$mu
sigma<-out0301$sigma

#Equation (4.2) Random numbers following the posterior distribution of the standard deviation and their squares
bunsan <-sigma^2
round( sigma[1:7],3)
round(bunsan[1:7],3)

#4.2.1 Variance (Equation 4.2)
max(bunsan)

#4.2.2 The square of the mean does not equal the mean of the squares (Equation 4.4)
round(mean(bunsan)  ,3)
round(mean(sigma)^2 ,3)

#4.2.3 The square of the 50% point is equal to the 50% point of the squares (Equation 4.5)
round(median(bunsan)  ,3)
round(median(sigma)^2 ,3)

#4.2.5 Coefficient of variation (Equations 4.6,4.7)
hendou <-sigma/mu

#Table 4.1 gqcal(x) is a numerical summary of generated quantities. x is a random number vector of generated quantities approximating the posterior distribution
#or a matrix of (random numbers x generated quantities). Here, it is illustrated with a matrix as it is the first appearance.
gqcal(cbind(bunsan,sigma,hendou)); #Variance,Standard Deviation and Coefficient of Variation

#4.3.1 Difference between the reference point and the mean value (Equations 4.8,4.9)
sa15<- 15-mu

#4.3.2 Difference between the standardized reference point and the mean (Equations 4.10,4.11)
delta15<-(mu-15.0)/sigma

#4.3.3 Quantile & Percentile (Equations 4.12,4.13)
quarter3<-mu+0.675*sigma;        #75% percentile

#4.3.4 Observation probability in a specific interval (Equations 4.14,4.15)
yosokbu<- 1-pnorm(15.0,mu,sigma)

#4.3.5 Ratio of reference point to mean (Equation 4.16,4.17)
hirit15mu<-mu/15

# Table 4.2 Numerical summary of generated quantities with a reference point (illustrated as a vector) 
gqcal(sa15)      #Difference between the 15-day reference point and the mean (Equation 4.8)
gqcal(delta15)   #Difference between the standardized reference point and the mean (Equation 4.10)
gqcal(quarter3)  #Third quartile                                    (Equations 4.12,4.13)
gqcal(yosokbu)   #Probability of turning negative after 15 days     (Equations 4.14,4.15)
gqcal(hirit15mu) #Ratio of mean values against the reference point of the 15th day (Equations 4.16,4.17)

#Table 4.2 A Table Compiled by Gathering Data
gqcal(cbind(sa15,delta15,quarter3,yosokbu,hirit15mu))

#Figure 4.1 Posterior distribution of variance
hist(bunsan ,breaks=200,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0401.pdf")

#Figure 4.2 Posterior distribution of the coefficient of variation
hist(hendou ,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0402.pdf")

#Figure 4.3 Posterior distribution of d15
hist(sa15 ,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0403.pdf")

#Figure 4.4 Posterior distribution of δ15
hist(delta15,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0404.pdf")

#Figure 4.5 Posterior distribution of the third quartile
hist(quarter3,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0405.pdf")

#Figure 4.6 Posterior distribution of observation probability in specific intervals
hist(yosokbu,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0406.pdf")

#Figure 4.7 Posterior distribution of the ratio of baseline to mean
hist(hirit15mu,breaks=200,xlim=c(0.8,1.0),
cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file="fig0407.pdf")

########################################################################
#Chapter 4 4.5 Practical assignment
#Table 1.3 Measurement Results (seconds) (Data entry for "Perception Time")
y<-c(31.43,31.09,33.38,30.49,29.62,35.40,32.58,28.96,29.43,28.52,
     25.39,32.68,30.51,30.15,32.33,30.43,32.50,32.07,32.35,31.57)

#Inference on the normal distribution
out0302 <-G1mean(y,prior=F)#Inference on the normal distribution
#save(out0302,  file="script/obje/out0302")
#load(file="script/obje/out0302"); #Load the results of a pre-run of MCMC

#Calculation of generated quantities
bunsan <-out0302$sigma^2
hendou <-out0302$sigma/out0302$mu
delta30<-(out0302$mu-30)/out0302$sigma
quarter<-out0302$mu-0.675*out0302$sigma
yosokbu<-pnorm(30.5,out0302$mu,out0302$sigma)-pnorm(29.5,out0302$mu,out0302$sigma)
hirit30<-out0302$xaste/30

# Estimation results of generated quantities (Summary of distribution)
gqcal(bunsan )#Variance
gqcal(hendou )#Coefficient of variation
gqcal(delta30)#Delta 30
gqcal(quarter)#First quartile
gqcal(yosokbu)#Probability of observing data 29.5 seconds < y < 30.5 seconds
gqcal(hirit30)#Ratio to 85g

#Various posterior distributions
hist(bunsan ,breaks=100)
hist(hendou ,breaks=100)
hist(delta30,breaks=100)
hist(quarter,breaks=100)
hist(yosokbu,breaks=100)
hist(hirit30,breaks=100)

#########################################################################
#Chapter 4.6 Practical Questions
#From the perspective of the following RQs,compare the properties of "New Drug A" and "New Drug B". 
# [RQ.8] Point estimate and interval estimate of the difference between the 15-day reference point and the mean value.
#[RQ.11] Point estimate and interval estimate of the probability of taking more than 15 days to become negative.

xB<-scan("dat/AntigenTestB.dat",sep="")
out0303 <-G1mean(xB,prior=F)
#save(out0303,  file="script/obje/out0303")
#load(file="script/obje/out0303"); #Load the results of a pre-run of MCMC

#Calculation of generated quantities
sa15B<-15-out0303$mu
yosokbuB<- 1-pnorm(15.0,out0303$mu,out0303$sigma)

# Estimation results of generated quantities (Summary of distribution) for New Drug B
gqcal(sa15B )    #Difference between the 15-day reference point and the mean value         Equation (4.8)             *[RQ.8]
gqcal(yosokbuB)  # Probability of becoming negative after the 15-day reference point       (Equations 4.14, 4.15) *[RQ.11]

#Estimation results of generated quantities (Summary of distribution) for New Drug A
gqcal(sa15 )     #Difference between the 15-day reference point and the mean value        Equation(4.8)
gqcal(yosokbu)   #Probability of becoming negative after the 15-day reference point      (Equations 4.14, 4.15)

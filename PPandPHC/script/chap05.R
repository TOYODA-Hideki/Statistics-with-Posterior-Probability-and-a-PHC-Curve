########################################################################
#Chapter 5
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)           #load and attach add-on package cmdstanr
library(posterior)          #load and attach add-on package posterior

# Table 1.1: Antigen Test Data - Measured Results (Days)
x <- c(
  11.5, 14.0, 15.0, 10.0, 14.0, 14.5, 12.5, 12.5, 12.5, 12.0, 
  13.0, 12.0, 12.5, 13.5, 15.0, 13.5, 12.5, 12.0, 12.0, 17.0
)

# Inference regarding the Normal Distribution
out0301 <- G1mean(x, prior = F)
# save(out0301,  file = "./script/obje/out0301")
# load(file = "./script/obje/out0301"); # Load the results of a pre-run of MCMC

# Figure 5.1: Area representing the probability attributed to the event μ <= 13.5
hist(out0301$mu[out0301$mu<13.5] ,breaks=70,col=gray(0.6),xlim=c(11,15),
     ylim=c(0,5500),cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="") 
par(new=T)
hist(out0301$mu[out0301$mu>=13.5],breaks=30,col=gray(1.0),xlim=c(11,15),
     ylim=c(0,5500),cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")  
# dev.copy2pdf(file = "../English_tex/fig/chap05/fig05_01.pdf")

# 5.1.4 Probability of research hypothesis μ < 13.5 being validated (Examples for equations 5.2, 5.3, and 5.4)
temp01 <- out0301$mu[11:20]
round(temp01, 2)
ifelse(temp01 < 13.5, 1, 0)

# The function PHC01() draws a PHC curve or creates a PHC table.
# seq01: Arithmetic sequence (or scalar) representing the reference point c (horizontal axis of the curve)
# a    : Numeric vector, parameter or generated value or random values from predictive distribution
# b    : Scalar or numeric vector, parameter or generated value or random values from predictive distribution, default is 0
# cc   : "gtc" implies a - b > c, "rope" implies abs(a-b) < c, "ltc" implies a - b < c
#        Since the default value of b is 0, if not specified, cc serves as an indicator to a.
# byoga: "yes" for drawing the PHC curve, "no" to create a PHC table
# dedits: Integer, number of digits after the decimal point in the output

# 5.2 "Probability that the Research Hypothesis is Correct" PHC (Equation 5.6)
PHC01(13.5, out0301$mu, 0, cc="ltc", byoga="no")

# 5.2.3 "Effectively the Same Range" ROPE (Equations 5.7, 5.8, 5.9)
PHC01(1.5, out0301$mu, 15.0, cc="rope", byoga="no")

# Probability that the hypothesis "Clearly Longer" is correct - p.67, last equation
PHC01(16.5, out0301$mu, 0, cc="gtc", byoga="no", dedits=5)

# 5.3.3 The reference point c doesn't necessarily have to be determined in advance
# Table 5.1 PHC(μ<=c) PHC table (New Drug A)
PHC01(seq(11.5,15.0,0.5), out0301$mu, 0, cc="ltc", byoga="no")

# Figure 5.2 PHC(μ<=c) PHC curve (New Drug A) n=20
PHC01(seq(11,16,0.2), out0301$mu, 0, cc="ltc", byoga="yes", cex.lab=2.2, cex.axis=2.0)
text(14.25, 0.8, "★", cex=4.0)
# dev.copy2pdf(file="../English_tex/fig/chap05/fig05_02.pdf", family="Japan1")

# 5.3.4 When n is small, claims favorable to the analyst are less likely to be supported
# Posterior distribution of parameters for Antigen Test B
xB <- scan("dat/AntigenTestB.dat", sep="")

out0303 <- G1mean(xB, prior=F) # Inference related to the normal distribution
#save(out0303,  file="./script/obje/out0303")
#load(file="./script/obje/out0303"); # Load the results of a pre-run of MCMC

# Figure 5.3 PHC(μ<=c) PHC curve (New Drug B) n=1000
PHC01(seq(11,16,0.02), out0303$mu, cc="ltc", byoga="yes", cex.lab=2.2, cex.axis=2.0)
arrows(14.25, 0.8, 14.25, 0.1, lwd=7.0)
arrows(15.25, 0.2, 15.25, 0.9, lwd=7.0)
#dev.copy2pdf(file="../English_tex/fig/chap05/fig05_03.pdf")

# Table 5.2 PHC(c<=c) PHC Table (New Drug B)
PHC01(seq(14.5,15.2,0.1), out0303$mu, 0, cc="ltc", byoga="no");

# 5.3.7 ROPE PHC Curve
# Table 5.3 PHC(|μ-15,|< c) PHC Table (Evaluation of ROPE)
PHC01(seq(0,4.0,0.5), out0301$mu, 15.0, cc="rope", byoga="no");

# Figure 5.4 PHC(|μ-15,|< c) PHC Curve
PHC01(seq(0,4,0.1), out0301$mu, 15, cc="rope", byoga="yes", cex.lab=2.2, cex.axis=2)
#dev.copy2pdf(file="../English_tex/fig/chap05/fig05_04.pdf")

# 5.4 PHC Curve for Generated Value
# Table 5.4 PHC Table for Parameters and Generated Value
sa15 <- 15 - out0301$mu                   # Difference between reference point of 15 days and mean value (Equation 4.8)
delta15 <- (out0301$mu - 15.0) / out0301$sigma; # Delta δ (Equation 4.10)
hirit15mu <- out0301$mu / 15            # Mean ratio to 15 days (Equations 4.16, 4.17)
PHC01(seq(1,3,0.25), out0301$sigma, 0, cc="ltc", byoga="no");
PHC01(round(seq(2/3,10/3,1/3),3), sa15, 0, cc="gtc", byoga="no");
PHC01(seq(-2.0,0.0,0.25), delta15, 0, cc="ltc", byoga="no");
PHC01(seq(0.8,1.0,0.025), hirit15mu, 0, cc="ltc", byoga="no");

# Figure 5.5 PHC Curve for Parameters and Generated Value
#dev.new()
par(mfrow=c(4,1))
PHC01(seq(1,3,0.1), out0301$sigma, 0, cc="ltc", byoga="yes", cex.lab=2.5, cex.axis=2.0, xlab="a. SD");
PHC01(seq(2/3,10/3,1/30), sa15, 0, cc="gtc", byoga="yes", cex.lab=2.5, cex.axis=2.0, xlab="b. Difference between Reference Point and Mean");
PHC01(seq(-2.0,0.0,0.05), delta15, 0, cc="ltc", byoga="yes", cex.lab=2.5, cex.axis=2.0, xlab="c. SD between Reference Point and Mean");
PHC01(seq(0.8,1.0,0.005), hirit15mu, 0, cc="ltc", byoga="yes", cex.lab=2.5, cex.axis=2.0, xlab="d. Ratio of Reference Point to Mean");
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap05/fig05_05.pdf")
#dev.off()
########################################################################
# Chapter 5, Practical assignment 5.6
# Table 1.3 Measurement Results (seconds) (Data entry for "Perception Time")
y <- c(31.43,31.09,33.38,30.49,29.62,35.40,32.58,28.96,29.43,28.52,
       25.39,32.68,30.51,30.15,32.33,30.43,32.50,32.07,32.35,31.57)

out0302 <- G1mean(y, prior=F) # Inference related to the normal distribution
#save(out0302, file="./script/obje/out0302")
#load(file="./script/obje/out0302"); # Load the results of a pre-run of MCMC

# Illustrative script is omitted

########################################################################
# Chapter 5, Practice Problem 5.7, p.75

xB <- scan("dat/AntigenTestB.dat", sep="")

out0303 <- G1mean(xB, prior=F) # Inference related to the normal distribution
#save(out0303,  file="./script/obje/out0303")
#load(file="./script/obje/out0303"); # Load the results of a pre-run of MCMC

# Scripts for the textbook's practical problems are omitted
# The results to address the following self-study are available in the textbook, so no new calculations are needed.

########################################################################

# Self-Study:
# 1) What can be inferred from Tables 5.1 and 5.2? From the perspective of the PHC curve, 
#    which drug, New Drug A or B, do you consider to be useful?
# 2) In the significance test for mu=15.0, 
#     New Drug A: p-value=0.00002188
#     New Drug B: p-value=0.0000009352
#    Clearly, the p-value for New Drug B is smaller. How does this align with 1)?








########################################################################
###The working directory should be 'PPandPHC'
#chapter15
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
library(cmdstanr)           #Load and attach add-on package cmdstanr
library(posterior)          #Load and attach add-on package posterior
library(psych)
library(rgl)
prob<-c(0.025, 0.05, 0.5, 0.95, 0.975) #Definition of probability points


################### 15.1 Analysis of Rectal Cancer Data
####### Table 15.1 Supply of food and rectal cancer corrected mortality rate
(can <- read.csv("dat/rectal_cancer.csv",header=T))

####### Table 15.2 Summary statistics of "Rectal Cancer Data"
(mean<-apply(can,2,"mean")); # Table 15.2
vari<-function(x){mean((x-mean(x))^2)} # Definition of function to calculate sample variance
(SD<-sqrt(apply(can,2,"vari"))) # Standard deviation

####### Figure 15.1 Multivariate scatter plot
pairs.panels(can,smooth=F,density=F,ellipses=F,pch=16,rug=F,cex.axis=1.0)
#dev.copy2pdf(file="fig15_1.pdf") 

####### Table 15.3 Covariance matrix of "Rectal Cancer Data"
round(var(can),0)

####### Table 15.4 Correlation matrix of "Rectal Cancer Data"
round(cor(can),2)

################### 15.2 Multiple Regression Model

############## Execution of multiple regression analysis using functions
out1<-Reg(can[,3],can[,c(2,1)]/100)
#save(out1,file="./script/obje/out1501")
#load(file="./script/obje/out1501");#Load object previously created with stan()

####### Table 15.5 Summary of posterior distribution
gqcal(out1$ext[,1:6])
pri1<-print(out1,3)

####### Figure 15.2 Posterior distribution of the coefficient of determination (left) and
####### Posterior distribution of the standardized partial regression coefficient of "Dairy Products" (right)
par(mfrow=c(1,2))
hist(out1$r2,breaks=50,cex.axis=2.0,cex.lab=2.0,main="",
     xlab="Coefficient of Determination", ylab="")
hist(out1$sb[,1],breaks=50,cex.axis=2.0,cex.lab=2.0,main="",
     xlab="Standardized Partial Regression Coefficient of 'Dairy Products'", ylab="")
par(mfrow=c(1,1))
#dev.copy2pdf(file="fig15_2.pdf")


####### Figure 15.3 ROPE for 'Dairy Products' and PHC curve for 'Total Calories'
par(mfrow=c(2,1))
PHC01(seq(0,0.5,0.01),out1$sb[,1],cc="rope",byoga="yes",
      xlab="ROPE of Standardized Partial Regression Coefficient of 'Dairy Products'")
PHC01(seq(0.65,1.2,0.01),out1$sb[,2],cc="gtc",byoga="yes",
      xlab="Standardized Partial Regression Coefficient of 'Total Calories'")
par(mfrow=c(1,1))
#dev.copy2pdf(file="fig15_3.pdf")


####### Table 15.6 PHC Table for Regression Coefficients and Coefficient of Determination
PHC01(seq(0,0.3,0.05),out1$sb[,1],cc="rope",byoga="no")
PHC01(seq(0.68,0.8,0.02),out1$sb[,2],cc="gtc",byoga="no")


######## Figure 15.4 3D Scatter Plot and Regression Plane (Unit: kcal)
plot3d(can[,"Total.Calories"],can[,"Dairy.Products"],can[,"Rectal.Cancer"],size=7,
       lwd=1.5,type="h",xlab="", ylab="", zlab="");                      
planes3d(a=0.00594,b=-0.00403,c=-1,d=-10.97805, alpha = 0.1)
#rgl.snapshot("fig15_4.png")


######## Figure 15.5 Residual Plot (Scatter plot for predicted values ^y and residuals e)
plot(pri1$yhat,pri1$resi,type="n",xlab="Predicted Values",ylab="Residual",
     cex.axis=2.0,cex.lab=1.5)
abline(h=0,lwd=2)
text(pri1$yhat,pri1$resi,rownames(can),cex=1.5)
#dev.copy2pdf(file="fig15_5.pdf")

######## Table 15.7 Predicted values, posterior standard deviation, 95th percentile of the posterior distribution, and 95th percentile of the posterior predictive distribution
Xnew <- matrix(c(240, 2800, 140, 2700, 340, 2900),,2,T); # Order of Dairy Products and Total Calories
XnewS <- Xnew/100;             # Dividing by 100 to match the state when estimated
pri3 <- print(out1, Xnew=XnewS)
(tab1406 <- data.frame(
  Dairy_Products=Xnew[,1],
  Total_Calories=Xnew[,2],
  Predicted_Value=round(pri3$yhatc[,"EAP"],1),
  Postsd=round(pri3$yhatc[,"post.sd"],1),
#  Confidence_025=round(pri3$yhatc[,"0.025"],1),
  Confidence_95=round(pri3$yhatc[,"0.95"],1),
  SD=round(pri3$yastc[,"sd"],1),
#  Predicted_025=round(pri3$yastc[,"0.025"],1),
  Predicted_95=round(pri3$yastc[,"0.95"],1)))

######## Self-Study Questions
#### 1. Draw a histogram of the posterior distribution of the standardized partial regression coefficient for 'Total Calories'.

#### 2. Draw the phc curve of the ROPE of the multiple correlation coefficient in the interval (0.0, 1.0).

#### 3. Create a phc table of the ROPE of the multiple correlation coefficient with seq(0.6,0.9,0.05).


########################################################################
# Chapter 13 Script (Inference of Multinomial Distribution)
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
   library(cmdstanr)           #load and attach add-on package cmdstanr
   library(posterior)          #load and attach add-on package posterior

########### 13.2   Estimation of Ratios (One Multinomial Distribution) #####################

######## Pet Problem 
# Table 13.1 First Pet Owned (Number of People)
  x1 <- c(32,29,18,15,10)


# Estimation of proportions for k categories
out1301<-Mu01(x1)
#save(out1301,  file="./script/obje/out1301")
#load(file="./script/obje/out1301"); # Load the results of a pre-run of MCM

# Table 13.2 Posterior distribution of proportions for the 'Pet Problem'
gqcal(out1301$pi)                 ; #

# Table 13.3 Probability that the proportion in row i is greater than the proportion in column j
PHC02(c=0.0 ,out1301$pi,cc="gtc") ; #

# Figure 13.1 PHC curves of the difference and ratio of proportions (comparison of dogs and birds)
par(mfrow=c(2,1))
PHC01(seq(0,0.25,0.01),out1301$pi[,1],out1301$pi[,4],cc="gtc",byoga="yes",xlab="Difference in proportions p1-p4");
PHC01(seq(0.0,1.0,0.02),out1301$pi[,4]/out1301$pi[,1],cc="ltc",byoga="yes",xlab="Ratio of proportions p4/p1");
par(mfrow=c(1,1))
#dev.copy2pdf(file="../English_tex/fig/chap13/z1301.pdf")

# Table 13.4 PHC table of the difference and ratio of proportions (comparison of dogs and birds)
PHC01(seq(0.02,0.08,0.01),out1301$pi[,1],out1301$pi[,4],cc="gtc",byoga="no");
PHC01(seq(0.7, 1.0 ,0.05),out1301$pi[,4]/out1301$pi[,1],cc="ltc",byoga="no");

################ 13.3 Inference of 2x2 Contingency Tables with Paired Data #########################

######## Brand Recognition Problem 2
# Table 13.5 Cross Table of Awareness of Two Brands (Number of People)
  x2 <- matrix(c(
       70,30,
       28,72),2,2,T)
apply(x2,1,"sum");apply(x2,2,"sum")

# Table 13.8 Sample Ratios of Brand Awareness
x2/200;
apply(x2,1,"sum")/200;apply(x2,2,"sum")/200

# Analysis of 2x2 Contingency Tables with Paired Data
out1302<-Mu02(x2)
#save(out1302,  file="./script/obje/out1302")
#load(file="./script/obje/out1302"); # Load the results of a pre-run of MCM

#Table13.9 Posterior Distribution of Parameters for'Brand Recognition Problem2'
# pim: Parental Ratio for row i, column j (a × b)
gqcal(out1302$pim)

#Table13.12 Summary Statistics of Posterior Distribution of eij·V (Upper Brand)
# res: Pearson Residual (a × b)   V: Cramer's Coefficient of Association
gqcal(out1302$res)
gqcal(out1302$V)

# Table 13.13 Summary Statistics of Posterior Distribution of 
#    Ratio Lij of Joint Probabilities and Product of Marginal Probabilities
gqcal(out1302$L)

# Table 13.14 PHC Table for L11 and L12
PHC01(seq(1.24,1.30,0.01),out1302$L[,1],cc="gtc",byoga="no");
PHC01(seq(0.69,0.75,0.01),out1302$L[,3],cc="ltc",byoga="no");


###### Playing Cards
# Table 13.10 Cross Table of Cards
  x3 <- matrix(c(
       20,20,
        6, 6),2,2,T)
apply(x3,1,"sum");apply(x3,2,"sum")

# Table 13.11 Probabilities of Cards Cross Table
round(x3/52,3)
round(apply(x3,1,"sum")/52,3)
round(apply(x3,2,"sum")/52,3)

# Analysis of 2x2 Contingency Tables with Paired Data
out1303<-Mu02(x3)
#save(out1303,  file="./script/obje/out1303")
#load(file="./script/obje/out1303"); # Load the results of a pre-run of MCM

#Table13.12 Summary Statistics of Posterior Distribution of eij·V (Lower Cards)
gqcal(out1303$res)
gqcal(out1303$V)

#Table 13.13 Summary Statistics of Posterior Distribution of 
#    Ratio Lij of Joint Probabilities and Product of Marginal Probabilities
gqcal(out1303$L)

###### 13.4 Inference of a×b Contingency Tables with Paired Data ###########
##### Pasta Problem
# Table 13.15 Toppings Chosen for Pasta
  x4 <- matrix(c(
       19,  9,  6,
       10, 19,  5,
       15, 14, 18),3,3,T)
colnames(x4)<-c('Basil','Truffle','None')
rownames(x4)<-c('Cold Tomato','Carbonara','Peperoncino')
x4
apply(x4,1,"sum");apply(x4,2,"sum");sum(x4)

# Analysis of a×b Contingency Tables with Paired Data
out1304<-Mu02(x4)
#save(out1304,  file="./script/obje/out1304")
#load(file="./script/obje/out1304"); # Load the results of a pre-run of MCM

# Table 13.18 Numerical Summary of Posterior Distribution of Joint 
#     Probabilities (parameters) for 'Pasta Problem'
gqcal(out1304$pim)

# Table 13.19 Numerical Summary of Posterior Distribution of Cramer's 
#    Coefficient of Association and Marginal Probabilities for 'Pasta Problem'
# pa:Marginal Probability of row i (a) pb:Marginal Probability of column j (b)
gqcal(out1304$V)
gqcal(out1304$pa)
gqcal(out1304$pb)

# Table 13.20 Probability that Joint Probabilities are Greater than 
#     the Product of Marginal Probabilities
# Up: Probability that Pearson Residual is Greater than 0 (a × b)
# Um: Probability that Pearson Residual is Less than or Equal to 0 (a × b)
gqcal(out1304$Up)
gqcal(out1304$Um)


# Table 13.21 Summary Statistics of Posterior Distribution of 
#     Ratio of Proportions (First Half)
gqcal(out1304$L)

# Table 13.21 PHC Table for Ratio of Proportions (Second Half)
PHC01(seq(1.10,1.40,0.05),out1304$L[,5],cc="gtc",byoga="no");
PHC01(seq(1.10,1.40,0.05),out1304$L[,9],cc="gtc",byoga="no");

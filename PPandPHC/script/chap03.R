########################################################################
###The working directory should be 'PPandPHC'
#Chapter 3
(n_wd<-getwd())                #Confirmation of working directory
source('myfunc/myfunc.R')      #Loading self-made functions
   library(cmdstanr)           #Load and attach add-on package cmdstanr
   library(posterior)          #Load and attach add-on package posterior
   library(bayesplot)          #Load and attach add-on package bayesplot

#Table 1.1: Measured Results (Days)
#antigen test data
x<-c(
11.5, 14.0, 15.0, 10.0, 14.0, 14.5, 12.5, 12.5, 12.5, 12.0, 
13.0, 12.0, 12.5, 13.5, 15.0, 13.5, 12.5, 12.0, 12.0, 17.0)

#G1mean is a inference regarding the Normal Distribution. x is data in vector form.
#Prior specifies the range of posterior distribution if T; if F, does not specify (leave it to stan)
out0301 <-G1mean(x,prior=F)

#save(out0301,  file="./script/obje/out0301")
#load(file="./script/obje/out0301"); #Load the results of a pre-run of MCMC

#figure 3.1 Trace plot of mean μ (upper panel) and standard deviation σ (lower panel)
mcmc_trace(out0301$draw,pars ="mu" )
mcmc_trace(out0301$draw,pars ="sigma" )
#dev.copy2pdf(file=”../English tex/fig/chap03/fig3_1.pdf”) 

#Table 3.1 Evaluation of random number sequences
#Table 3.2 Estimated results for the parameters of 'Antigen Testing.'
gqcal(out0301$ext[,1:3])
summary(out0301$draw)

#Figure 3.3 Posterior distribution of mean μ
hist(out0301$mu,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file=”../English tex/fig/chap03/fig3_3.pdf”) 

#Figure 3.4 Posterior distribution of standard deviation σ
hist(out0301$sigma,breaks=150,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file=”../English tex/fig/chap03/fig3_4.pdf”) 

#Figure 3.5 Posterior predictive distribution
hist(out0301$xaste,breaks=100,cex.axis=2.0,cex.lab=2.0,main="",xlab="",ylab="")
#dev.copy2pdf(file=”../English tex/fig/chap03/fig3_5.pdf”) 

#Upper and lower bounds for conditional predictive distribution with EAP
round(13.07- 1.96*1.65, 2)
round(13.07+ 1.96*1.65, 2)


########################################################################
#Chapter 3 3.6 Practical assignment
#Table 1.3 Measurement Results ( sec.) (Data input for 'Perceived Time')
y<-c(31.43,31.09,33.38,30.49,29.62,35.40,32.58,28.96,29.43,28.52,
     25.39,32.68,30.51,30.15,32.33,30.43,32.50,32.07,32.35,31.57)

#Inference regarding the normal distribution
out0302 <-G1mean(y,prior=F) #Inference regarding the normal distribution
#save(out0302,  file="./script/obje/out0302")
#load(file="./script/obje/out0302"); #Load the results of a pre-run of MCMC

#Equivalent to Figure 3.1 and Table 3.2
mcmc_trace(out0302$draw,pars ="mu" )
mcmc_trace(out0302$draw,pars ="sigma" )
gqcal(out0302$ext[,1:3])
summary(out0302$draw)

hist(out0302$mu,breaks=100)       #Posterior distribution of means equivalent to Figure 3.3
hist(out0302$sigma,breaks=100)    #Posterior distribution of sd equivalent to Figure 3.4
hist(out0302$xaste,breaks=100)    #posterior predictive distribution equivalent to Figure 3.5

########################################################################
#Self-study Chapter 3 3.7 Practice Questions
#1) Calculate and show below the results of the estimation of the parameters equivalent to Table 3.2 with data B.
#2) Compare the µEAPs of data A and B, compare the effects, and evaluate and interpret them.
#3) Compare the σEAP of data A and B, compare the properties, and evaluate and interpret them.
#4) Calculate the probability that the number of days to negative result will be longer for new drug A and new drug B
#   than the average of 15 days without medication, compare the probabilities, and determine which is superior.

xB<-scan("dat/AntigenTestB.dat",sep="")

out0303 <-G1mean(xB,prior=F) #Inference regarding the normal distribution
#save(out0303,  file="./script/obje/out0303")
#load(file="./script/obje/out0303"); #Load the results of a pre-run of MCMC

# 1), 2), 3)　Equivalent to Table 3.2
gqcal(out0303$ext[,1:3])
summary(out0303$draw)

# 4) Hint: Use the point estimates of µ and σ to use the script from textbook 2.6, Practical assignment 2).
#    Evaluate 1-F(15|µ,σ).　What was the prefix of the distribution function?
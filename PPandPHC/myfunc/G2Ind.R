##########################################################################
#G2Ind, Inferring the difference between two independent groups
###Arguments
# x1: vector. Group 1 data (Group with large mean value)
# x2: vector. Group 2 data (Group with small mean value)
# EQU: logical. If 1, two variances are treated as equal.
#               If 0, two variances are treated separately
# prior: logical. If it is TRUE, specify the range of the prior distribution.
#                 If it is FALSE, stan default is used.
# mL, mH, sL, sH: Prior distribution parameters
#    Upper and lower bounds of the mean (sd) value of a uniform distribution.
#prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points reported
#see=1234,cha=5,war=1000,ite=21000 : Parameter that specifies MCMC
###Value
# fit: output of stan, par: parameter list, prob: probability vector,
# mu1: 1st group average, mu2: 2nd group average, 
# sigma1: 1st group standard deviation, sigma2: 2nd group standard deviation,
# xaste1: Predicted distribution 1, xaste2: Predicted distribution 2, 
# log_lik: Log likelihood
##########################################################################
G2Ind <- function(x1,x2,EQU=1,prior=F,mL=-1000, mH=1000, sL=0, sH=100,
         prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
         see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)  
   library(posterior) 
   if (prior) {
     dat <- list(n1=length(x1),n2=length(x2),x1=x1,x2=x2,EQU=EQU,prior=prior,mL=mL,mH=mH,sL=sL,sH=sH)
     scr <- "stan/G2IndPT.stan" 
   } else {
     dat <- list(n1=length(x1),n2=length(x2),x1=x1,x2=x2,EQU=EQU)
     scr <- "stan/G2IndPF.stan" 
   }
   par<-c("mu","sigma1","sigma2","xaste","log_lik")
      mod <- cmdstan_model(scr) 
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see) 
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   mu1<-ext$mu1; mu2<-ext$mu2; sigma1<-ext$sigma1; sigma2<-ext$sigma2; 
   xaste1<-ext$xaste1; xaste2<-ext$xaste2; log_lik<-ext$log_lik
   out<-list(EQU=EQU,prior=prior,ext=ext,draw=draw,par=par,
             prob=prob,mu1=mu1,mu2=mu2,sigma1=sigma1,sigma2=sigma2,
             xaste1=xaste1,xaste2=xaste2,log_lik=log_lik)
   class(out)<-'G2Ind'
   return(invisible(out))
}

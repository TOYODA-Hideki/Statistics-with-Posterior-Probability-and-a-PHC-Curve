##########################################################################
#Estimating the difference between two matched groups.
#■Input
#x:Data in a matrix format with n rows and 2 columns.
#  (The group with the larger average value is designated in the first column)
#EQU=1 : Boolean value. 1 for equal variance. 0 for different variances.
#prior=F : Boolean value.  T to specify the range of the prior distribution. F not to specify.
#mL=-1000, mH=1000, sL=0, sH=100 : Parameters of the prior distribution
#prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to be reported in the posterior distribution
#see=1234,cha=5,war=1000,ite=21000 : Parameters related to MCMC
#■Output
# fit:Output of stan, par:list of parameters, prob:probability vector,
# mu1:Average of the first group, mu2:Average of the second group,
# sigma1:Standard deviation of the first group, sigma2:Standard deviation of the second group,, 
# rho:Correlation coefficient, xaste1:Predictive distribution 1, 
# xaste2:predictive distribution 2, log_lik:log likelihood
##########################################################################
G2pair <- function(x,EQU=1,prior=F,mL=-1000, mH=1000, sL=0, sH=100,
            prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
            see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              #Load and attach add-on package cmdstanr
   library(posterior)             #Load and attach add-on package posterior
   if (prior) {
     dat <- list(n=nrow(x),x=x,EQU=EQU,mL=mL,mH=mH,sL=sL,sH=sH)
     scr <- "stan/G2pairPT.stan"                           # Stan's filename
   } else {
     dat <- list(n=nrow(x),x=x,EQU=EQU)
     scr <- "stan/G2pairPF.stan"                           # Stan's filename
   }
   par<-c("mu","sigma1","sigma2","rho","xaste","log_lik")# sampling variable
      mod <- cmdstan_model(scr)                   #Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see)    #MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   mu1<-ext$mu1; mu2<-ext$mu2; 
   sigma1<-ext$sigma1; sigma2<-ext$sigma2; rho<-ext$rho;
   xaste1<-ext$xaste1; xaste2<-ext$xaste2; log_lik<-ext$log_lik
   out<-list(EQU=EQU,prior=prior,ext=ext,draw=draw,par=par,
       prob=prob,mu1=mu1,mu2=mu2,sigma1=sigma1,sigma2=sigma2,rho=rho,
       xaste1=xaste1,xaste2=xaste2,log_lik=log_lik)
   class(out)<-'G2pair'
   return(invisible(out))
}

##########################################################################
#Inference regarding the normal distribution of one group
#■input
#x: Data in vector form
#prior=F : Logical value. if T, specifies range of prior distribution; if F, not specified.
#mL, mH, sL, sH : parameters of the prior distribution
#prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points reported in the posterior distribution
#see=1234,cha=5,war=1000,ite=21000 : Parameters related to MCMC
#■output
#fit:stan output, par: parameter list, prob: probability vector 
#mu: mean, sigma: standard deviation, xaste: predictive distribution, log_lik: log likelihood
##########################################################################
G1mean <- function(x,prior=F, mL=-1000, mH=1000, sL=0, sH=100,
prob=c(0.025, 0.05, 0.5, 0.95, 0.975),see=1234,cha=5,war=1000,ite=21000)
{
   library(cmdstanr)              #load and attach add-on package cmdstanr
   library(posterior)             #load and attach add-on package posterior
   if (prior) {
     dat <- list(n=length(x),x=x,mL=mL,mH=mH,sL=sL,sH=sH)
     scr <- "stan/G1meanPT.stan"                           # Stan File Name
   } else {
     dat <- list(n=length(x),x=x)
     scr <- "stan/G1meanPF.stan"                           # Stan File Name
   }
   par<-c("mu","sigma","xaste","log_lik")        # sampling variable
      mod <- cmdstan_model(scr)            # compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see)    # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   mu<-ext$mu; sigma<-ext$sigma; xaste<-ext$xaste; log_lik<-ext$log_lik
   out<-list(ext=ext,draw=draw,par=par,prob=prob,mu=mu,
             sigma=sigma,xaste=xaste,log_lik=log_lik)
   class(out)<-'G1mean'
   return(invisible(out))
}
##########################################################################
#Methods for printing
#x: object of class 'G1mean'
#degits=3 : Rounding of decimals
#cr1=F : base point of effect size (specified by measured value)
#cr2a=F : distribution function upper bound (specified by measured value)
#cr2b=F : lower limit of distribution function (specified by measured value)
#cr3=F : Reference point for ratio to measured value (specified by measured value)
#cr4=F : Probability that the mean is less than cr4 (specified by measured value)
#cr5=F : Probability that the effect size by cr1 is smaller than cr5 (specified by effect size)
#pr1=F : Base point of % point (specified by probability from the bottom)
#pr2=F : Probability that the area of cr2 is observed is greater than pr2 (specified by probability)
##########################################################################
print.G1mean<-function(x,degits=3,cr1=F,cr2a=F,cr2b=F,cr3=F,cr4=F,cr5=F,pr1=F,pr2=F)
{
   mu<-x$mu; sigma<-x$sigma; xaste<-x$xaste; log_lik<-x$log_lik
   prob<-x$prob
   G<-matrix(0,length(mu),6)
   G[,1]<- sigma^2;                  #variance
   G[,2]<- sigma/mu;                 #coefficient of variation
   co<-2; lab<-c("variance","coefficient of variation")
   if(is.numeric(cr1)){co<- co+1;G[,co]<- (mu-cr1)/sigma;
      lab<-c(lab,paste("effect size(",cr1,")",sep=""))}
   if(is.numeric(pr1)){co<- co+1;G[,co]<- mu+qnorm(pr1)*sigma;
      lab<-c(lab,paste("percentile(",pr1,")",sep=""))}   
   if((is.numeric(cr2a))&(!is.numeric(cr2b)))
                      {co<- co+1;G[,co]<- pnorm(cr2a,mu,sigma);
      lab<-c(lab,paste("Probability that measured value is less than (",cr2a,")",sep=""))}
   if((is.numeric(cr2b))&(!is.numeric(cr2a)))
                      {co<- co+1;G[,co]<- 1-pnorm(cr2b,mu,sigma);
      lab<-c(lab,paste("Probability that the measured value is greater than (",cr2b,")",sep=""))}
   if((is.numeric(cr2a))&(is.numeric(cr2b)))
                      {co<- co+1;G[,co]<- pnorm(cr2a,mu,sigma)-pnorm(cr2b,mu,sigma);
      lab<-c(lab,paste("Probability that the measured value is greater than (",cr2b,") and less than (",cr2a,")",sep=""))}
   if(is.numeric(cr3)){co<- co+1;G[,co]<- xaste/cr3;
      lab<-c(lab,paste("Ratio to reference point (",cr3,")",sep=""))}
   Gc<-cbind(
      apply(G[,1:co],2,mean),
      apply(G[,1:co],2,sd),
      t(apply(G[,1:co],2,quantile, probs=prob))
   )
   colnames(Gc)<-c("EAP","post.sd",prob)
   rownames(Gc)<-lab
   U<-matrix(0,length(mu),4)
   co<-0;lab<-NULL
   if(is.numeric(cr4)){co<- co+1;U[,co]<- ifelse(mu< cr4,    1, 0);
         lab<-c(lab,paste("Probability that mean is less than (",cr4,")",sep=""))}
   if((is.numeric(cr2a))&(!is.numeric(cr2b)))
                      {co<- co+1;U[,co]<- ifelse(xaste< cr2a,1,0);
      lab<-c(lab,paste("Probability that measured value is less than (",cr2a,") (point estimate)",sep=""))}
   if((is.numeric(cr2b))&(!is.numeric(cr2a)))
                      {co<- co+1;U[,co]<- ifelse(xaste> cr2b,1,0);
      lab<-c(lab,paste("Probability that measured value is greater than (",cr2b,") (point estimate)",sep=""))}
   if((is.numeric(cr2a))&(is.numeric(cr2b)))
                      {co<- co+1;U[,co]<- ifelse((xaste< cr2a)&(xaste> cr2b),1,0);
      lab<-c(lab,paste("Probability that measured value is greater than (",cr2b,") and less than (",cr2a,") (point estimate)",sep=""))}
   if(is.numeric(cr5)){co<- co+1;U[,co]<- ifelse(((mu-cr1)/sigma)< cr5, 1, 0);
         lab<-c(lab,paste("Probability that the effect size due to (",cr1,") is less than (",cr5,")",sep=""))}
   if((is.numeric(cr2a))&(!is.numeric(cr2b)))
                      {co<- co+1;U[,co]<- ifelse(pnorm(cr2a,mu,sigma)>pr2,1,0);
      lab<-c(lab,paste("Probability that measured value is less than (",cr2a,") but greater than (",pr2,")",sep=""))}
   if((is.numeric(cr2b))&(!is.numeric(cr2a)))
                      {co<- co+1;U[,co]<- ifelse(1-pnorm(cr2b,mu,sigma)>pr2,1,0);
      lab<-c(lab,paste("Probability that the measured value is greater than (",cr2b,") is greater than (",pr2,")",sep=""))}
   if((is.numeric(cr2a))&(is.numeric(cr2b)))
                      {co<- co+1;U[,co]<- ifelse(pnorm(cr2a,mu,sigma)-pnorm(cr2b,mu,sigma)>pr2,1,0);
      lab<-c(lab,paste("Probability that the measured value is greater than (",cr2b,") and less than (",cr2a,") is greater than (",pr2,")",sep=""))}
   if(co>0.9){
      Uc<-matrix(colMeans(as.matrix(U[,1:co])),co,1);
      rownames(Uc)<-lab
      }else{Uc<-0.0}
   print(x$fit,pars=x$par,digits_summary=degits,probs=x$prob)
   print(round(Gc,degits))
   print(round(Uc,degits))
   out<-list(G=G,Gc=Gc,Uc=Uc)
   return(invisible(out))
}

##########################################################################
# Estimation of Ratios (One Binomial Distribution)
#■ Input
# x: Number of positive responses
# n: Number of Bernoulli trials
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : MCMC related parameters
#■ Output
# fit: Output from Stan, par: Parameter list, prob: Probability vector,
# theta: Parent ratio, xaste: Predictive distribution, Odds: Odds
# log_lik: Log likelihood
##########################################################################
Bi01 <- function(x,n,prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   dat <- list(x=x,n=n)
   scr <- "stan/Bi01.stan"                                # Stan file name
   par<-c("theta","xaste","Odds","log_lik")               # Sampling variables
   if (pack == "rstan") {
      fit<-stan(file=scr,data=dat,pars=par,
                seed=see,chains=cha,warmup=war,iter=ite)
   }else{
      mod <- cmdstan_model(scr)                   # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see)    # MCMC
      fit <- rstan::read_stan_csv(csrfit$output_files()) # Convert to Stan format
   }
   ext<-extract(fit,par);
   out<-list(fit=fit,par=par,prob=prob,
        theta=ext$theta,xaste=ext$xaste, Odds=ext$Odds,log_lik=ext$log_lik);
   class(out)<-'Bi01'
   return(invisible(out))
}
##########################################################################
# Print method
#■ Input
# x: Object of class 'Bi01'
# degits=3 : Rounding of decimals
# pr1=F : Minimum value for theta
# cr1=F : Minimum value for x*
# cr2=F : Minimum value for Odds
##########################################################################
print.Bi01<-function(x,degits=3,pr1=F,cr1=F,cr2=F)
{
   log_lik<-x$log_lik;theta<-x$theta;xaste<-x$xaste;Odds<-x$Odds;
   waic<- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
   print(x$fit,pars=x$par,digits_summary=degits,probs=x$prob)
   U<-matrix(0,length(theta),3)
   co<-0;lab<-NULL
   if(is.numeric(pr1)){co<- co+1;U[,co]<- ifelse(theta> pr1, 1, 0);
         lab<-c(lab,paste("Proportion greater than (",pr1,")",sep=""))}
   if(is.numeric(cr1)){co<- co+1;U[,co]<- ifelse(xaste> cr1, 1, 0);
         lab<-c(lab,paste("x* greater than (",cr1,")",sep=""))}
   if(is.numeric(cr2)){co<- co+1;U[,co]<- ifelse(Odds>  cr2, 1, 0);
         lab<-c(lab,paste("Odds greater than (",cr2,")",sep=""))}
   if(co>0.9){
      Uc<-matrix(colMeans(as.matrix(U[,1:co])),co,1);
      rownames(Uc)<-lab
      }else{Uc<-0.0}
   print(round(Uc,degits))
   print(paste("waic=",round(waic,degits),sep=""))
   out<-list(waic=waic,Uc=Uc)
   return(invisible(out))
}

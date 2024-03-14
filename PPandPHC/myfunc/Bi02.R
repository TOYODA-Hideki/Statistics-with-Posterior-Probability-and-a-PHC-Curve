##########################################################################
# Statistical Inference for Ratios in Two Independent Groups (Two Binomial Distributions)
#■ Input
# x: Number of positive responses (vector of length 2, with the larger sample ratio first)
# n: Number of Bernoulli trials (vector of length 2)
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : Parameters related to MCMC
#■ Output
# fit: Output from Stan, par: List of parameters, prob: Probability vector, 
# p: Parent ratio (2)
# xaste: Predictive distribution of x* (2)
# p_sa: Difference in proportions
# p_hi: Ratio of proportions
# Odds: Odds (2)
# Odds_hi: Odds ratio
# log_lik: Log likelihood
##########################################################################
Bi02 <- function(x,n,prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   dat <- list(x=x,n=n)
   scr <- "stan/Bi02.stan"                                # Stan file name
   par<-c("p","xaste","p_sa","p_hi","Odds","Odds_hi","log_lik") # Sampling variables
   if (pack == "rstan") {
      fit<-stan(file=scr,data=dat,pars=par,
                seed=see,chains=cha,warmup=war,iter=ite)
   } else {
      mod <- cmdstan_model(scr)                           # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
                           iter_warmup=war,parallel_chains=cha,seed=see) # MCMC
      fit <- rstan::read_stan_csv(csrfit$output_files())  # Convert to Stan format
   }
   ext<-extract(fit,par);
   out<-list(fit=fit,par=par,prob=prob,
             p=ext$p, xaste=ext$xaste, log_lik=ext$log_lik, p_sa=ext$p_sa, 
             p_hi=ext$p_hi,Odds=ext$Odds, Odds_hi=ext$Odds_hi);
   class(out)<-'Bi02'
   return(invisible(out))
}

##########################################################################
# Print method
#■ Input
# x: Object of class 'Bi02'
# digits=3 : Rounding of decimals
# cr1=F : Minimum value for difference in proportions
# cr2=F : Minimum value for ratio of proportions
# cr3=F : Minimum value for odds ratio
##########################################################################
print.Bi02 <- function(x,digits=3,cr1=F,cr2=F,cr3=F)
{
   log_lik <- x$log_lik; p <- x$p; xaste <- x$xaste; p_sa <- x$p_sa; 
   p_hi <- x$p_hi; Odds <- x$Odds;  Odds_hi <- x$Odds_hi;
   waic <- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
   print(x$fit,pars=x$par,digits_summary=digits,probs=x$prob)
   U <- matrix(0,length(p),3)
   co <- 0; lab <- NULL
   if (is.numeric(cr1)) {co <- co+1; U[,co] <- ifelse(p_sa > cr1, 1, 0);
                         lab <- c(lab, paste("Proportion difference greater than (",cr1,")",sep=""))}
   if (is.numeric(cr2)) {co <- co+1; U[,co] <- ifelse(p_hi >  cr2, 1, 0);
                         lab <- c(lab, paste("Proportion ratio greater than (",cr2,")",sep=""))}
   if (is.numeric(cr3)) {co <- co+1; U[,co] <- ifelse(Odds_hi > cr3, 1, 0);
                         lab <- c(lab, paste("Odds ratio greater than (",cr3,")",sep=""))}
   if (co > 0.9) {
      Uc <- matrix(colMeans(as.matrix(U[,1:co])), co, 1);
      rownames(Uc) <- lab
   } else {
      Uc <- 0.0
   }
   print(round(Uc, digits))
   print(paste("waic=", round(waic, digits), sep=""))
   out <- list(waic=waic, Uc=Uc)
   return(invisible(out))
}

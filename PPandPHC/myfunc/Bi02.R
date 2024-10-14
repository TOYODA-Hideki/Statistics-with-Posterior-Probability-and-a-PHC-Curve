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
                 see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   dat <- list(x=x,n=n)
   scr <- "stan/Bi02.stan"                                # Stan file name
   par<-c("p","xaste","p_sa","p_hi","Odds","Odds_hi","log_lik") # Sampling variables
      mod <- cmdstan_model(scr)                           # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
                           iter_warmup=war,parallel_chains=cha,seed=see) # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   out<-list(ext=ext,draw=draw,par=par,prob=prob,
         p=as.matrix(ext[,1:2]), xaste=as.matrix(ext[,3:4]), 
         p_sa=as.matrix(ext$p_sa), p_hi=as.matrix(ext$p_hi),
         Odds=as.matrix(ext[,7:8]), Odds_hi=as.matrix(ext$Odds_hi),
         log_lik=ext$log_lik);  
   class(out)<-'Bi02'
   return(invisible(out))
}

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
                 see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   dat <- list(x=x,n=n)
   scr <- "stan/Bi01.stan"                                # Stan file name
   par<-c("theta","xaste","Odds","log_lik")          # sampling variables
      mod <- cmdstan_model(scr)                   # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see)    # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   out<-list(ext=ext,par=par,prob=prob,
        theta=ext$theta,xaste=ext$xaste, Odds=ext$Odds,log_lik=ext$log_lik);
   class(out)<-'Bi01'
   return(invisible(out))
}

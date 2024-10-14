##########################################################################
# Statistical Inference for Ratios in g Independent Groups (g Binomial Distributions)
#■ Input
# x: Number of positive responses (vector of length g)
# n: Number of Bernoulli trials (vector of length g)
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : Parameters related to MCMC
#■ Output
# fit: Output from Stan, par: List of parameters, prob: Probability vector, 
# p: Parent ratios (g)
# xaste: Predictive distribution of x* (g)
# U2: Matrix indicating if the row's parent ratio is greater than the column's (g*g)
# log_lik: Log likelihood
##########################################################################
Bi03 <- function(x,n,prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   g <- length(x)
   dat <- list(g=g,x=x,n=n)
   scr <- "stan/Bi03.stan"                                # Stan file name
   par<-c("p","xaste","U2","log_lik")                     # Sampling variables
      mod <- cmdstan_model(scr)                           # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
                      iter_warmup=war,parallel_chains=cha,seed=see) # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   out<-list(ext=ext,draw=draw,par=par,prob=prob,
             p=ext[,1:g], xaste=ext[,c((g+1):(2*g))], 
             log_lik=ext$log_lik);
   class(out)<-'Bi03'
   return(invisible(out))
}

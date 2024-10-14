##########################################################################
# Statistical Inference for Proportions of k Categories
#■ Input
# x: Number of positive responses (vector of length k)
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : Parameters related to MCMC
#■ Output
# fit: Output from Stan, par: List of parameters, prob: Probability vector, 
# pi: Parent proportions (k)
# U2: Matrix indicating if the row's parent proportion is greater than the column's (k*k)
# log_lik: Log likelihood
##########################################################################
Mu01 <- function(x,prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   k <- length(x)
   dat <- list(k=k,x=x)
   scr <- "stan/Mu01.stan"                                # Stan file name
   par<-c("pi","U2","log_lik")                            # Sampling variables
      mod <- cmdstan_model(scr)                           # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
                           iter_warmup=war,parallel_chains=cha,seed=see) # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   out<-list(ext=ext,draw=draw,par=par,prob=prob,
        pi=as.matrix(ext[,1:k]), U2=as.matrix(ext[,c((k+1):(k+k*k))]),
       log_lik=ext$log_lik);
   class(out)<-'Mu01'
   return(invisible(out))
}

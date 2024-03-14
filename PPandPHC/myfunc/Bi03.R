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
                 see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   g <- length(x)
   dat <- list(g=g,x=x,n=n)
   scr <- "stan/Bi03.stan"                                # Stan file name
   par<-c("p","xaste","U2","log_lik")                     # Sampling variables
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
             p=ext$p, xaste=ext$xaste, U2=ext$U2, log_lik=ext$log_lik);
   class(out)<-'Bi03'
   return(invisible(out))
}

##########################################################################
# Print method
#■ Input
# x: Object of class 'Bi03'
# digits=3 : Decimal rounding
##########################################################################
print.Bi03 <- function(x,digits=3)
{
   log_lik <- x$log_lik; p <- x$p; xaste <- x$xaste; U2 <- x$U2; 
   waic <- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
   print(x$fit,pars=x$par,digits_summary=digits,probs=x$prob)
   print(paste("waic=", round(waic, digits), sep=""))
   out <- list(waic=waic)
   return(invisible(out))
}

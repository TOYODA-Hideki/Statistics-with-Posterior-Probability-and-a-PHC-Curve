##########################################################################
# Statistical Inference for a×b Contingency Tables with Paired Data
#■ Input
# x: Contingency table (specified in a × b matrix format)
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : Parameters related to MCMC
#■ Output
# fit: Output from Stan, par: List of parameters, prob: Probability vector, 
# pim: Parent proportion for row i, column j (a × b)
# xaste: Predictive distribution of x* (a × b)
# res: Pearson residuals (a × b)
# V: Cramer's V coefficient
# pa: Marginal probability for row i (a)
# pb: Marginal probability for column j (b)
# Up: Probability that Pearson residual is greater than 0 (a × b)
# Um: Probability that Pearson residual is less than or equal to 0 (a × b)
# L: Ratio of joint probabilities to the product of marginal probabilities (a × b)
# log_lik: Log likelihood
##########################################################################
Mu02 <- function(x,prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                 see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   library(rstan)
   dat <- list(a=nrow(x),b=ncol(x),x=x)
   scr <- "stan/Mu02.stan"                                    # Stan file name
   par<-c("pim","xaste","log_lik","res","V","pa","pb","Up","Um","L") # Sampling variables
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
             pim=ext$pim, xaste=ext$xaste,log_lik=ext$log_lik,res=ext$res,
             V=ext$V,pa=ext$pa,pb=ext$pb,Up=ext$Up,Um=ext$Um,L=ext$L);
   class(out)<-'Mu02'
   return(invisible(out))
}

##########################################################################
# Print method
#■ Input
# x: Object of class 'Mu02'
# digits=3 : Decimal rounding
##########################################################################
print.Mu02 <- function(x,digits=3)
{
   log_lik <- x$log_lik; 
   waic <- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
   print(x$fit,pars=x$par,digits_summary=digits,probs=x$prob)
   print(paste("waic=", round(waic, digits), sep=""))
   out <- list(waic=waic)
   return(invisible(out))
}

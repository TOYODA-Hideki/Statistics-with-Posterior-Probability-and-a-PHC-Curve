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
                 see=1234, cha=5, war=1000, ite=21000)
{
   library(cmdstanr)              # Load the cmdstanr package
   library(posterior)             # Load the posterior package
   a<-nrow(x); b<-ncol(x);  dat <- list(a=a, b=b, x=x); ab<-a*b;
   scr <- "stan/Mu02.stan"                              # Stan file name
   par<-c("pim","xaste","res","Up","Um","L","pa","pb","V","log_lik") 
      mod <- cmdstan_model(scr)                           # Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
                           iter_warmup=war,parallel_chains=cha,seed=see) # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
   out<-list(ext=ext,draw=draw,par=par,prob=prob,
      pim=as.matrix(ext[,1:ab]), xaste=as.matrix(ext[,c((ab+1):(2*ab))]),
      res=as.matrix(ext[,c((2*ab+1):(3*ab))]),
      Up =as.matrix(ext[,c((3*ab+1):(4*ab))]),
      Um =as.matrix(ext[,c((4*ab+1):(5*ab))]),
      L  =as.matrix(ext[,c((5*ab+1):(6*ab))]),
      pa =as.matrix(ext[,c((6*ab+1):(6*ab+a))]),
      pb =as.matrix(ext[,c((6*ab+a+1):(6*ab+a+b))]),
      V=ext$V,log_lik=ext$log_lik);
   class(out)<-'Mu02'
   return(invisible(out))
}

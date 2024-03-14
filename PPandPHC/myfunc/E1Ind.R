##########################################################################
# Analyzing data from an independent single-factor experiment
# ■ Input
# y: Characteristic values in vector form
# A: Levels in vector form (specified as integers from 1 to a, without gaps)
# prior=F : Logical value. If T, specify the range of the prior distribution. If F, do not specify.
# mL=-1000, mH=1000, sL=0, sH=100: Parameters of the prior distribution
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000: Parameters related to MCMC
# ■ Output
# fit: output from stan, par: list of parameters, prob: probability vector, 
# muA: mean of levels, sigmaE: within-level sd, sigmaA: between-level sd, eta2: explanatory rate
# delta: effect size, mu: overall mean, aj: effect of levels, Ubig: effect of levels greater than 0,
# Usma: effect of levels less than or equal to 0, U2: probability that row i is greater than column j, log_lik: log likelihood
##########################################################################
E1Ind <- function(y, A, prior=F, mL=-1000, mH=1000, sL=0, sH=100,
                  prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                  see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
  library(cmdstanr)              # Load and attach add-on package cmdstanr
  library(posterior)             # Load and attach add-on package posterior
  library(rstan)                 # Load and attach add-on package rstan
  if (prior) {
    dat <- list(n=length(y), a=max(A), y=y, A=A, mL=mL, mH=mH, sL=sL, sH=sH)
    scr <- "stan/E1IndPT.stan"    # Stan file name
  } else {
    dat <- list(n=length(y), a=max(A), y=y, A=A)
    scr <- "stan/E1IndPF.stan"    # Stan file name
  }
  par<-c("muA", "sigmaE", "sigmaA", "eta2", "delta",
         "mu", "aj", "Ubig", "Usma", "U2", "log_lik")  # Sampling variables
  if (pack == "rstan") {
    fit<-stan(file=scr, data=dat, pars=par, seed=see, chains=cha, warmup=war, iter=ite)
  }else{
    mod <- cmdstan_model(scr)                   # Compile
    csrfit <- mod$sample(data=dat, chains=cha, iter_sampling=ite,
                         iter_warmup=war, parallel_chains=cha, seed=see)    # MCMC
    fit <- rstan::read_stan_csv(csrfit$output_files()) # Convert to stan format
  }
  ext<-extract(fit, par); muA<-ext$muA; sigmaE<-ext$sigmaE; sigmaA<-ext$sigmaA; 
  eta2<-ext$eta2; delta<-ext$delta; mu<-ext$mu; aj<-ext$aj; Ubig<-ext$Ubig;
  Usma<-ext$Usma; U2<-ext$U2; log_lik<-ext$log_lik;
  out<-list(prior=prior, fit=fit, par=par, prob=prob, muA=muA,
            sigmaE=sigmaE, sigmaA=sigmaA, eta2=eta2, delta=delta, mu=mu, aj=aj,
            Ubig=Ubig, Usma=Usma, U2=U2, log_lik=log_lik)
  class(out)<-'E1Ind'
  return(invisible(out))
}
##########################################################################
# Print method
# ■ Input
# x: An object of class 'E1Ind'
# digits=3 : Rounding of decimals
##########################################################################
print.E1Ind<-function(x, digits=3)
{
  prior<-x$prior; log_lik<-x$log_lik;
  if (prior) {print("******************Prior Distribution Specified******************")}
  else {print("******************Default Prior Distribution*********************")}
  waic<- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
  print(x$fit, pars=x$par, digits_summary=digits, probs=x$prob)
  print(paste("waic=", round(waic, digits), sep=""))
  out<-list(waic=waic)
  return(invisible(out))
}  
##########################################################################
# Probability that a conjunctive proposition holds (common function for E1Ind and E2Ind)
# ■ Input
# x: Class U2[,i,j] where 1 if i is greater than j, 0 otherwise, or U2A, U2B
# digits=3 : Rounding of decimals
# IJ : A matrix with 2 columns (the first column for rows, the second for columns)
##########################################################################
printIJ<-function(x, digits=3, IJ)
{
  pro<-1
  for (i in 1:nrow(IJ)){
    pro<-pro*x[,IJ[i,1],IJ[i,2]]
  }
  print(round(mean(pro), digits))
}
##########################################################################
# Comparison of two levels (outputs superiority rate, threshold rate with only one reference value)
# ■ Input
# x: An object of class 'E1Ind'
# digits=3 : Rounding of decimals
# I: Integer, level with the higher mean value
# J: Integer, level with the lower mean value
# cr1: Reference value for threshold rate
##########################################################################
E1betw_level<-function(x, digits=3, I, J, cr1=F)
{
  G<-matrix(0, length(x$mu), 5)
  G[,1] <- x$muA[,I] - x$muA[,J];
  G[,2] <- G[,1] / x$sigmaE;
  G[,3] <- pnorm(x$muA[,I], x$muA[,J], x$sigmaE);
  G[,4] <- pnorm((G[,2] / sqrt(2)), 0.0, 1.0);
  lab<-c("Difference in Means", "Effect Size", "Measure of Nonoverlap", "Superiority Rate")
  co<-4
  if(is.numeric(cr1)){  co<-co+1;
  G[,5] <- pnorm(((G[,1]-cr1)/(sqrt(2)*x$sigmaE)), 0.0, 1.0);
  lab<-c(lab, paste("Threshold Rate(", cr1, ")", sep=""));}
  Gc<-cbind(
    apply(G[,1:co], 2, mean),
    apply(G[,1:co], 2, sd),
    t(apply(G[,1:co], 2, quantile, probs=c(0.025, 0.05, 0.5, 0.95, 0.975)))
  )
  colnames(Gc)<-c("EAP", "post.sd", "2.5%", "5%", "50%", "95%", "97.5%")
  rownames(Gc)<-lab
  round(Gc, digits)
}
##########################################################################
# Comparison of two levels (does not output superiority rate, threshold rate specified with a vector)
# ■ Input
# x: An object of class 'E1Ind'
# digits=3 : Rounding of decimals
# I: Integer, level with the higher mean value
# J: Integer, level with the lower mean value
# cr1: Vector specifying reference points for threshold rate
##########################################################################
E1between<-function(x, digits=3, I, J, cr1=0)
{
  G<-matrix(0, length(x$mu), 3)
  G[,1] <- x$muA[,I] - x$muA[,J];
  G[,2] <- G[,1] / x$sigmaE;
  G[,3] <- pnorm(x$muA[,I], x$muA[,J], x$sigmaE);
  lab<-c("Mean Difference", "Delta", "Measure of Nonoverlap")
  H<-matrix(0, length(x$mu), length(cr1))
  for (i in 1:length(cr1)){
    H[,i] <- pnorm(((G[,1]-cr1[i])/(sqrt(2)*x$sigmaE)), 0.0, 1.0);
    lab<-c(lab, paste("Threshold Rate(", cr1[i], ")", sep=""));}
  Gc<-matrix(0, 3+length(cr1), 7)
  for (i in 1:3)          {Gc[i,]<-gqcal(G[,i])}
  for (i in 1:length(cr1)){Gc[i+3,]<-gqcal(H[,i])}
  colnames(Gc)<-c("EAP", "post.sd", "2.5%", "5%", "50%", "95%", "97.5%")
  rownames(Gc)<-lab
  round(Gc, digits)
}
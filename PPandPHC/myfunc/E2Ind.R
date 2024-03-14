##########################################################################
# Analyzing data from an independent two-factor experiment
# ■ Input
# y: Characteristic values in vector form
# A: Levels in vector form (specified as integers from 1 to a, without gaps)
# B: Levels in vector form (specified as integers from 1 to b, without gaps)
# prior=F : Logical value. If T, specify the range of the prior distribution. If F, do not specify.
# mL=-1000, mH=1000, sL=0, sH=100 : Parameters of the prior distribution
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution
# see=1234, cha=5, war=1000, ite=21000 : Parameters related to MCMC
# ■ Output
# fit: output from stan, par: list of parameters, prob: probability vector, 
# mu: overall mean, muA: mean of A, muB: mean of B, muAB: interaction, 
# cellmean: cell means, sigmaA: standard deviation of A, sigmaB: standard deviation of B, sigmaAB: standard deviation of AB,
# sigmaE: standard deviation of E, eta2A: explanatory rate of A, eta2B: explanatory rate of B, eta2AB: explanatory rate of AB, eta2T: total explanatory rate,
# deltaA: effect size of A, deltaB: effect size of B, deltaAB: effect size of AB, UbigA[a]: probability of A being greater than or equal to 0, 
# UsmaA[a]: probability of A being less than or equal to 0, UbigB[b]: probability of B being greater than or equal to 0, UsmaB[b]: probability of B being less than or equal to 0,
# UbigAB[a,b]: probability of AB being greater than or equal to 0, UsmaAB[a,b]: probability of AB being less than or equal to 0, 
# U2A[a,a]: comparison of two levels of A, U2B[b,b]: comparison of two levels of B, log_lik: log likelihood
##########################################################################
E2Ind <- function(y, A, B, prior=F, mL=-1000, mH=1000, sL=0, sH=100,
                  prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                  see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
  library(cmdstanr)              # Load and attach add-on package cmdstanr
  library(posterior)             # Load and attach add-on package posterior
  library(rstan)
  if (prior) {
    dat <- list(n=length(y), a=length(unique(A)), b=length(unique(B)),
                y=y, A=A, B=B, mL=mL, mH=mH, sL=sL, sH=sH)
    scr <- "stan/E2IndPT.stan"    # Stan file name
  } else {
    dat <- list(n=length(y), a=length(unique(A)), b=length(unique(B)),
                y=y, A=A, B=B)
    scr <- "stan/E2IndPF.stan"    # Stan file name
  }
  par<-c("mu", "muA", "muB", "muAB", "cellmean", "sigmaA", "sigmaB", "sigmaAB",
         "sigmaE", "eta2A", "eta2B", "eta2AB", "eta2T", "deltaA", "deltaB", "deltaAB",
         "UbigA", "UsmaA", "UbigB", "UsmaB", "UbigAB", "UsmaAB", "U2A", "U2B", "log_lik")
  if (pack == "rstan") {
    fit<-stan(file=scr, data=dat, pars=par, seed=see, chains=cha, warmup=war, iter=ite)
  }else{
    mod <- cmdstan_model(scr)                   # Compile
    csrfit <- mod$sample(data=dat, chains=cha, iter_sampling=ite,
                         iter_warmup=war, parallel_chains=cha, seed=see)    # MCMC
    fit <- rstan::read_stan_csv(csrfit$output_files())  # Convert to stan format
  }
  ext<-extract(fit, par);
  out<-list(prior=prior, fit=fit, par=par, prob=prob,
            mu=ext$mu, muA=ext$muA, muB=ext$muB, muAB=ext$muAB, cellmean=ext$cellmean,
            sigmaA=ext$sigmaA, sigmaB=ext$sigmaB, sigmaAB=ext$sigmaAB, sigmaE=ext$sigmaE,
            eta2A=ext$eta2A, eta2B=ext$eta2B, eta2AB=ext$eta2AB, eta2T=ext$eta2T,
            deltaA=ext$deltaA, deltaB=ext$deltaB, deltaAB=ext$deltaAB, UbigA=ext$UbigA,
            UsmaA=ext$UsmaA, UbigB=ext$UbigB, UsmaB=ext$UsmaB, UbigAB=ext$UbigAB,
            UsmaAB=ext$UsmaAB, U2A=ext$U2A, U2B=ext$U2B, log_lik=ext$log_lik)
  class(out)<-'E2Ind'
  return(invisible(out))
}
##########################################################################
# Print method
# ■ Input
# x: An object of class 'E2Ind'
# digits=3 : Rounding of decimals
##########################################################################
print.E2Ind<-function(x, digits=3)
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
# Comparison of two cells
# ■ Input
# x: An object of class 'E2Ind'
# digits=3 : Rounding of decimals
# H="A": The factor to fix, either A or B
# F=1: Integer, the level of the factor to fix
# I=1: Integer, level with the higher mean value of the factor to compare
# J=2: Integer, level with the lower mean value of the factor to compare
# cr1: Reference value for threshold rate
##########################################################################
E2betw_level<-function(x, digits=3, H="A", F=1, I=1, J=2, cr1=F)
{
  if (H=='A') { big<-x$cellmean[,F,I]; sma<-x$cellmean[,F,J]; comparison<-"B";}
  else { big<-x$cellmean[,I,F]; sma<-x$cellmean[,J,F]; comparison<-"A";}
  G<-matrix(0, length(x$mu), 5)
  G[,1] <- big-sma;
  G[,2] <- G[,1]/x$sigmaE;
  G[,3] <- pnorm(big, sma, x$sigmaE);
  G[,4] <- pnorm((G[,2]/sqrt(2)), 0.0, 1.0);
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
  cat("Fixing level", F, "of factor", H, "and comparing levels", I, J, "of factor", comparison, "\n")
  round(Gc, digits)
}
##########################################################################

# Comparison of two cells (does not output superiority rate, threshold rate specified with a vector)
# x: An object of class 'E2Ind'
# H="A": The factor to fix, either A or B
# F=1: Integer, the level of the factor to fix
# I=1: Integer, level with the higher mean value of the factor to compare
# J=2: Integer, level with the lower mean value of the factor to compare
# cr1: Vector specifying reference points for threshold rate
# digits=3 : Rounding of decimals
##########################################################################
E2between<-function(x, H="A", F=1, I=1, J=2, cr1=0, digits=3)
{
  if (H=='A') { big<-x$cellmean[,F,I]; sma<-x$cellmean[,F,J]; comparison<-"B";}
  else { big<-x$cellmean[,I,F]; sma<-x$cellmean[,J,F]; comparison<-"A";}
  G<-matrix(0, length(x$mu), 3)
  G[,1] <- big-sma;
  G[,2] <- G[,1]/x$sigmaE;
  G[,3] <- pnorm(big, sma, x$sigmaE);
  lab<-c("Mean Difference", "Delta", "Measure of Nonoverlap")
  K<-matrix(0, length(x$mu), length(cr1))
  for (i in 1:length(cr1)){
    K[,i] <- pnorm(((G[,1]-cr1[i])/(sqrt(2)*x$sigmaE)), 0.0, 1.0);
    lab<-c(lab, paste("Threshold Rate(", cr1[i], ")", sep=""));}
  Gc<-matrix(0, 3+length(cr1), 7)
  for (i in 1:3)          {Gc[i,]<-gqcal(G[,i])}
  for (i in 1:length(cr1)){Gc[i+3,]<-gqcal(K[,i])}
  colnames(Gc)<-c("EAP", "post.sd", "2.5%", "5%", "50%", "95%", "97.5%")
  rownames(Gc)<-lab
  cat("Fixing level", F, "of factor", H, "and comparing levels", I, J, "of factor", comparison, "\n")
  round(Gc, digits)
}

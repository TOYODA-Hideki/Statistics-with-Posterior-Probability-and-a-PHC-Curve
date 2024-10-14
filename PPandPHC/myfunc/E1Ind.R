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
# Usma: effect of levels less than or equal to 0, log_lik: log likelihood
##########################################################################
E1Ind <- function(y, A, prior=F, mL=-1000, mH=1000, sL=0, sH=100,
                  prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
                  see=1234, cha=5, war=1000, ite=21000)
{
  library(cmdstanr)              # Load and attach add-on package cmdstanr
  library(posterior)             # Load and attach add-on package posterior
  if (prior) {
    dat <- list(n=length(y), a=max(A), y=y, A=A, mL=mL, mH=mH, sL=sL, sH=sH)
    scr <- "stan/E1IndPT.stan"    # Stan file name
  } else {
    dat <- list(n=length(y), a=max(A), y=y, A=A)
    scr <- "stan/E1IndPF.stan"    # Stan file name
  }
  par<-c("muA", "sigmaE", "sigmaA", "eta2", "delta",
         "mu", "aj", "Ubig", "Usma", "U2", "log_lik")  # Sampling variables
    mod <- cmdstan_model(scr)                   # Compile
    csrfit <- mod$sample(data=dat, chains=cha, iter_sampling=ite,
                  iter_warmup=war, parallel_chains=cha, seed=see)    # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
      a<-max(A)
   muA<-as.matrix(ext[,1:a]);  sigmaE<-as.matrix(ext$sigmaE);
   sigmaA<-as.matrix(ext$sigmaA); eta2<-as.matrix(ext$eta2);
   delta<-as.matrix(ext$delta); mu<-as.matrix(ext$mu); 
   aj<-as.matrix(ext[,c((a+6):(2*a+5))]); 
   Ubig<-as.matrix(ext[,c((2*a+6):(3*a+5))]); 
   Usma<-as.matrix(ext[,c((3*a+6):(4*a+5))]);
   log_lik<-as.matrix(ext$log_lik);
   out<-list(prior=prior,ext=ext,draw=draw,par=par,prob=prob,muA=muA,
            sigmaE=sigmaE, sigmaA=sigmaA, eta2=eta2, delta=delta, mu=mu, aj=aj,
            Ubig=Ubig, Usma=Usma, log_lik=log_lik)
  class(out)<-'E1Ind'
  return(invisible(out))
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
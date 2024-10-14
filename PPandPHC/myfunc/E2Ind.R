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
                  see=1234, cha=5, war=1000, ite=21000)
{
  library(cmdstanr)              # Load and attach add-on package cmdstanr
  library(posterior)             # Load and attach add-on package posterior
   a<-length(unique(A)); b<-length(unique(B));
   if (prior) {
     dat <- list(n=length(y),a=a,b=b,y=y,A=A,B=B,mL=mL,mH=mH,sL=sL,sH=sH)
     scr <- "stan/E2IndPT.stan" 
   } else {
     dat <- list(n=length(y),a=a,b=b,y=y,A=A,B=B)
     scr <- "stan/E2IndPF.stan" 
   }
  par<-c("mu", "muA", "muB", "muAB", "cellmean", "sigmaA", "sigmaB", "sigmaAB",
       "sigmaE","eta2A","eta2B","eta2AB","eta2T","deltaA","deltaB","deltaAB",
       "UbigA","UsmaA","UbigB","UsmaB","UbigAB","UsmaAB","U2A","U2B","log_lik")
    mod <- cmdstan_model(scr)                   # Compile
    csrfit <- mod$sample(data=dat, chains=cha, iter_sampling=ite,
             iter_warmup=war, parallel_chains=cha, seed=see)    # MCMC
      draw<-csrfit$draws(par)
      ext<-as_draws_df(draw)
      colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext))
   out<-list(prior=prior,ext=ext,draw=draw,par=par,prob=prob,
      mu=as.matrix(ext$mu),muA=as.matrix(ext[,c(2:(a+1))]),
      muB=as.matrix(ext[,c((a+2):(a+b+1))]),
      muAB=as.matrix(ext[,c((a+b+2):(a+b+a*b+1))]),
      cellmean=as.matrix(ext[,c((a+b+a*b+2):(a+b+2*a*b+1))]),
      sigmaA=as.matrix(ext$sigmaA), sigmaB=as.matrix(ext$sigmaB),
      sigmaAB=as.matrix(ext$sigmaAB), sigmaE=as.matrix(ext$sigmaE),
      eta2A=as.matrix(ext$eta2A), eta2B=as.matrix(ext$eta2B),
      eta2AB=as.matrix(ext$eta2AB) ,eta2T=as.matrix(ext$eta2T),
      log_lik=ext$log_lik)
  class(out)<-'E2Ind'
  return(invisible(out))
}
##########################################################################
# Comparison of two cells (does not output superiority rate, threshold rate specified with a vector)
# x: An object of class 'E2Ind'
# I=1: Integer, cell with the higher mean value
# J=2: Integer, cell with the lower mean value 
# cr1: Vector specifying reference points for threshold rate
# digits=3 : Rounding of decimals
##########################################################################
E2between<-function(x,I=2,J=3,cr1=0,digits=3)
{
  big<-x$cellmean[,I]; sma<-x$cellmean[,J];
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
  cat("Compare the cell in column", I, "with the cell in column", J, "\n")
  round(Gc, digits)
}

##########################################################################
#Estimating the difference between two matched groups.
#■Input
#x:Data in a matrix format with n rows and 2 columns.
#  (The group with the larger average value is designated in the first column)
#EQU=1 : Boolean value. 1 for equal variance. 0 for different variances.
#prior=F : Boolean value.  T to specify the range of the prior distribution. F not to specify.
#mL=-1000, mH=1000, sL=0, sH=100 : Parameters of the prior distribution
#prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to be reported in the posterior distribution
#see=1234,cha=5,war=1000,ite=21000 : Parameters related to MCMC
#■Output
# fit:Output of stan, par:list of parameters, prob:probability vector,
# mu1:Average of the first group, mu2:Average of the second group,
# sigma1:Standard deviation of the first group, sigma2:Standard deviation of the second group,, 
# rho:Correlation coefficient, xaste1:Predictive distribution 1, 
# xaste2:predictive distribution 2, log_lik:log likelihood
##########################################################################
G2pair <- function(x,EQU=1,prior=F,mL=-1000, mH=1000, sL=0, sH=100,
            prob=c(0.025, 0.05, 0.5, 0.95, 0.975),
            see=1234, cha=5, war=1000, ite=21000, pack="cmds")
{
   library(cmdstanr)              #Load and attach add-on package cmdstanr
   library(posterior)             #Load and attach add-on package posterior
   library(rstan)                 #Load and attach add-on package rstan
   if (prior) {
     dat <- list(n=nrow(x),x=x,EQU=EQU,mL=mL,mH=mH,sL=sL,sH=sH)
     scr <- "stan/G2pairPT.stan"                           # Stan's filename
   } else {
     dat <- list(n=nrow(x),x=x,EQU=EQU)
     scr <- "stan/G2pairPF.stan"                           # Stan's filename
   }
   par<-c("mu","sigma1","sigma2","rho","xaste","log_lik")# sampling variable
   if (pack == "rstan") {
      fit<-stan(file=scr,data=dat,pars=par,seed=see,chains=cha,warmup=war,iter=ite)
   }else{
      mod <- cmdstan_model(scr)                   #Compile
      csrfit <- mod$sample(data=dat,chains=cha,iter_sampling=ite,
         iter_warmup=war,parallel_chains=cha,seed=see)    #MCMC
      fit <- rstan::read_stan_csv(csrfit$output_files())#Conversion to stan format
   }
   ext<-extract(fit, par);    mu1<-ext$mu[,1]; mu2<-ext$mu[,2]; 
   sigma1<-ext$sigma1; sigma2<-ext$sigma2; rho<-ext$rho;
   xaste1<-ext$xaste[,1]; xaste2<-ext$xaste[,2]; log_lik<-ext$log_lik
   out<-list(EQU=EQU,prior=prior,fit=fit,par=par,prob=prob,mu1=mu1,mu2=mu2,
       sigma1=sigma1,sigma2=sigma2,rho=rho,
       xaste1=xaste1,xaste2=xaste2,log_lik=log_lik)
   class(out)<-'G2pair'
   return(invisible(out))
}
##########################################################################
#Printing Methods
#■Input
#x:Object of class 'G2pair
#onlydiff=T Only outputs related to difference scores, or not.
#degits=3 : Rounding decimals
#cr1=F : Reference point for the difference of means
#cr2=F : Reference point for probability beyond threshold  
#cr3=F : Reference point of effect size 
#cr4=F : Reference point for standard deviation of difference scores 
#ra= 1.0:Correlation upper bound
#rb=-1.0:Correlation lower bound
#pr1=F : Base probability of measure of nonoverlap
#pr2=F : Base probability of predominance ratio
#pr3=F : Base probability of probability beyond threshold 
#pr4=F : Base probability of homothetic ratio
##########################################################################
print.G2pair<-function(x,onlydiff=T,degits=3,cr1=F,cr2=F,cr3=F,cr4=F,ra= 1.0,rb= -1.0,pr1=F,pr2=F,pr3=F,pr4=F)
{
   EQU<-x$EQU; prior<-x$prior;
   if (EQU>0.5) {print("******************equidistribution model******************")}
           else {print("******************heteroscedasticity model******************")}
   if (prior) {print("******************Specify prior distribution******************")}
             else {print("******************Default prior distribution******************")}
   mu1<-x$mu1; mu2<-x$mu2; sigma1<-x$sigma1; sigma2<-x$sigma2; 
   rho<-x$rho; xaste1<-x$xaste1; xaste2<-x$xaste2; log_lik<-x$log_lik
   prob<-x$prob
   G<-matrix(0,length(mu1),12)
   G0<-sqrt(sigma1^2 + sigma2^2)
   G[,1] <- mu1 - mu2;                         #mean difference
   G[,2] <- G[,1]/sigma1;                      #Effect size1
   G[,3] <- G[,1]/sigma2;                      #Effect size2
   G[,4] <- pnorm(mu1,mu2,sigma2);             #measure of nonoverlap 1
   G[,5] <- 1-pnorm(mu2,mu1,sigma1);           #measure of nonoverlap 2
   G[,6] <- pnorm((G[,1]/G0), 0.0, 1.0);       #predominance ratio
   #Below is an analysis using rho
   G[,7] <- sqrt(sigma1^2 + sigma2^2 -2*rho*sigma1*sigma2);  #Difference score sd
   G[,8] <- G[,1]/G[,7];                       #Effect size of difference score
   G[,9] <- pnorm(G[,8], 0.0, 1.0);            #Predominance ratio of difference score
   G[,10] <- 0.5+asin(rho)/pi                  #Homothetic ratio
   co<-10; lab<-c("Mean Difference","Effect size1","Effect size2","Measure of nonoverlap group 1","Measure of nonoverlap group 2","Predominance ratio","Difference score sd","Effect size of difference score","predominance ratio of difference score","homothetic ratio")
   if(is.numeric(cr2)){
         co<- co+1; G[,co]<- pnorm((G[,1]-cr2)/G0, 0.0, 1.0);
         lab<-c(lab,paste("Probability beyond threshold (",cr2,")",sep=""))
         co<- co+1; G[,co]<- pnorm((G[,1]-cr2)/G[,7], 0.0, 1.0);
         lab<-c(lab,paste("Difference score probability beyond threshold (",cr2,")",sep=""))
   }
   if (onlydiff) {if (is.numeric(cr2)) {Grow<-c(1,7:10,12)} else {Grow<-c(1,7:10)}}
   else {if (is.numeric(cr2)) {Grow<-c(1:12)} else {Grow<-c(1:10)}};
   Gc<-cbind(
      apply(G[,Grow],2,mean),
      apply(G[,Grow],2,sd),
      t(apply(G[,Grow],2,quantile, probs=prob))
   )
   colnames(Gc)<-c("EAP","post.sd",prob)
   rownames(Gc)<-lab[Grow]

   U<-matrix(0,length(mu1),21)
   U[,1] <- ifelse(G[,1]>0.0,  1, 0);
   U[,2] <- ifelse(G[,1]<=0.0, 1, 0);
   U[,3] <- ifelse(xaste1-xaste2>0.0, 1, 0);
   U[,4] <- ifelse((rb<=rho)&(rho<=ra), 1, 0);
   co<-4;lab<-c("Probability that μ1-μ2 is greater than 0","Probability that μ1-μ2 is less than or equal to 0"
                ,"Predominance ratio of difference scores(direct comparison)",paste("Probability that the correlation is between",rb,"and",ra,sep=""))
   if(is.numeric(cr4)){co<- co+1;U[,co]<- ifelse(G[,7]< cr4,1, 0);
         lab<-c(lab,paste("Probability that sd of difference scores is less than(",cr4,")",sep=""))}
   if(is.numeric(cr2)){co<- co+1;U[,co]<- ifelse(xaste1-xaste2> cr2,1, 0);
         lab<-c(lab,paste("Probability of difference scores beyond threshold (",cr2,")(direct comparison)",sep=""))}
   if(is.numeric(cr1)){co<- co+1;U[,co]<- ifelse(G[,1]> cr1,    1, 0);
         lab<-c(lab,paste("Probability that μ1-μ2 is greater than (",cr1,")",sep=""))}
   if(is.numeric(cr3)){co<- co+1;U[,co]<- ifelse(G[,8]> cr3, 1, 0);
         lab<-c(lab,paste("Probability that Effect size of difference score is greater than (",cr3,")",sep=""));}
   if(is.numeric(pr2)){co<- co+1;U[,co]<- ifelse(G[,9]> pr2,    1, 0);
         lab<-c(lab,paste("Probability that the difference score's predominance ratio is greater than (",pr2,")",sep=""))}
   if((is.numeric(cr2))&(is.numeric(pr3)))
                      {co<- co+1;U[,co]<- ifelse(G[,12]> pr3,    1, 0);
         lab<-c(lab,paste("Probability beyond threshold (",cr2,") of difference score being greater than (",pr3,")",sep=""))}
   if(is.numeric(pr4)){co<- co+1;U[,co]<- ifelse(G[,10]> pr4,    1, 0);
         lab<-c(lab,paste("Probability that the same order rate is greater than (",pr4,")",sep=""))}
   if(!onlydiff & is.numeric(cr3)){co<- co+1;U[,co]<- ifelse(G[,2]> cr3, 1, 0);
         lab<-c(lab,paste("Effect Probability that size1 is greater than (",cr3,")",sep=""));
                       co<- co+1;U[,co]<- ifelse(G[,3]> cr3,    1, 0);
         lab<-c(lab,paste("Effect Probability that size2 is greater than (",cr3,")",sep=""))}
   if(!onlydiff & is.numeric(pr1)){co<- co+1;U[,co]<- ifelse(G[,4]> pr1, 1, 0);
         lab<-c(lab,paste("Probability that the measure of nonoverlap in the first group is greater than (",pr1,")",sep=""));
                       co<- co+1;U[,co]<- ifelse(G[,5]> pr1,    1, 0);
         lab<-c(lab,paste("Probability that the measure of nonoverlap in the second group is greater than (",pr1,")",sep=""))}
   if(!onlydiff & is.numeric(pr2)){co<- co+1;U[,co]<- ifelse(G[,6]> pr2, 1, 0);
         lab<-c(lab,paste("Probability that predominance ratio is greater than (",pr2,")",sep=""))}
   if(!onlydiff & (is.numeric(cr2))&(is.numeric(pr3)))
                      {co<- co+1;U[,co]<- ifelse(G[,11]> pr3,    1, 0);
         lab<-c(lab,paste("Probability beyond threshold (",cr2,") is greater than (",pr3,")",sep=""))}
   if(co>0.9){
      Uc<-matrix(colMeans(as.matrix(U[,1:co])),co,1);
      rownames(Uc)<-lab
      }else{Uc<-0.0}
   waic<- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
   print(x$fit,pars=x$par,digits_summary=degits,probs=x$prob)
   print(round(Gc,degits))
   print(round(Uc,degits))
   print(paste("waic=",round(waic,degits),sep=""))
   out<-list(G=G,Gc=Gc,Uc=Uc,waic=waic)
   return(invisible(out))
}  

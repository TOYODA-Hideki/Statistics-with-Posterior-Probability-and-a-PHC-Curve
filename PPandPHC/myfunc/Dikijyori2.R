# Dikijyori2(　) : Function to return a curve/table of probability beyond threshold of difference score
#seq01: Vector (or scalar) of reference point c
#out  : MCMC object containing a random vector of mu1,mu2,sigma1,sigma2 and rho
#If byoga="yes", the PROBABILITY BEYOND THRESHOLD curve,
# in other cases, the probability beyond threshold table is output.
#xlab="",ylab="",cex.lab=2.2,cex.axis=1.5,lwd=2.0
#Returns are table
Dikijyori2<-function(seq01,out,byoga="yes", 
      dedits=3,lwd=2,lty=1, xlab="" ,ylab="",cex.lab=2.0, cex.axis=2.0){
  mu1<-out$mu1;mu2<-out$mu2;sigma1<-out$sigma1;sigma2<-out$sigma2;rho<-out$rho
  noc<-length(seq01)
  Dikijyo<-matrix(0,length(mu1),noc)
  colnames(Dikijyo)<-paste("Dikijyo(",seq01,")",sep="")
  for (i in 1:noc){
    Dikijyo[,i]<-pnorm((mu1-mu2-seq01[i])/sqrt(sigma1^2+
        sigma2^2-2*rho*sigma1*sigma2),0,1)}
  DikijyoT<-gqcal(Dikijyo,dedits)
  if (byoga=="yes") {
    a<-m<-b<-numeric(noc)
    a<-DikijyoT[,7]
    m<-DikijyoT[,1]
    b<-DikijyoT[,3]
    plot(seq01,a,type="l",xlim=c(min(seq01),max(seq01)),ylim=c(0,1.0),
      cex.axis=cex.axis,cex.lab=cex.lab,main="",xlab="",ylab="");  par(new=T)
    plot(seq01,m,type="l",xlim=c(min(seq01),max(seq01)),ylim=c(0,1.0),lwd=lwd,
      cex.axis=cex.axis,cex.lab=cex.lab,main="",xlab="",ylab="");  par(new=T)
    plot(seq01,b,type="l",xlim=c(min(seq01),max(seq01)),ylim=c(0,1.0),
      cex.axis=cex.axis,cex.lab=cex.lab,main="",xlab="",ylab="")
    return(invisible(DikijyoT))
  }
  return(DikijyoT)
}  

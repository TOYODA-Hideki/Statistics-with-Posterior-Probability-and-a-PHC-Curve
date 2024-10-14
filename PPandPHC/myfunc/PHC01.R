# Drawing of a PHC curve indicating that a is larger (or smaller) than b by c
# seq01 A vector representing the horizontal axis (reference point c of the difference). For example, seq01=seq(0,1.2,0.05)
# a, b, Random number vectors approximating the posterior distribution, calculating PHC for a - b > c or a - b < c
# To handle a alone, set b=0.
# If cc="gtc", then a - b > c, if "rope" then abs(a-b) < c,
#              otherwise ("ltc" is assumed) a - b < c
# byoga="yes" means draw the PHC curve, otherwise, output the PHC
# xlab="", ylab="", cex.lab=2.2, cex.axis=1.5, lwd=1.5
# The return is the PHC vector
PHC01<-function(seq01, a, b=0
                , cc="gtc", byoga="yes", dedits=3, lwd=2, lty=1,
                xlab="" , ylab="", cex.lab=2.2, cex.axis=1.5){
  av<-numeric(length(seq01))
  j<-0
  if (cc == "gtc") {
    for(c in seq01){ j<-j+1; av[j]<-mean(a-b>c);}
  }else if(cc == "rope"){
    for(c in seq01){ j<-j+1; av[j]<-mean(abs(a-b)<c);}
  }else{
    for(c in seq01){ j<-j+1; av[j]<-mean(a-b<c);}
  }
  names(av)<-seq01
  if (byoga=="yes") {
    plot(seq01, av, lwd=lwd, lty=lty, type="l", ylab="", ylim=c(0,1),
         xlab=xlab, xlim=range(seq01), cex.lab=cex.lab, cex.axis=cex.axis)
    grid(lwd=1.5)
    return(invisible(round(av, dedits)))
  }
  return(round(av, dedits))
}

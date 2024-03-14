#z：乱数
#幅：階級幅
con_hist<-function(z,幅=0.5,横軸="",digits=3,上空=300,横空=100)
{
   q0<-min(z);q1<-quantile(z,0.025);q2<-quantile(z,0.05);q3<-quantile(z,0.95)
   q4<-quantile(z,0.975);q5<-max(z);e<-(q5-q0)/横空
   ymax<-max(hist(z,breaks=seq(q0-e,q5+e,幅))$counts)
   hist(z,breaks=seq(q0-e,q5+e,幅),xlab=横軸,cex.lab=1.5,cex.axis=1.5,main="",ylab="",ylim=c(0,ymax*1.1))
   eap<-mean(z)
   segments(eap,0,eap,ymax,lwd=2.0);text(eap,ymax+上空,paste("EAP",round(eap,digits)),cex=1.5)
   par(new=T)
   hist(z[((z<q2)|(q3<z))],breaks=seq(q0-e,q5+e,幅),xlab='',cex.lab=1.5,cex.axis=1.5,main="",ylab="",ylim=c(0,ymax*1.1),col = gray(0.8))
   par(new=T)
   hist(z[((z<q1)|(q4<z))],breaks=seq(q0-e,q5+e,幅),xlab='',cex.lab=1.5,cex.axis=1.5,main="",ylab="",ylim=c(0,ymax*1.1),col = gray(0.4))
}

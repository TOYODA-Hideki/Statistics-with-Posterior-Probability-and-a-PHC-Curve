#Calculate the summary statistics of the generated quantities
#x is a vector of generated quantities (random number × 1) or
#a matrix of (random numbers × generated quantities)
gqcal<-function(x,digits =3,probs=c(0.025,0.05,0.5,0.95,0.975) ){
   xx<-as.matrix(x)
   y1<-rbind(apply(xx,2,mean),apply(xx,2,sd),apply(xx,2,quantile,probs))
   y<-round(t(y1),digits)
   rownames(y)<-colnames(xx)
   colnames(y)<-c("EAP","post.sd",probs)
   return(y)
}

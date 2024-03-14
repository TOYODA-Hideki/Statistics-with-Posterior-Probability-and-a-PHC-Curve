#Correlation between test and items
disc.co<-function(x){
  m<-ncol(x)
  test.sco <-  apply(x, 1,'sum')
  test.sco.a<-matrix(test.sco,nrow(x),ncol(x),byrow=F)
  test.sco.b<-test.sco.a - x
  dis<-c(1:m)
  for (i in 1:m){
    dis[i]<-cor(test.sco.b[,i], x[,i])
  }
  names(dis)<-colnames(x)
  return(dis)
}

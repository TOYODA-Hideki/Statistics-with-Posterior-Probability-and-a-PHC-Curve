#alpha coefficient
###Arguments x (Subject, item) matrix
alpha.co <- function(x){
   m <- ncol(x)
   alpha.co<-(m/(m-1))*(1-(sum(apply(x, 2,var))/var(apply(x, 1, sum))))
   return(round(alpha.co,3))
}

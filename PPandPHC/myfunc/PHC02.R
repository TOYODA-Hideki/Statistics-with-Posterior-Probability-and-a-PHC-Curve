# Matrix of probabilities that the statement "there is a difference greater than c between multiple levels" is correct
# c   Scalar, reference point
# ext A matrix with random numbers in rows and parameters in columns
# If cc="gtc", then row parameter - column parameter > c,
#  If "rope", then abs(row parameter - column parameter) < c,
#      otherwise, row parameter - column parameter < c
PHC02<-function(c=0, ext, cc="gtc", digits=3){
  J<-ncol(ext)
  pro_matrix<-matrix(0,J,J)
  colnames(pro_matrix)<-colnames(ext)
  rownames(pro_matrix)<-colnames(ext)
  if (cc == "gtc") {
    for(i in 1:J){for(j in 1:J){pro_matrix[i,j]<-mean(ext[,i]-ext[,j]>c)}}
  }else if(cc == "rope"){
    for(i in 1:J){for(j in 1:J){pro_matrix[i,j]<-mean(abs(ext[,i]-ext[,j])<c)}}
  }else{
    for(i in 1:J){for(j in 1:J){pro_matrix[i,j]<-mean(ext[,i]-ext[,j]<c)}}
  }
  round(pro_matrix,digits)
}

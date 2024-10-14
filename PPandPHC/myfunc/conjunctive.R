# Probability that two research hypotheses are simultaneously true under the conjunctive criteria of reference point vectors c1, c2
# x: Matrix constituted by random numbers in column a, not (0,1) but random numbers
# IJ: A 2*2 matrix, evaluates the probability that IJ[1,1]-IJ[1,2] and IJ[2,1]-IJ[2,2] both hold
# c1, c2: Vectors of each reference point
conjunctive<-function(x, IJ, c1=0, c2=0, dedits=3){
  phcT<-matrix(0, length(c1), length(c2))
  ii<-0
  for (i in c1){
    ii<-ii+1; jj<-0;
    for (j in c2){
      jj<-jj+1
      phcT[ii,jj]<-mean((x[,IJ[1,1]]-x[,IJ[1,2]]>i)*(x[,IJ[2,1]]-x[,IJ[2,2]]>j))
    }
  }
  colnames(phcT)<-c1
  rownames(phcT)<-c2
  return(round(phcT, dedits))
}

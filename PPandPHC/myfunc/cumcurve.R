cumcurve<-function(x,ラベル='')
{
   plot(as.numeric(rownames(table(x))),
      cumsum(table(x))/length(x),type='l',
      ylab='確率', xlab=ラベル)
}
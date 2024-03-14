# mean_sd_plot( ) Plotting cell mean, standard deviation, and mean for a 2-factor experiment
# y: Measurement value vector
# A: Level specification (more than 2 levels)
# B: Level specification (2 levels)
mean_sd_plot<-function(y, A, B, ylab="", xlab="",
                       cex.lab=1.5, cex.axis=1.5, lwd=2.0, ylim=c(ymin, ymax)){
  cat("Cell Mean\n")
  print(meany<-round(tapply(y, list(A, B), mean), 2))
  cat("Cell Standard Deviation\n")
  print(sdny <-round(tapply(y, list(A, B), sdn), 2))   # Standard deviation with denominator n
  ymax<-max(meany) + 0.5*max(sdny)
  ymin<-min(meany) - 0.5*max(sdny)
  plot(rownames(meany), meany[,1], type="b", xlim=c(0, max(A)+1), ylim=ylim,
       ylab="", xlab="", cex.lab=cex.lab, cex.axis=cex.axis, lwd=lwd)
  par(new=T)
  plot(rownames(meany), meany[,2], type="b", xlim=c(0, max(A)+1), ylim=ylim,
       ylab=ylab, xlab=xlab, cex.lab=cex.lab, cex.axis=cex.axis, lwd=lwd)
}

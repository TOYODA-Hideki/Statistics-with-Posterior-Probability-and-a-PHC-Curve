##########################################################################
# Multiple Regression Model
# ■ Input
# y: Reference variable (n)
# X: Predictor variable matrix (n*p)
# see=1234, cha=5, war=1000, ite=21000: MCMC related parameters

# ■ Output
# ext: Stan output, par: parameter list,
# a: Intercept
# b: Regression coefficient vector (p)
# sigma: Error standard deviation
# r2: Coefficient of determination
# r: Multiple correlation coefficient
# sb: Standardized regression coefficients
# vyhat: Variance of yhat
# log_lik: Log likelihood
##########################################################################
Reg <- function(y, X, see=1234, cha=5, war=1000, ite=20000)
{
  library(cmdstanr)              # Load and attach add-on package cmdstanr
  library(posterior)             # Load and attach add-on package posterior
  y <- as.vector(y); X <- as.matrix(X);
  n <- length(y); p <- ncol(X)
  X0 <- cbind(X, 1); bi <- as.vector(solve(t(X0) %*% X0) %*% t(X0) %*% y);
  sigmai <- as.vector(sqrt((sum(y^2) - t(bi) %*% t(X0) %*% y) / n));
  initi <- function() { list(bb = bi, sigma = sigmai) }
  dat <- list(n = n, p = p, y = y, X = X0)
  scr <- "stan/Reg.stan"                             # Stan file name
  par <- c("a", "b", "sigma", "r2", "r", "log_lik")  # Sampling variables
    mod <- cmdstan_model(scr)                       # Compile
    csrfit <- mod$sample(data = dat, chains = cha, iter_sampling = ite,
       init=initi,iter_warmup= war,parallel_chains=cha, seed = see)  # MCMC
    ext<-as_draws_df(csrfit$draws(par))
    (colnames(ext) <- gsub("\\[|\\]|,", "", colnames(ext)))
  ifelse(p==1, b<-as.matrix(ext[,2]), b<-as.matrix(ext[,c(2,2+p-1)]));
  vari <- function(xx) { mean((xx - mean(xx))^2) };
  sdX <- sqrt(apply(X, 2, vari))
  sb <- matrix(0, nrow(b), ncol(b))
  for (i in 1:ncol(b)) { sb[, i] <- b[, i] * sdX[i] / sqrt(vari(y)); }
  out<-list(ext=ext,par=par,a=ext$a, b=b, sigma = ext$sigma, r2 = ext$r2,
              r = ext$r, sb = sb, X = X, y = y, log_lik = ext$log_lik);
  class(out) <- 'Reg'
  return(invisible(out))
}
##########################################################################
# Print method
# ■ Input
# x: Object of class 'Reg'
# degits=3: Decimal rounding
# prob=c(0.025, 0.05, 0.5, 0.95, 0.975): Probability points to report in the posterior distribution,
# Xnew=F: Matrix of size (n*P), vector of size n for simple regression, Data for predictive distribution
# ■ Output
# waic=waic, y=y, X=X
# yhat=Predicted values, resi=Residuals (based on X)
# Gc: Standardized regression coefficients
# yhata, yhatc: Posterior distribution of regression equation (based on Xnew)
# yasta, yastc: Predictive distribution (based on Xnew)
##########################################################################
print.Reg <- function(x, degits=3, prob=c(0.025, 0.05, 0.5, 0.95, 0.975), Xnew=F)
{
  Gc=NULL; yhata=NULL; yhatc=NULL; yasta=NULL; yastc=NULL;
  a <- x$a; b <- x$b; sigma <- x$sigma; r2 <- x$r2; r <- x$r; sb <- x$sb; X <- x$X;
  y <- x$y; log_lik <- x$log_lik;
  waic <- (-2)*(log(mean(exp(log_lik)))) + 2*(var(log_lik))
#  print(x$ext, pars=x$par, digits_summary=degits, probs=prob)
  cat("Standardized regression coefficients\n")
  Gc <- cbind(apply(sb, 2, mean), apply(sb, 2, sd),
              t(apply(sb, 2, quantile, probs=prob)))
  colnames(Gc) <- c("EAP", "post.sd", prob)
  rownames(Gc) <- paste("b", 1:ncol(b), sep="")
  print(round(Gc, degits))
  # Predicted values and residuals based on X
  yhat <- apply((b %*% t(X)) + a %*% matrix(1, 1, nrow(X)), 2, mean);
  names(yhat) <- rownames(X)
  resi <- y - yhat;
  # Commented out: Print predicted values and residuals
  # cat("Predicted values\n")
  # print(round(yhat, degits))
  # cat("Residuals\n")
  # print(round(resi, degits))
  cat("Information criterion\n")
  print(paste("waic=", round(waic, degits), sep=""))
  # Predicted distribution based on Xnew
  if (!(is.logical(Xnew))) {
    if (length(dim(Xnew)) < 2) { Xnew <- matrix(Xnew, , 1) }
    cat("Posterior distribution of regression formula\n")
    yhata <- matrix(0, nrow(b), nrow(Xnew))
    yhata <- (b %*% t(Xnew)) + a %*% matrix(1, 1, nrow(Xnew));
    yhatc <- cbind(apply(yhata, 2, mean), apply(yhata, 2, sd),
                   t(apply(yhata, 2, quantile, probs=prob)))
    colnames(yhatc) <- c("EAP", "post.sd", prob)
    rownames(yhatc) <- rownames(Xnew)
    print(round(yhatc, degits))
    cat("Predictive distribution for reference variable\n")
    yasta <- matrix(0, nrow(b), nrow(Xnew))
    for (j in 1:nrow(Xnew)) {
      yasta[, j] <- rnorm(nrow(b), yhata[, j], sigma)
    }
    yastc <- cbind(apply(yasta, 2, mean), apply(yasta, 2, sd),
                   t(apply(yasta, 2, quantile, probs=prob)))
    colnames(yastc) <- c("EAP", "sd", prob)
    rownames(yastc) <- rownames(Xnew)
    print(round(yastc, degits))
  }
  out <- list(waic=waic, y=y, X=X, yhat=yhat, resi=resi, sb=sb, Gc=Gc, yhata=yhata,
              yhatc=yhatc, yasta=yasta, yastc=yastc)
  return(invisible(out))
}

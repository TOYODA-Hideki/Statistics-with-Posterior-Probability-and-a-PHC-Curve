data {
  int<lower=0> n; //Number of data   
  array[n] vector[2] x; //Data     
  real EQU; //SDÅA1 - CommonÅA0 - Different
  real mL;
  real mH;
  real sL;
  real sH; //Prior distribution 
}
parameters {
  vector<lower=mL, upper=mH>[2] mu; //Mean(Range specification)
  real<lower=sL, upper=sH> sigma1; //Standard deviation (Range specification)
  real<lower=sL, upper=sH> dummy; //SD for dummies(Range specification)
  real<lower=-1, upper=1> rho; //Correlation
}
transformed parameters {
  real<lower=0> sigma2; //Standard deviation 2
  cov_matrix[2] Sigma;
  sigma2 = EQU > 0.5 ? sigma1 : dummy;
  Sigma[1, 1] = pow(sigma1, 2); //Covariance matrix 
  Sigma[2, 2] = pow(sigma2, 2);
  Sigma[1, 2] = sigma1 * sigma2 * rho;
  Sigma[2, 1] = Sigma[1, 2];
}
model {
  for (i in 1 : n) {
    x[i] ~ multi_normal(mu, Sigma);
  } //Bivariate normal distribution
}
generated quantities {
  vector[2] xaste;
  real log_lik;
  xaste = multi_normal_rng(mu, Sigma); //Predictive distribution 
  log_lik = multi_normal_lpdf(x | mu, Sigma); //Log-likelihood
}

data {
  int<lower=0> n1;
  int<lower=0> n2; //The number of data
  array[n1] real x1;
  array[n2] real x2; //data
  real EQU; //sd, 1 - equal, 0 - not equal
  real mL;
  real mH;
  real sL;
  real sH; //prior distribution
}
parameters {
  vector<lower=mL, upper=mH>[2] mu; //mean
  real<lower=sL, upper=sH> sigma1; //standard deviation
  real<lower=sL, upper=sH> dummy; //standard deviation of dummy
}
transformed parameters {
  real<lower=0> sigma2; //standard deviation 2
  sigma2 = EQU > 0.5 ? sigma1 : dummy;
}
model {
  x1 ~ normal(mu[1], sigma1); //normal distribution 1
  x2 ~ normal(mu[2], sigma2); //normal distribution 2
}
generated quantities {
  array[2] real xaste;
  real log_lik;
  xaste[1] = normal_rng(mu[1], sigma1); //predictive distribution 1
  xaste[2] = normal_rng(mu[2], sigma2); //predictive distribution 2
  log_lik = normal_lpdf(x1 | mu[1], sigma1) + 
            normal_lpdf(x2 | mu[2], sigma2); //Log likelihood
}

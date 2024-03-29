//ÎĖČĒgQĖäĶÉÖ·évIŠ
data {
  int<lower=0> g; //QĖ 
  array[g] int<lower=0> x; //ģ― 
  array[g] int<lower=0> n; //f[^ 
}
parameters {
  array[g] real<lower=0, upper=1> p; //ęäĶ
}
transformed parameters {
  
}
model {
  for (i in 1 : g) {
    x[i] ~ binomial(n[i], p[i]);
  }
}
generated quantities {
  array[g] int<lower=0> xaste;
  real log_lik;
  array[g, g] real<lower=0, upper=1> U2; //2är
  for (i in 1 : g) {
    U2[i, i] = 0;
  }
  for (i in 1 : (g - 1)) {
    for (j in (i + 1) : g) {
      U2[i, j] = p[i] - p[j] > 0 ? 1 : 0;
      U2[j, i] = !(U2[i, j] != 0.0);
    }
  }
  log_lik = 0.0;
  for (i in 1 : g) {
    xaste[i] = binomial_rng(n[i], p[i]); //\ŠŠz
    log_lik = log_lik + binomial_lpmf(x[i] | n[i], p[i]);
  } //ÎmĶ
}

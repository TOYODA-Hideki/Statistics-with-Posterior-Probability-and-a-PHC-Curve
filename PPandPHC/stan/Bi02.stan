//対応のない2群の比率に関する統計的推測
data {
  array[2] int<lower=0> x; //正反応数 
  array[2] int<lower=0> n; //データ数 
}
parameters {
  array[2] real<lower=0, upper=1> p; //母比率
}
transformed parameters {
  
}
model {
  for (i in 1 : 2) {
    x[i] ~ binomial(n[i], p[i]);
  }
}
generated quantities {
  array[2] int<lower=0> xaste;
  real p_sa;
  real p_hi;
  real Odds_hi;
  array[2] real Odds;
  real log_lik;
  p_sa = p[1] - p[2];
  p_hi = p[1] / p[2];
  Odds[1] = p[1] / (1 - p[1]);
  Odds[2] = p[2] / (1 - p[2]);
  Odds_hi = Odds[1] / Odds[2];
  log_lik = 0.0;
  for (i in 1 : 2) {
    xaste[i] = binomial_rng(n[i], p[i]); //予測分布
    log_lik = log_lik + binomial_lpmf(x[i] | n[i], p[i]);
  } //対数確率
}

//1群の比率に関する統計的推測
data {
  int<lower=0> x; //正反応数 
  int<lower=0> n; //データ数 
}
parameters {
  real<lower=0, upper=1> theta; //母比率
}
transformed parameters {
  
}
model {
  x ~ binomial(n, theta);
}
generated quantities {
  int<lower=0> xaste;
  real<lower=0> Odds;
  real log_lik;
  xaste = binomial_rng(n, theta); //予測分布
  Odds = theta / (1 - theta); //オッズ
  log_lik = binomial_lpmf(x | n, theta); //対数確率
}


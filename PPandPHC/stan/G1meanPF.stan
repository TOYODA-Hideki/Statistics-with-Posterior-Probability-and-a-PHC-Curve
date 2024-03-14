//1群の正規分布の分析
data {
  int<lower=0> n; //データ数
  array[n] real<lower=0> x; //データ
}
parameters {
  real mu; //平均(範囲指定)
  real<lower=0> sigma; //標準偏差(範囲指定)
}
transformed parameters {
  
}
model {
  x ~ normal(mu, sigma); //正規分布
}
generated quantities {
  real xaste;
  real log_lik;
  xaste = normal_rng(mu, sigma); //予測分布
  log_lik = normal_lpdf(x | mu, sigma); //対数確率
}

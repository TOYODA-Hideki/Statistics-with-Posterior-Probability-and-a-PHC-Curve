//重回帰モデル
data {
  int<lower=0> n; //データ数 
  int<lower=0> p; //予測変数数   
  vector[n] y; //基準変数   
  matrix[n, p] X; //予測変数行列   
}
parameters {
  real a; //切片
  vector[p] b; //偏回帰係数ベクトル  
  real<lower=0> sigma; //誤差標準偏差
}
transformed parameters {
  vector[n] yhat; //yhat  
  yhat = a + X * b;
}
model {
  y ~ normal(yhat, sigma); //正規分布モデル
}
generated quantities {
  real log_lik; //対数尤度
  real<lower=0> vyhat; //yhatの分散
  real<lower=0, upper=1> r2; //決定係数
  vyhat = (variance(yhat) * (n - 1)) / n;
  r2 = vyhat / (vyhat + sigma ^ 2);
  log_lik = normal_lpdf(y | a + X * b, sigma);
}

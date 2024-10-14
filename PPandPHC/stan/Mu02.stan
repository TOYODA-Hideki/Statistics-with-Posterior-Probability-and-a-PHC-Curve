//対応のあるa×bのクロス表に関する統計的推測
data {
  int<lower=0> a; //行数 
  int<lower=0> b; //列数 
  array[a, b] int<lower=0> x; //反応数の行列形式 
}
transformed data {
  int<lower=0> N; //合計反応数 
  int<lower=0> ab; //セル合計数 
  array[a * b] int<lower=0> xv; //反応数のベクトル形式
  ab = a * b;
  for (i in 1 : a) {
    for (j in 1 : b) {
      xv[(i - 1) * b + j] = x[i, j];
    }
  }
  N = sum(xv);
}
parameters {
  simplex[ab] pi; //和が1の母比率のベクトル形式
}
transformed parameters {
  array[a, b] real<lower=0, upper=1> pim; //和が1の母比率の行列形式
  for (i in 1 : a) {
    for (j in 1 : b) {
      pim[i, j] = pi[(i - 1) * b + j];
    }
  }
}
model {
  xv ~ multinomial(pi); //多項分布
}
generated quantities {
  array[ab] int<lower=0> xastev; //予測分布ベクトル形式
  array[a, b] int<lower=0> xaste; //予測分布行列形式
  real log_lik;
  array[a] real pa;
  array[b] real pb;
  real V;
  array[a, b] real res;
  array[a, b] real<lower=0, upper=1> Up; //ピアソン残差が＋
  array[a, b] real<lower=0, upper=1> Um; //ピアソン残差がー
  array[a, b] real<lower=0> L; //同時確率と周辺確率の積との比
  xastev = multinomial_rng(pi, N); //予測分布
  V = 0;
  for (i in 1 : a) {
    pa[i] = 0;
  } //周辺確率
  for (j in 1 : b) {
    pb[j] = 0;
  }
  for (i in 1 : a) {
    for (j in 1 : b) {
      pa[i] = pa[i] + pim[i, j];
      pb[j] = pb[j] + pim[i, j];
    }
  }
  V = 0;
  for (i in 1 : a) {
    for (j in 1 : b) {
      xaste[i, j] = xastev[(i - 1) * b + j];
      res[i, j] = (pim[i, j] - pa[i] * pb[j]) / sqrt(pa[i] * pb[j]);
      L[i, j] = pim[i, j] / (pa[i] * pb[j]);
      V = V + pow(res[i, j], 2);
      Up[i, j] = res[i, j] > 0 ? 1 : 0;
      Um[i, j] = !(Up[i, j] != 0.0);
    }
  }
  V = sqrt(V / (min(a, b) - 1));
  log_lik = multinomial_lpmf(xv | pi); //対数確率
}

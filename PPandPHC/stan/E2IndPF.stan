//2要因対応なし アンバランスでも計算可能
data {
  int<lower=1> n; //データ数 
  int<lower=2> a; //A水準数  
  int<lower=2> b; //B水準数  
  vector[n] y; //特性値   
  array[n] int<lower=1> A; //A水準   
  array[n] int<lower=1> B; //B水準   
}
parameters {
  real mu; //全平均
  vector[a - 1] m1A; //A平均1つ少ない
  vector[b - 1] m1B; //B平均1つ少ない
  matrix[a - 1, b - 1] m1AB; //交互作用1つ少ない
  real<lower=0> sigmaE; //E標準偏差
}
transformed parameters {
  vector[a] muA; //A平均
  vector[b] muB; //B平均
  matrix[a, b] muAB; //交互作用
  vector[a - 1] m1a; //途中計算用
  vector[b - 1] m1b; //途中計算用
  for (i in 1 : (a - 1)) {
    muA[i] = m1A[i];
    m1a[i] = 0.0;
  }
  muA[a] = -sum(m1A);
  for (j in 1 : (b - 1)) {
    muB[j] = m1B[j];
    m1b[j] = 0.0;
  }
  muB[b] = -sum(m1B);
  for (i in 1 : (a - 1)) {
    for (j in 1 : (b - 1)) {
      muAB[i, j] = m1AB[i, j];
      m1a[i] = m1a[i] + m1AB[i, j];
    }
  }
  for (j in 1 : (b - 1)) {
    for (i in 1 : (a - 1)) {
      m1b[j] = m1b[j] + m1AB[i, j];
    }
  }
  for (i in 1 : (a - 1)) {
    muAB[i, b] = (-1) * m1a[i];
  }
  for (j in 1 : (b - 1)) {
    muAB[a, j] = (-1) * m1b[j];
  }
  muAB[a, b] = sum(m1a); //ここはマイナスではない
}
model {
  for (i in 1 : n) {
    y[i] ~ normal(mu + muA[A[i]] + muB[B[i]] + muAB[A[i], B[i]], sigmaE);
  }
}
generated quantities {
  real<lower=0> sigmaA; //A標準偏差
  real<lower=0> sigmaB; //B標準偏差
  real<lower=0> sigmaAB; //AB標準偏差
  matrix[a, b] cellmean; //セル平均値
  real log_lik; //対数尤度
  real<lower=0, upper=1> eta2A; //A説明率
  real<lower=0, upper=1> eta2B; //B説明率
  real<lower=0, upper=1> eta2AB; //AB説明率
  real<lower=0, upper=1> eta2T; //全説明率
  real<lower=0> ABAB; //中間変数A+B+AB分散
  real<lower=0> deltaA; //A効果量
  real<lower=0> deltaB; //B効果量
  real<lower=0> deltaAB; //AB効果量
  array[a] real<lower=0, upper=1> UbigA; //A論理値0以上確率
  array[a] real<lower=0, upper=1> UsmaA; //A論理値0以下確率
  array[b] real<lower=0, upper=1> UbigB; //B論理値0以上確率
  array[b] real<lower=0, upper=1> UsmaB; //B論理値0以下確率
  array[a, b] real<lower=0, upper=1> UbigAB; //AB論理値0以上確率
  array[a, b] real<lower=0, upper=1> UsmaAB; //AB論理値0以下確率
  array[a, a] real<lower=0, upper=1> U2A; //Aの2水準比較
  array[b, b] real<lower=0, upper=1> U2B; //Bの2水準比較
  sigmaA = sqrt(variance(muA) * (a - 1) / a);
  sigmaB = sqrt(variance(muB) * (b - 1) / b);
  sigmaAB = sqrt(variance(muAB) * (a * b - 1) / (a * b));
  ABAB = pow(sigmaA, 2) + pow(sigmaB, 2) + pow(sigmaAB, 2);
  eta2A = pow(sigmaA, 2) / (ABAB + pow(sigmaE, 2));
  eta2B = pow(sigmaB, 2) / (ABAB + pow(sigmaE, 2));
  eta2AB = pow(sigmaAB, 2) / (ABAB + pow(sigmaE, 2));
  eta2T = ABAB / (ABAB + pow(sigmaE, 2));
  deltaA = sigmaA / sigmaE;
  deltaB = sigmaB / sigmaE;
  deltaAB = sigmaAB / sigmaE;
  for (i in 1 : a) {
    for (j in 1 : b) {
      cellmean[i, j] = mu + muA[i] + muB[j] + muAB[i, j];
    }
  }
  for (i in 1 : a) {
    UbigA[i] = muA[i] > 0 ? 1 : 0;
    UsmaA[i] = muA[i] < 0 ? 1 : 0;
    U2A[i, i] = 0;
  }
  for (i in 1 : b) {
    UbigB[i] = muB[i] > 0 ? 1 : 0;
    UsmaB[i] = muB[i] < 0 ? 1 : 0;
    U2B[i, i] = 0;
  }
  for (i in 1 : a) {
    for (j in 1 : b) {
      UbigAB[i, j] = muAB[i, j] > 0 ? 1 : 0;
      UsmaAB[i, j] = muAB[i, j] < 0 ? 1 : 0;
    }
  }
  for (i in 1 : (a - 1)) {
    for (j in (i + 1) : a) {
      U2A[i, j] = muA[i] - muA[j] > 0 ? 1 : 0;
      U2A[j, i] = !(U2A[i, j] != 0.0);
    }
  }
  for (i in 1 : (b - 1)) {
    for (j in (i + 1) : b) {
      U2B[i, j] = muB[i] - muB[j] > 0 ? 1 : 0;
      U2B[j, i] = !(U2B[i, j] != 0.0);
    }
  }
  log_lik = 0.0;
  for (i in 1 : n) {
    log_lik = log_lik
              + normal_lpdf(y[i] | mu + muA[A[i]] + muB[B[i]]
                                   + muAB[A[i], B[i]], sigmaE);
  }
}



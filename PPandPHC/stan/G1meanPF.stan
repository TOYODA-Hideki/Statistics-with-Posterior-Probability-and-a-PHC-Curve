//1�Q�̐��K���z�̕���
data {
  int<lower=0> n; //�f�[�^��
  array[n] real<lower=0> x; //�f�[�^
}
parameters {
  real mu; //����(�͈͎w��)
  real<lower=0> sigma; //�W���΍�(�͈͎w��)
}
transformed parameters {
  
}
model {
  x ~ normal(mu, sigma); //���K���z
}
generated quantities {
  real xaste;
  real log_lik;
  xaste = normal_rng(mu, sigma); //�\�����z
  log_lik = normal_lpdf(x | mu, sigma); //�ΐ��m��
}

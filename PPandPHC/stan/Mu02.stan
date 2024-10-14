//�Ή��̂���a�~b�̃N���X�\�Ɋւ��铝�v�I����
data {
  int<lower=0> a; //�s�� 
  int<lower=0> b; //�� 
  array[a, b] int<lower=0> x; //�������̍s��`�� 
}
transformed data {
  int<lower=0> N; //���v������ 
  int<lower=0> ab; //�Z�����v�� 
  array[a * b] int<lower=0> xv; //�������̃x�N�g���`��
  ab = a * b;
  for (i in 1 : a) {
    for (j in 1 : b) {
      xv[(i - 1) * b + j] = x[i, j];
    }
  }
  N = sum(xv);
}
parameters {
  simplex[ab] pi; //�a��1�̕�䗦�̃x�N�g���`��
}
transformed parameters {
  array[a, b] real<lower=0, upper=1> pim; //�a��1�̕�䗦�̍s��`��
  for (i in 1 : a) {
    for (j in 1 : b) {
      pim[i, j] = pi[(i - 1) * b + j];
    }
  }
}
model {
  xv ~ multinomial(pi); //�������z
}
generated quantities {
  array[ab] int<lower=0> xastev; //�\�����z�x�N�g���`��
  array[a, b] int<lower=0> xaste; //�\�����z�s��`��
  real log_lik;
  array[a] real pa;
  array[b] real pb;
  real V;
  array[a, b] real res;
  array[a, b] real<lower=0, upper=1> Up; //�s�A�\���c�����{
  array[a, b] real<lower=0, upper=1> Um; //�s�A�\���c�����[
  array[a, b] real<lower=0> L; //�����m���Ǝ��ӊm���̐ςƂ̔�
  xastev = multinomial_rng(pi, N); //�\�����z
  V = 0;
  for (i in 1 : a) {
    pa[i] = 0;
  } //���ӊm��
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
  log_lik = multinomial_lpmf(xv | pi); //�ΐ��m��
}

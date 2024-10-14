//‘Î‰‚Ì‚È‚¢gŒQ‚Ì”ä—¦‚ÉŠÖ‚·‚é“Œv“I„‘ª
data {
  int<lower=0> g; //ŒQ‚Ì” 
  array[g] int<lower=0> x; //³”½‰” 
  array[g] int<lower=0> n; //ƒf[ƒ^” 
}
parameters {
  array[g] real<lower=0, upper=1> p; //•ê”ä—¦
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
  array[g, g] real<lower=0, upper=1> U2; //2…€”äŠr
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
    xaste[i] = binomial_rng(n[i], p[i]); //—\‘ª•ª•z
    log_lik = log_lik + binomial_lpmf(x[i] | n[i], p[i]);
  } //‘Î”Šm—¦
}

// Multiple Linear Regression Model
data {
  int<lower=0> n; // Number of data points
  int<lower=0> p; // Number of predictor variables
  vector[n] y; // Dependent variable
  matrix[n, p + 1] X; // Predictor variable matrix
}
parameters {
  vector[p + 1] bb; // Regression coefficients vector
  real<lower=0> sigma; // Error standard deviation
}
transformed parameters {
  vector[n] yhat; // Predicted values yhat
  yhat = X * bb;
}
model {
  y ~ normal(yhat, sigma); // Normal distribution model
}
generated quantities {
  real a; // Intercept
  vector[p] b; // Regression coefficients (without intercept)
  real log_lik; // Log-likelihood
  real<lower=0> vyhat; // Variance of yhat
  real<lower=0, upper=1> r2; // Coefficient of determination (R-squared)
  real<lower=0, upper=1> r; // Multiple correlation coefficient
  
  a = bb[p + 1];
  for (i in 1 : p) {
    b[i] = bb[i];
  }
  vyhat = (variance(yhat) * (n - 1)) / n;
  r2 = vyhat / (vyhat + sigma ^ 2);
  r = sqrt(r2);
  log_lik = 0;
  for (i in 1 : n) {
    log_lik = log_lik + normal_lpdf(y[i] | yhat[i], sigma);
  }
}

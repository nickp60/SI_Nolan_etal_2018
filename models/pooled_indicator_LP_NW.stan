data {
  int<lower=0> N;           // total number of datapoints
  int<lower=0> Ntrt;        // total number of treatments
  vector[N] y;              // bacterial counts
  vector[N] t;              // timepoints
  int<lower=1, upper=Ntrt> trt[N];            // treatments
}
parameters {
  vector[Ntrt] m1;                  // slope - early stage
  vector[Ntrt] c1;                  // intercept - early stage
  vector[Ntrt] m2;                  // slope - late stage
  vector[Ntrt] c2;                  // intercept - late stage
  real<lower=0> sigma;      // measurement error
  real mu_m1;                       // pooled parameter - early stage slope
  real<lower=0> sigma_m1;           // pooled parameter - early stage slope
  real mu_m2;                       // pooled parameter - late stage slope
  real<lower=0> sigma_m2;           // pooled parameter - late stage slope
  real mu_c1;                       // pooled parameter - early stage intercept
  real<lower=0> sigma_c1;           // pooled parameter - early stage intercept
  real mu_c2;                       // pooled parameter - late stage intercept
  real<lower=0> sigma_c2;           // pooled parameter - late stage intercept
}
transformed parameters {
  vector[N] y_hat;         // the *real* measured value, without measurement error
  
  for (i in 1:N) {
    if (t[i] <= 7) {                      // hardcode transition timepoint
      y_hat[i] = m1[trt[i]] * t[i]+ c1[trt[i]];           // simple linear fit, early stage
    } else {
      y_hat[i] = m2[trt[i]] * t[i] +c2[trt[i]];           // simple linear fit, late stage
    }
  }
}
model {
//  mu_m1 ~ normal(-0.4, 0.1);
  m1 ~ cauchy(mu_m1, sigma_m1);
  c1 ~ cauchy(mu_c1, sigma_c1);
//  mu_m2 ~ normal(0, 0.1);
  m2 ~ cauchy(mu_m2, sigma_m2);
  c2 ~ cauchy(mu_c2, sigma_c2);
  
  y ~ normal(y_hat, sigma);               // fit with error
}

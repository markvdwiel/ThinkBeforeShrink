#------------------------------------------------------------------------------#
# Bayesian logistic regression with global and local shrinkage priors
#------------------------------------------------------------------------------#
# global shrinkage with half-Cauchy prior.  Sampling from the half-Cauchy as in 
# Makalic and Schmidt (2016), IEEE Signal Processing Letters:
# tau ~ C+(0,sc) if tau^2 | gam^2 ~ IG(1/2,1/gam^2) and gam^2 ~ IG(1/2,1/sc^2)
code_bin_brg_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> n;
    matrix[n,p] X;
    int<lower=0, upper=1> y[n];
    real<lower=0> a;
    real<lower=0> b;
    real<lower=0> sc;
  }
  parameters {
    real alpha;
    vector[p] beta;
    real<lower=0> gam2;
    real<lower=0> tau2;
  }
  model {
      alpha ~ normal(0, sqrt(100));
      gam2 ~ inv_gamma(a, 1/sc^2);
      tau2 ~ inv_gamma(b, 1/gam2);
      for(i in 1:p){
        beta[i] ~ normal(0, sqrt(tau2));
      }
      y ~ bernoulli_logit(X*beta + alpha);
  }
"
# local shrinkage with beta-prime prior
code_bin_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> n;
    matrix[n,p] X;
    int<lower=0, upper=1> y[n];
    real<lower=0> a;
    real<lower=0> b;
    real<lower=0> sc;
    int g[p];
    int<lower=0> K;
  }
  parameters {
    real alpha;
    vector[p] beta;
    real<lower=0>  gam2[K];
    real<lower=0>  tau2[K];
  }
  model {
    alpha ~ normal(0, sqrt(100));
    for (k in 1:K) {
      gam2[k] ~ inv_gamma(a, 1/sc^2);
      tau2[k] ~ inv_gamma(b, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(tau2[g[i]]));
    }
    y ~ bernoulli_logit(X*beta + alpha);
  }
"

#------------------------------------------------------------------------------#
# Bayesian logistic regression with global and local shrinkage priors
#------------------------------------------------------------------------------#
# global shrinkage with beta-prime prior. Sampling from the half-Cauchy as in 
# Makalic and Schmidt (2016), IEEE Signal Processing Letters:
# tau ~ C+(0,sc) if (tau^2 | gam^2) ~ IG(1/2,1/gam^2) and gam^2 ~ IG(1/2,1/sc^2)
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

simandfit <- function(iter,n,k,S,beta0,beta1,Xtest,prob.test){
  mse = matrix(0,iter,6)
  time = matrix(0,iter,6)
  allbs <- list()
  for (i in 1:iter){
    printi(i,5)
    #i<-2
    x = mvrnorm(n,mu=rep(0,k),Sigma=S)
    lp = beta0 + x %*% beta1
    prob = 1/(1+exp(-lp))
    y = rbinom(n,1,prob=prob)
    
    #-------------
    # maximum likelihood, with Firth in case of separation
    pt1 <- proc.time()
    fit_glm = glm(y~x,family="binomial")
    b_glm = coef(fit_glm)
    pt_glm <- proc.time() - pt1
    
    pt_f <- proc.time()
    fit_f = hush(logistf(y~x))
    b_f <- coef(fit_f)
    pt_f <- proc.time() - pt_f
    
    
    #-------------
    # glmnet + fixed tau2 (default prior)
    # shrinkage with ridge with fixed prior N(0,tau^2)
    pt1 <- proc.time()
    tau2 = 0.5               # Sullivan and Greenland (2013) p. 311
    suppressWarnings({
      fit_glmnet = hush(glmnet(x, y, alpha = 0, family="binomial", lambda=1/(n*tau2)))
    }
    )
    b_glmnet = coef(fit_glmnet,s=1/(n*tau2))
    b_glmnet = drop(b_glmnet)
    pt_glmnet <- proc.time() - pt1
    
    
    #-------------
    # glmnet + cv
    pt1 <- proc.time()
    tryCatch({
      fit_glmnet_cv = hush(cv.glmnet(x, y, alpha = 0, family="binomial"))
    }, error = function(e) {NULL}
    )
    b_glmnet_cv = coef(fit_glmnet_cv, s="lambda.min")   # or choose s="lambda.1se"
    b_glmnet_cv = drop(b_glmnet_cv)
    pt_glmnetcv <- proc.time() - pt1
    
    # bayesian local shrinkage, scale=0.5
    
    # data for stan
    dat <- list(y = y, X = x, g = 1:ncol(x), n = nrow(x), p = ncol(x), K = ncol(x), a = 0.5, b = 0.5, sc=sqrt(0.5))
    
    pt1b <- proc.time()
    
    # sample
    fit_rstanb <- hush(sampling(stan_bin_brl_HalfCauchy, data = dat, iter = 25000, warmup = 5000, thin = 10, verbose = FALSE, chains = 1))
    
    # posterior summary betas
    b_bayes_locb <- c(mean(extract(fit_rstanb, "alpha")$alpha), colMeans(extract(fit_rstanb, "beta")$beta))
    pt_bayes_locb <- proc.time() - pt1b
    
    #-------------
    # bayesian global shrinkage
    
    # data for stan
    dat <- list(y = y, X = x, n = nrow(x), p = ncol(x), a = 0.5, b = 0.5, sc=sqrt(0.5))
    
    pt1b <- proc.time()
    
    # sample
    fit_rstan2b <-  hush(sampling(stan_bin_brg_HalfCauchy, data = dat, iter = 25000, warmup = 5000, thin = 10, verbose = FALSE, chains = 1))
  
    
    # posterior summary betas
    b_bayes_glob <- c(mean(extract(fit_rstan2b, "alpha")$alpha), colMeans(extract(fit_rstan2b, "beta")$beta))
    pt_bayes_glob <- proc.time() - pt1
    
    #-------------
    # results
    bs <- cbind(b_glm,b_glmnet,b_glmnet_cv, b_bayes_locb,b_bayes_glob,b_f)
    
    # compute estimated probabilities in test set
    prob_glm = 1/(1+exp(-Xtest %*% b_glm))
    prob_glmnet = 1/(1+exp(-Xtest %*% b_glmnet)) 
    prob_glmnet_cv = 1/(1+exp(-Xtest %*% b_glmnet_cv))
    prob_brlb = 1/(1+exp(-Xtest %*% b_bayes_locb))
    prob_brgb = 1/(1+exp(-Xtest %*% b_bayes_glob))
    prob_f = 1/(1+exp(-Xtest %*% b_f))
    
    
    # compute MSE of probabilities over test set
    mse[i,1] = sqrt(mean((prob_glm - prob.test)^2))
    mse[i,2] = sqrt(mean((prob_glmnet - prob.test)^2))
    mse[i,3] = sqrt(mean((prob_glmnet_cv - prob.test)^2))
    mse[i,4] = sqrt(mean((prob_brlb - prob.test)^2))
    mse[i,5] = sqrt(mean((prob_brgb - prob.test)^2))
    mse[i,6] = sqrt(mean((prob_f - prob.test)^2))
    
    # time
    time[i,1] = pt_glm[3]
    time[i,2] = pt_glmnet[3]
    time[i,3] = pt_glmnetcv[3]
    time[i,4] = pt_bayes_locb[3]
    time[i,5] = pt_bayes_glob[3]
    time[i,6] = pt_f[3]
    allbs <- c(allbs,list(bs))
  }
  colnames(mse) <- colnames(time) <- c("ML", "glmnet05","glmnet_cv","bayes_loc05", "bayes_glo05", "firth")
  return(list(allbs = allbs,mse = mse, time = time))
}



#auxiliary function
printi <- function(i,printn=50){
  if(i%%printn == 0) print(i)
}

hush=function(code){
  sink(nullfile()) # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}
#' # Think before you shrink, Example script
#' 
#' Script: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl;
#' Package 'shrinkage', version 1.2, written by Gwenaël Leday, g.g.r.leday@wur.nl
#' 
#' This script illustrates usage of R-packages 'mgcv', 'rstan' and 'shrinkage' 
#' for the purpose of regression-based prediction of continuous and binary response. 
#' We focus on the assets of these methods as compared to default alternatives 
#' like simple, one-penalty ridge regression. In particular:
#' 
#' 1) Covariate-group specific shrinkage
#' 2) Confidence intervals that propagate uncertainty of the penalty parameter(s)
#' 3) Local shrinkage using a hierarchial prior (Bayesian setting: rstan, shrinkage) 
#' 
#' The manuscript "Think before you shrink: Alternatives to default shrinkage methods can 
#' improve prediction accuracy, calibration and coverage" illustrates the potential
#' benefits of 1), 2) and 3) for improving prediction accuracy, coverage and 
#' calibration, respectively.
#' 
#' Note that scripts for reproducing the results presented in the manuscript are 
#' available from the same repository. These contain more extensive coding,
#' also for other (more default) shrinkage methods, such as ridge and Firth, as
#' discussed in the manuscript.
#'  
#' 
#' ## Preliminaries
#' loading libraries; using suppressWarnings() to avoid clutter
library(mgcv) #available from CRAN

#' To install 'shrinkage' package from github run:
#' library(devtools); install_github("gleday/shrinkage").
#' If you haven't installed devtools yet, first run: install.packages("devtools"). 
#' Alternatively (without devtools), download the entire folder of source material 
#' from GitHub repo gleday/shrinkage and install from source using
#' install.packages("[path]/shrinkage-master/", repos=NULL, type="source")
library(shrinkage)

#' ## Linear regression
#' Load example training and test data sets (train and test objects). 
#' These are subsets of the large synthetic Helius data set, also available 
#' on the repository. ADJUST PATH, OR USE setwd()!
load("dsynth_example.Rdata") 

dim(train)
head(train)
colnames(train)

#' Retrieve response (sbp) and covariates
Y <- train[,1]
X <- train[,-1]

dim(test)
head(test)
Yte <- test[,1]
Xte <- test[,-1]

#' ## Fit models
#' 
#' ### mgcv

#' nr of covariates
ncov <- ncol(X)
ncov

#' names of covariates
colnames(X)

#' Covariate group 1: age (1), gender (2), bmi (10). Covariate group 2: 
#' all others.
whknown <- c(1,2,10)

#' Create two vectors with 1's at position of the covariates that are in 
#' the group and 0's elsewhere
diag1 <- rep(0,ncov);diag1[whknown] <- 1;diag2 <- rep(0,ncov);diag2[-whknown] <- 1;
diag1
diag2

#' Specification of penalty matrix for gam. Argument 'sp=-1' indicates that 
#' penalties need to be estimated. 
PP2 <- list(X=list(diag(diag1),diag(diag2)),sp=-1)

#' Fit ridge regression model and estimate multiple penalties using maximum 
#' (marginal) likelihood. Argument 'scale=-1' indicates that error variance 
#' should be estimated as well
pred2 <- gam(Y ~ 1+X,family="gaussian",paraPen=PP2,method="ML",scale=-1)

#' Coefficients and their uncertainties
summary(pred2)

#' Retrieve estimated penalties for both groups
pred2$sp

#' Prediction for test samples; Argument 'se.fit=T' indicates that standard errors 
#' should be returned as well; Argument 'unconditional = T' indicates that these
#' are retrieved from the Bayesian posterior covariance matrix
pred2test <- predict.gam(pred2,list(X=Xte),se.fit=T,unconditional=T)
pred2test$fit[1:10]
pred2test$se.fit[1:10]

#' Confidence intervals of predictions
cis <- rbind(pred2test$fit - 1.96*pred2test$se.fit, pred2test$fit + 1.96*pred2test$se.fit)
cis[,1:10]

#' ### shrinkage
#' 
#' The shrinkage package fits Bayesian ridge variations with prior beta ~ N(0,tau^2 sigma^2),
#' where sigma^2 is the error variance. The default shrinkage prior for tau^2 
#' is a BetaPrime(a=0.5, b=0.5) prior, which is the same as a half-Cauchy(0,1) on tau; see Polson & Scott, 
#' Bayesian Analysis, 2012). 

#' Add intercept to design matrix
Xint <- cbind(intercept=rep(1,nrow(X)),X)

#' Fit Bayesian ridge with global shrinkage using mcmc. 
bay <- brg(y = Y, X = Xint)

#' Fit Bayesian ridge with local prior, meaning all covariates (including the intercept) 
#' are separate groups
ncov1 <- ncol(Xint)
bay_loc <- brl(y = Y, X = Xint, g = 1:ncov1)

#' Now fit with two groups of covariates as above (hence allowing different shrinkage); 
#' As the intercept is added to the design matrix, it is added as a separate group
groups <- rep(3,ncol(Xint));groups[1] <- 1; groups[whknown+1] <- 2; 
cbind(colnames(Xint),groups)

#' fit Bay2 (2 covariate groups as above)
bay2 <- brl(y = Y, X = Xint, g = groups)

#' Summary of posteriors of coefficients
bay2sum <- bay2$betas_summary
bay2sum

#' Predictions test set, obtain from posterior means of betas 
pred_bay2 <- bay2sum[1,1] +  Xte %*% bay2sum[-1,1,drop=FALSE]

#' Obtain credible intervals for predictions test set. Three steps: 1. 
#' retrieve posterior samples for beta's; 2: compute posteriors for test sample 
#' predictions (including the intercept); 3: compute 95 perc cred interval from the quantiles
postbay2 <- bay2$betas
postbay2_test <- postbay2[1,] + Xte %*% postbay2[-1,]
cibay2_test <- apply(postbay2_test,1,quantile,probs=c(0.025,0.975))
cibay2_test[,1:10]

#' ### Comparison ridge2 and Bay2

#' Plot 20 test samples (of 100) to compare results of ridge2 (mgcv; black) with Bay2 
#' (shrinkage; red); point predictions and intervals
par(mar=c(2,2,2,1))
ord <- order(pred2test$fit)[seq(1,100,length.out=20)]
plot(1:20,pred2test$fit[ord], pch=16,cex=0.7, ylim=c(-2.3,2), xlab="",ylab="")
points(1:20+0.2,pred_bay2[ord], pch=17,cex=0.7,col="red")
segments(1:20,cis[1,ord],1:20,cis[2,ord])
segments(1:20+0.2,cibay2_test[1,ord],1:20+0.2,cibay2_test[2,ord], col="red")
legend(13,-1.6, legend=c("ridge2 (mgcv)","Bay2 (shrinkage)"), col=c("black","red"),pch=c(16,17),
       lty=c(1,1))

#' ## Logistic regression 
#' loading libraries; using suppressWarnings() to avoid clutter
suppressWarnings(library(rstan)) #available from CRAN
library(glmnet)

#' Load example training and test data set. 
#' These are simulated as in the manuscript, n=50, p=5 (+int), moderate setting 
#' ADJUST PATH, OR USE setwd()!
load("binary_example.Rdata") 

dim(train_bin)
head(train_bin)
colnames(train_bin)

#' Retrieve binary response and covariates
Y <- train_bin[,1]
X <- train_bin[,-1]

dim(test_bin)
head(test_bin)
Yte <- test_bin[,1]
Xte <- test_bin[,-1]


#' Glmnet with fixed ridge penalty, lambda = 1/tau^2 = 1/0.5 (ridge05), as advocated in 
#' Sullivan & Greenland (2013), Int J Epid. Note: glmnet scales loglik with n 
#' (see ?glmnet -> Details), so penalty parameter needs to be scaled accordingly
n <- length(Y)
ridge05 <- glmnet(X, Y, alpha = 0, family="binomial", lambda=1/(n*0.5)) #alpha = 0 -> ridge regression

#' Coefficients, including intercept
c(ridge05$a0,as.numeric(ridge05$beta))

#' Predicted probabilities from glmnet on test set. Uncertainties are not provided.
predridge05 <- predict(ridge05, newx=Xte, type="response")
predridge05[1:10]

#' Stan models need to be compiled once before fitting
source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Github/ExampleThinkShrink_functions.R')

#' Compiles stan code for Bayesian fits; takes some time.
pmt <- proc.time()
stan_bin_brg_HalfCauchy <- stan_model(model_name = "stan_bin_brg_HalfCauchy", model_code = code_bin_brg_HalfCauchy)
stan_bin_brl_HalfCauchy <- stan_model(model_name = "stan_bin_brl_HalfCauchy", model_code = code_bin_brl_HalfCauchy)
proc.time()-pmt


#' Bayesian global shrinkage: all beta ~ N(0,tau^2), tau ~ half-Cauchy(0,sqrt(0.5)). The 
#' latter implies prior mean of tau^2 equals 0.5, which matches with the frequentist 
#' penalty above. Sampling from half-Cauchy as in  Makalic and Schmidt (2016), IEEE Signal 
#' Processing Letters: tau ~ C+(0,sc) if (tau^2 | gam^2) ~ IG(1/2,1/gam^2) and gam^2 ~ IG(1/2,1/sc^2)
dat <- list(y = Y, X = X, n = nrow(X), p = ncol(X), a = 0.5, b = 0.5, sc=sqrt(0.5))
Bay_glob05 <- sampling(stan_bin_brg_HalfCauchy, data = dat, iter = 25000, 
                           warmup = 5000, thin = 10, chain = 1)
#' Extensive summary
summary(Bay_glob05)$summary

#' Posteriors of coefficients
intpost <-  extract(Bay_glob05, "alpha")$alpha #posterior mcmc samples intercept
betapost <- extract(Bay_glob05, "beta")$beta #posterior mcmc samples five betas
dim(betapost)

#' Point estimates of coefficients (posterior means)
intest <- mean(intpost)
betaest <- apply(betapost,2,mean)
c(intest,betaest)

#' Compute posteriors of linear predictor and of predicted probabilities for test set
etapost <- cbind(rep(1,nrow(Xte)),Xte) %*% t(cbind(intpost,betapost)) 
expit <- function(x) exp(x)/(1+exp(x))
probtestpost <- expit(etapost)
dim(probtestpost)

#' Predicted probabilities on test samples (posterior mean)
predBay_glob05 <- apply(probtestpost,1,mean)
predBay_glob05[1:10]

#' Credible intervals for test set predictions
ciBay_glob05  <- apply(probtestpost,1,quantile,probs=c(0.025,0.975))
ciBay_glob05[,1:10]

#' ### Comparison ridge05 and Bay_glob05
#' 
#' Plot 20 test samples (of 50) to compare results of ridge05 (glmnet; black) with 
#' Bay_glob05 (stan; red); point predictions (both) and intervals (Bay_glob05 only)
par(mar=c(2,2,2,1))
ord <- order(predridge05)[seq(1,50,length.out=20)]
plot(1:20,predridge05[ord,1], pch=16,cex=0.7, ylim=c(0,1), xlab="",ylab="")
points(1:20+0.2,predBay_glob05[ord], pch=17,cex=0.7,col="red")
segments(1:20+0.2,ciBay_glob05[1,ord],1:20+0.2,ciBay_glob05[2,ord], col="red")
legend(13,0.15, legend=c("ridge05 (glmnet)","Bay_glob05 (stan)"), col=c("black","red"),pch=c(16,17),
       lty=c(1,1))

#' ### Bayesian local (or grouped) shrinkage for logistic regression
#' 
#' Bayesian local shrinkage, beta[k] ~ N(0,tau2[k]), tau[k] ~ half-Cauchy(0,sqrt(0.5))
#' g indicates the groups, in this case each covariate is its own group; K is the number 
#' of covariate groups
dat <- list(y = Y, X = X, g = 1:ncol(X), n = nrow(X), p = ncol(X), K = ncol(X), 
            a = 0.5, b = 0.5, sc=sqrt(0.5))
#'
#' Note: stan may give warnings for the sampling, as the local shrinkage factors (tau2[k]) 
#' are hard to estimate. Manuscript shows however that prediction results are 
#' are solid, and superior to those of max lik (which is comparable in the sense
#' that it does not group covariates either)
Bay_loc05 <- sampling(stan_bin_brl_HalfCauchy, data = dat, iter = 25000, 
                      warmup = 5000, thin = 10, chain = 1)
#' Extensive summary
summary(Bay_loc05)$summary

#' Covariate groups; e.g. global shrinkage for cov 1,2,3,5 and local for cov 4.
#' Set number of covariate groups K=2 
covgroups <- c(1,1,1,2,1) 
dat <- list(y = Y, X = X, g = covgroups, n = nrow(X), p = ncol(X), K = 2, 
            a = 0.5, b = 0.5, sc=sqrt(0.5))
Bay_205 <- sampling(stan_bin_brl_HalfCauchy, data = dat, iter = 25000, 
                      warmup = 5000, thin = 10, chain = 1)

#' Extensive summary (should shrink beta4 estimate less than Bay_glob05)
summary(Bay_205)$summary
summary(Bay_glob05)$summary

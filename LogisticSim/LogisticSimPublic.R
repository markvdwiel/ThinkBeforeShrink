#' # Think before you shrink
#' 
#' Script: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl
#' 
#' Script allows for repeating calculations as performed in the Manuscript:
#' "Think before you shrink: Alternatives to default shrinkage methods can 
#' improve prediction accuracy, calibration and coverage". This script presents 
#' results for logistic regression on simulated data. Results for linear 
#' regression are presented in the HeliusSubsetsPublic.R script.
#' Note: results may deviate slightly from those in the manuscript, due to use
#' of different simulated sets.
#'  
#' Script evaluates predictive performance in terms of accuracy and calibration 
#' for several shrinkage methods, including some Bayesian ones. 
#' Focus lies on variations of ridge-type shrinkage.
#' 
#' ## Preliminaries
#'
#' Load libraries; Using suppressWarnings() to avoid clutter
suppressWarnings(library(logistf)) #for Firth
library(glmnet) #for logistic ridge
library(MASS) #for sampling covariates from mvt normal
library(rstan) #for Bayesian ridge variations


#' Source code that includes the fitting functions. 
#setwd() #to change working directory
source('LogisticSimPublic_functions.R')

#' Compiles stan code for Bayesian fits; takes some time
pmt <- proc.time()
stan_bin_brl_HalfCauchy <- stan_model(model_name = "stan_bin_brl_HalfCauchy", model_code = code_bin_brl_HalfCauchy)
stan_bin_brg_HalfCauchy <- stan_model(model_name = "stan_bin_brg_HalfCauchy", model_code = code_bin_brg_HalfCauchy)
proc.time()-pmt

#' ## Simulation
#' 
#' Simulation set-up.
#' Manuscript presents three (x3) scenarios. All can be reproduced here.
#' 1. Sample n=50, beta1 = c(0.2,0.2,0.2,0.5,0.8) [Main document]
#' 2. As 1, but with n=100 [Supp Mat]
#' 3. As 1. but with 5 zero's added to the beta vector. [Supp Mat]
#' For all these scenarios, three settings of beta are considered:
#' weak: beta = beta1/3, normal: beta=beta1, strong: beta = 3*beta1.`
#' 
#' Sample size
n = 50
#n =100

#' Add 5 zeros to beta vector?
addzeros <- FALSE
#addzeros <- TRUE

#' weak, normal or strong betas?; Correspond to left, middle and right panels 
#' in manuscript
setting <- "normal"
#setting < - "weak"
#setting <- "strong"

#' Event probability
PE = 0.5           

#' Nr of simulations.  Equals 50 for manuscript, but you may want to set it
#' to a lower number to try the fitting first
iter = 50   

#' Correlation in X
rho = 0.5    

#' Regression coefficients
beta1 = c(0.2,0.2,0.2,0.5,0.8) 
if(setting == "weak") beta1 = beta1/3
if(setting == "strong") beta1 = beta1*3
if(addzeros) beta1=c(rep(0,5),beta1)           # add 5 "noise" predictors
k = length(beta1)

#' File for storing results; state working directory for storing results
wd <- "C:\\ExternData\\Helius\\"
fileout <- paste(wd,"logisticsim","_n",n,"_rho",rho,"_PE",PE,"_addzeros",addzeros,"_betas",setting,".Rdata",sep="")
fileout

#' Simulate test covariates
set.seed(123)
S = matrix(rho,k,k)                  
diag(S) = 1
x.test = mvrnorm(10^5,mu=rep(0,k),Sigma=S)

#' Compute beta0 to meet desired PE
tmp = function(beta0){
  lp = beta0 + x.test %*% beta1
  mean(1/(1+exp(-lp))) - PE
}
beta0 = uniroot(tmp,interval=c(-10,10))$"root"
beta = c(beta0,beta1)
beta

#' Results test set using true coefficients
Xtest = cbind(1, x.test)
lp.test <- Xtest %*% beta
prob.test = 1/(1+exp(-lp.test))   
ntest <- length(lp.test)
ytest <- rbinom(ntest,1,prob.test)

#' Simulate training data, perform fits, and compute MSEs. 
#' Takes considerable time; consider running with a low value of iter first.
#' Note: you may ignore any error message of the type: 
#' 'possible Error in x$.self$finalize() : attempt to apply non-function'.
#' This is harmless. suppressWarnings() is used to avoid clutter.
pt0 <- proc.time()
simres <- suppressWarnings(simandfit(iter,n,k,S,beta0,beta1,Xtest,prob.test))
totaltime <- proc.time()-pt0
totaltime
allbs <- simres$allbs;mse <- simres$mse; time <- simres$time
save(setting, addzeros, totaltime, allbs, mse, time, beta1, Xtest, lp.test, file = fileout)

#' ## Analysis and plot simulation results
#' 
#' Load simulation results if stored in a file
#load("logisticsim_n50_rho0.5_PE0.5_addzerosFALSE_betasnormal.Rdata")

#' Computes cslope
pmt <- proc.time()
slope = matrix(0,length(allbs),6)
for(i in 1:length(allbs)){
bs <- allbs[[i]]
lps <- Xtest %*% bs
  for(k in 1:6){
  slope[i,k] = coef(lm(lps[,k] ~ lp.test))[2]
  }
}
proc.time()-pmt

#' Computes classical calibration slopes; takes some time 
pmt <- proc.time()
calslope = matrix(0,length(allbs),6)
for(i in 1:length(allbs)){
  bs <- allbs[[i]]
  lps <- Xtest %*% bs
  for(k in 1:6){
    calslope[i,k] = coef(glm(ytest ~ lps[,k],family="binomial"))[2]
  }
}
proc.time()-pmt

#' Plots cslope and MSE
  dfslopemat0 <- apply(apply(slope,c(1,2),min, y=3), c(1,2),max, y=1/3)
  cols <- c(1,6,3,2,5,4)
  dfslopemat <- dfslopemat0[,cols]
  colnames(dfslopemat) <- c("ML","Firth","ridgeCV","ridge05","Bay_glo05","Bay_loc05")
  if(setting == "weak") la <- "Slope; Weak signal"
  if(setting == "normal") la <- "Slope; Moderate signal"
  if(setting == "strong") la <- "Slope; Strong signal"
  par(yaxt="n")
  boxplot(log(dfslopemat,3),xlab="",main=la, ylim=c(-1,1))
  abline(h=0)
  par(yaxt="s")
  axis(2,at=c(-1,log(1/2,3),0,log(2,3),1),labels=c("1/3","1/2","1","2","3"))
  
  dfMSEmat0 <- mse
  dfMSEmat <- dfMSEmat0[,cols]
  colnames(dfMSEmat) <- c("ML","Firth","ridgeCV","ridge05","Bay_glo05","Bay_loc05")
  if(setting == "weak") la <- "MSEp; Weak signal"
  if(setting == "normal") la <- "MSEp; Moderate signal"
  if(setting == "strong") la <- "MSEp; Strong signal"
  boxplot(dfMSEmat,xlab="", main=la)
  


#' Plots inverted cslope and classical calibration slope 
  dfslopemat0 <- apply(apply(1/slope,c(1,2),min, y=10), c(1,2),max, y=1/10)
  c("ML", "glmnet05","glmnet_cv","bayes_loc05", "bayes_glo05", "firth")
  cols <- c(1,6,3,2,5,4)
  dfslopemat <- dfslopemat0[,cols]
  colnames(dfslopemat) <- c("ML","Firth","ridgeCV","ridge05","Bay_glo05","Bay_loc05")
  if(setting == "weak") la <- "Inverted cslope; Weak signal"
  if(setting == "normal") la <- "Inverted cslope; Moderate signal"
  if(setting == "strong") la <- "Inverted cslope; Strong signal"
  par(yaxt="n")
  boxplot(log(dfslopemat,3),xlab="",main=la, ylim=c(-1,1))
  abline(h=0)
  par(yaxt="s")
  axis(2,at=c(log(1/10,3),-1,log(1/2,3),0,log(2,3),1,log(10,3)),labels=c("1/10", "1/3","1/2","1","2","3","10"))

  dfslopemat0 <- apply(apply(calslope,c(1,2),min, y=10), c(1,2),max, y=1/10)
  c("ML", "glmnet05","glmnet_cv","bayes_loc05", "bayes_glo05", "firth")
  cols <- c(1,6,3,2,5,4)
  dfslopemat <- dfslopemat0[,cols]
  colnames(dfslopemat) <- c("ML","Firth","ridgeCV","ridge05","Bay_glo05","Bay_loc05")
  if(setting == "weak") la <- "Calibration slope; Weak signal"
  if(setting == "normal") la <- "Calibration slope; Moderate signal"
  if(setting == "strong") la <- "Calibration slope; Strong signal"
  par(yaxt="n")
  boxplot(log(dfslopemat,3),xlab="",main=la, ylim=c(-1,1))
  abline(h=0)
  par(yaxt="s")
  axis(2,at=c(log(1/10,3),-1,log(1/2,3),0,log(2,3),1,log(10,3)),labels=c("1/10", "1/3","1/2","1","2","3","10"))
  

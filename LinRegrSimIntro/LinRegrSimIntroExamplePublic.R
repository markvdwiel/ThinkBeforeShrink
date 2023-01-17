#' # Think before you shrink
#' 
#' Script: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl
#' Package 'shrinkage', version 1.2, written by Gwenaël Leday, g.g.r.leday@wur.nl
#' 
#' Script allows for repeating simulations presented in the Intro of the Manuscript:
#' "Think before you shrink: Alternatives to default shrinkage methods can 
#' improve prediction accuracy, calibration and coverage" on a synthetic 
#' copy of the data provided with this script. Note: results may deviate 
#' somewhat from those in the manuscript due to difference in simulated data. 
#' This script presents results for a linear regression example on simulated data. 
#' Results for logistic regression on simulated data are presented in the 
#' LogisticSimPublic.R script. Results on subsets of Helius data are presented in 
#' HeliusSubsetsPublic.R script.
#'  
#' Script evaluates predictive performance in terms of accuracy and coverage.
#' Focus lies on variations of ridge-type shrinkage. Main prupose is to show
#' the potential benefit of grouped shrinkage.
#' 
#' ## Preliminaries
library(mgcv) #available from CRAN
library(glmnet) #available from CRAN

#' install shrinkage package from github run:
#' library(devtools); install_github("gleday/shrinkage") 
#' OR if you haven't installed devtools yet, run:
#' install.packages("devtools"); library(devtools); 
library(shrinkage) 


#' Source functions
#setwd() #to change working directory
source('LinRegrSimIntroExamplePublic_functions.R')

#' One may skip the most time-consuming part, fitting of the models, if
#' these are already available 
dofit <- T

#' State working directory for storing results
wd <- "C:\\ExternData\\Helius\\"
set.seed(34623)

#' ## Simulation
#' 
#' Set simulation parameters
ngen <- 1000 #number of training sets
n <- 50 #sample size training sets
p <- 6 #nr of covariates besides Treatment
CorX <- 0.2 #correlation between non-Treatment covariates
prop1 <- 0.5 #proportion Treated (binary covariate)
noisesd <- 1 #noise sd

#' Simulation setting (see Manuscript)
setting1 <- T
#' Regression coefficients. First one: intercept; second one: Treatment 
#' coefficient; other: coefficients of remaining p covariates 
if(setting1) { #Setting 1: stronger treatment effect
betavec <- c(0,-0.25,-0.05,-0.05,-0.05,0.05,0.05,0.05) 
fileout <- paste(wd,"penCsim1000.Rdata",sep="")
} else { #Setting 2: equally strong effects
betavec <- c(0,-0.1,-0.1,-0.1,-0.1,0.1,0.1,0.1) 
fileout <-  paste(wd,"penCsimequalbeta.Rdata", sep="")
}

#' Generate ngen small training sets
traindata <- lapply(1:ngen,simdatfun, n=n, p=p, CorX=CorX, prop1=prop1,
                    noisesd=noisesd, betavec=betavec)

#' Generate one large test set
ntest <- 1000
testdata <- lapply(c(1),simdatfun,n=ntest,p=p,CorX=CorX, prop1=prop1,
                   noisesd=noisesd, betavec=betavec)[[1]]

#' Remove test response
testdataX <- testdata[,-1]

pmt <- proc.time()
if(dofit){
penCsim <- lapply(1:ngen,penCompareSim, traindat=traindata, testdatX=testdataX)
#save(penCsim,file=fileout)
  } else {
  load(fileout)
  }
proc.time() - pmt

#' ## Analyse simulation results
#'
#' Retrieve results for ridge, ridge2 and Bay2
#' Latter two methods penalize Treatment coefficient separately from the others
penCmgcv1 <- lapply(penCsim,function(ex) ex[1]) #ridge
penCmgcv2 <- lapply(penCsim,function(ex) ex[2]) #ridge2
penCbayes <- lapply(penCsim,function(ex) ex[3]) #Bay2

#' Retrieve estimated coefficients for ridge, ridge2 and Bay2
coefsmgcv1 <- lapply(penCmgcv1,coeffun,method="mgcv") #ridge
coefsmgcv2 <- lapply(penCmgcv2,coeffun,method="mgcv") #ridge2
coefsbayes <- lapply(penCbayes,coeffun,method="shrinkage") #Bay2

#' ### Prediction accuracy
#' 
#' Obtain true prediction after adding the intercept to the design matrix X
predtrue <- cbind(rep(1,nrow(testdataX)),testdataX) %*% matrix(betavec,ncol=1)

#' Mean squared error of predictions (MSEp); mean and 90 percent quantile 
#' across training sets
psemgcv1 <- PSE(coefsmgcv1, predtesttrue=predtrue, testX=testdataX)
psemgcv2 <- PSE(coefsmgcv2, predtesttrue=predtrue, testX=testdataX)
psebayes<- PSE(coefsbayes, predtesttrue=predtrue, testX=testdataX)

allpse <- cbind(psemgcv1,psemgcv2,psebayes)
colnames(allpse) <- c("ridge", "ridge2", "Bay2")
rownames(allpse) <- c("mean", "q_90")
allpse

#' ### Coverage
#' 
#' Coverage of confidence intervals of predictions for ridge, ridge2, Bay2, 
#' averaged across training sets and test samples
cover1 <- covermgcv_onetest(fit=penCsim, wh =1, testdatX = testdataX, uncond = T, bench=predtrue)
cover2 <- covermgcv_onetest(fit=penCsim, wh =2, testdatX = testdataX, uncond = T, bench=predtrue)

#' Intervals for Bayesian fitter are already available in penCsim object
coverbayes <- coverbayes_onetest(fit=penCsim, wh = 3, bench=predtrue)
covermean <- apply(data.frame(ridge=cover1,ridge2=cover2,Bay2 = coverbayes),2,mean)
covermean

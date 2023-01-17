#generates p covariates with equicorrelation CorX
simdatfun <- function(k,n,p, CorX, prop1, noisesd, betavec){
  X6 <-  matrix(rep(rnorm(n,sd=sqrt(CorX/(1-CorX))),times=p),n,p)+ matrix(rnorm(n*p),n,p)
  Xtr <- 2*rbinom(n,1,prop1)-1
  X <- cbind(Xtr,X6)
  Y <- betavec[1] + X %*% matrix(betavec[-1],ncol=1) + rnorm(n,mean=0,sd=noisesd)
  simdat <- cbind(Y,X)
  return(simdat)
}

#Note that BetaPrime(0.5,0.5) prior is the same as the half-Cauchy(0,1) prior
penCompareSim <- function(setind,traindat,testdatX){
  #setind <- 2
  printi(setind)
  data1 <- traindat[[setind]]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  X1int <- cbind(Intercept=rep(1,nrow(X1)),X1)
  testdatXint <- cbind(Intercept=rep(1,nrow(testdatX)),testdatX)
  ncov <- ncol(X1)
  whknown <- c(1)
  PP1 <- list(X1=list(diag(1,ncov),sp=-1))
  pred1 <- gam(Y1 ~ 1 + X1,family="gaussian",paraPen=PP1,method="ML",scale=-1)
  diag1 <- rep(0,ncov);diag1[whknown] <- 1;diag2 <- rep(0,ncov);diag2[-whknown] <- 1;
  PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)
  pred2 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
  
  groups <- rep(3,ncol(X1int));groups[1] <- 1; groups[whknown+1] <- 2; #2+1 groups: intercept, known covs and unknowns
  res_local_bp_group <- brl(y = Y1, X = X1int, g = groups, prior = "BetaPrime", a = 0.5, b = 0.5,
                            mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  postbeta <- res_local_bp_group$betas
  postetatest <- testdatXint %*% postbeta
  cipredtest <- t(apply(postetatest,1,quantile,probs=c(0.025,0.975)))
  whrem <- which(names(res_local_bp_group )%in% c("betas", "tau2s","sigma2s"))
  res_local_bp_group <-  res_local_bp_group[-whrem]
  res_local_bp_group <- c(res_local_bp_group,list(cipred=cipredtest))
  allres <- list(pred1=pred1,pred2=pred2,res_bay_glo2=res_local_bp_group)
  return(allres)
}

#output="mgcv", "glmnet", "shrinkage"
coeffun <- function(penC1, method="mgcv"){
  #penC1 <- penClasso[[1]];method="glmnet"
  nmod <- length(penC1)
  if(method =="mgcv") {
    coefs <- sapply(1:nmod,function(mod) coef(penC1[[mod]]))
  }
  if(method=="glmnet"){
    coefs <- sapply(1:nmod,function(mod) coef(penC1[[mod]]$fit,s=penC1[[mod]]$optlam))[[1]]
  }
  if(method=="shrinkage"){
    coefs <- sapply(1:nmod,function(mod){res <- penC1[[mod]]; res$betas_summary[, c("Mean")]})
  }
  if(nmod==1) coefs <- matrix(coefs,ncol=1)
  #if(!is.na(whichreorder[1])) coefs[,whichreorder] <- coefs[reorder,whichreorder]
  colnames(coefs) <- names(penC1)
  return(coefs)
}

#prediction on test set from estimated coefficients
predfromcoef <- function(i,coefs, testX){
  #i <- 1; coefs <- coefsmgcv
  coefi <- coefs[[i]]
  preds <- testX %*% coefi
  return(preds)
}


#prediction squared error, mean and 90 perc quantile across training sets
PSE <- function(coefs, predtesttrue, testX){
  #coefs <- coefsmgcv 
  testXint <- cbind(rep(1,nrow(testX)),testX)
  sumPSE <- 0;
  nfit <- length(coefs)
  pses <- c()
  for(i in 1: nfit){
    preds <- predfromcoef(i, coefs, testXint)
    predtesttrue <- predtrue
    pse <- apply(preds,2,function(predsi) mean((predsi-predtesttrue)^2))
    pses <- rbind(pses,pse)
  }
  pse <- mean(pses)
  pseq01 <- apply(pses,2,quantile,probs=c(0.9))
  return(c(pse,pseq01))
}

#computes coverage for ridge methods, fit by mgcv
covermgcv_onetest <- function(fit,wh, sets,uncond=T,bench, testdatX){
  nfit <- length(fit)
  allcov <- sapply(1:nfit,function(i) { 
    printi(i)
    #i<-1
    penC1 <- fit[[i]][[wh]]
    pg <- predict.gam(penC1,list(X1=testdatX),se.fit=T,unconditional=uncond)
    cis <- cbind(pg$fit - 1.96*pg$se.fit, pg$fit + 1.96*pg$se.fit)
    covs <- (bench <= cis[,2]) * (bench >= cis[,1])
    return(covs)
  })
  return(apply(allcov,1,mean))
}

#computes coverage for Bayesian methods, fit by shrinkage
coverbayes_onetest <- function(fit,wh,bench){
  nfit <- length(fit)
  allcov <- sapply(1:nfit,function(i) { 
    penC1 <- fit[[i]][[wh]]
    cis <- penC1$cipred
    covs <- (bench <= cis[,2]) * (bench >= cis[,1])
    return(covs)
  })
  return(apply(allcov,1,mean))}

#auxiliary function
printi <- function(i,printn=50){
  if(i%%printn == 0) print(i)
}


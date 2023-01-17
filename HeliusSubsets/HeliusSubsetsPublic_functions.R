#penCompare fits multiple ridge models using mgcv's gam function. 
penCompare <- function(setind,trainind){
  printi(setind)
  data1 <- helius4[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  ncov <- ncol(X1)
  
  #group 1 covariates
  whknown <- c(1,2,10)
  
  #ethnicity covariates
  whetn <- 3:6
  
  #OLS; checked agaist lm() results
  pred0 <- gam(Y1 ~ 1 + X1,family="gaussian",method="ML",scale=-1)
  
  #step
  predstep <- step(lm(Y1 ~.,data=data.frame(Y1,X1)),trace=0)
  pscoef <- predstep$coefficients; cnX1 <- c("(Intercept)",colnames(X1))
  ma <- match(names(pscoef), cnX1); pscoefnew <- rep(0,ncov+1); pscoefnew[ma] <- pscoef
  names(pscoefnew) <- cnX1
  predstep$coefficients <- pscoefnew
  
  #ridge
  PP1 <- list(X1=list(diag(1,ncov),sp=-1))
  pred1 <- gam(Y1 ~ 1 + X1,family="gaussian",paraPen=PP1,method="ML",scale=-1)
  
  #ridge2
  diag1 <- rep(0,ncov);diag1[whknown] <- 1;diag2 <- rep(0,ncov);diag2[-whknown] <- 1;
  PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)
  pred2 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
  
  #ridge2un
  PP2un <- list(X1=list(diag(diag2)),sp=-1) #do not penalize known covariates
  pred2un <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2un,method="ML",scale=-1)
  
  #ridge2r (random grouping of covariates)
  set.seed(2322 + setind) #for reproducibility 
  whr <- sort(sample(1:ncov,3))
  diag1 <- rep(0,ncov);diag1[whr] <- 1;diag2 <- rep(0,ncov);diag2[-whr] <- 1;
  
  PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)
  pred2r <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
  
  #ridge2unr
  PP2un <- list(X1=list(diag(diag2)),sp=-1)
  pred2run <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2un,method="ML",scale=-1)
  
  
  #ridge3 (ethnicity as extra covariate group)
  diag1 <- rep(0,ncov);diag1[whknown] <- 1;diag2 <- rep(0,ncov);diag2[whetn] <- 1;
  diag3 <- rep(0,ncov);diag3[-c(whknown,whetn)] <- 1
  PP3 <- list(X1=list(diag(diag1),diag(diag2),diag(diag3)),sp=-1)
  pred3 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP3,method="ML",scale=-1)
  
  allres <- list(pred0=pred0,predstep=predstep,pred1=pred1,pred2=pred2, pred2r=pred2r,
                 pred2un=pred2un,pred2run=pred2run,pred3 =pred3)
  return(allres)
}

#fits the lasso with glmnet, including cv
penComparelasso <- function(setind,trainind){
  #setind <- 159
  printi(setind)
  data1 <- helius4[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  n <- length(Y1)
  set.seed(345234)
  cvlas <- cv.glmnet(X1,Y1,standardize=FALSE)
  lam <- cvlas$lambda.min
  pred1las0 <- glmnet(X1,Y1,standardize=FALSE)
  pred1las <- list(fit=pred1las0,optlam = lam)
  return(list(pred1las = pred1las))
}

#fit Bayesian models; Note that BetaPrime(0.5,0.5) prior is the same as the half-Cauchy(0,1) prior
penCompareBayes <- function(setind,trainind,testind){
  printi(setind)
  data1 <- helius4[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  testdata <- as.matrix(helius4[testind,-1])
  ncov <- ncol(X1)
  colnames(X1)
  whknown <- c(1,2,10)
  whetn <- 3:6
  #Xkn1 = X1[,whknown];Xu1=X1[,-whknown]
  X1int <- cbind(Intercept=rep(1,nrow(X1)),X1)
  #Bayes ridge with MML tau2
  res_global_mml <- brg(y = Y1, X = X1int, prior = "ml", mcmc = 5000, verbose=FALSE)
  # fit Bayesian ridge with BetaPrime(0.5, 0.5) on global variance (= half-Cauchy on standard deviation)
  res_global_bp <- brg(y = Y1, X = X1int, prior = "BetaPrime", a = 0.5, b = 0.5,mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  # fit Bayesian ridge with invGamma(1e-05, 1e-05) on global variance (= half-Cauchy on standard deviation)
  res_global_ig <- brg(y = Y1, X = X1int, prior = "invGamma", a = 1e-05, b = 1e-05, mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  # fit Bayesian ridge with local prior
  res_local_bp <- brl(y = Y1, X = X1int, g = 1:(ncov+1), prior = "BetaPrime", a = 0.5, b = 0.5,
                      mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  groups <- rep(2,ncol(X1int));groups[1] <- 1; groups[whknown+1] <- 3; #2+1 groups: intercept, known covs and unknowns
  #fit Bay2 (2 covariate groups as above)
  res_local_bp_group <- brl(y = Y1, X = X1int, g = groups, prior = "BetaPrime", a = 0.5, b = 0.5,
                            mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  #fit Bay2 with 2 random groups
  set.seed(2322 + setind) #for reproducibility 
  whr <- sort(sample(1:ncov,3)) #random groups
  groupsr <- rep(2,ncol(X1int));groupsr[1] <- 1; groupsr[whr+1] <- 3; #2+1 groups: intercept, known covs and unknowns
  res_local_bp_groupr <- brl(y = Y1, X = X1int, g = groupsr, prior = "BetaPrime", a = 0.5, b = 0.5,
                             mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  #fit Bay3 (3 covariate groups as above)
  groups3 <- rep(2,ncol(X1int));groups3[1] <- 1; groups3[whknown+1] <- 3;groups3[whetn + 1] <- 4 #3+1 groups: intercept, 2 known covs and unknowns
  res_local_bp_group3 <- brl(y = Y1, X = X1int, g = groups3, prior = "BetaPrime", a = 0.5, b = 0.5,
                             mcmc = 5000L, burnin = 5000L, thin = 10L, verbose=FALSE)
  allres <- list(pred_mml = res_global_mml,predglo_ig= res_global_ig, 
                 predglo_bp = res_global_bp,predloc_bp_group=res_local_bp_group,predloc_bp_groupr=res_local_bp_groupr,predloc_bp = res_local_bp,
                 predloc_bp_group3=res_local_bp_group3)  
  #95 percent c.i. for predictions on test set
   addci <- function(res){
    postbeta <- res$betas
    postetatest <- postbeta[1,] + testdata %*% postbeta[-1,]
    cipredtest <- t(apply(postetatest,1,quantile,probs=c(0.025,0.975)))
    res <- c(res,list(cipred=cipredtest))
    return(res)
  }
  allres <- lapply(allres,addci)
  allres <- lapply(allres,function(res){whrem <- which(names(res) %in% c("betas", "tau2s","sigma2s"));return(res[-whrem])})
  
  return(allres)
}

#Same as penCompareBayes but for only a few methods, and with more mcmc samples 
#(for obtaining a more fine-grained posterior)
penCompareBayes2post <- function(setind,trainind,testind,nmcmc=5000){
  #setind <- 1;trainind <- td50_400;testind <- testind50
  printi(setind)
  data1 <- helius4[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  testdata <- as.matrix(helius4[testind,-1])
  ncov <- ncol(X1)
  colnames(X1)
  whknown <- c(1,2,10)
  whetn <- 3:6
  #Xkn1 = X1[,whknown];Xu1=X1[,-whknown]
  X1int <- cbind(Intercept=rep(1,nrow(X1)),X1)
  #Bayes ridge with MML tau2
  res_global_ig <- brg(y = Y1, X = X1int, prior = "invGamma", a = 1e-05, b = 1e-05, mcmc = nmcmc, burnin = 5000L, thin = 10L, verbose=FALSE)
  # fit Bayesian ridge with local prior
  res_local_bp <- brl(y = Y1, X = X1int, g = 1:(ncov+1), prior = "BetaPrime", a = 0.5, b = 0.5,
                      mcmc = nmcmc, burnin = 5000L, thin = 10L, verbose=FALSE)
  allres <- list(predglo_ig= res_global_ig,predloc_bp = res_local_bp)  
  addci <- function(res){
    #res <-  res_global_bp
    postbeta <- res$betas
    postetatest <- postbeta[1,] + testdata %*% postbeta[-1,]
    cipredtest <- t(apply(postetatest,1,quantile,probs=c(0.025,0.975)))
    res <- c(res,list(cipred=cipredtest))
    return(res)
  }
  allres <- lapply(allres,addci)
  # allres <- lapply(allres,function(res){whrem <- which(names(res) %in% c("betas", "tau2s","sigma2s"));return(res[-whrem])})
  return(allres)
}

#Only OLS, ridge and ridge2 (runs faster)
penCompare3 <- function(setind,trainind){
  #setind <- 1
  printi(setind)
  data1 <- helius4[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  ncov <- ncol(X1)
  whknown <- c(1,2,10)
  whetn <- 3:6
  #with mgcv; ML for hyperparameters, penlties + sigma
  
  PP1 <- list(X1=list(diag(1,ncov),sp=-1))
  pred0 <- gam(Y1 ~ 1 + X1,family="gaussian",method="ML",scale=-1)
  pred1 <- gam(Y1 ~ 1 + X1,family="gaussian",paraPen=PP1,method="ML",scale=-1)
  
  diag1 <- rep(0,ncov);diag1[whknown] <- 1;diag2 <- rep(0,ncov);diag2[-whknown] <- 1;
  PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)
  pred2 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
  allres <- list(pred0=pred0, pred1=pred1,pred2=pred2)
  #allres <- lapply(allres,function(res){whrem <- which(names(res) %in% c("family", "model","paraPen","formula","control","terms"));return(res[-whrem])})
  return(allres)
}

#Fit OLS, ridge and ridge2 for various training sample sizes
pencompare3n <- function(en,nfit){
  print(paste("Training sample size:",en))
  sets <- makeTraining(en,nfit=nfit) 
  resn <- lapply(1:nfit, penCompare3,trainind=sets)
  return(resn)
}

#retrieve coefficients from fits
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

#compute predictions on test set i (note: assume X and sets are known objects)
predfromcoef <- function(i,coefs){
  coefi <- coefs[[i]]
  Xtesti <- as.matrix(X[-sets[[i]],])
  Xtesti <- cbind(rep(1,nrow(Xtesti)),Xtesti)
  preds <- Xtesti %*% coefi 
  return(preds)
}

#compute prediction squared error on test sets (note: assume X, sets and predtrue 
#are known objects)
RPSE <- function(coefs,sets){
  pses <- c();
  nfit <- length(coefs)
  for(i in 1:nfit){
    #i <- 1; coefs <- coefsmgcv
    preds <- predfromcoef(i,coefs)
    predtesttrue <- predtrue[-sets[[i]]]
    pse <- apply(preds,2,function(predsi) mean((predsi-predtesttrue)^2))
    pses <- cbind(pses,pse)
  }
  return(pses)
}

#computes cslope by regressing estimated prediction on true ones for test sets
calibrateperset <- function(coefs,sets, wh, reverse=F){
  #nfit <- 5; coefs <- coefsmgcv;wh<- c(1,3)
  nfit <- length(coefs)
  allslopes <- c()
  for(i in 1:nfit){
    #i <- 1;
    printi(i)
    preds <- predfromcoef(i,coefs)
    predtesttrue <- predtrue[-sets[[i]]]
    if(reverse) slopes <- sapply(wh,function(k) return(lm(x~0+y,data=data.frame(y=preds[,k],x=predtesttrue))$coef)) else slopes <- sapply(wh,function(k) return(lm(y~0+x,data=data.frame(y=preds[,k],x=predtesttrue))$coef))
    allslopes <- rbind(allslopes,slopes)
  }
  colnames(allslopes) <- colnames(coefs[[1]][,wh])
  return(allslopes)
}

#computes coverage for ridge methods, fit by mgcv
covermgcv <- function(fit,wh, sets,uncond=T,bench,testind){
  nfit <- length(fit)
  allcov <- sapply(1:nfit,function(i) { 
    printi(i)
    #i<-1
    penC1 <- fit[[i]][[wh]]
    testdat <- as.matrix(helius4[testind,-1])
    benchi <- bench[testind]
     pg <- predict.gam(penC1,list(X1=testdat),se.fit=T,unconditional=uncond)
    cis <- cbind(pg$fit - 1.96*pg$se.fit, pg$fit + 1.96*pg$se.fit)
    covs <- (benchi <= cis[,2]) * (benchi >= cis[,1])
    return(covs)
  })
  return(apply(allcov,1,mean))
}

#computes coverage for Bayesian methods, fit by shrinkage
coverbayes <- function(fit,wh,testind,bench){
  #fit <- penCbayes; wh <- 3;sets <- testind50;uncond<- T;bench=predtrue
  nfit <- length(fit)
  allcov <- sapply(1:nfit,function(i) { 
    #i<-1
    penC1 <- fit[[i]][[wh]]
    benchi <- bench[testind]
    cis <- penC1$cipred
    covs <- (benchi <= cis[,2]) * (benchi >= cis[,1])
    return(covs)
  })
  return(apply(allcov,1,mean))}

#computes widths of intervals for ridge methods, fit by mgcv
widthmgcv <- function(fit,wh,sets,uncond=T,bench,testind){
  #fit <- penCmgcv; wh <- 1;sets <- testind50;uncond<- T;bench=predtrue
  nfit <- length(fit)
  allwid <- sapply(1:nfit,function(i) { 
    printi(i)
    penC1 <- fit[[i]][[wh]]
    testdat <- as.matrix(helius4[testind,-1])
    benchi <- bench[testind]
    pg <- predict.gam(penC1,list(X1=testdat),se.fit=T,unconditional=uncond)
    cis <- cbind(pg$fit - 1.96*pg$se.fit, pg$fit + 1.96*pg$se.fit)
    wids <- cis[,2] -  cis[,1]
    return(wids)
  })
  return(apply(allwid,1,mean))
}

#computes widths of intervals for Bayesian methods, fit by shrinkage
widthbayes <- function(fit,wh,testind,bench){
  #fit <- penCbayes; wh <- 3;sets <- testind50;uncond<- T;bench=predtrue
  nfit <- length(fit)
  allwid <- sapply(1:nfit,function(i) { 
    #i<-1
    penC1 <- fit[[i]][[wh]]
    benchi <- bench[testind]
    cis <- penC1$cipred
    wids <- cis[,2] -  cis[,1]
    return(wids)
  })
  return(apply(allwid,1,mean))}

#auxiliary functions
printi <- function(i,printn=50){
  if(i%%printn == 0) print(i)
}


#function for creating random test sets 
CVfolds <- function(Y,model=NULL,balance=TRUE,kfold=10,fixedfolds=TRUE,nrepeat=1){ 
  #Y: response vector, length n
  #model: "linear", "logistic", "cox", etc
  #balance: should the splits balance levels of the response?
  #kfold: scalar, the number of folds
  #fixedfolds: should the folds be fixed? (for reproducibility)
  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold
  
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  response <- Y
  if(model=="linear") balance <- FALSE
  CVfoldsrep <- function(rep){
    nsam <- length(response)
    if (fixedfolds) set.seed(3534+rep-1) #else set.seed(NULL)  #changed 19/4
    if (!balance) {
      rand <- sample(1:nsam)
      grs1 <- floor(nsam/kfold)
      grs2 <- grs1 + 1
      ngr1 <- kfold * grs2 - nsam
      folds <- lapply(1:kfold, function(xg) {
        if (xg <= ngr1)
          els <- rand[(1 + (xg - 1) * grs1):(xg * grs1)]
        else els <- rand[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                               1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
        return(sort(els))
      })
    }
else {
  if (model == "logistic")
    if (class(response) == "factor")
      nev <- which((as.numeric(response) - 1) == 1)
    else nev <- which(response == 1)
    if (model == "cox") nev <- which(response[, 2] == 1)
    nsamev <- length(nev)
    randev <- sample(nev)
    grs1 <- floor(nsamev/kfold)
    grs2 <- grs1 + 1
    ngr1 <- kfold * grs2 - nsamev
    foldsev <- lapply(1:kfold, function(xg) {
      if (xg <= ngr1)
        els <- randev[(1 + (xg - 1) * grs1):(xg * grs1)]
      else els <- randev[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                               1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
      return(els)
    })
    nonev <- setdiff(1:nsam, nev)
    nsamnonev <- length(nonev)
    randnonev <- sample(nonev)
    grs1 <- floor(nsamnonev/kfold)
    grs2 <- grs1 + 1
    ngr1 <- kfold * grs2 - nsamnonev
    foldsnonev <- lapply(1:kfold, function(xg) {
      if (xg <= ngr1)
        els <- randnonev[(1 + (xg - 1) * grs1):(xg *
                                                  grs1)]
      else els <- randnonev[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                                  1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
      return(els)
    })
    folds <- lapply(1:kfold, function(i) sort(c(foldsev[[i]],
                                                foldsnonev[[i]])))
}
return(folds)
}
return(unlist(lapply(1:nrepeat,CVfoldsrep),recursive=FALSE))
}
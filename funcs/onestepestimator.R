#-------------------------------------------------------------------------------#
# Last updated date: 7/10/2020
# Author: Duzhe Wang (dwang282@wisc.edu)
#--------------------------------------------------------------------------------#


library(glmnet)   ## lasso algorithm 
library(flare)    ## solve CLIME 
#library(fastclime)  ## solve CLIME




# one step estimator using flare R package in CLIME method for Model 1

onestepestimator1.flare=function(data, lassolambda, climemu, L2){
  Atrain=data$X  
  btrain=data$Y
  N=nrow(Atrain)
  # step 1: Lasso 
  lasso=glmnet(Atrain, btrain, alpha=1, lambda=lassolambda, intercept = FALSE) # alpha=1 is lasso
  betahat=as.vector(coef(lasso)[-1])
  
  # step 2: set W hat and Sigma hat based on betahat and each model
  What=(1/25)*exp(0.5*abs(Atrain%*%betahat))
  trimmedWhat=rep(0, N)
  for (i in 1:N){
    if (What[i]<=L2){
      trimmedWhat[i]=What[i]
    } else{
      trimmedWhat[i]=L2
    }
  }
  trimmedWhatDiagInverse=diag(1/trimmedWhat)
  Sigmahat=(1/N)*t(Atrain)%*%trimmedWhatDiagInverse%*%Atrain
  
  # step 3: CLIME using flare R package
  climeout=sugm(Sigmahat, lambda=climemu, method="clime")
  Thetahat=climeout$icov1[[1]]
  
  # step 4: one-step estimator
  betatilde=betahat+(1/N)*Thetahat%*%t(Atrain)%*%trimmedWhatDiagInverse%*%(btrain-Atrain%*%betahat)
  # output 
  return(list(betahat=betahat, betatilde=betatilde, trimmedWhat=trimmedWhat, Thetahat=Thetahat))
}




# one step estimator using flare R package in CLIME method for Model 2

onestepestimator2.flare=function(data, lassolambda, climemu, L2){
  Atrain=data$X  
  btrain=data$Y
  N=nrow(Atrain)
  # step 1: Lasso 
  lasso=glmnet(Atrain, btrain, alpha=1, lambda=lassolambda, intercept = FALSE) # alpha=1 is lasso
  betahat=as.vector(coef(lasso)[-1])
  
  # step 2: set W hat and Sigma hat based on betahat and each model
  What=(1/50)*exp(0.05*(Atrain%*%betahat)^2)
  trimmedWhat=rep(0, N)
  for (i in 1:N){
    if (What[i]<=L2){
      trimmedWhat[i]=What[i]
    } else{
      trimmedWhat[i]=L2
    }
  }
  trimmedWhatDiagInverse=diag(1/trimmedWhat)
  Sigmahat=(1/N)*t(Atrain)%*%trimmedWhatDiagInverse%*%Atrain
  
  # step 3: CLIME using flare R package
  climeout=sugm(Sigmahat, lambda=climemu, method="clime")
  Thetahat=climeout$icov1[[1]]
  
  # step 4: one-step estimator
  betatilde=betahat+(1/N)*Thetahat%*%t(Atrain)%*%trimmedWhatDiagInverse%*%(btrain-Atrain%*%betahat)
  # output 
  return(list(betahat=betahat, betatilde=betatilde, trimmedWhat=trimmedWhat, Thetahat=Thetahat))
}




# one step estimator using flare R package in CLIME method for Model 3
onestepestimator3.flare=function(data, lassolambda, climemu, L1, L2){
  Atrain=data$X  
  btrain=data$Y
  N=nrow(Atrain)
  # step 1: Lasso 
  lasso=glmnet(Atrain, btrain, alpha=1, lambda=lassolambda, intercept = FALSE) # alpha=1 is lasso
  betahat=as.vector(coef(lasso)[-1])
  
  # step 2: set W hat and Sigma hat based on betahat and each model
  What=0.25*abs(Atrain%*%betahat)
  trimmedWhat=rep(0, N)
  for (i in 1:N){
    if (What[i]<=L1){
      trimmedWhat[i]=L1
    } else if( (What[i]>L1) & (What[i]<=L2) ){
      trimmedWhat[i]=What[i]
    } else{
      trimmedWhat[i]=L2
    }
  }
  trimmedWhatDiagInverse=diag(1/trimmedWhat)
  Sigmahat=(1/N)*t(Atrain)%*%trimmedWhatDiagInverse%*%Atrain
  
  # step 3: CLIME using flare R package
  climeout=sugm(Sigmahat, lambda=climemu, method="clime")
  Thetahat=climeout$icov1[[1]]
  
  # step 4: one-step estimator
  betatilde=betahat+(1/N)*Thetahat%*%t(Atrain)%*%trimmedWhatDiagInverse%*%(btrain-Atrain%*%betahat)
  # output 
  return(list(betahat=betahat, betatilde=betatilde, trimmedWhat=trimmedWhat, Thetahat=Thetahat))
}




# one step estimator using flare R package in CLIME method for Model 4
onestepestimator4.flare=function(data, lassolambda, climemu, L1, L2){
  Atrain=data$X  
  btrain=data$Y
  N=nrow(Atrain)
  # step 1: Lasso 
  lasso=glmnet(Atrain, btrain, alpha=1, lambda=lassolambda, intercept = FALSE) # alpha=1 is lasso
  betahat=as.vector(coef(lasso)[-1])
  
  # step 2: set W hat and Sigma hat based on betahat and each model
  What=(1/16)*((Atrain%*%betahat)^2)
  trimmedWhat=rep(0, N)
  for (i in 1:N){
    if (What[i]<=L1){
      trimmedWhat[i]=L1
    } else if( (What[i]>L1) & (What[i]<=L2) ){
      trimmedWhat[i]=What[i]
    } else{
      trimmedWhat[i]=L2
    }
  }
  trimmedWhatDiagInverse=diag(1/trimmedWhat)
  Sigmahat=(1/N)*t(Atrain)%*%trimmedWhatDiagInverse%*%Atrain
  
  # step 3: CLIME using flare R package
  climeout=sugm(Sigmahat, lambda=climemu, method="clime")
  Thetahat=climeout$icov1[[1]]
  
  # step 4: one-step estimator
  betatilde=betahat+(1/N)*Thetahat%*%t(Atrain)%*%trimmedWhatDiagInverse%*%(btrain-Atrain%*%betahat)
  # output 
  return(list(betahat=betahat, betatilde=betatilde, trimmedWhat=trimmedWhat, Thetahat=Thetahat))
}









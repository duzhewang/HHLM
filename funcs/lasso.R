#-------------------------------------------------------------------------------#
# Last updated date: 7/10/2020
# Author: Duzhe Wang (dwang282@wisc.edu)
#--------------------------------------------------------------------------------#


library(glmnet)   ## lasso algorithm 

#---------tuning parameter selection--------------#

# Criteria: choose lambda with the smallest l2 estimation error 
lassotuning.l2=function(data, truebeta, lambdaseq){
  # data: tuning dataset  
  # truebeta: true regression coefficients 
  # lambdaseq: a sequence of lambda values 
  
  l2estimationerror=rep(NA, length(lambdaseq))   # l2 estimation error for each tuning parameter 
  coefsupplength=rep(NA, length(lambdaseq))    # length of support of estimated beta 
  coefmat=matrix(NA, nrow=length(lambdaseq), ncol=ncol(data$X))
  A=data$X
  b=data$Y
  for (i in 1:length(lambdaseq)){
      lasso=glmnet(A, b, alpha=1, lambda=lambdaseq[i], intercept = FALSE) # alpha=1 is lasso
      lassocoef=as.vector(coef(lasso)[-1])
      coefmat[i, ]=lassocoef
      coefsupplength[i]=length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
  }
  bestlambda=lambdaseq[which(l2estimationerror== min(l2estimationerror), arr.ind = TRUE)]  # best lambda
  bestl2estimationerror=min(l2estimationerror)
  return(list(bestlambda=bestlambda, bestl2estimationerror=bestl2estimationerror,
              l2estimationerror=l2estimationerror, coefsupplength=coefsupplength, 
              coefmat=coefmat)) 
}


# Criteria: we choose a lambda, such that 
# 1. the sparsity of the corresponding estimated beta is equal to the true sparsity and if 
#    there are multiple such lambdas, then choose the one with the smallest l2 estimation error
# 2. if there is no lambda such that the corresponding estimated beat is equal to the true sparsity, then we choose
#    the one which sparisty is the closest to the true sparsity and has the smallest l2 estimation error

lassotuning.vs=function(data, truebeta, lambdaseq){
  # data: tuning dataset  
  # truebeta: true regression coefficients 
  # lambdaseq: a sequence of lambda values 
  
  s=length(which(truebeta!=0))
  l2estimationerror=rep(NA, length(lambdaseq))   # l2 estimation error for each tuning parameter 
  coefsupplength=rep(NA, length(lambdaseq))
  coefmat=matrix(NA, nrow=length(lambdaseq), ncol=ncol(data$X))
  for (i in 1:length(lambdaseq)){
    A=data$X  
    b=data$Y  
    lasso=glmnet(A, b, alpha=1, lambda=lambdaseq[i], intercept = FALSE) # alpha=1 is lasso
    lassocoef=as.vector(coef(lasso)[-1])
    coefmat[i, ]=lassocoef
    coefsupplength[i]=length(which(lassocoef!=0))
    l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
  }
  
  # choose bestlambda
  if (any(coefsupplength==s)=="TRUE"){
    index=which(coefsupplength==s, arr.ind=TRUE)
    bestlambda=lambdaseq[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
    bestlambda.shat=coefsupplength[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
    bestlambda.l2estimationerror=l2estimationerror[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
  }else{
    index=which(abs(coefsupplength-s)==min(abs(coefsupplength-s)), arr.ind=TRUE)
    bestlambda=lambdaseq[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
    bestlambda.shat=coefsupplength[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
    bestlambda.l2estimationerror=l2estimationerror[which(l2estimationerror==min(l2estimationerror[index]), arr.ind = TRUE)]
  }
  
  return(list(bestlambda=bestlambda, bestlambda.shat=bestlambda.shat, bestlambda.l2estimationerror=bestlambda.l2estimationerror, 
              l2estimationerror=l2estimationerror, coefsupplength=coefsupplength, coefmat=coefmat)) 
}


# Criteria: choose lambda with the smallest mean squared prediction error using K-fold CV
lassotuning.cv=function(data, K, truebeta, lambdaseq){
  # data: dataset which is used to do cross validation
  # K: number of folds
  # lambdaseq: a sequence of lambda values

  N=nrow(data$X)                                   # sample size
  len=N/K                                          # length of each fold
  index=sample(1:N)                                # shuffle the dataset
  SoS=rep(NA, K)                                   # sum of squares
  MSE=rep(NA, length(lambdaseq))                   # mean squared error 
  for (i in 1:length(lambdaseq)){
    for (k in 1:K){
      print(paste("Current lambda:", lambdaseq[i], "Current fold:", k))
      Atrain=data$X[-index[((k-1)*len+1):(k*len)], ]  # training set
      btrain=data$Y[-index[((k-1)*len+1):(k*len)]]    # training set
      Avali=data$X[index[((k-1)*len+1):(k*len)], ]    # validation set
      bvali=data$Y[index[((k-1)*len+1):(k*len)]]      # validation set
      lasso=glmnet(Atrain, btrain, alpha=1, lambda=lambdaseq[i], intercept = FALSE) # alpha=1 is lasso
      lassocoef=as.vector(coef(lasso)[-1])
      SoS[k]=sum((bvali-Avali%*%lassocoef)^2)  ## sum of squares
    }
    MSE[i]=sum(SoS)/N
  }
  bestlambda=lambdaseq[which(MSE== min(MSE), arr.ind = TRUE)]  # best lambda
  A=data$X
  b=data$Y
  wholelasso=glmnet(A, b, alpha=1, lambda=bestlambda, intercept = FALSE) # alpha=1 is lasso
  wholelassocoef=as.vector(coef(wholelasso)[-1])
  wholelassocoefsupplength=length(which(wholelassocoef!=0))
  wholelassol2estimationerror=sqrt(sum((wholelassocoef-truebeta)^2))
  return(list(bestlambda=bestlambda,
              wholelassocoefsupplength=wholelassocoefsupplength,
              wholelassol2estimationerror=wholelassol2estimationerror,
              wholelassocoef=wholelassocoef,
              MSE=MSE))
}





# Criteria: AIC 
lassotuning.AIC=function(data, truebeta, lambdaseq, model, L1, L2){
  # data: tuning dataset  
  # truebeta: true regression coefficients 
  # lambdaseq: a sequence of lambda values 
  # model: indicate which model (model 1, 2, 3 and 4)
  # L1
  # L2
  
  l2estimationerror=rep(NA, length(lambdaseq))   # l2 estimation error for each tuning parameter 
  coefsupplength=rep(NA, length(lambdaseq))    # length of support of estimated beta 
  coefmat=matrix(NA, nrow=length(lambdaseq), ncol=ncol(data$X))
  AIC=rep(NA, length(lambdaseq))
  A=data$X
  b=data$Y
  N=nrow(data$X)  
  for (i in 1:length(lambdaseq)){
    lasso=glmnet(A, b, alpha=1, lambda=lambdaseq[i], intercept = FALSE) # alpha=1 is lasso
    lassocoef=as.vector(coef(lasso)[-1])
    if (model==1){
      W=(1/25)*exp(0.5*abs(A%*%lassocoef))
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L2){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      AIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+2*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else if (model==2){
      W=(1/50)*exp(0.05*(A%*%lassocoef)^2)
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L2){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      AIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+2*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else if (model==3){
      W=0.25*abs(A%*%lassocoef)
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L1){
          trimmedW[i]=L1
        } else if( (W[m]>L1) & (W[m]<=L2) ){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      AIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+2*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else{
      W=(1/16)*(A%*%lassocoef)^2
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L1){
          trimmedW[m]=L1
        } else if( (W[m]>L1) & (W[m]<=L2) ){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      AIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+2*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    }
  }
  bestlambda=lambdaseq[which(AIC== min(AIC), arr.ind = TRUE)]  # best lambda
  bestAIC=min(AIC)
  bestlambda.supplength=coefsupplength[which(AIC== min(AIC), arr.ind = TRUE)]
  bestlambda.l2estimationerror=l2estimationerror[which(AIC== min(AIC), arr.ind = TRUE)]
  return(list(bestlambda=bestlambda, bestAIC=bestAIC, bestlambda.supplength=bestlambda.supplength,
              bestlambda.l2estimationerror=bestlambda.l2estimationerror,
              AIC=AIC,l2estimationerror=l2estimationerror, coefsupplength=coefsupplength)) 
}



# Criteria: BIC 
lassotuning.BIC=function(data, truebeta, lambdaseq, model, L1, L2){
  # data: tuning dataset  
  # truebeta: true regression coefficients 
  # lambdaseq: a sequence of lambda values 
  # model: indicate which model (model 1, 2, 3 and 4)
  # L1
  # L2
  
  l2estimationerror=rep(NA, length(lambdaseq))   # l2 estimation error for each tuning parameter 
  coefsupplength=rep(NA, length(lambdaseq))    # length of support of estimated beta 
  coefmat=matrix(NA, nrow=length(lambdaseq), ncol=ncol(data$X))
  BIC=rep(NA, length(lambdaseq))
  A=data$X
  b=data$Y
  N=nrow(data$X)  
  for (i in 1:length(lambdaseq)){
    lasso=glmnet(A, b, alpha=1, lambda=lambdaseq[i], intercept = FALSE) # alpha=1 is lasso
    lassocoef=as.vector(coef(lasso)[-1])
    if (model==1){
      W=(1/25)*exp(0.5*abs(A%*%lassocoef))
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L2){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      BIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+log(N)*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else if (model==2){
      W=(1/50)*exp(0.05*(A%*%lassocoef)^2)
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L2){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      BIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+log(N)*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else if (model==3){
      W=0.25*abs(A%*%lassocoef)
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L1){
          trimmedW[i]=L1
        } else if( (W[m]>L1) & (W[m]<=L2) ){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      BIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+log(N)*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    } else{
      W=(1/16)*(A%*%lassocoef)^2
      trimmedW=rep(0, N)
      for (m in 1:N){
        if (W[m]<=L1){
          trimmedW[m]=L1
        } else if( (W[m]>L1) & (W[m]<=L2) ){
          trimmedW[m]=W[m]
        } else{
          trimmedW[m]=L2
        }
      }
      BIC[i]=t(b-A%*%lassocoef)%*%diag(1/trimmedW)%*%(b-A%*%lassocoef)+sum(log(trimmedW))+log(N)*length(which(lassocoef!=0))
      l2estimationerror[i]=sqrt(sum((lassocoef-truebeta)^2))
      coefsupplength[i]=length(which(lassocoef!=0))
    }
  }
  bestlambda=lambdaseq[which(BIC== min(BIC), arr.ind = TRUE)]  # best lambda
  bestBIC=min(BIC)
  bestlambda.supplength=coefsupplength[which(BIC== min(BIC), arr.ind = TRUE)]
  bestlambda.l2estimationerror=l2estimationerror[which(BIC== min(BIC), arr.ind = TRUE)]
  return(list(bestlambda=bestlambda, bestBIC=bestBIC, bestlambda.supplength=bestlambda.supplength,
              bestlambda.l2estimationerror=bestlambda.l2estimationerror,
              BIC=BIC,l2estimationerror=l2estimationerror, coefsupplength=coefsupplength)) 
}






#-----------Lasso estimation---------------------------#

lassoest=function(data, lambda, truebeta){
  Atrain=data$X  
  btrain=data$Y
  lasso=glmnet(Atrain, btrain, alpha=1, lambda=lambda, intercept = FALSE) # alpha=1 is lasso
  lassocoef=as.vector(coef(lasso)[-1])
  l2error=sqrt(sum((lassocoef-truebeta)^2))
  l1error=sum(abs(lassocoef-truebeta))
  Shatlength=length(which(lassocoef!=0))
  MSE=sum((btrain-Atrain%*%lassocoef)^2)/nrow(Atrain)
  return(list(lassocoef=lassocoef, l2error=l2error, l1error=l1error, Shatlength=Shatlength, MSE=MSE)) 
}





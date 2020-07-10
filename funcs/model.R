#-------------------------------------------------------------------------------#
# Last updated date: 7/10/2020
# Author: Duzhe Wang (dwang282@wisc.edu)
#--------------------------------------------------------------------------------#

library(MASS)  # mvrnorm function



#-------------------------------
#  Generate Sigma_x
#-------------------------------

Sigmax=function(rho, n, s){
  Sigma=matrix(0, nrow=n, ncol=n)
  for (i in 1:s){
    for (j in 1:s){
      Sigma[i, j]=rho^{abs(i-j)}
    }
  }
  for (i in (s+1):n){
    for (j in (s+1):n){
      Sigma[i, j]=rho^{abs(i-j)}
    }
  }
  
  return(Sigma)
}



#------------------------------------------------------------
#     Model 0: homoscedastic linear model 
#------------------------------------------------------------
HHR0<-function(N, n, Sigmax, truebeta, L2){
  # n: dimension of parameter, i.e., number of columns for design matrix X. 
  # N: sample size, i.e., number of rows for design matrix.
  # Sigmax: covariance matrix
  # truebeta: true regression coefficients
  # L2: sd 
  
  # generate X
  X<-mvrnorm(N, rep(0, n), Sigmax) # mean 0 and covariance Sigmax
  
  # generate epsilon
  epsilon=rnorm(N, mean=0, sd=L2)
  
  # generate Y
  Y<-X%*%truebeta+epsilon
  
  # output 
  data<-list(X=X, Y=Y, epsilon=epsilon)
  return(data)
}



#------------------------------------------------------------
#     Model 1
#------------------------------------------------------------
HHR1<-function(N, n, Sigmax, truebeta, L2){
  # n: dimension of parameter, i.e., number of columns for design matrix X. 
  # N: sample size, i.e., number of rows for design matrix.
  # Sigmax: covariance matrix
  # truebeta: true regression coefficients
  # L2: upper bound of variance 
  

  # generate X
  X<-mvrnorm(N, rep(0, n), Sigmax) # mean 0 and covariance Sigmax
  
  # generate variance
  W=(1/25)*exp(0.5*abs(X%*%truebeta))
  
  # generate epsilon
  epsilon=rep(0, N)
  trimmedW=rep(0, N)
  for (i in 1:N){
    if (W[i]<=L2){
      trimmedW[i]=W[i]
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else{
      trimmedW[i]=L2
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    }
  }
  
  # generate Y
  Y<-X%*%truebeta+epsilon
  
  # output 
  data<-list(X=X, Y=Y, W=W, trimmedW=trimmedW, epsilon=epsilon)
  return(data)
}


#------------------------------------------------------------
#     Model 2
#------------------------------------------------------------
HHR2<-function(N, n, Sigmax, truebeta, L2){
  # n: dimension of parameter, i.e., number of columns for design matrix X. 
  # N: sample size, i.e., number of rows for design matrix.
  # Sigmax: covariance matrix
  # truebeta: true regression coefficients
  # L2: upper bound of variance 
  
  
  # generate X
  X<-mvrnorm(N, rep(0, n), Sigmax) # mean 0 and covariance Sigmax
  
  # generate variance
  W=(1/50)*exp(0.05*(X%*%truebeta)^2)
  
  # generate epsilon
  epsilon=rep(0, N)
  trimmedW=rep(0, N)
  for (i in 1:N){
    if (W[i]<=L2){
      trimmedW[i]=W[i]
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else{
      trimmedW[i]=L2
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    }
  }
  
  # generate Y
  Y<-X%*%truebeta+epsilon
  
  # output 
  data<-list(X=X, Y=Y, W=W, trimmedW=trimmedW, epsilon=epsilon)
  return(data)
}


#------------------------------------------------------------
#     Model 3
#------------------------------------------------------------
HHR3<-function(N, n, Sigmax, truebeta, L1, L2){
  # n: dimension of parameter, i.e., number of columns for design matrix X. 
  # N: sample size, i.e., number of rows for design matrix.
  # Sigmax: covariance matrix
  # truebeta: true regression coefficients
  # L1: lower bound of variance 
  # L2: upper bound of variance 
  
  
  # generate X
  X<-mvrnorm(N, rep(0, n), Sigmax) # mean 0 and covariance Sigmax
  
  # generate variance
  W=0.25*abs(X%*%truebeta)
  
  # generate epsilon
  epsilon=rep(0, N)
  trimmedW=rep(0, N)
  for (i in 1:N){
    if (W[i]<=L1){
      trimmedW[i]=L1
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else if( (W[i]>L1) & (W[i]<=L2) ){
      trimmedW[i]=W[i]
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else{
      trimmedW[i]=L2
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    }
  }
  
  # generate Y
  Y<-X%*%truebeta+epsilon
  
  # output 
  data<-list(X=X, Y=Y, W=W, trimmedW=trimmedW, epsilon=epsilon)
  return(data)
}


#------------------------------------------------------------
#     Model 4
#------------------------------------------------------------
HHR4<-function(N, n, Sigmax, truebeta, L1, L2){
  # n: dimension of parameter, i.e., number of columns for design matrix X. 
  # N: sample size, i.e., number of rows for design matrix.
  # Sigmax: covariance matrix
  # truebeta: true regression coefficients
  # L1: lower bound of variance 
  # L2: upper bound of variance 
  
  
  # generate X
  X<-mvrnorm(N, rep(0, n), Sigmax) # mean 0 and covariance Sigmax
  
  # generate variance
  W=(1/16)*(X%*%truebeta)^2
  
  # generate epsilon
  epsilon=rep(0, N)
  trimmedW=rep(0, N)
  for (i in 1:N){
    if (W[i]<=L1){
      trimmedW[i]=L1
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else if( (W[i]>L1) & (W[i]<=L2) ){
      trimmedW[i]=W[i]
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    } else{
      trimmedW[i]=L2
      epsilon[i]=rnorm(1, mean=0, sd=sqrt(trimmedW[i]))
    }
  }
  
  # generate Y
  Y<-X%*%truebeta+epsilon
  
  # output 
  data<-list(X=X, Y=Y, W=W, trimmedW=trimmedW, epsilon=epsilon)
  return(data)
}































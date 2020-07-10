#-------------------------------------------------------------------------------#
# Last updated date: 7/10/2020
# Author: Duzhe Wang (dwang282@wisc.edu)
#--------------------------------------------------------------------------------#

rm(list=ls())
setwd("/Users/peterwang/Desktop/HHLM/model3")
source("/Users/peterwang/Desktop/HHLM/funcs/model.R")
source("/Users/peterwang/Desktop/HHLM/funcs/lasso.R")
set.seed(2020)

#---------------Step 1: Set parameters--------------------------------#
Nvec=c(30, 60, 90, 120, 150, 180, 210, 240)                                        # sample size
n=150                                                                              # ambient dimension 
nsim=200                                                                           # number of simulation times 
rho=0.5                 
lambdaseqconstant=seq(from=0.700, to=0.900, 0.05)
L1=0.05
L2=5

# true beta
s=6
Supp=c(1,2,3, 4, 5, 6)
truebeta=rep(0, n)
truebeta[1:s]=c(3, 4, 3, 1.5, 2, 1.5)

# generate covariance matrix 
Sigmax=Sigmax(rho, n, s)
maxeigen=max(eigen(Sigmax)$values)   # maxeigen is around 3.
#-------------------------------------------------------------------#

#-------------Step 2: SIMULATION-------------------------#

# each row records N, lambda, simu iter, l2 estimation error, l1 estimation error, sparsity of estimaed beta and MSE
lassoestiresults=matrix(NA, nrow=length(Nvec)*nsim, ncol=7*length(lambdaseqconstant))

# each row records estimated beta in each iteration
lassocoefmat=matrix(NA, nrow=length(Nvec)*nsim, ncol=n*length(lambdaseqconstant))


for (k in 1:length(lambdaseqconstant)){
    for (i in 1:length(Nvec)){
        lambda=lambdaseqconstant[k]*sqrt(maxeigen*L2*log(n)/Nvec[i])
        for (iter in 1:nsim){
        print(paste("lambda constant:",lambdaseqconstant[k],"N:", Nvec[i], "iteration:", iter))
        trainingdata=HHR3(Nvec[i], n, Sigmax, truebeta, L1, L2)
        lassosimu=lassoest(trainingdata, lambda=lambda, truebeta)
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+1]=Nvec[i]
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+2]=lambda
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+3]=iter
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+4]=lassosimu$l2error
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+5]=lassosimu$l1error
        lassoestiresults[(i-1)*nsim+iter,(k-1)*7+6]=lassosimu$Shatlength
        lassoestiresults[(i-1)*nsim+iter,k*7]=lassosimu$MSE
        lassocoefmat[(i-1)*nsim+iter, ((k-1)*n+1):(k*n)]=lassosimu$lassocoef
        }
      }
}

#--------------------------------------------------------------#

#---------------Step 3: Analysis of results--------------------#

## Estimation 
avgl2=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
avgl1=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
avgMSE=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
sdl2=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
sdl1=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
sdMSE=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))


for (k in 1:length(lambdaseqconstant)){
    lassoestiresults.df=as.data.frame(lassoestiresults[,((k-1)*7+1):(k*7)])
    colnames(lassoestiresults.df)=c("N", "lambda", "iter", "l2error", "l1error", "betahatsparsity", "MSE")
    for (i in 1:length(Nvec)){
         avgl2[i,k]=mean(lassoestiresults.df$l2error[lassoestiresults.df$N==Nvec[i]])
         avgl1[i,k]=mean(lassoestiresults.df$l1error[lassoestiresults.df$N==Nvec[i]])
         avgMSE[i,k]=mean(lassoestiresults.df$MSE[lassoestiresults.df$N==Nvec[i]])
         sdl2[i,k]=sd(lassoestiresults.df$l2error[lassoestiresults.df$N==Nvec[i]])
         sdl1[i,k]=sd(lassoestiresults.df$l1error[lassoestiresults.df$N==Nvec[i]])
         sdMSE[i,k]=sd(lassoestiresults.df$MSE[lassoestiresults.df$N==Nvec[i]])
    }
}


## Variable selection
Supphat=matrix(NA, nrow=length(Nvec)*nsim, ncol=n*length(lambdaseqconstant)) # support set: 1 represents nonzero  
shat=matrix(NA, nrow=length(Nvec)*nsim, ncol=length(lambdaseqconstant))  # size of support set 
successindex=matrix(NA, nrow=length(Nvec)*nsim, ncol=length(lambdaseqconstant))
recall=matrix(NA, nrow=length(Nvec)*nsim, ncol=length(lambdaseqconstant))
precision=matrix(NA, nrow=length(Nvec)*nsim, ncol=length(lambdaseqconstant))
FDR=matrix(NA, nrow=length(Nvec)*nsim, ncol=length(lambdaseqconstant))
successprob=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
avgrecall=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
avgprecision=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))
avgFDR=matrix(NA, nrow=length(Nvec), ncol=length(lambdaseqconstant))

for (k in 1:length(lambdaseqconstant)){
    for(i in 1:(length(Nvec)*nsim)){
       Supphat[i, ((k-1)*n+1):(k*n)]=as.integer(lassocoefmat[i,((k-1)*n+1):(k*n)]!=0)
       shat[i,k]=sum(Supphat[i, ((k-1)*n+1):(k*n)])
       supportset=which(lassocoefmat[i,((k-1)*n+1):(k*n)]!=0)
       recall[i,k]=length(intersect(supportset, Supp))/s
       precision[i,k]=length(intersect(supportset, Supp))/shat[i,k]
       FDR[i,k]=1-precision[i,k]
       if (isTRUE(shat[i,k]!=s)){
       successindex[i,k]=0  # variable selection fail
       }else if (all(supportset%in%Supp)=="FALSE"){
       successindex[i,k]=0  # variable selection fail
       } else{
       successindex[i,k]=1
       }
    }
}

 

for (k in 1:length(lambdaseqconstant)){ 
    for(i in 1:length(Nvec)){
    successprob[i, k]=sum(successindex[((i-1)*nsim+1):(i*nsim), k])/nsim
    avgrecall[i, k]=mean(recall[((i-1)*nsim+1):(i*nsim), k])
    avgprecision[i, k]=mean(precision[((i-1)*nsim+1):(i*nsim), k]) 
    avgFDR[i, k]=mean(FDR[((i-1)*nsim+1):(i*nsim), k])
    }
}



#-------------------------------------------------------------#


#-------------Print and save results--------------------------#
finalresults=list(avgl2=avgl2, avgl1=avgl1, avgMSE=avgMSE,
                  sdl2=sdl2, sdl1=sdl1, sdMSE=sdMSE,
                  successprob=successprob, avgrecall=avgrecall, 
                  avgprecision=avgprecision, avgFDR=avgFDR)

#-------------------------------------------------------------#


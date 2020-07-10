#-------------------------------------------------------------------------------#
# Last updated date: 7/10/2020
# Author: Duzhe Wang (dwang282@wisc.edu)
#--------------------------------------------------------------------------------#

rm(list=ls())
setwd("/Users/peterwang/Desktop/HHLM/model1")
source("/Users/peterwang/Desktop/HHLM/funcs/model.R")
source("/Users/peterwang/Desktop/HHLM/funcs/lasso.R")
source("/Users/peterwang/Desktop/HHLM/funcs/onestepestimator.R")
library(ggplot2)
library(latex2exp)
set.seed(2020)

#---------------Step 1: Set parameters--------------------------------#
Nvec=120
n=150                                                                              # ambient dimension 
nsim=200                                                                           # number of simulation times 
rho=0.5                 
L1=1/25
L2=5

# tuning parameter of CLIME
museqconstant=0.0001
  
# true beta
s=6
Supp=c(1, 2, 3, 4, 5, 6)
truebeta=rep(0, n)
truebeta[1:s]=c(3, 4, 3, 1.5, 2, 1.5)

# generate covariance matrix 
Sigmax=Sigmax(rho, n, s)
maxeigen=max(eigen(Sigmax)$values)   # maxeigen is around 3.
#-------------------------------------------------------------------#

#-------------Step 2: SIMULATION------------------------------------#
betahatmat=matrix(NA, nrow=length(Nvec)*nsim, ncol=n*length(museqconstant))
betatildemat=matrix(NA, nrow=length(Nvec)*nsim, ncol=n*length(museqconstant))
Thetahatlist=vector(mode = "list", length = nsim*length(Nvec)*length(museqconstant))


for (k in 1:length(museqconstant)){
    for (i in 1:length(Nvec)){
        lassolambda=0.75*sqrt(maxeigen*L2*log(n)/Nvec[i])
        climemu=museqconstant[k]*((1/L1)*sqrt(log(n)/Nvec[i])+s*sqrt(L2*log(n)/Nvec[i]))
        for (iter in 1:nsim){
            print(paste("muconstant:", museqconstant[k], "N:", Nvec[i], "iteration:", iter))
            trainingdata=HHR1(Nvec[i], n, Sigmax, truebeta, L2)
            onestepsimu=onestepestimator1.flare(trainingdata, lassolambda, climemu, L2)
            betahatmat[(i-1)*nsim+iter,((k-1)*n+1):(k*n)]=onestepsimu$betahat
            betatildemat[(i-1)*nsim+iter,((k-1)*n+1):(k*n)]=onestepsimu$betatilde
            Thetahatlist[[((k-1)*nsim*length(Nvec)+(i-1)*nsim+iter)]]=onestepsimu$Thetahat
        }
    }
}
#------------------------------------------------------------------#


#---------------------------Analyze results-------------------------#

# Visualization for the paper: choose beta5 and beta7 as an example

# betatilde5
data=sqrt(Nvec)*(betatildemat[, 5]-rep(truebeta[5], nsim))
pvaluetilde5=shapiro.test(data)$p.value
df5=data.frame(data=data)

ggplot(df5, aes(sample=data))+stat_qq(size=1.5)+stat_qq_line(col="black", size=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", size=0.5, col="red")+
  geom_hline(yintercept = 0, linetype="dashed", size=0.5, col="red")+
  theme(plot.title = element_text(hjust=0.5))+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  ggtitle("Normal Q-Q plot")


# betatilde7
data=sqrt(Nvec)*(betatildemat[, 7]-rep(truebeta[7], nsim))
pvaluetilde7=shapiro.test(data)$p.value
df7=data.frame(data=data)

ggplot(df7, aes(sample=data))+stat_qq(size=1.5)+stat_qq_line(col="black", size=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", size=0.5, col="red")+
  geom_hline(yintercept = 0, linetype="dashed", size=0.5, col="red")+
  theme(plot.title = element_text(hjust=0.5))+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  ggtitle("Normal Q-Q plot")



# betahat5
data=sqrt(Nvec)*(betahatmat[, 5]-rep(truebeta[5], nsim))
pvaluehat5=shapiro.test(data)$p.value
df5=data.frame(data=data)

ggplot(df5, aes(sample=data))+stat_qq(size=1.5)+stat_qq_line(col="black", size=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", size=0.5, col="red")+
  geom_hline(yintercept = 0, linetype="dashed", size=0.5, col="red")+
  theme(plot.title = element_text(hjust=0.5))+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  ggtitle("Normal Q-Q plot")



# betahat7
data=sqrt(Nvec)*(betahatmat[, 7]-rep(truebeta[7], nsim))
pvaluehat7=shapiro.test(data)$p.value
df7=data.frame(data=data)

ggplot(df7, aes(sample=data))+stat_qq(size=1.5)+stat_qq_line(col="black", size=0.25)+
  geom_vline(xintercept = 0, linetype="dashed", size=0.5, col="red")+
  geom_hline(yintercept = 0, linetype="dashed", size=0.5, col="red")+
  theme(plot.title = element_text(hjust=0.5))+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  ggtitle("Normal Q-Q plot")

#------------------------------------------------------------#





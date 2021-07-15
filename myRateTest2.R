# myRateTestWTCCC.R
#Written by Emily Winn
# Last updated 6/30/2021

# This file we are going to use to test the new RATE functions on WTccc data. We will call on 
# RATE2.R to make these. First I need to import data.

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(corpcor)
library(doParallel)
library(MASS)
library(Matrix)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(svd)

### Load in the RATE R functions ### (Path set by user in both)
#source("C:/Users/etwin/git_repos/RATE/Software/RATE.R") #Changing path for etwin PC.
source("RATE.R")

### Load in the C++ BAKR functions ###
#sourceCpp("C:/Users/etwin/git_repos/BAKR-master/BAKR-master/Rcpp/BAKRGibbs.cpp")
sourceCpp("BAKRGibbs.cpp")

# Importing the data for the WTCCC. It's a matrix called "X.select"
#load("C:/Users/etwin/Downloads/WTCCC_Sub.RData")
#load("Z:/RATE/RATE-master/RATE-master/Software/WTCCC_Sub.RData")
load("WTCCC_Sub.RData")

n = dim(X.select)[1] #Sample size
p = dim(X.select)[2] #Number of markers or genes
pve=0.6; rho=1;

#Now need to simulate the y vector, which is the same code as before.
#Number of causal snps
set.seed(11151990)
ncausal = 3
s=c(1:ncausal)

#Marginal Effects Only
Xmarginal=X.select[,s]
beta1=rep(1,ncausal)
y_marginal=c(Xmarginal%*%beta1)
beta1=beta1*sqrt(pve*rho/var(y_marginal))
y_marginal=Xmarginal%*%beta1

#Pairwise Epistatic Effects
Xepi=cbind(X.select[,s[1]]*X.select[,s[3]],X.select[,s[2]]*X.select[,s[3]])
beta2=c(1,1)
y_epi=c(Xepi%*%beta2)
beta2=beta2*sqrt(pve*(1-rho)/var(y_epi))
y_epi=Xepi%*%beta2

#Error Terms
y_err=rnorm(n)
y_err=y_err*sqrt((1-pve)/var(y_err))

y=c(y_marginal+y_epi+y_err); #Full Model
colnames(X.select) = paste("SNP",1:ncol(X.select),sep="")

######################################################################################
######################################################################################
######################################################################################

### Running the Bayesian Gaussian Process (GP) Regression Model ###

### Create the Nonlinear Covariance Matrix with the Gaussian Kernel ###
#This function takes on two arguments:
#(1) The Genotype matrix X. This matrix should be fed in as a pxn matrix. That is, predictor
#are the rows and subjects/patients/cell lines are the columns.
#(2) The bandwidth (also known as a smoothing parameter or lengthscale) h. For example, the 
#Gaussian kernel can be specified as k(u,v) = exp{||uâv||^2/2h^2}.

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
Kn = GaussKernel(t(X.select)); diag(Kn)=1 # 

### Center and Scale K_tilde ###
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Kn=M%*%Kn%*%M
Kn=Kn/mean(diag(Kn))

### Gibbs Sampler ###
sigma2 = 1e-3
fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
fhat.rep = rmvnorm(1e4,fhat,Kn - Kn %*% solve(Kn+diag(sigma2,n),Kn))

### Compute the KL Divergence to find Marginal Importance ###
cores = detectCores()
registerDoParallel(cores=cores)
print("Starting RATE")
### Run the RATE Function ###
nl = NULL
start = Sys.time()
res = RATE(X=X.select,f.draws=fhat.rep,snp.nms = colnames(X.select))#,cores = cores)
end = Sys.time()
rate1_time = end-start
print(rate1_time)

save(res, file = "rate1_WTCCC.Rdata")

#Now load next RATE function
source("RATE2.R")
### Run the RATE Function ###
print("Starting RATE2")
nl = NULL
start = Sys.time()
res2 = RATE(X=X.select,f.draws=fhat.rep,snp.nms = colnames(X.select),cores = cores)
end = Sys.time()
rate2_time = end-start
print(rate2_time)

save(res2, file =  "rate2_WTCCC.Rdata")
save.image(file="rate_WTCCC_test_Take2.Rdata")
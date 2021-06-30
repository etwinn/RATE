# newRate.R
# Written by Emily Winn
# Last updated 6/30/2021

# This file we are going to use to simulate and test the new RATE functions. We will call on 
# RATE2.R to make these. First I need to simulate data.

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

### Load in the RATE R functions ### (Path set by user in both)
source("C:/Users/etwin/git_repos/RATE/Software/RATE2.R") #Changing path for etwin PC.

### Load in the C++ BAKR functions ###
sourceCpp("C:/Users/etwin/git_repos/BAKR-master/BAKR-master/Rcpp/BAKRGibbs.cpp")

# Data simulation - want 10 variables, 3 pairwise interactions. Using Lorin's code from 
# Centrality_Tutorial.R with some changes to make it our scale. Still going to use Gaussian process.
### Set the random seed to reproduce research ###
set.seed(11151990)

n = 2e3; p = 10; pve=0.6; rho=1;

### Define the Number of Causal SNPs
ncausal = 3

### Simulate Synthetic Data ###
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
s=c(8:10)

#Marginal Effects Only
Xmarginal=X[,s]
beta1=rep(1,ncausal)
y_marginal=c(Xmarginal%*%beta1)
beta1=beta1*sqrt(pve*rho/var(y_marginal))
y_marginal=Xmarginal%*%beta1

#Pairwise Epistatic Effects
Xepi=cbind(X[,s[1]]*X[,s[3]],X[,s[2]]*X[,s[3]])
beta2=c(1,1)
y_epi=c(Xepi%*%beta2)
beta2=beta2*sqrt(pve*(1-rho)/var(y_epi))
y_epi=Xepi%*%beta2

#Error Terms
y_err=rnorm(n)
y_err=y_err*sqrt((1-pve)/var(y_err))

y=c(y_marginal+y_epi+y_err); #Full Model
colnames(X) = paste("SNP",1:ncol(X),sep="")

######################################################################################
######################################################################################
######################################################################################

### Running the Bayesian Gaussian Process (GP) Regression Model ###

### Create the Nonlinear Covariance Matrix with the Gaussian Kernel ###
#This function takes on two arguments:
#(1) The Genotype matrix X. This matrix should be fed in as a pxn matrix. That is, predictor
#are the rows and subjects/patients/cell lines are the columns.
#(2) The bandwidth (also known as a smoothing parameter or lengthscale) h. For example, the 
#Gaussian kernel can be specified as k(u,v) = exp{||u−v||^2/2h^2}.

n = dim(X)[1] #Sample size
p = dim(X)[2] #Number of markers or genes

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
Kn = GaussKernel(t(X)); diag(Kn)=1 # 

### Center and Scale K_tilde ###
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Kn=M%*%Kn%*%M
Kn=Kn/mean(diag(Kn))

### Gibbs Sampler ###
sigma2 = 1e-3
fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
fhat.rep = rmvnorm(1e4,fhat,Kn - Kn %*% solve(Kn+diag(sigma2,n),Kn))


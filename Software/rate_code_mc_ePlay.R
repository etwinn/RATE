# THIS IS MARYCLARE'S CODE FOR STUFF, LOOK THROUGH AND MAKE IT WORK.
# Now renaming it ePlay
# Emily is modifying on Oct 19 2021.
# rate_code_mc.R
# Same as the tutorial but with added sampling to get the g value as needed for
# running MaryClare's suggestion for RATE.


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
#sourceCpp("C:/Users/etwin/git_repos/BAKR-master/BAKR-master/Rcpp/BAKRGibbs.cpp")
sourceCpp("C:/Users/etwin/git_repos/RATE/Software/BAKRGibbs.cpp")

# Data simulation - want 10 variables, 3 pairwise interactions. Using Lorin's code from 
# Centrality_Tutorial.R with some changes to make it our scale. Still going to use Gaussian process.
### Set the random seed to reproduce research ###
set.seed(11151990)

n = 200 # 2000; 
p = 10; pve=0.6; rho=0.7; #rho=0 makes all epistatic, 1 is all marginal

### Define the Number of Causal SNPs
ncausal = 3

### Simulate Synthetic Data ###
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)
#Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
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

### Running the Bayesian Gaussian Process (GP) Regression Model, sampling for delta_j = (g(j)-f) ###

### Create the Nonlinear Covariance Matrix with the Gaussian Kernel ###
#This function takes on two arguments:
#(1) The Genotype matrix X. This matrix should be fed in as a pxn matrix. That is, predictor
#are the rows and subjects/patients/cell lines are the columns.
#(2) The bandwidth (also known as a smoothing parameter or lengthscale) h. For example, the 
#Gaussian kernel can be specified as k(u,v) = exp{||uâv||^2/2h^2}.
#Same with cokernel for the other matrices.

### Gibbs Sampler ###
sample_size=100 # 1e4
sigma2 = 1e-3

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
B = GaussKernel(t(X)); diag(B)=1
A <- B + sigma2*diag(1,nrow=n,ncol=n)
A_svd <- svd(A)
Ainv = nearPD(A_svd$v%*%diag(1/A_svd$d, nrow=n, ncol=n)%*%t(A_svd$u))$mat
#Ainv = nearPD(solve(A))$mat
Aiy <- Ainv%*%y
BAiy <- B %*% Aiy
IAiB <- (diag(1,nrow=n, ncol=n)-Ainv%*%B)
BIAiB <- B%*%IAiB

### Find basis and kernel matrix, center and scale, Gibbs Sampler for g ###
cores = detectCores()
registerDoParallel(cores=cores)
start = Sys.time()
#Calculate Delta, which ends up being p x sample_size matrix (add dopar instead of do once off windows)
delta = ComputeESA(as.matrix(X), as.vector(Aiy), as.vector(BAiy))
#This (the parallel loop) does not work in windows but does in Linux.
#delta = foreach(k = 1:p, .combine='c', .packages=c("Rcpp", "RcppArmadillo"), .noexport=c("GaussKernel", "GaussCoKernel", "mvrnormArma"))%dopar%{
#delta = matrix(0,nrow=n, ncol=p) #nrow=sample_size
#for(k in 1:p){
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  #sourceCpp("C:/Users/etwin/git_repos/RATE/Software/RateParFunc.cpp")
 # new_X = X 
#  new_X[,k] <- new_X[,k]+1
  # MCG: Need to be careful here - the predictor only takes on values 0-2, may want to be careful
 # Cj = GaussCoKernel(t(X), t(new_X))
  
  #CtAiy <- t(Cj) %*% Aiy
  #AiC <- Ainv%*%Cj
  #del_cov = nearPD(B+BIAiB-t(Cj)%*%(AiC +2*IAiB))
  
  #del = mvrnormArma(sample_size, as.matrix(BAiy-CtAiy), as.matrix(del_cov$mat))
  #delta[, k] <- colMeans(del)
  #delt = colMeans(del)
  #delt
#}
cat("Done with delta! \n")
Sys.time()-start

### Run the RATE Function ###
rate_choice = "RATE MC"
nl = NULL
#start = Sys.time()
res = RATE_MC(X=X,beta.draws=delta, snp.nms = colnames(X),cores = cores)
#end = Sys.time()
#print(end-start)

### Get the Results ###
rates = res$RATE
DELTA = res$Delta
ESS = res$ESS



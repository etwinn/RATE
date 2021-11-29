### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Change the Library Path ###
.libPaths("/data/mukherjeelab/RLibs")

### Load in the R libraries ###
library(coda)
library(MASS)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(BGLR)
library(monomvn)
library(kernlab)
library(CompQuadForm)
library(adegenet)

### Load in the C++ functions ###
sourceCpp("/home/lac55/BAKRGibbs.cpp")

######################################################################################
######################################################################################
######################################################################################

### Set the working directory ###
setwd("/data/mukherjeelab/BAKR")

### Set the Seed for the analysis ###
set.seed(11151990)

### Number of Datasets ###
n.datasets = 100

### Call the available cores accross the computer ###
registerDoParallel(cores = detectCores())

######################################################################################
######################################################################################
######################################################################################

# LINEAR REGRESSION EXAMPLE
# -------------------------
# Data are 2,000 uncorrelated ("unlinked") single nucleotide
# polymorphisms (SNPs) with simulated genotypes. (Strategy by Peter Carbonetto)

### Set the Seed for the analysis ###
set.seed(11151990)

### Set the Results Directory ###
setwd("/home/lcrawfo1/Results/KLD/")
fn = "Add"

n.datasets = 100

### Load in the Genotypes ###
load("/home/lcrawfo1/Data/WTCCC/WTCCC_Sub.RData")
unique.snps = apply(X.select,1,function(x) length(unique(x)))
X.select = X.select[unique.snps>1,]
maf=apply(X.select, 1, mean)/2
X = t(X.select[maf>0.05,]);
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
ind = dim(X)[1]; nsnp = dim(X)[2]

### Compute the Top PCs ###
PCs = ComputePCs(X,10)

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
Kn = GaussKernel(t(X)); diag(Kn) = 1

v=matrix(1, ind, 1)
M=diag(ind)-v%*%t(v)/ind
Kn=M%*%Kn%*%M
Kn=Kn/mean(diag(Kn))

# simulation parameters
pve=0.3; rho=1; pc.var = 0; ncausal = 30
ncausal1= ncausal/6 #Set 1 of causal SNPs 
ncausal2 = ncausal-ncausal1 #Set 2 of Causal SNPs

### Create a matrix to save the results ###
Final = matrix(nrow = n.datasets,ncol = 11,byrow=TRUE)
head(Final)
colnames(Final) = c("BRR","BL","BLMM","SVM","BAKR_h5","BAKR_h2","BAKR_h1","BAKR_h.5","BAKR_h.25","BAKR_h.05","Optimal")

### Run the splits for prediction ###
for(j in 1:n.datasets){
  
  s=sample(1:nsnp,ncausal,replace = FALSE)
  
  #Select Causal SNPs
  s1=sample(s, ncausal1, replace=F)
  s2=sample(s[s%in%s1==FALSE], ncausal2, replace=F)
  
  # Generate the ground-truth regression coefficients for the variables
  # (X). Adjust the effects so that
  # the variables (SNPs) explain x percent of the variance in the
  # outcome.
  Xcausal1=X[,s1]; Xcausal2=X[,s2];
  Xepi=c()
  for(i in 1:ncausal1){
    Xepi=cbind(Xepi,Xcausal1[,i]*Xcausal2)
  }
  dim(Xepi)
  
  # Marginal Effects Only
  Xmarginal=cbind(X[,s])
  beta=rnorm(dim(X[,s])[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
  
  #Pairwise Epistatic Effects
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta
  
  ### Define the effects of the PCs ###
  beta=rnorm(dim(PCs)[2])
  y_pcs=c(PCs%*%beta)
  beta=beta*sqrt(pc.var/var(y_pcs))
  y_pcs=PCs%*%beta
  
  # error
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve-pc.var)/var(y_err))
  
  ### Simulate the Response ###
  y=c(y_marginal+y_epi+y_pcs+y_err) #Full Model
  y=(y-mean(y))/(sd(y))
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Create the training and test sets ###
  train.idc = sample(1:length(y), size=0.8*length(y), replace=FALSE)
  
  X_train = X[train.idc,]; y_train = as.numeric(y[train.idc])
  X_test = X[-train.idc,]; y_test = as.numeric(y[-train.idc])
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,p)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h1 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,2*p)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h2 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,p/2)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h3 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,p/4)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h4 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,p/20)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h5 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = ApproxGaussKernel(t(X_train),p,p*5)
  
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Find the Eigenvalue Decomposition of K ###
  evd = EigDecomp(Kn)
  
  ### Truncate the data based on the desired cumulative variance explained ###
  explained_var = cumsum(evd$lambda/sum(evd$lambda))
  q = 1:min(which(explained_var >= 0.99))
  Lambda = diag(sort(evd$lambda,decreasing = TRUE)[q]^(-1)) # Matrix of Eigenvalues
  U = evd$U[,q] # Unitary Matrix of Eigenvectors
  
  ### Define Inverse Mapping ###
  B = InverseMap(t(X_train),U)
  
  ### Set up the number of MCMC samples and burn-ins ###
  mcmc.iter = 2e3
  mcmc.burn = 1e3
  
  ### Run BAKR ### 
  Gibbs = BAKRGibbs(U,y_train,Lambda,mcmc.iter,mcmc.burn)
  
  ### Look at the Posterior Summaries ###
  theta.out = PostMean(Gibbs$theta)
  beta.out = PostBeta(B,theta.out); names(beta.out) = colnames(X)
  BAKR_pred = X_test%*%beta.out
  
  ### Get Diagnostics ###
  MSPE_BAKR_h6 = mean((y_test-BAKR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Bayesian Ridge Regression ###
  ETA = list(list(X = X_train, model="BRR"))
  
  ### Run the Gibbs Sampler ###
  reg.BRR = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)
  
  ### Get the posterior of the missing variables ###
  BRR_pred = X_test%*%reg.BRR$ETA[[1]]$b
  
  ### Get Diagnostics ###
  MSPE_BRR = mean((y_test-BRR_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Bayesian BLUP ###
  K = GetLinearKernel(t(X_train))
  ETA = list(list(K = K, model="RKHS"))
  
  ### Run the Gibbs Sampler ###
  reg.BBLUP = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)
  reg.BBLUP_b = ginv(X_train)%*%reg.BBLUP$ETA[[1]]$u
  
  ### Get the posterior of the missing variables ###
  BBLUP_pred = X_test%*%reg.BBLUP_b
  
  ### Get Diagnostics ###
  MSPE_BBLUP = mean((y_test-BBLUP_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Set up List ###
  ETA = list(list(X = X_train, model="BL"))
  
  ### Run the Gibbs Sampler ###
  reg.BL = BGLR(y=y_train, ETA=ETA, nIter=mcmc.iter, burnIn=mcmc.burn, verbose=FALSE)
  
  ### Get the posterior of the missing variables ###
  BL_pred = X_test%*%reg.BL$ETA[[1]]$b
  
  ### Find MSE and Correlations ###
  MSPE_BL = mean((y_test-BL_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### SVM Model ###
  reg.svm  = ksvm(y=y_train,x=X_train,type="nu-svr",kernel="rbfdot")
  SVM_pred = predict(reg.svm, X_test, type="response")
  
  ### Get Diagnostics ###
  MSPE_SVM = mean((y_test-SVM_pred)^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  scores = c(MSPE_BRR,MSPE_BL,MSPE_BBLUP,MSPE_SVM,MSPE_BAKR_h6,MSPE_BAKR_h2,MSPE_BAKR_h1,MSPE_BAKR_h3,MSPE_BAKR_h4,MSPE_BAKR_h5)
  names(scores) = c("BRR","BL","BLMM","SVM","BAKR_h5","BAKR_h2","BAKR_h1","BAKR_h.5","BAKR_h.25","BAKR_h.05")
  Final[j,] = c(scores,names(scores)[which(scores == min(scores))])
  
  ### Report Status ###
  cat("Completed Dataset", j, "\n", sep = " ")
}

save(Final, file = "Pred_02_Results.RData")
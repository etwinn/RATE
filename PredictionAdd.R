### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

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
library(glmnet)
library(svd)

### Load in the C++ functions ###
sourceCpp("BAKRGibbs.cpp")

######################################################################################
######################################################################################
######################################################################################

### Set the Seed for the analysis ###
set.seed(11151990)

### Number of Datasets ###
n.datasets = 100

### Call the available cores accross the computer ###
cores=detectCores()
registerDoParallel(cores = cores)

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
#setwd("/home/lcrawfo1/Results/KLD/")
fn = "Add"

n.datasets = 100

### Load in the Genotypes ###
load("~/data/ewinn/WTCCC_Sub.RData")
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

# simulation parameters
pve=0.3; rho=1; pc.var = 0; ncausal = 30
ncausal1= ncausal/6 #Set 1 of causal SNPs 
ncausal2 = ncausal-ncausal1 #Set 2 of Causal SNPs

### Create a matrix to save the results ###
FinalRMSE = matrix(nrow = n.datasets,ncol = 7,byrow=TRUE)
head(FinalRMSE)
colnames(FinalRMSE) = c("GP_NL","GPLIN","BRR","LASSO","ENET","SVM","Optimal")

FinalR2 = matrix(nrow = n.datasets,ncol = 7,byrow=TRUE)
head(FinalR2)
colnames(FinalR2) = c("GP_NL","GPLIN","BRR","LASSO", "ENET","SVM","Optimal")

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
  tss = sum((y_test-mean(y_test))^2)
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Defined extra parameters needed to run the analysis ###
  n = dim(X_train)[1] #Sample size
  p = dim(X_train)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  B = GaussKernel(t(X_train)) #deleting p,p in the GaussKernel function (unused arguments)
  
  ### Gibbs Sampler ###
  sigma2 = 1e-3
  sample_size = 5e3
  fhat = B %*% solve(B + diag(sigma2,n), y_train)

  ### Calculate Delta ###
  A <- B + diag(sigma2,nrow=n,ncol=n)
  A_svd <- propack.svd(A)
  Ainv = nearPD(A_svd$v%*%diag(1/A_svd$d, nrow=n, ncol=n)%*%t(A_svd$u))$mat
  Aiy <- Ainv%*%y_train
  BAiy <- B %*% Aiy
  
  #Calculate Delta, which ends up being p x sample_size matrix (add dopar instead of do once off windows)
  delta = ComputeESAFast(as.matrix(X_train), as.matrix(B), as.vector(Aiy), as.vector(BAiy), cores=cores)
  cat("calculated delta \n")
  GP_NL_pred = X_test%*%colMeans(delta)
  
  ### Get Diagnostics ###
  RMSE_GP_NL = sqrt(mean((y_test-GP_NL_pred)^2))
  rss = sum((y_test-GP_NL_pred)^2)
  R2_GP_NL = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  cat("GPlin \n")
  ### GP linear effect sizes ####
  beta.draws = ginv(X_train)%*%fhat
  
  GPLIN_pred = X_test%*%beta.draws
  
  ### Get Diagnostics ###
  RMSE_GPLIN = sqrt(mean((y_test-GPLIN_pred)^2))
  rss = sum((y_test-GPLIN_pred)^2)
  R2_GPLIN = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  cat("BRR \n")
  ### Bayesian Ridge Regression ###
  fit= cv.glmnet(X_train, y_train,intercept=FALSE,alpha=0)
  brr = as.matrix(coef(fit,s = fit$lambda.1se))
  brr = c(brr[-1,])
  BRR_pred = X_test%*%brr
  
  ### Get Diagnostics ###
  RMSE_BRR = sqrt(mean((y_test-BRR_pred)^2))
  rss = sum((y_test-BRR_pred)^2)
  R2_BRR = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  cat("Lasso \n")
  ### LASSO ###
  fit= cv.glmnet(X_train, y_train,intercept=FALSE,alpha=1)
  lasso = as.matrix(coef(fit,s = fit$lambda.1se))
  lasso = c(lasso[-1,])
  LASSO_pred = X_test%*%lasso
  
  ### Get Diagnostics ###
  RMSE_LASSO = sqrt(mean((y_test-LASSO_pred)^2))
  rss = sum((y_test-LASSO_pred)^2)
  R2_LASSO = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  cat("enet \n")
  ### Elastic Net ###
  fit= cv.glmnet(X_train, y_train,intercept=FALSE,alpha=0.5)
  enet = as.matrix(coef(fit,s = fit$lambda.1se))
  enet = c(enet[-1,])
  
  ### Get the posterior of the missing variables ###
  ENET_pred = X_test%*%enet
  
  ### Find MSE and Correlations ###
  RMSE_ENET = sqrt(mean((y_test-ENET_pred)^2))
  rss = sum((y_test-ENET_pred)^2)
  R2_ENET = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### SVM Model ###
  reg.svm  = ksvm(y=y_train,x=X_train,type="nu-svr",kernel="rbfdot")
  SVM_pred = predict(reg.svm, X_test, type="response")
  
  ### Get Diagnostics ###
  RMSE_SVM = sqrt(mean((y_test-SVM_pred)^2))
  rss = sum((y_test-SVM_pred)^2)
  R2_SVM = 1-rss/tss
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  scoresRMSE = c(RMSE_GP_NL, RMSE_GPLIN, RMSE_BRR,RMSE_ENET,RMSE_LASSO,RMSE_SVM)
  scoresR2 = c(R2_GP_NL, R2_GPLIN, R2_BRR,R2_ENET,R2_LASSO,R2_SVM)
  names(scoresRMSE) = c("GP_NL", "GPLIN", "BRR","ENET","LASSO","SVM")
  names(scoresR2) = c("GP_NL", "GPLIN", "BRR","ENET","LASSO","SVM")
  FinalRMSE[j,] = c(scoresRMSE,names(scoresRMSE)[which(scoresRMSE == min(scoresRMSE))])
  FinalR2[j,] = c(scoresR2,names(scoresR2)[which(scoresR2 == min(scoresR2))])
  
  ### Report Status ###
  cat("Completed Dataset", j, "\n", sep = " ")
}

file = paste("~/data/ewinn/RATEMC/",fn,"Pred_02_Results.RData",sep="")
save(FinalRMSE, FinalR2, file = file)
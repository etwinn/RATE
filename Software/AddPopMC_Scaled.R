### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(Matrix)
library(mvtnorm)
library(LaplacesDemon)
library(MASS)
library(corpcor)
library(varbvs)
library(glmnet)
#library(gdmp)
library(BGLR)
library(svd)

### Load in the RATE R functions ###
source("RATE (1).R")

### Load in the C++ BAKR functions ###
sourceCpp("BAKRGibbs.cpp")

cat("RATE and BAKR locked and loaded. \n")
### Define the Compute Power Function ###
compute.power <- function(pvals,SNPs){
  nsnps = length(pvals)
  Pos = SNPs #True Positive Set
  Negs = names(pvals)[which(names(pvals)%in%SNPs==FALSE)] #True Negative Set
  x = foreach(i = 1:nsnps)%dopar%{
    v = sort(pvals,decreasing = TRUE)[1:i] #Test Positives
    z = pvals[which(names(pvals)%in%names(v)==FALSE)] #Test Negatives
    
    TP = length(which(names(v)%in%Pos==TRUE))
    FP = length(which(names(v)%in%Pos==FALSE))
    TN = length(which(names(z)%in%Negs==TRUE))
    FN = length(which(names(z)%in%Negs==FALSE))
    
    TPR = TP/(TP+FN); FPR = FP/(FP+TN); FDR = FP/(FP+TP)
    c(TPR,FPR,FDR)
  }
  return(matrix(unlist(x),ncol = 3,byrow = TRUE))
}

######################################################################################
######################################################################################
######################################################################################

### Set the Seed for the analysis ###
set.seed(11151990)

### Set the Results Directory ###
fn = "AddPopMC_Scaled"

n.datasets = 50

### Load in the Genotypes ###
load("~/data/ewinn/WTCCC_Sub.RData")
#load("C:/Users/etwin/Downloads/WTCCC_Sub.RData")
unique.snps = apply(X.select,1,function(x) length(unique(x)))
X.select = X.select[unique.snps>1,]
maf=apply(X.select, 1, mean)/2
X = t(X.select[maf>0.05,]);
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); Xcenter=t((t(X)-Xmean)/Xsd)
ind = dim(X)[1]; nsnp = dim(X)[2]
cat("data locked and loaded \n")
### Compute the Top PCs ###
#PCs = ComputePCs(X,10)
PCs_center = ComputePCs(Xcenter, 5)

### Defined extra parameters needed to run the analysis ###
n = dim(X)[1] #Sample size
p = dim(X)[2] #Number of markers or genes

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
B = GaussKernel(t(X))
Bcenter = GaussKernel(t(Xcenter))
### To get the kernel for the GP_Lin, we just have to center B, as centering B 
### and centering the kernel calculated from Xcenter gives the same matrix.
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n
Bcenter = M%*%Bcenter%*%M
Bcenter = Bcenter/mean(diag(Bcenter))

# simulation parameters
pve=0.3; rho=1; pc.var = 0.1; ncausal = 30
ncausal1= ncausal/6 #Set 1 of causal SNPs 
ncausal2 = ncausal-ncausal1 #Set 2 of Causal SNPs

######################################################################################
######################################################################################
######################################################################################
cat("Entering for loop \n")
### Set the Saving Mechanisms ###
Final = list();
cores = detectCores()
registerDoParallel(cores=cores)
for(o in 1:n.datasets){
  start = Sys.time()
  s=sample(1:nsnp,ncausal,replace = FALSE)
  
  #Select Causal SNPs
  s1=sample(s, ncausal1, replace=F)
  s2=sample(s[s%in%s1==FALSE], ncausal2, replace=F)
  
  # Generate the ground-truth regression coefficients for the variables
  # (Xcenter). Adjust the effects so that
  # the variables (SNPs) explain x percent of the variance in the
  # outcome.
  Xcausal1=Xcenter[,s1]; Xcausal2=Xcenter[,s2];
  Xepi=c()
  for(i in 1:ncausal1){
    Xepi=cbind(Xepi,Xcausal1[,i]*Xcausal2)
  }
  dim(Xepi)
  
  # Marginal Effects Only
  Xmarginal=cbind(Xcenter[,s])
  beta=rnorm(dim(Xcenter[,s])[2])
  y_marginal=c(Xmarginal%*%beta)
  beta=beta*sqrt(pve*rho/var(y_marginal))
  y_marginal=Xmarginal%*%beta
  
  #Pairwise Epistatic Effects
  beta=rnorm(dim(Xepi)[2])
  y_epi=c(Xepi%*%beta)
  beta=beta*sqrt(pve*(1-rho)/var(y_epi))
  y_epi=Xepi%*%beta
  
  ### Define the effects of the PCs ###
  beta=rnorm(dim(PCs_center)[2])
  y_pcs=c(PCs_center%*%beta)
  beta=beta*sqrt(pc.var/var(y_pcs))
  y_pcs=PCs_center%*%beta
  
  # error
  y_err=rnorm(ind)
  y_err=y_err*sqrt((1-pve-pc.var)/var(y_err))
  
  
  ### Simulate the Response ###
  y=c(y_marginal+y_epi+y_pcs+y_err) #Full Model
  y=(y-mean(y))/(sd(y))
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Gibbs Sampler ###
  sigma2 = 1e-3
  sample_size = 5e3
  fhat = Bcenter %*% solve(Bcenter + diag(sigma2,n), y)
  fhat.rep = mvrnormArma(sample_size,fhat,Bcenter - Bcenter %*% solve(Bcenter+diag(sigma2,n),Bcenter))
  
  ### Calculate Delta ###
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  A <- B + diag(sigma2,nrow=n,ncol=n)
  A_svd <- propack.svd(A)
  Ainv = nearPD(A_svd$v%*%diag(1/A_svd$d, nrow=n, ncol=n)%*%t(A_svd$u))$mat
  #Ainv = nearPD(solve(A))$mat
  Aiy <- Ainv%*%y
  BAiy <- B %*% Aiy
  #IAiB <- (diag(1,nrow=n, ncol=n)-Ainv%*%B)
  #BIAiB <- B%*%IAiB
  
  #Calculate Delta, which ends up being p x sample_size matrix (add dopar instead of do once off windows)
  cat("Calculating delta \n")
  delta = ComputeESAFast(as.matrix(X), as.matrix(B), as.vector(Aiy), as.vector(BAiy))
  cat("Calculated delta! \n")
  
  ### Run the RATE_MC function ###
  ptm <- proc.time() #Start clock
  res = RATE(X=X,beta.draws=delta,snp.nms = colnames(X),cores = cores)
  ratesMC = res$RATE
  proc.time() - ptm #Stop clock
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Run the RATE Function ###
  beta.tilde = t(apply(fhat.rep,1,function(x) return(cov(Xcenter,x))))
  
  ptm <- proc.time() #Start clock
  res = RATE(X=Xcenter,beta.draws=beta.tilde,snp.nms = colnames(X),cores = cores)
  rates = res$RATE
  proc.time() - ptm #Stop clock
  
  ### LASSO ###
  fit= cv.glmnet(Xcenter, y,intercept=FALSE,alpha=1)
  lasso = as.matrix(coef(fit,s = fit$lambda.1se))
  lasso = c(lasso[-1,])
  
  ### Elastic Net ###
  fit= cv.glmnet(Xcenter, y,intercept=FALSE,alpha=0.5)
  enet = as.matrix(coef(fit,s = fit$lambda.1se))
  enet = c(enet[-1,])
  
  ### Scan One ###
  lp_scanone = sapply(1:ncol(Xcenter),function(i) -log10(summary(lm(y~Xcenter[,i]))$coef[2,4]))
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Compute the Power of Each Metric ###
  b = ratesMC
  power.kld.mc = compute.power(b,colnames(X)[s])
  
  b = rates
  power.kld = compute.power(b,colnames(X)[s])
  
  b = abs(lasso); names(b) = colnames(X)
  power.lasso= compute.power(b,colnames(X)[s])
  
  b = abs(enet); names(b) = colnames(X)
  power.enet= compute.power(b,colnames(X)[s])
  
  b = lp_scanone; names(b) = colnames(X)
  power.scanone = compute.power(b,colnames(X)[s])
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  
  ### Save the Results ###
  Final[[o]] = cbind(power.kld.mc, power.kld, power.lasso,power.enet,power.scanone)
  
  ### Report Status ###
  cat("Completed Dataset", o, "\n", sep = " ")
  cat(Sys.time()-start, "\n")
}

### Save the Results ###
file = paste("~/data/ewinn/RATEMC/",fn,"50.RData",sep="")
save(Final, file = file)
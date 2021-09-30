#RATE experiments for OSCAR (no plots).

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
#source("C:/Users/etwin/git_repos/RATE/Software/RATE2.R") #Changing path for etwin PC.
source("RATE2.R")

### Load in the C++ BAKR functions ###
#sourceCpp("C:/Users/etwin/git_repos/BAKR-master/BAKR-master/Rcpp/BAKRGibbs.cpp")
sourceCpp("BAKRGibbs.cpp")

#Register cores
cores = detectCores()
registerDoParallel(cores=cores)

######################################################################################
######################################################################################
######################################################################################
#Run z rounds
z=10
#Save variables via a list, need RATE_OG, RATE_MC, RATE_quad, RATE_quad2, RATE_combo
# Then need 4 rounds each. So let's make a list of each run, each run has
# each function, function, each function has a list of values from each iteration of RATE, which is a list
# Also need to append the data for each round (n, p, pve, rho, X, seed?)
RATE_MA_same = list()

#Want to parallelize, but may just not because it's giving linux a hard time
# foreach (k = 1:z) %dopar% {
for (k in 1:z) {
  RATES = list()
  ### Set the random seed to reproduce research ###
  set.seed(11151990+k)
  
  n = 2e3; p = 25; pve=0.6; rho=0.5;
  
  ### Define the Number of Causal SNPs
  ncausal = 3
  
  ### Simulate Synthetic Data ###
  maf <- 0.05 + 0.45*runif(p)
  X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
  X   <- matrix(as.double(X),n,p,byrow = TRUE)
  #Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
  s=c(23:25)
  
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
  #Gaussian kernel can be specified as k(u,v) = exp{||uâv||^2/2h^2}.
  
  n = dim(X)[1] #Sample size
  p = dim(X)[2] #Number of markers or genes
  
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  Kn = GaussKernel(t(X)); diag(Kn)=1 # 
  #Kn = X%*%t(X)/p
  ### Center and Scale K_tilde ###
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  Kn=M%*%Kn%*%M
  Kn=Kn/mean(diag(Kn))
  
  ### Gibbs Sampler ###
  sigma2 = 1e-3
  sample_size = 1e4
  fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
  fhat.rep = rmvnorm(sample_size,fhat,Kn - Kn %*% solve(Kn+diag(sigma2,n),Kn))
  
  
  #######################
  # Now run different versions of RATE.
  
  ##############################################
  ### Run the RATE_OG Function ###
  ##############################################
  
  nl = NULL
  res = RATE(X=X,f.draws=fhat.rep,snp.nms = colnames(X),cores = cores)
  
  ### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###
  
  ### Run the RATE Function ###
  top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))  
  res2 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res3 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ###Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res4 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  RATES = append(RATES, list("RATE_OG" = list("run1" = res, "run2"=res2, "run3"=res3, "run4"= res4)))
  
  ##############################################
  ### Run the RATE_quad Function ###
  ##############################################
  
  nl = NULL
  res = RATE_quad(X=X,f.draws=fhat.rep,snp.nms = colnames(X),cores = cores)
  
  ### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###
  
  ### Run the RATE Function ###
  top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))  
  res2 = RATE_quad(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res3 = RATE_quad(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ###Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res4 = RATE_quad(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  RATES = append(RATES, list("RATE_quad" = list("run1" = res, "run2"=res2, "run3"=res3, "run4"= res4)))
  
  ##############################################
  ### Run the RATE_quad2 Function ###
  ##############################################
  
  nl = NULL
  res = RATE_quad2(X=X,f.draws=fhat.rep,snp.nms = colnames(X),cores = cores)
  
  ### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###
  
  ### Run the RATE Function ###
  top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))  
  res2 = RATE_quad2(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res3 = RATE_quad2(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ###Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res4 = RATE_quad2(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  RATES = append(RATES, list("RATE_quad2" = list("run1" = res, "run2"=res2, "run3"=res3, "run4"= res4)))
  
  ##############################################
  ### Run the RATE_combo Function ###
  ##############################################
  
  nl = NULL
  res = RATE_combo(X=X,f.draws=fhat.rep,snp.nms = colnames(X),cores = cores)
  
  ### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###
  
  ### Run the RATE Function ###
  top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))  
  res2 = RATE_combo(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res3 = RATE_combo(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ###Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res4 = RATE_combo(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  RATES = append(RATES, list("RATE_combo" = list("run1" = res, "run2"=res2, "run3"=res3, "run4"= res4)))
  
  ##############################################
  ### Run the RATE_MC Function ###
  ##############################################
  start = Sys.time()
  delta = matrix(0,nrow=sample_size, ncol=p)
  #foreach(j = 1:p, combine = 'rbind')%dopar%{
  for (j in 1:p){
    #g = matrix(0,2000,25) #fhat 1 by 2000, g_j is 1 by 2000... but fhat rep is 10000 by 2000
    #for(j in 1:p)
    ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
    new_X = X+cbind(matrix(0,n,j-1),matrix(1,n,1),matrix(0,n,p-j))
    Kn_g = GaussKernel(t(new_X)); diag(Kn_g)=1 # 
    
    ### Center and Scale K_tilde ###
    #v=matrix(1, n, 1)
    #M=diag(n)-v%*%t(v)/n
    #Kn_g=M%*%Kn_g%*%M
    #Kn_g=Kn_g/mean(diag(Kn_g))
    #g #Don't need to sample, just get the expected value.
    g = Kn_g %*% solve(Kn_g + diag(sigma2,n), y)
    #g.rep is a 10000 by 2000 matrix
    g.rep = rmvnorm(sample_size,g,Kn_g - Kn_g %*% solve(Kn_g+diag(sigma2,n),Kn_g))
    delta[,j] = t(rowMeans(g.rep-fhat.rep))
  }

  end = Sys.time()
  
  nl = NULL
  res = RATE_MC(X=X,beta.draws = delta,snp.nms = colnames(X),cores = cores)
  res$Time = res$Time+(end-start)
  
  ### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###
  
  ### Run the RATE Function ###
  top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))  
  res2 = RATE_MC(X=X,beta.draws = delta,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res3 = RATE_MC(X=X,beta.draws = delta,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  ###Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###
  
  ### Run the RATE Function ###
  top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
  nl = c(nl,as.numeric(top))
  res4 = RATE_MC(X=X,beta.draws=delta,nullify = nl,snp.nms = colnames(X),cores = cores)
  
  RATES = append(RATES, list("RATE_MC" = list("run1" = res, "run2"=res2, "run3"=res3, "run4"= res4)))
  
  ##### SAVE ALL VARIABLES AND THE DATA ###
  data = list("X"=X, "n"=n, "p"=p, "rho"= rho, "pve"=pve)
  RATE_MA_same = append(RATE_MA_same, list("data" = data, "RATES" = RATES))
  assign(paste0("Trial_",k), RATE_MA_same[[k]])
}

##Export all saved variables

save(RATE_MA_same, file="~/scratch/data/RATE_MA_same_noXmean.Rdata")
# THIS IS MARYCLARE'S CODE FOR STUFF, LOOK THROUGH AND MAKE IT WORK.
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
sourceCpp("C:/Users/etwin/git_repos/BAKR-master/BAKR-master/Rcpp/BAKRGibbs.cpp")

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

### Running the Bayesian Gaussian Process (GP) Regression Model ###

### Create the Nonlinear Covariance Matrix with the Gaussian Kernel ###
#This function takes on two arguments:
#(1) The Genotype matrix X. This matrix should be fed in as a pxn matrix. That is, predictor
#are the rows and subjects/patients/cell lines are the columns.
#(2) The bandwidth (also known as a smoothing parameter or lengthscale) h. For example, the 
#Gaussian kernel can be specified as k(u,v) = exp{||uâv||^2/2h^2}.

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
B = GaussKernel(t(X)); diag(B)=1 # 

### Center and Scale K_tilde ###
# v=matrix(1, n, 1)
# M=diag(n)-v%*%t(v)/n
# Kn=M%*%Kn%*%M
# Kn=Kn/mean(diag(Kn))

### Gibbs Sampler ###
sample_size=100 # 1e4
sigma2 = 1e-3
A <- B + sigma2*diag(1, nrow = n, ncol = n)
Aiy <- solve(A, y)
BAiy <- B %*% Aiy
BAiB <- B %*% solve(A, B)
BmBAiB <- B - BAiB
# fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
ef <- BAiy
fhat.rep = rmvnorm(sample_size,BAiy,BmBAiB)

### Find basis and kernel matrix, center and scale, Gibbs Sampler for g ###
cores = detectCores()
registerDoParallel(cores=cores)

#Calculate Delta, which ends up being p x sample_size matrix (add dopar instead of do once off windows)
# f is sample_size x n matrix.... how solve this...
#Need to figure out list for this.

#delta = foreach(j = 1:p, .combine='rbind')%do%{
#g = matrix(0,2000,25) #fhat 1 by 2000, g_j is 1 by 2000... but fhat rep is 10000 by 2000
delta = matrix(0,nrow=sample_size, ncol=p)
for(k in 1:p){
  #cat("k=", k, "\n")
  ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
  new_X = X 
  new_X[,k] <- new_X[, k]+1
  # MCG: Need to be careful here - the predictor only takes on values 0-2, may want to be careful
  # Kn_g = GaussKernel(t(new_X)); diag(Kn_g)=1 # 
  Cj = GaussCoKernel(t(X), t(new_X)); diag(Cj)=1
  Dj = B
  
  CtAiy <- t(Cj) %*% Aiy
  CtAiC <- t(Cj) %*% solve(A, Cj)
  CtAiB <- t(Cj) %*% solve(A, B)
  BAiC <- B %*% solve(A, Cj)
  DmCtAiC <- Dj - CtAiC
  CtmCtAiB <- t(Cj) - CtAiB
  CmBAiC <- Cj - BAiC
  # fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
  eg <- CtAiy
  fghat.rep = rmvnorm(sample_size,c(BAiy, CtAiy),
                      rbind(cbind(BmBAiB, CmBAiC),
                            cbind(CtmCtAiB, DmCtAiC)))
  del <- t(apply(fghat.rep, 1, function(fg) {
    fg[length(fg)/2 + 1:(length(fg)/2)] - fg[1:(length(fg)/2)]
  }))
  delta[, k] <- rowMeans(del)
  
  # ### Center and Scale K_tilde ###
  #v=matrix(1, n, 1)
  #M=diag(n)-v%*%t(v)/n
  # # Kn_g=M%*%Kn_g%*%M
  # # Kn_g=Kn_g/mean(diag(Kn_g))
  # #g #Don't need to sample, just get the expected value.
  # g = Kn_g %*% solve(Kn_g + diag(sigma2,n), y)
  # #g.rep is a 10000 by 2000 matrix
  # g.rep = rmvnorm(sample_size,g,Kn_g - Kn_g %*% solve(Kn_g+diag(sigma2,n),Kn_g))
  # diff <- g.rep-fhat.rep
  # # this is a t by n matrix (t is number of draws), average over n
  # delta[,j] = t(rowMeans(diff))
  # #return t by 1 matrix, thus delta is a t by p matrix
}

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

### Plot the results with the uniformity line ###
par(mar=c(5,5,4,2))
barplot(rates,xlab = "Covariates",ylab=expression(RATE(tilde(beta)[j])),names.arg ="",col = ifelse(c(1:p)%in%s,"blue","grey80"),border=NA,cex.names = 0.6,ylim=c(0,0.6),cex.lab=1.25,cex.axis = 1.25)
lines(x = 0:length(rates)*1.5,y = rep(1/(p-length(nl)),length(rates)+1),col = "red",lty=2,lwd=2)
legend("topleft",legend=c(as.expression(bquote(DELTA~"="~.(round(DELTA,3)))),as.expression(bquote("ESS ="~.(round(ESS,2))*"%")), as.expression(bquote("RATE fn ="~.(rate_choice)))),bty = "n",pch = 19,cex = 1.25,col = "red")

######################################################################################
######################################################################################
######################################################################################

### Find Second Order Centrality by Nullifying the Top Associated Predictor Variable ###

### Run the RATE Function ###
top = substring(names(res$KLD)[order(res$KLD,decreasing=TRUE)[1]],first = 4)
nl = c(nl,as.numeric(top))  
res2 = RATE_MC(X=X,beta.draws=delta, nullify=nl, snp.nms = colnames(X),cores = cores)

### Get the Results ###
rates = res2$RATE
DELTA = res2$Delta
ESS = res2$ESS

### Plot the results with the uniformity line ###
par(mar=c(5,5,4,2))
barplot(rates,xlab = "Covariates",ylab=bquote(RATE(tilde(beta)[j]~"|"~tilde(beta)[.(as.integer(nl))]=="0")),names.arg = "",col = ifelse(c(1:p)[-nl]%in%s,"blue","grey80"),border=NA,cex.names = 0.6,ylim=c(0,0.6),cex.lab=1.25,cex.axis = 1.25)
lines(x = 0:length(rates)*1.5,y = rep(1/(p-length(nl)),length(rates)+1),col = "red",lty=2,lwd=2)
legend("topleft",legend=c(as.expression(bquote(DELTA~"="~.(round(DELTA,3)))),as.expression(bquote("ESS ="~.(round(ESS,2))*"%")), as.expression(bquote("RATE fn ="~.(rate_choice)))),bty = "n",pch = 19,cex = 1.25,col = "red")

######################################################################################
######################################################################################
######################################################################################

### Find Third Order Centrality by Nullifying the Top 2 Associated Predictor Variables ###

### Run the RATE Function ###
top = substring(names(res2$KLD)[order(res2$KLD,decreasing=TRUE)[1]],first = 4)
nl = c(nl,as.numeric(top))
res3 = RATE_MC(X=X,beta.draws=delta, nullify = nl, snp.nms = colnames(X),cores = cores)

### Get the Results ###
rates = res3$RATE
DELTA = res3$Delta
ESS = res3$ESS

### Plot the results with the uniformity line ###
par(mar=c(5,5,4,2))
barplot(rates,xlab = "Covariates",ylab=bquote(RATE(tilde(beta)[j]~"|"~tilde(beta)[.(as.integer(nl[1]))]==~tilde(beta)[.(as.integer(nl[2]))]=="0")),names.arg = "",col = ifelse(c(1:p)[-nl]%in%s,"blue","grey80"),border=NA,cex.names = 0.6,ylim=c(0,0.6),cex.lab=1.25,cex.axis = 1.25)
lines(x = 0:length(rates)*1.5,y = rep(1/(p-length(nl)),length(rates)+1),col = "red",lty=2,lwd=2)
legend("topleft",legend=c(as.expression(bquote(DELTA~"="~.(round(DELTA,3)))),as.expression(bquote("ESS ="~.(round(ESS,2))*"%")),as.expression(bquote("RATE fn ="~.(rate_choice)))),bty = "n",pch = 19,cex = 1.25,col = "red")

######################################################################################
######################################################################################
######################################################################################

### Find Fourth Order Centrality by Nullifying the Top 3 Associated Predictor Variables ###

### Run the RATE Function ###
top = substring(names(res3$KLD)[order(res3$KLD,decreasing=TRUE)[1]],first = 4)
nl = c(nl,as.numeric(top))
res4 = RATE_MC(X=X,beta.draws=delta, nullify=nl, snp.nms = colnames(X),cores = cores)

### Get the Results ###
rates = res4$RATE
DELTA = res4$Delta
ESS = res4$ESS

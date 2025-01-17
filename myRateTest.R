# myRateTest.R
# Written by Emily Winn
# Last updated 9/1/2021

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

n = 2e3; p = 25; pve=0.1; rho=0.5; #rho=0 makes all epistatic, 1 is all marginal

### Define the Number of Causal SNPs
ncausal = 3

### Simulate Synthetic Data ###
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
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
#Gaussian kernel can be specified as k(u,v) = exp{||u−v||^2/2h^2}.

### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
Kn = GaussKernel(t(X)); diag(Kn)=1 # 

### Center and Scale K_tilde ###
v=matrix(1, n, 1)
M=diag(n)-v%*%t(v)/n

Kn=Kn/mean(diag(Kn))

### Gibbs Sampler ###
sigma2 = 1e-3
fhat = Kn %*% solve(Kn + diag(sigma2,n), y)
fhat.rep = rmvnorm(1e4,fhat,Kn - Kn %*% solve(Kn+diag(sigma2,n),Kn))

### Compute the KL Divergence to find Marginal Importance ###
cores = detectCores()
registerDoParallel(cores=cores)

### Run the RATE Function ###
rate_choice = "RATE"
nl = NULL
#start = Sys.time()
res = RATE(X=X,f.draws=fhat.rep,snp.nms = colnames(X),cores = cores)
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
res2 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)

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
res3 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)

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
res4 = RATE(X=X,f.draws=fhat.rep,nullify = nl,snp.nms = colnames(X),cores = cores)

### Get the Results ###
rates = res4$RATE
DELTA = res4$Delta
ESS = res4$ESS

### Plot the results with the uniformity line ###
par(mar=c(5,5,4,2))
barplot(rates,xlab = "Covariates",ylab=bquote(RATE(tilde(beta)[j]~"|"~tilde(beta)[.(as.integer(nl[1]))]==~tilde(beta)[.(as.integer(nl[2]))]==~tilde(beta)[.(as.integer(nl[3]))]=="0")),names.arg = "",col = ifelse(c(1:p)[-nl]%in%s,"blue","grey80"),border=NA,cex.names = 0.6,ylim=c(0,0.6),cex.lab=1.25,cex.axis = 1.25)
lines(x = 0:length(rates)*1.5,y = rep(1/(p-length(nl)),length(rates)+1),col = "red",lty=2,lwd=2)
legend("topleft",legend=c(as.expression(bquote(DELTA~"="~.(round(DELTA,3)))),as.expression(bquote("ESS ="~.(round(ESS,2))*"%")),as.expression(bquote("RATE fn ="~.(rate_choice)))),bty = "n",pch = 19,cex = 1.25,col = "red")
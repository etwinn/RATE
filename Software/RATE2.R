### Variable Prioritization via RelATive cEntrality (RATE) centrality measures ###

### NOTE: This function assumes that one has already obtained (posterior) draws/estimates of a nonparametric or nonlinear function as suggested in Crawford et al. (2018) ###

### Review of the function parameters ###
#'X' is the nxp design matrix (e.g. genotypes) where n is the number of samples and p is the number of dimensions. This is the original input data
#'f.draws' is the Bxn matrix of the nonparametric model estimates (i.e. f.hat) with B being the number of sampled (posterior) draws;
#'prop.var' is the desired proportion of variance that the user wants to explain when applying singular value decomposition (SVD) to the design matrix X (this is preset to 1);
#'low.rank' is a boolean variable detailing if the function will use low rank matrix approximations to compute the RATE values --- note that this highly recommended in the case that the number of covariates (e.g. SNPs, genetic markers) is large; 
#'rank.r' is a parameter detailing which matrix rank order to compute the RATE measures --- the default is min(n,p); however, this can be specified to a smaller value at the possible cost of accuracy and power;
#'nullify' is an optional vector specifying a given predictor variable effect that user wishes to "remove" from the data set;
#'snp.nms' is an optional vector specifying the names of the predictor variables;
#'cores' is a parameter detailing the number of cores to parallelize over. If too many are assigned, RATE will set this variable to maximum number of cores that are available on the operating machine.

######################################################################################
######################################################################################
######################################################################################

RATE = function(X = X, f.draws = f.draws,prop.var = 1, low.rank = FALSE, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = snp.nms, cores = 1, kl = "OG"){
  
  ### Install the necessary libraries ###
  #usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  usePackage("matrixcalc")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  registerDoParallel(cores=cores)
  
  ### First Run the Matrix Factorizations ###  
  svd_X = propack.svd(X,rank.r); 
  dx = svd_X$d > 1e-10
  px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
  r_X = dx&px 
  u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
  v = svd_X$v[,r_X]  
  
  if(low.rank==TRUE){
  # Now, calculate Sigma_star
  SigmaFhat = cov(f.draws)
  Sigma_star = u %*% SigmaFhat %*% t(u)
  
  # Now, calculate U st Lambda = U %*% t(U)
  svd_Sigma_star = propack.svd(Sigma_star,rank.r)
  r = svd_Sigma_star$d > 1e-10
  U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
  
  V = v%*%Sigma_star%*%t(v) #Variances
  
  mu = v%*%u%*%colMeans(f.draws) #Effect Size Analogues 
  }else{
    beta.draws = t(ginv(X)%*%t(f.draws))
    V = cov(beta.draws); #V = as.matrix(nearPD(V)$mat)
    D = ginv(V)
    svd_D = svd(D)
    r = sum(svd_D$d>1e-10)
    U = with(svd_D,t(sqrt(d[1:r])*t(u[,1:r])))
    
    mu = colMeans(beta.draws)
  }
  
  ### Create Lambda ###
  Lambda = tcrossprod(U)
  
  ### Compute the Kullback-Leibler divergence (KLD) quadratic alternative for Each Predictor ###
  int = 1:length(mu); l = nullify
  
 if(length(l)>0){int = int[-l]}

 if(kl == "OG"){
   KLD = foreach(j = int, .combine='c', .export='sherman_r')%dopar%{ #Adding the .export part for windows
     q = unique(c(j,l))
     m = abs(mu[q])
     
     U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
     #U_Lambda_sub = U_Lambda_sub[-q,-q]
     alpha = crossprod(U_Lambda_sub[-q,q],crossprod(U_Lambda_sub[-q,-q],U_Lambda_sub[-q,q]))
     kld = (t(m)%*%alpha%*%m)/2
     names(kld) = snp.nms[j]
     kld
   } 
 }else if(kl == "MC"){
   KLD = foreach(j=int, .combine='c', .export = 'hadamard.prod')%dopar%{
    q = unique(c(j,l))
    m = abs(mu[q])
  
    #U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
    #U_Lambda_sub = U_Lambda_sub[-q,-q]
    theta = crossprod(-Lambda[-q,-q],Lambda[-q,q])
    alpha = t(theta)%*%Lambda[-q,-q]%*%theta
    kld = sum(hadamard.prod(Lambda[-q,-q],V[-q,-q])) + ((beta.draws[q]-m)^2)*alpha
    names(kld) = snp.nms[j]
    kld
   }
 }else{
   print('KLD specification not among the options.')
   KLD = 1
 }

  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}

RATE_MC = function(X = X, f.draws = f.draws,prop.var = 1, low.rank = FALSE, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = snp.nms, cores = 1, kl = "OG"){
  
  ### Install the necessary libraries ###
  #usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  usePackage("matrixcalc")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  registerDoParallel(cores=cores)
  
  ### First Run the Matrix Factorizations ###  
  svd_X = propack.svd(X,rank.r); 
  dx = svd_X$d > 1e-10
  px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
  r_X = dx&px 
  u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
  v = svd_X$v[,r_X]  
  
  ### Calculate g(j) if necessary (note: dopar loop does not work on windows machines):
    n = dim(X)[1]
    p = dim(X)[2]
    g = foreach(j = 1:p, .combine='c')%dopar%{
    #g = matrix(0,2000,25) #fhat 1 by 2000, g_j is 1 by 2000... but fhat rep is 10000 by 2000
    #for(j in 1:p){
    E_j = matrix(0, nrow=n, ncol=p)
    E_j[,j] = 1
    new_X = X + E_j
    ### Find the Approximate Basis and Kernel Matrix; Choose N <= D <= P ###
    Kn_g = GaussKernel(t(new_X)); diag(Kn_g)=1 # 
      
    ### Center and Scale K_tilde ###
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    Kn_g=M%*%Kn_g%*%M
    Kn_g=Kn_g/mean(diag(Kn_g))
    }
      
  ### Gibbs Sampler ###
  sigma2 = 1e-3
  g[,j] = Kn_g %*% solve(Kn_g + diag(sigma2,n), y)
  #g #Don't need to sample, just get the expected value.
  
  print("Dim g is ", dim(g))
  
  if(low.rank==TRUE){
    # Now, calculate Sigma_star
    SigmaFhat = cov(f.draws)
    Sigma_star = u %*% SigmaFhat %*% t(u)
    
    # Now, calculate U st Lambda = U %*% t(U)
    svd_Sigma_star = propack.svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    V = v%*%Sigma_star%*%t(v) #Variances
    
    mu = v%*%u%*%colMeans(f.draws) #Effect Size Analogues 
  }else{
    #Beta draws j needs to be the mean of g_j - f.
    fmean = colMeans(f.draws)
    beta.draws = foreach(j=1:p)%dopar%{
      beta.draws = g[,j]-fmean
      beta.draws
    }
    V = cov(beta.draws); #V = as.matrix(nearPD(V)$mat)
    D = ginv(V)
    svd_D = svd(D)
    r = sum(svd_D$d>1e-10)
    U = with(svd_D,t(sqrt(d[1:r])*t(u[,1:r])))
    
    mu = colMeans(beta.draws)
  }
  
  ### Create Lambda ###
  Lambda = tcrossprod(U)
  
  ### Compute the Kullback-Leibler divergence (KLD) quadratic alternative for Each Predictor ###
  int = 1:length(mu); l = nullify
  
  if(length(l)>0){int = int[-l]}
  
  if(kl == "OG"){
    KLD = foreach(j = int, .combine='c', .export='sherman_r')%dopar%{ #Adding the .export part for windows
      q = unique(c(j,l))
      m = abs(mu[q])
      
      U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
      #U_Lambda_sub = U_Lambda_sub[-q,-q]
      alpha = crossprod(U_Lambda_sub[-q,q],crossprod(U_Lambda_sub[-q,-q],U_Lambda_sub[-q,q]))
      kld = (t(m)%*%alpha%*%m)/2
      names(kld) = snp.nms[j]
      kld
    } 
  }else if(kl == "MC"){
    KLD = foreach(j=int, .combine='c', .export = 'hadamard.prod')%dopar%{
      q = unique(c(j,l))
      m = abs(mu[q])
      
      #U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
      #U_Lambda_sub = U_Lambda_sub[-q,-q]
      theta = crossprod(-Lambda[-q,-q],Lambda[-q,q])
      alpha = t(theta)%*%Lambda[-q,-q]%*%theta
      kld = sum(hadamard.prod(Lambda[-q,-q],V[-q,-q])) + ((beta.draws[q]-m)^2)*alpha
      names(kld) = snp.nms[j]
      kld
    }
  }else{
    print('KLD specification not among the options.')
    KLD = 1
  }
  
  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}


#Calculate the Beta's using just the linear model coefficients (super naive way)
RATE_lin <- function(X = X, f.draws = f.draws,prop.var = 1, low.rank = FALSE, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = snp.nms, cores = 1){
  ### Install the necessary libraries ###
  #usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  #registerDoParallel(cores=cores)
  
  ### First Run the Matrix Factorizations ###  
  svd_X = propack.svd(X,rank.r); 
  dx = svd_X$d > 1e-10
  px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
  r_X = dx&px 
  u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
  v = svd_X$v[,r_X]  
  
  if(low.rank==TRUE){
    # Now, calculate Sigma_star
    SigmaFhat = cov(f.draws)
    Sigma_star = u %*% SigmaFhat %*% t(u)
    
    # Now, calculate U st Lambda = U %*% t(U)
    svd_Sigma_star = propack.svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    mu = v%*%u%*%colMeans(f.draws)
  }else{
    beta.draws = t(apply(f.draws, 1, function(y) solve(crossprod(X))%*%t(X)%*%y))
    V = cov(beta.draws); #V = as.matrix(nearPD(V)$mat)
    D = ginv(V)
    svd_D = svd(D)
    r = sum(svd_D$d>1e-10)
    U = with(svd_D,t(sqrt(d[1:r])*t(u[,1:r])))
    
    mu = colMeans(beta.draws)
  }
  
  ### Create Lambda ###
  Lambda = tcrossprod(U)
  
  ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
  int = 1:length(mu); l=nullify;
  
  if(length(l)>0){int = int[-l]}
  
  KLD = foreach(j = int, .combine='c', .export='sherman_r')%dopar%{ #Adding the .export part for windows
    q = unique(c(j,l))
    m = abs(mu[q])
    
    U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
    #U_Lambda_sub = U_Lambda_sub[-q,-q]
    alpha = crossprod(U_Lambda_sub[-q,q],crossprod(U_Lambda_sub[-q,-q],U_Lambda_sub[-q,q]))
    kld = (t(m)%*%alpha%*%m)/2
    names(kld) = snp.nms[j]
    kld
  }
  
  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}

#Calculate the Beta's using just the quadratic model coefficients (super naive way)
RATE_quad <- function(X = X, f.draws = f.draws,prop.var = 1, low.rank = FALSE, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = snp.nms, cores = 1){
  ### Install the necessary libraries ###
  usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  #registerDoParallel(cores=cores)
  
  ### First Run the Matrix Factorizations ###  
  svd_X = propack.svd(X,rank.r); 
  dx = svd_X$d > 1e-10
  px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
  r_X = dx&px 
  u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
  v = svd_X$v[,r_X]  
  
  if(low.rank==TRUE){
    # Now, calculate Sigma_star
    SigmaFhat = cov(f.draws)
    Sigma_star = u %*% SigmaFhat %*% t(u)
    
    # Now, calculate U st Lambda = U %*% t(U)
    svd_Sigma_star = propack.svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    mu = v%*%u%*%colMeans(f.draws)
  }else{
    Z = X*X
    beta.draws = t(apply(f.draws, 1, function (y) ginv(Z)%*%(diag(n)-ginv(X))%*%y))
    V = cov(beta.draws); #V = as.matrix(nearPD(V)$mat)
    D = ginv(V)
    svd_D = svd(D)
    r = sum(svd_D$d>1e-10)
    U = with(svd_D,t(sqrt(d[1:r])*t(u[,1:r])))
    
    mu = colMeans(beta.draws)
  }
  
  ### Create Lambda ###
  Lambda = tcrossprod(U)
  
  ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
  int = 1:length(mu); l=nullify;
  
  if(length(l)>0){int = int[-l]}
  
  KLD = foreach(j = int, .combine='c', .export='sherman_r')%dopar%{ #Adding the .export part for windows
    q = unique(c(j,l))
    m = abs(mu[q])
    
    U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
    #U_Lambda_sub = U_Lambda_sub[-q,-q]
    alpha = crossprod(U_Lambda_sub[-q,q],crossprod(U_Lambda_sub[-q,-q],U_Lambda_sub[-q,q]))
    kld = (t(m)%*%alpha%*%m)/2
    names(kld) = snp.nms[j]
    kld
  }
  
  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}

sherman_r <- function(Ap, u, v) {
  Ap - (Ap %*% u %*% t(v) %*% Ap)/drop(1 + u %*% t(v) %*% Ap)
}

### Define the Package ###
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

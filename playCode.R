#Code to play with from Lorin

X =matrix(rnorm(1e4),nrow = 1e3)
beta1 = rnorm(1) #Additive Effect
beta2 = rnorm(1) #Additive Effect
alpha = rnorm(1) #Interaction Effect

y = X[,1]*beta1 + X[,2]*beta2 + (X[,1]*X[,2])*alpha + rnorm(nrow(X))

### Models ###
mod = lm(y~.-1,data = data.frame(X)); summary(mod) ### Linear Model
mod2 = lm(y~(.)^2-1,data = data.frame(X)); summary(mod2) ### Linear Model with Effects

### Effect Sizes 
beta_hat = solve(crossprod(X))%*%t(X)%*%y

### Check Accuracy of Model Coefficients ###
cbind(c(beta1,beta2,rep(0,ncol(X)-2)),beta_hat, mod$coefficients,mod2$coefficients[1:ncol(X)])
rm(list=ls())
set.seed(123)

library(MASS)
library(MCMCpack)
library(bayesm)

source("Function-LM.R")
set.seed(1234)


## settings
out.X <- F     # T or F
aa <- 10    # outlier location: 10 or 20
R <- 20000   # number of Monte Carlo replications
n <- 50   # the number of observations
p <- 3    # the number of predictor variables

mc <- 1500    # MCMC length
bn <- 500    # length of burn-in 



om.set <- c(0, 0.05, 0.1, 0.2)    # contamination ratio
if( aa!=10 ){ om.set <- om.set[-1] }
J <- length(om.set)


## parameters
sigma <- 0.5
Beta <- c(0.3, 0, 1)   # coefficients
int <- 0.5    # intercept
Para <- c(int, Beta)


## functions
quant <- function(x){ quantile(x, prob=c(0.025, 0.975)) } 
IF.compute <- function(x){ as.vector(numEff(x)$f) }



## results box
meth <- c("LPM", "aLPM","LPT1", "LPT2", "C", "T3", "aT", "MT", "N")
L <- length(meth)
MSE <- array(NA, c(R, L, J))
MSE.sigma <- array(NA, c(R, L, J))
CP <- array(NA, c(R, L, J))
AL <- array(NA, c(R, L, J))
IF <- array(NA, c(R, L, J))
dimnames(MSE)[[2]] <- dimnames(MSE.sigma)[[2]] <- dimnames(CP)[[2]] <- dimnames(AL)[[2]] <- meth




###   Monte Carlo replications   ###
for(j in 1:J){
  om <- om.set[j]
  print(paste0("omega=", om))
  
  for(r in 1:R){
    # covariates
    rho <- 0.2
    mat <- matrix(NA, p, p)
    for(k in 1:p){
      for(l in 1:p){ mat[k,l] <- rho^(abs(k-l)) }
    }
    X <- mvrnorm(n, rep(0, p), mat) 
    
    # high leverage points 
    if(out.X){
      omX <- 0.1
      for(i in 1:n){
        chX <- rbinom(1, 1, omX)
        if(chX==1){
          sel <- sample(1:p, 1)
          X[i, sel] <- X[i, sel] + 5
        }
      }
    }
    
    # data generation 
    ch <- rbinom(n, 1, om)
    noise <- (1-ch)*rnorm(n, 0, 1) + ch*rnorm(n, aa, 1)
    Y <- int + as.vector(X%*%Beta) + sigma*noise
    
    plot(X[,1], Y)
  
    # fitting 
    Est.Beta <- matrix(NA, p+1, L)
    Est.Sig <- c()
    CI.Beta <- list()
    for(l in 1:L){
      if(l==1){ fit <- EHE(Y, X, gam.est=F, mc=mc, burn=bn) }
      if(l==2){ fit <- EHE(Y, X, gam.est=T, mc=mc, burn=bn) }
      if(l==3){ fit <- LPTN(Y, X, rho=0.95, mc=mc, burn=bn) }
      if(l==4){ fit <- LPTN(Y, X, rho=0.8, mc=mc, burn=bn) }
      if(l==5){ fit <- TR(Y, X, mc=mc, burn=bn, nu=1) }
      if(l==6){ fit <- TR(Y, X, mc=mc, burn=bn, nu=3) }
      if(l==7){ fit <- TR(Y, X, mc=mc, burn=bn, nu=3, estimation=T) }
      if(l==8){ fit <- MT(Y, X, mc=mc, burn=bn, nu=1/2) }
      if(l==9){ fit <- LR(Y, X, mc=mc, burn=bn) }
      
      Est.Beta[,l] <- apply(fit$Beta, 2, mean)
      Est.Sig[l] <- mean(fit$Sig)
      CI.Beta[[l]]  <- apply(fit$Beta, 2, quant)
      IF[r,l,j] <- mean(apply(fit$Beta, 2, IF.compute))
    }
    
    # performance measure 
    MSE[r,,j] <- apply((Est.Beta-Para)^2, 2, mean)    # MSE (Beta)
    MSE.sigma[r,,j] <- (Est.Sig-sigma)^2      # MSE (Sigma)
    for(l in 1:L){
      CP[r,l,j] <- mean(CI.Beta[[l]][1,]<Para & CI.Beta[[l]][2,]>Para)   # CP 
      AL[r,l,j] <- mean(CI.Beta[[l]][2,] - CI.Beta[[l]][1,])   # AL
    }
    
    # print
    if(round(r/100)==(r/100)){ print(r) }
  }
}


100*sqrt(apply(MSE.sigma, c(2,3), mean))
100*sqrt(apply(MSE, c(2,3), mean))
apply(CP, c(2,3), mean)
apply(AL, c(2,3), mean)


if(out.X){
  save(list=ls(), file=paste0("sim-a", aa, "-out.RData"))
}else{
  save(list=ls(), file=paste0("sim-a", aa, ".RData"))
}





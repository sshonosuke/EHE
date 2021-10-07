library(glmnet)
library(MCMCpack)
library(statmod)
library(GIGrvg)
library(robust)


###  EHE distribution  ###
EHE <- function(Y, X, gam.est=F, gam=1, mc=2000, burn=500, gam.hp=c(100,100)){
  delta <- 1   # hyperparameter for sigma
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  pp <- dim(XX)[[2]]
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  om.pos <- c()
  if(gam.est){ Gam.pos <- c() }
  
  # initial value
  init.fit <- lmRob(Y~X)
  Beta <- coef(init.fit)
  res <- init.fit$residuals
  Sig <- init.fit$scale
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  W <- rep(1, n)
  U1 <- rep(1, n)   # non-outlier (fixed)
  U2 <- rep(10, n)   # outlier (sampled)
  Gam <- gam
  
  # MCMC
  for(k in 1:mc){
    UU <- (U2^Z) * (U1^(1-Z))
    # Beta and Sigma
    mat <- solve(t(XX)%*%(XX/UU))
    Beta <- mvrnorm(1, mat%*%t(XX)%*%(Y/UU), Sig^2*mat)
    Beta.pos[k,] <- Beta
    mu <- as.vector(XX%*%Beta)
    resid <- Y - mu
    Sig <- sqrt( rinvgamma(1, delta/2+n/2, delta/2+sum(resid^2/UU)/2 ) )
    Sig.pos[k] <- Sig
    
    # Z
    mu <- as.vector(XX%*%Beta)
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    # omega
    om <- rbeta(1, 1/2+sum(Z), 1/2+n-sum(Z))
    om.pos[k] <- om 
    
    # Gamma
    if(gam.est){
      Gam <- rgamma(1, gam.hp[1]+n, gam.hp[2]+sum(log(1+log(1+U2))) )
      Gam.pos[k] <- Gam
    }
    
    # V and W
    W <- rgamma(n, 1+Gam, 1+log(1+U2))
    V <- rgamma(n, 1+W, 1+U2)
    
    # UU 
    Chi <- resid^2/Sig^2
    GIG <- c()
    for(i in 1:n){ 
      GIG[i] <- rgig(1, lambda=1/2, chi=Chi[i], psi=2*V[i])
    }
    U2 <- Z*GIG + (1-Z)*rgamma(n, 1, V)
  }
  
  # Summary
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  Sig.pos <- Sig.pos[-omit]
  om.pos <- om.pos[-omit]
  Res <- list(Beta=Beta.pos, Sig=Sig.pos, Om=om.pos)
  return(Res)
}




###  LPTN distribution  ###
LPTN <- function(Y, X, rho=0.8, mc=3000, burn=1000){
  ## preparation 
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  pp <- dim(XX)[[2]]
  tau <- qnorm( (1+rho)/2 )
  lam <- 2*(1-rho)^(-1)*dnorm(tau)*tau*log(tau)
  
  LL <- function(beta, sig){
    z <- (Y - as.vector(XX%*%beta))/sig
    z1 <- z[ abs(z)<=tau ]
    z2 <- z[ abs(z)>tau ]
    dens <- rep(0, n)
    dens[ abs(z)<=tau ] <- dnorm(z1)
    dens[ abs(z)>tau ] <- dnorm(tau)*tau/abs(z2)*(log(tau)/log(abs(z2)))^(lam+1)
    return( log(dens/sig) )
  }
  
  ## MCMC
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  init.fit <- lmRob(Y~X)
  Beta <- coef(init.fit)
  Sig <- init.fit$scale
  
  for(k in 1:mc){
    Beta <- Beta
    Sig0 <- Sig
    # Beta
    new.Beta <- Beta + rnorm(pp, 0, 0.05)
    log.prob <- sum( LL(new.Beta, Sig) - LL(Beta, Sig) ) - sum(new.Beta^2)/(2*1000^2) + sum(Beta^2)/(2*1000^2) 
    prob <- min(1, exp(log.prob))
    ch <- rbinom(1, 1, prob)
    Beta <- Beta + ch*(new.Beta - Beta)
    # Sig
    new.Sig <- Sig + rnorm(1, 0, 0.05)
    if(new.Sig < 0.001){ new.Sig <- 0.001 }
    if(new.Sig > 1000){ new.Sig <- 1000 }
    log.prob <- sum( LL(Beta, new.Sig) - LL(Beta, Sig) ) 
    prob <- min(1, exp(log.prob))
    ch <- rbinom(1, 1, prob)
    Sig <- Sig + ch*(new.Sig - Sig)
    
    Beta.pos[k, ] <- Beta
    Sig.pos[k] <- Sig
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om, ]
  Sig.pos <- Sig.pos[-om]
  return(list(Beta=Beta.pos, Sig=Sig.pos))
}






###  Regression with t-distribution   ###
TR <- function(Y, X, mc=2000, burn=500, nu=3, estimation=F){
  ## preparation 
  cc <- 1   # prior for sigma^2
  nu.grid <- c(1:5, 8, 10, 15, 20, 30, 50)
  L <- length(nu.grid)
  nu.pos <- c()
  Int.pos <- c()
  Beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  
  ## initial value
  init.fit <- lmRob(Y~X)
  Beta <- coef(init.fit)
  sig <- init.fit$scale
  Int <- Beta[1]
  Beta <- Beta[-1]
  sig <- 1
  V <- rep(1, n)  
  
  ## MCMC
  for(k in 1:mc){
    # Int
    resid <- as.vector(Y-X%*%Beta)
    Int <- rnorm(1,sum(V*resid)/sum(V),sig/sqrt(sum(V)))
    Int.pos[k] <- Int
    # Beta
    Yt <- Y-Int
    W <- diag(V)
    invA <- ginv(t(X)%*%W%*%X/sig^2)
    Beta <- mvrnorm(1,invA%*%t(X)%*%W%*%Yt/sig^2,invA)
    Beta.pos[k,] <- Beta
    # sigma
    resid <- as.vector(Y-Int-X%*%Beta)
    sig2 <- rinvgamma(1,n/2+cc,sum(V*resid^2)/2+cc)
    sig <- sqrt(sig2)
    sig.pos[k] <- sig
    # V
    V <- rgamma(n,nu/2+1/2,nu/2+resid^2/(2*sig^2))
    # nu (if estimated)
    if(estimation){
      logdens <- matrix(NA, n, L)
      for(l in 1:L){
        logdens[,l] <- dgamma(V, nu.grid[l]/2+1/2, nu.grid[l]/2+resid^2/(2*sig^2), log=T)
      }
      logdens <- apply(logdens, 2, mean)
      pp <- exp(logdens)/sum(exp(logdens))
      nu <- sample(nu.grid, 1, prob=pp)
      nu.pos[k] <- nu
    }
  }
  ## Summary
  om <- 1:burn
  Beta.pos <- cbind(Int.pos, Beta.pos)[-om,]
  Sig.pos <- sig.pos[-om]
  if(estimation){ nu.pos <- nu.pos[-om] }
  Res <- list(Beta=Beta.pos, Sig=Sig.pos, Nu=nu.pos)
  return(Res)
}







###  Regression with t-distribution (two-component mixture)   ###
MT <- function(Y, X, nu=0.5, mc=2000, burn=500){
  delta <- 1    # hyperparameter for sigma
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  pp <- dim(XX)[[2]]
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  om.pos <- c()
  
  # initial value
  init.fit <- lmRob(Y~X)
  Beta <- coef(init.fit)
  Sig <- init.fit$scale
  res <- init.fit$residuals
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  U1 <- rep(1, n)   # non-outlier (fixed)
  U2 <- rep(5, n)   # outlier (sampled)
  
  # MCMC
  for(k in 1:mc){
    UU <- (U2^Z) * (U1^(1-Z))
    # Beta and Sigma
    mat <- solve(t(XX)%*%(XX*UU))
    Beta <- mvrnorm(1, mat%*%t(XX)%*%(Y*UU), Sig^2*mat)
    Beta.pos[k,] <- Beta
    mu <- as.vector(XX%*%Beta)
    resid <- Y - mu
    Sig <- sqrt( rinvgamma(1, delta/2+n/2, delta/2+sum(resid^2*UU)/2 ) )
    Sig.pos[k] <- Sig
    
    # Z
    mu <- as.vector(XX%*%Beta)
    logdens1 <- dnorm(Y, mu, Sig/sqrt(U1), log=T)
    logdens2 <- dnorm(Y, mu, Sig/sqrt(U2), log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    # omega
    om <- rbeta(1, 1/2+sum(Z), 1/2+n-sum(Z))
    om.pos[k] <- om 
    
    # UU 
    resid <- Y - mu
    Chi <- resid^2/Sig^2
    U2 <- Z*rgamma(n, nu/2+1/2, nu/2+Chi/2) + (1-Z)*rgamma(n, nu/2, nu/2)
  }
  
  # Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Sig.pos <- Sig.pos[-om]
  return(list(Beta=Beta.pos, Sig=Sig.pos))
}








###  Regression with normal error distribution   ###
LR <- function(Y, X, mc=2000, burn=500){
  Int.pos <- c()
  Beta.pos <- matrix(NA,mc,p)
  sig.pos <- c()
  
  Int <- 1
  Beta <- rep(0,p)
  sig <- 1
  V <- rep(1,n)  
  cc <- 1   # prior for sigma^2
  
  for(k in 1:mc){
    # Int
    resid <- as.vector(Y-X%*%Beta)
    Int <- rnorm(1,sum(V*resid)/sum(V),sig/sqrt(sum(V)))
    Int.pos[k] <- Int
    # Beta
    Yt <- Y-Int
    W <- diag(V)
    invA <- ginv(t(X)%*%W%*%X/sig^2)
    Beta <- mvrnorm(1,invA%*%t(X)%*%W%*%Yt/sig^2,invA)
    Beta.pos[k,] <- Beta
    # sigma
    resid <- as.vector(Y-Int-X%*%Beta)
    sig2 <- rinvgamma(1,n/2+cc,sum(V*resid^2)/2+cc)
    sig <- sqrt(sig2)
    sig.pos[k] <- sig
  }
  
  # Summary
  om <- 1:burn
  Beta.pos <- cbind(Int.pos, Beta.pos)[-om,]
  Sig.pos <- sig.pos[-om]
  return(list(Beta=Beta.pos, Sig=Sig.pos))
}








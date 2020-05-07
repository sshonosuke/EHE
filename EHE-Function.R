library(glmnet)
library(MCMCpack)
library(statmod)
library(GIGrvg)



##-------------------------------------------------------------##
##     (Proposed) Robust regression with EHE distribution      ##
##-------------------------------------------------------------##
EHE <- function(Y, X, gam.est=F, mc=2000, burn=500){
  a <- 10^8    # tuning parameter for non-outlying error 
  delta <- 1   # hyperparameter for sigma
  om.a <- 1   # hyperparameter for omega
  om.b <- 1   # hyperparameter for omega
  XX <- cbind(1, X)
  n <- dim(XX)[[1]]
  pp <- dim(XX)[[2]]
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  om.pos <- c()
  if(gam.est){ Gam.pos <- c() }
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef(fit )
  res <- fit$residuals
  Sig <- sqrt( mean(res^2) )
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  W <- rep(1, n)
  U1 <- rep(1, n)
  U2 <- rep(5, n)
  Gam <- 1
  
  # MCMC
  for(k in 1:mc){
    UU <- U2*Z + U1*(1-Z)
    # Beta and Sigma
    mat <- solve(t(XX)%*%(XX/UU))
    Beta <- mvrnorm(1, mat%*%t(XX)%*%(Y/UU), Sig^2*mat)
    Beta.pos[k,] <- Beta
    mu <- as.vector(XX%*%Beta)
    resid <- Y - mu
    Sig <- sqrt( rinvgamma(1, delta/2+n/2, delta/2+sum(resid^2/UU)/2 ) )
    Sig.pos[k] <- Sig
    
    # U 
    Chi <- resid^2/Sig^2
    for(i in 1:n){ 
      U1[i] <- (1-Z[i])*rgig(1, lambda=a, chi=Chi[i], psi=2*a) + Z[i]*rgamma(1, a, a)
      U2[i] <- Z[i]*rgig(1, lambda=1/2, chi=Chi[i], psi=2*W[i]) + (1-Z[i])*rgamma(1, 1, W[i]) 
    }
    UU <- U2*Z + U1*(1-Z)
    # V
    V <- rgamma(n, 1+Gam, 1+log(1+U2))
    # W
    W  <- rgamma(n, 1+V, 1+U2)
    # omega
    om <- rbeta(1, om.a+sum(Z), om.b+n-sum(Z))
    om.pos[k] <- om 
    
    # Z
    mu <- as.vector(XX%*%Beta)
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    # Gamma
    if(gam.est){
      Gam <- rgamma(1, delta+n, delta+sum(log(1+log(1+U2))) )
      Gam.pos[k] <- Gam
    }
  }
  
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  Sig.pos <- Sig.pos[-omit]
  om.pos <- om.pos[-omit]
  Res <- list(Beta=Beta.pos, Sig=Sig.pos, Om=om.pos)
  return(Res)
}






##--------------------------------------------------------------------------##
##     (Proposed) Robust regression with EHE distribution and HS prior      ##
##--------------------------------------------------------------------------##
EHE.HS <- function(Y, X, gam.est=F, mc=2000, burn=500, om.hp=c(1,1)){
  a <- 10^8    # tuning parameter for non-outlying error 
  delta <- 1   # hyperparameter for sigma
  n <- dim(X)[1]
  p <- dim(X)[2]
  Beta.pos <- matrix(NA, mc, p)
  Alpha.pos <- c()
  Sig.pos <- c()
  om.pos <- c()
  if(gam.est){ Gam.pos <- c() }
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )[-1]
  Alpha <- coef( fit )[1]
  res <- fit$residuals
  Sig <- sqrt( mean(res^2) )
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  W <- rep(1, n)
  U1 <- rep(1, n)
  U2 <- rep(5, n)
  Gam <- 1
  Lam2 <- rep(1, p)   # local parameter in HS
  Nu <- rep(1, p)    # latent for Lam
  Tau2 <- 1    # scale parameter in HS
  Xi <- 1    # latent for Tau
  
  # MCMC
  for(k in 1:mc){
    UU <- U2*Z + U1*(1-Z)
    # Alpha
    Mu <- as.vector(X%*%Beta)
    A <- sum(1/UU)/Sig^2
    B <- sum((Y-Mu)/UU)/Sig^2
    Alpha <- rnorm(1, B/A, sqrt(1/A))
    Alpha.pos[k] <- Alpha
    # Beta 
    sY <- Y - Alpha
    invLmat <- diag(1/Lam2)/Tau2
    invAA <- solve( t(X)%*%(X/UU) + invLmat )
    Beta <- mvrnorm(1, invAA%*%t(X)%*%(sY/UU), Sig^2*invAA)
    Beta.pos[k,] <- Beta
    # Sigma
    Mu <- as.vector(X%*%Beta)
    resid <- sY - Mu
    Sig <- sqrt( rinvgamma(1, delta/2+(n+p)/2, delta/2+sum(resid^2/UU)/2+sum(Beta^2*diag(invLmat))/2 ) )
    Sig.pos[k] <- Sig
    # latent variables
    Lam2 <- rinvgamma(p, 1, 1/Nu + Beta^2/(2*Tau2*Sig^2))
    Tau2 <- rinvgamma(1, (p+1)/2, 1/Xi + sum(Beta^2/Lam2)/(2*Sig^2))
    Nu <- rinvgamma(p, 1, 1+1/Lam2)
    Xi <- rinvgamma(1, 1, 1+1/Tau2)
    # U 
    Chi <- resid^2/Sig^2
    for(i in 1:n){ 
      U1[i] <- (1-Z[i])*rgig(1, lambda=a, chi=Chi[i], psi=2*a) + Z[i]*rgamma(1, a, a)
      U2[i] <- Z[i]*rgig(1, lambda=1/2, chi=Chi[i], psi=2*W[i]) + (1-Z[i])*rgamma(1, 1, W[i]) 
    }
    UU <- U2*Z + U1*(1-Z)
    # V
    V <- rgamma(n, 1+Gam, 1+log(1+U2))
    # W
    W <- rgamma(n, 1+V, 1+U2)
    # omega
    om <- rbeta(1, om.hp[1]+sum(Z), om.hp[2]+n-sum(Z))
    om.pos[k] <- om 
    
    # Z
    mu <- as.vector(Alpha + X%*%Beta)
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    # Gamma
    if(gam.est){
      Gam <- rgamma(1, delta+n, delta+sum(log(1+log(1+U2))) )
      Gam.pos[k] <- Gam
    }
  }
  
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  Sig.pos <- Sig.pos[-omit]
  om.pos <- om.pos[-omit]
  Res <- list(Beta=Beta.pos, Alpha=Alpha.pos, Sig=Sig.pos, Om=om.pos)
  return(Res)
}





##--------------------------------------------------------------##
##      Robust regression with LPTN distribution                ##
##--------------------------------------------------------------##
LPTN <- function(Y, X, rho=0.8, mc=3000, burn=1000){
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
  
  # MCMC
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  Beta <- coef( lm(Y~X) )
  Sig <- sqrt( mean((Y-as.vector(XX%*%Beta))^2) )
  
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
  
  om <- 1:burn
  Beta.pos <- Beta.pos[-om, ]
  Sig.pos <- Sig.pos[-om]
  return(Beta.pos)
}










##--------------------------------------------------------------##
##       Robust regression with t distribution                  ##
##--------------------------------------------------------------##
tBR <- function(Y, X, mc=2000, burn=500, nu=3, estimation=F){
  Int.pos <- c()
  Beta.pos <- matrix(NA,mc,p)
  sig.pos <- c()
  
  cc=1   # prior for sigma^2 
  nu.grid <- c(1:5, 8, 10, 15, 20, 30, 50)    # grid for degrees of freedom 
  L <- length(nu.grid)
  nu.pos <- c()
  
  # initial value
  Int <- 1
  Beta <- rep(0,p)
  sig <- 1
  V <- rep(1,n)  
  
  # MCMC
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
  
  # Summary
  om <- 1:burn
  Res <- cbind(Int.pos, Beta.pos)[-om,]
  if(estimation){ attr(Res, "nu") <- nu.pos }
  return(Res)
}






##--------------------------------------------------------------##
##       Standard regression with t distribution                ##
##--------------------------------------------------------------##
BR <- function(Y, X, mc=2000, burn=500){
  Int.pos <- c()
  Beta.pos <- matrix(NA,mc,p)
  sig.pos <- c()
  
  # initial value
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
  Res <- cbind(Int.pos, Beta.pos)[-om,]
  return(Res)
}









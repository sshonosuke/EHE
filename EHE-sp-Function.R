library(glmnet)
library(MCMCpack)
library(statmod)
library(GIGrvg)



##-----------------------------------------------------##
##     Robust spatial linear regression with           ##     
##     EHE distribution and horseshoe prior            ##
##-----------------------------------------------------##
EHE.HS.sp <- function(Y, X, Sp, gam.est=F, mc=2000, burn=500, gam=1, om.hp=c(1,1), gam.hp=c(100,100)){
  ## settings
  delta <- 1   # hyperparameter for sigma and tau
  n <- dim(X)[1]
  p <- dim(X)[2]
  Distance <- as.matrix(dist(Sp))
  band <- quantile(Distance, prob=0.05)/10
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, p)
  Alpha.pos <- c()
  Sig.pos <- c()
  om.pos <- c()
  Gam.pos <- c() 
  Eta.pos <- matrix(NA, mc, n)    # spatial effect 
  h.pos <- c()   # bandwidth 
  Psi.pos <- c()    # variance for spatial effect
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )[-1]
  Alpha <- coef( fit )[1]
  res <- fit$residuals
  Sig <- 1
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  Eta <- rep(0, n)    # spatial effect
  Psi <- 1
  h <- band/3
  om <- mean(Z)   # contamination 
  W <- rep(1, n)
  U1 <- rep(1, n)   # fixed
  U2 <- rep(5, n)
  Gam <- gam
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
    sY <- Y - Mu - Eta
    B <- sum(sY/UU)/Sig^2
    Alpha <- rnorm(1, B/A, sqrt(1/A))
    Alpha.pos[k] <- Alpha
    # Beta 
    sY <- Y - Alpha - Eta
    invLmat <- diag(1/Lam2)/Tau2
    invAA <- solve( t(X)%*%(X/UU) + invLmat )
    Beta <- mvrnorm(1, invAA%*%t(X)%*%(sY/UU), Sig^2*invAA)
    Beta.pos[k,] <- Beta
    # Eta (spatial effect)
    sY <- Y - Alpha - as.vector(X%*%Beta)
    H <- exp(-Distance^2/(2*h^2)) 
    IH <- solve(H)
    mat <- solve( diag(1/UU)/Sig^2 + IH/Psi )
    mm <- as.vector(mat%*%(sY/UU))/Sig^2
    Eta <- mvrnorm(1, mm, mat)
    Eta.pos[k,] <- Eta
    # Psi 
    Psi <- rinvgamma(1, delta+n/2, delta+0.5*t(Eta)%*%IH%*%Eta)
    Psi.pos[k] <- Psi
    
    # h (random-walk MH)
    bb <- -sum(log(eigen(H)$values))
    new.h <- h + band*rnorm(1)
    new.h[new.h < 10^(-8)] <- 10^(-8)
    new.h[new.h > 0.01] <- 0.01
    new.H <- exp(-Distance^2/(2*new.h^2)) 
    new.IH <- solve(new.H)
    new.bb <- -sum(log(eigen(new.H)$values))
    val1 <- 0.5*bb - 0.5*t(Eta)%*%IH%*%Eta/Psi
    val2 <- 0.5*new.bb - 0.5*t(Eta)%*%new.IH%*%Eta/Psi
    prob <- min(1, exp(val2-val1))
    ch <- rbinom(1, 1, prob)
    h <- h + ch*(new.h-h)
    h.pos[k] <- h
    
    # Sigma
    Mu <- as.vector(X%*%Beta)
    resid <- Y - Alpha - Mu - Eta
    Sig <- sqrt( rinvgamma(1, delta+(n+p)/2, delta+sum(resid^2/UU)/2+sum(Beta^2*diag(invLmat))/2 ) )
    Sig.pos[k] <- Sig
    
    # latent variables (for HS)
    Lam2 <- rinvgamma(p, 1, 1/Nu + Beta^2/(2*Tau2*Sig^2))
    Tau2 <- rinvgamma(1, (p+1)/2, 1/Xi + sum(Beta^2/Lam2)/(2*Sig^2))
    Nu <- rinvgamma(p, 1, 1+1/Lam2)
    Xi <- rinvgamma(1, 1, 1+1/Tau2)
    
    # U 
    Chi <- resid^2/Sig^2
    for(i in 1:n){ 
      U1[i] <- 1
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
    mu <- as.vector(Alpha + X%*%Beta + Eta)
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    # Gamma
    if(gam.est){
      Gam <- rgamma(1, gam.hp[1]+n, gam.hp[2]+sum(log(1+log(1+U2))) )
    }
    Gam.pos[k] <- Gam
    
    if(round(k/100)==(k/100)){ print(k) }
  }
  
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Alpha.pos <- Alpha.pos[-om]
  Sig.pos <- Sig.pos[-om]
  Eta.pos <- Eta.pos[-om,]
  Psi.pos <- Psi.pos[-om]
  om.pos <- om.pos[-om]
  h.pos <- h.pos[-om]
  Gam.pos <- Gam.pos[-omit]
  Res <- list(Beta=Beta.pos, Alpha=Alpha.pos, Sig=Sig.pos, Om=om.pos, 
              Eta=Eta.pos, Psi=Psi.pos, h=h.pos, Gam=Gam.pos)
  return(Res)
}






##-----------------------------------------------------##
##     Robust spatial linear regression with           ##     
##     two-component mixture of t-distribution         ##
##                and horseshoe prior                  ##
##-----------------------------------------------------##
MT.HS.sp <- function(Y, X, Sp, nu=0.5, mc=2000, burn=500, om.hp=c(1,1)){
  ## settings
  delta <- 1   # hyperparameter for sigma and tau
  n <- dim(X)[1]
  p <- dim(X)[2]
  Distance <- as.matrix(dist(Sp))
  band <- quantile(Distance, prob=0.05)/10
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, p)
  Alpha.pos <- c()
  Sig.pos <- c()
  om.pos <- c()
  if(gam.est){ Gam.pos <- c() }
  Eta.pos <- matrix(NA, mc, n)    # spatial effect 
  h.pos <- c()   # bandwidth 
  Psi.pos <- c()    # variance for spatial effect
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )[-1]
  Alpha <- coef( fit )[1]
  res <- fit$residuals
  Sig <- 1
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  Eta <- rep(0, n)    # spatial effect
  Psi <- 1
  h <- band/3
  om <- mean(Z)   # contamination 
  U1 <- rep(1, n)   # non-outlier (fixed)
  U2 <- rep(5, n)   # outlier (sampled)
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
    sY <- Y - Mu - Eta
    B <- sum(sY/UU)/Sig^2
    Alpha <- rnorm(1, B/A, sqrt(1/A))
    Alpha.pos[k] <- Alpha
    # Beta 
    sY <- Y - Alpha - Eta
    invLmat <- diag(1/Lam2)/Tau2
    invAA <- solve( t(X)%*%(X/UU) + invLmat )
    Beta <- mvrnorm(1, invAA%*%t(X)%*%(sY/UU), Sig^2*invAA)
    Beta.pos[k,] <- Beta
    # Eta (spatial effect)
    sY <- Y - Alpha - as.vector(X%*%Beta)
    H <- exp(-Distance^2/(2*h^2)) 
    IH <- solve(H)
    mat <- solve( diag(1/UU)/Sig^2 + IH/Psi )
    mm <- as.vector(mat%*%(sY/UU))/Sig^2
    Eta <- mvrnorm(1, mm, mat)
    Eta.pos[k,] <- Eta
    # Psi 
    Psi <- rinvgamma(1, delta+n/2, delta+0.5*t(Eta)%*%IH%*%Eta)
    Psi.pos[k] <- Psi
    
    # h (random-walk MH)
    bb <- -sum(log(eigen(H)$values))
    new.h <- h + band*rnorm(1)
    new.h[new.h < 10^(-8)] <- 10^(-8)
    new.h[new.h > 0.01] <- 0.01
    new.H <- exp(-Distance^2/(2*new.h^2)) 
    new.IH <- solve(new.H)
    new.bb <- -sum(log(eigen(new.H)$values))
    val1 <- 0.5*bb - 0.5*t(Eta)%*%IH%*%Eta/Psi
    val2 <- 0.5*new.bb - 0.5*t(Eta)%*%new.IH%*%Eta/Psi
    prob <- min(1, exp(val2-val1))
    ch <- rbinom(1, 1, prob)
    h <- h + ch*(new.h-h)
    h.pos[k] <- h
    
    # Sigma
    Mu <- as.vector(X%*%Beta)
    resid <- Y - Alpha - Mu - Eta
    Sig <- sqrt( rinvgamma(1, delta/2+(n+p)/2, delta/2+sum(resid^2/UU)/2+sum(Beta^2*diag(invLmat))/2 ) )
    Sig.pos[k] <- Sig
    # latent variables (for HS)
    Lam2 <- rinvgamma(p, 1, 1/Nu + Beta^2/(2*Tau2*Sig^2))
    Tau2 <- rinvgamma(1, (p+1)/2, 1/Xi + sum(Beta^2/Lam2)/(2*Sig^2))
    Nu <- rinvgamma(p, 1, 1+1/Lam2)
    Xi <- rinvgamma(1, 1, 1+1/Tau2)
    # U 
    Chi <- resid^2/Sig^2
    U1 <- rep(1, n)
    U2 <- Z*rinvgamma(n, nu/2+1/2, nu/2+Chi/2) + (1-Z)*rinvgamma(n, nu/2, nu/2)
    UU <- U2*Z + U1*(1-Z)
    # omega
    om <- rbeta(1, om.hp[1]+sum(Z), om.hp[2]+n-sum(Z))
    om.pos[k] <- om 
    
    # Z
    mu <- as.vector(Alpha + X%*%Beta + Eta)
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(n, 1, pr)    # contamination indicator 
    
    if(round(k/100)==(k/100)){ print(k) }
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Alpha.pos <- Alpha.pos[-om]
  Sig.pos <- Sig.pos[-om]
  Eta.pos <- Eta.pos[-om,]
  Psi.pos <- Psi.pos[-om]
  om.pos <- om.pos[-om]
  h.pos <- h.pos[-om]
  Res <- list(Beta=Beta.pos, Alpha=Alpha.pos, Sig=Sig.pos, Om=om.pos, 
              Eta=Eta.pos, Psi=Psi.pos, h=h.pos)
  return(Res)
}














##-----------------------------------------------------##
##       Spatial linear regression with                ##     
##  normal error distribution and horseshoe prior      ##
##-----------------------------------------------------##
HS.sp <- function(Y, X, Sp, mc=3000, burn=1000){
  ## settings
  delta <- 1   # hyperparameter for sigma and tau
  n <- dim(X)[1]
  p <- dim(X)[2]
  Distance <- as.matrix(dist(Sp))
  band <- quantile(Distance, prob=0.05)/10
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, p)
  Alpha.pos <- c()
  Sig.pos <- c()
  Eta.pos <- matrix(NA, mc, n)    # spatial effect 
  h.pos <- c()   # bandwidth 
  Psi.pos <- c()    # variance for spatial effect
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )[-1]
  Alpha <- coef( fit )[1]
  Sig <- 1
  Eta <- rep(0, n)    # spatial effect
  Psi <- 1
  h <- band/3
  Lam2 <- rep(1, p)   # local parameter in HS
  Nu <- rep(1, p)    # latent for Lam
  Tau2 <- 1    # scale parameter in HS
  Xi <- 1    # latent for Tau
  
  # MCMC
  for(k in 1:mc){
    # Alpha
    Mu <- as.vector(X%*%Beta)
    A <- n/Sig^2
    sY <- Y - Mu - Eta
    B <- sum(sY)/Sig^2
    Alpha <- rnorm(1, B/A, sqrt(1/A))
    Alpha.pos[k] <- Alpha
    # Beta 
    sY <- Y - Alpha - Eta
    invLmat <- diag(1/Lam2)/Tau2
    invAA <- solve( t(X)%*%X + invLmat )
    Beta <- mvrnorm(1, invAA%*%t(X)%*%sY, Sig^2*invAA)
    Beta.pos[k,] <- Beta
    # Eta (spatial effect)
    sY <- Y - Alpha - as.vector(X%*%Beta)
    H <- exp(-Distance^2/(2*h^2)) 
    IH <- solve(H)
    mat <- solve( diag(n)/Sig^2 + IH/Psi )
    mm <- as.vector(mat%*%sY)/Sig^2
    Eta <- mvrnorm(1, mm, mat)
    Eta.pos[k,] <- Eta
    # Psi 
    Psi <- rinvgamma(1, delta+n/2, delta+0.5*t(Eta)%*%IH%*%Eta)
    Psi.pos[k] <- Psi
    
    # h (random-walk MH)
    bb <- -sum(log(eigen(H)$values))
    new.h <- h + band*rnorm(1)
    new.h[new.h < 10^(-8)] <- 10^(-8)
    new.h[new.h > 0.01] <- 0.01
    new.H <- exp(-Distance^2/(2*new.h^2)) 
    new.IH <- solve(new.H)
    new.bb <- -sum(log(eigen(new.H)$values))
    val1 <- 0.5*bb - 0.5*t(Eta)%*%IH%*%Eta/Psi
    val2 <- 0.5*new.bb - 0.5*t(Eta)%*%new.IH%*%Eta/Psi
    prob <- min(1, exp(val2-val1))
    ch <- rbinom(1, 1, prob)
    h <- h + ch*(new.h-h)
    h.pos[k] <- h
    
    # Sigma
    Mu <- as.vector(X%*%Beta)
    resid <- Y - Alpha - Mu - Eta
    Sig <- sqrt( rinvgamma(1, delta/2+(n+p)/2, delta/2+sum(resid^2)/2+sum(Beta^2*diag(invLmat))/2 ) )
    Sig.pos[k] <- Sig
    # latent variables (for HS)
    Lam2 <- rinvgamma(p, 1, 1/Nu + Beta^2/(2*Tau2*Sig^2))
    Tau2 <- rinvgamma(1, (p+1)/2, 1/Xi + sum(Beta^2/Lam2)/(2*Sig^2))
    Nu <- rinvgamma(p, 1, 1+1/Lam2)
    Xi <- rinvgamma(1, 1, 1+1/Tau2)
    
    if(round(k/100)==(k/100)){ print(k) }
  }
  
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Alpha.pos <- Alpha.pos[-om]
  Sig.pos <- Sig.pos[-om]
  Eta.pos <- Eta.pos[-om,]
  Psi.pos <- Psi.pos[-om]
  h.pos <- h.pos[-om]
  Res <- list(Beta=Beta.pos, Alpha=Alpha.pos, Sig=Sig.pos, Eta=Eta.pos, Psi=Psi.pos, h=h.pos)
  return(Res)
}





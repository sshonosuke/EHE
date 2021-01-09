library(glmnet)
library(MCMCpack)
library(statmod)
library(GIGrvg)




##-------------------------------------------------------------##
##    Robust random intercept model with EHE distribution      ##
##-------------------------------------------------------------##
EHE.RI <- function(Y, X, ID, gam.est=F, mc=2000, burn=500, gam=1, om.hp=c(1,1), gam.hp=c(100,100)){
  delta <- 1   # hyperparameter for variance parameter
  XX <- cbind(1, X)
  m <- length(unique(ID))
  n <- sum(ID==1)
  N <- dim(XX)[1]
  pp <- dim(XX)[[2]]
  Beta.pos <- matrix(NA, mc, pp)
  RE.pos <- matrix(NA, mc, m)
  Sig.pos <- c()
  Tau.pos <- c()
  om.pos <- c()
  Gam.pos <- c() 
  
  # initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )
  res <- fit$residuals
  Sig <- sqrt( mean(res^2) )
  Tau <- 0.5*Sig
  RE <- rep(0, m)
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  W <- rep(1, N)
  U1 <- rep(1, N)   # non-outlier (fixed)
  U2 <- rep(5, N)   # outlier (sampled)
  Gam <- gam
  
  C <- matrix(0, N, N)
  for(i in 1:m){
    ind <- (1+(i-1)*n):(n*i)
    C[ind, ind] <- 1
  }
  C <- as(C, "sparseMatrix")
  
  # MCMC
  for(k in 1:mc){
    UU <- (U2^Z) * (U1^(1-Z))
    DU <- as(diag(UU), "sparseMatrix")
    invSig <- solve(Sig^2*DU+Tau^2*C)
    mat <- solve(t(XX)%*%invSig%*%XX)
    mm <- as.vector( mat%*%t(XX)%*%invSig%*%Y )
    Beta <- mvrnorm(1, mm, mat)
    Beta.pos[k,] <- Beta
    mu <- as.vector(XX%*%Beta+RE[ID])
    resid <- Y - mu
    Sig <- sqrt( rinvgamma(1, delta+N/2, delta+sum(resid^2/UU)/2 ) )
    Sig.pos[k] <- Sig
    
    # RE and tau
    mU <- matrix(UU, n, m)
    A <- diag( 1/tau^2 + 1/sigma^2*apply(1/mU, 2, sum) )
    mY <- matrix((Y-as.vector(XX%*%Beta))/UU, n, m)
    B <- apply(mY, 2, sum)/Sig^2
    IA <- solve(A)
    RE <- rnorm(m, IA%*%B, IA)
    RE.pos[k,] <- RE
    Tau <- sqrt( rinvgamma(1, delta+m/2, delta+sum(RE^2)/2 ) )
    Tau.pos[k] <- Tau
    
    # Z
    mu <- as.vector(XX%*%Beta + RE[ID])
    logdens1 <- dnorm(Y, mu, sqrt(U1)*Sig, log=T)
    logdens2 <- dnorm(Y, mu, sqrt(U2)*Sig, log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(N, 1, pr)    # contamination indicator 
    
    # omega
    om <- rbeta(1, om.hp[1]+sum(Z), om.hp[2]+N-sum(Z))
    om.pos[k] <- om 
    
    # Gamma
    if(gam.est){
      Gam <- rgamma(1, gam.hp[1]+N, gam.hp[2]+sum(log(1+log(1+U2))) )
    }
    Gam.pos[k] <- Gam
    
    # V and W
    W <- rgamma(N, 1+Gam, 1+log(1+U2))
    V <- rgamma(N, 1+W, 1+U2)
    
    # UU 
    Chi <- resid^2/Sig^2
    for(i in 1:N){ 
      U2[i] <- Z[i]*rgig(1, lambda=1/2, chi=Chi[i], psi=2*V[i]) + (1-Z[i])*rgamma(1, 1, V[i])
    }
  }
  
  ## Summary
  omit <- 1:burn
  Beta.pos <- Beta.pos[-omit,]
  RE.pos <- RE.pos[-omit,]
  Tau.pos <- Tau.pos[-omit]
  Sig.pos <- Sig.pos[-omit]
  om.pos <- om.pos[-omit]
  Gam.pos <- Gam.pos[-omit] 
  Result <- list(Beta=Beta.pos, RE=RE.pos, Tau=Tau.pos, Sig=Sig.pos, Om=om.pos, Gam=Gam.pos)
  return(Result)
}





##--------------------------------------------------------------##
##    Robust random intercept model with t-distribution         ##
##--------------------------------------------------------------##
TBR.RI <- function(Y, X, ID, mc=2000, burn=500, nu=3, estimation=F){
  ## settings
  m <- length(unique(ID))
  n <- sum(ID==1)
  N <- length(Y)
  
  ## MCMC box
  Int.pos <- c()
  Beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  tau.pos <- c()
  RE.pos <- matrix(NA, mc, m)
  
  # initial values
  Int <- 1
  Beta <- rep(0,p)
  sig <- 1
  tau <- 0.5
  V <- rep(1,N) 
  RE <- rep(0, m)
  
  cc <- 1   # prior for sigma^2
  nu.grid <- c(1:5, 8, 10, 15, 20, 30, 50)
  L <- length(nu.grid)
  nu.pos <- c()
  
  ## MCMC
  for(k in 1:mc){
    # Int
    resid <- as.vector(Y-X%*%Beta-RE[ID])
    Int <- rnorm(1,sum(V*resid)/sum(V),sig/sqrt(sum(V)))
    Int.pos[k] <- Int
    
    # Beta
    Yt <- Y-Int-RE[ID]
    W <- diag(V)
    invA <- ginv(t(X)%*%W%*%X/sig^2)
    Beta <- mvrnorm(1,invA%*%t(X)%*%W%*%Yt/sig^2,invA)
    Beta.pos[k,] <- Beta
    
    # sigma
    resid <- as.vector(Y-Int-X%*%Beta-RE[ID])
    sig2 <- rinvgamma(1,N/2+cc,sum(V*resid^2)/2+cc)
    sig <- sqrt(sig2)
    sig.pos[k] <- sig
    
    # V
    V <- rgamma(N,nu/2+1/2,nu/2+resid^2/(2*sig^2))
    
    # RE
    mV <- matrix(V, n, m)
    A <- diag( 1/tau^2 + apply(mV, 2, sum)/sig^2 )
    mY <- matrix((Y-Int-as.vector(X%*%Beta))*V, n, m)
    B <- apply(mY, 2, sum)/sig^2
    IA <- solve(A)
    RE <- rnorm(m, IA%*%B, IA)
    RE.pos[k,] <- RE
    
    # tau
    tau <- sqrt( rinvgamma(1, cc+m/2, cc+sum(RE^2)/2 ) )
    tau.pos[k] <- tau
    
    # nu (if estimated)
    if(estimation){
      logdens <- matrix(NA, N, L)
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
  sig.pos <- sig.pos[-om]
  tau.pos <- tau.pos[-om]
  RE.pos <- RE.pos[-om,]
  Res <- list(Beta=Beta.pos, RE=RE.pos, tau=tau.pos, sig=sig.pos)
  return(Res)
}





##--------------------------------------------------------------##
##    Robust random intercept model with two-component          ##
##         mixture of normal and t-distribution                 ##
##--------------------------------------------------------------##
MTBR.RI <- function(Y, X, ID, nu=0.5, mc=3000, burn=1000){
  ## settings
  delta <- 1   # hyperparameter for sigma
  XX <- cbind(1, X)
  m <- length(unique(ID))
  n <- sum(ID==1)
  N <- dim(XX)[1]
  pp <- dim(XX)[[2]]
  
  ## MCMC box
  Beta.pos <- matrix(NA, mc, pp)
  Sig.pos <- c()
  om.pos <- c()
  RE.pos  <- matrix(NA, mc, m)
  Tau.pos <- c()
  
  ## initial value
  fit <- lm(Y~X)
  Beta <- coef( fit )
  res <- fit$residuals
  Sig <- sqrt( mean(res^2) )
  Tau <- 0.5*Sig
  RE <- rep(0, m)
  st.res <- res/Sig
  Z <- ifelse(abs(st.res)>2, 1, 0)
  om <- mean(Z)   # contamination 
  U1 <- rep(1, N)   # non-outlier (fixed)
  U2 <- rep(5, N)   # outlier (sampled)
  
  ## MCMC
  for(k in 1:mc){
    UU <- (U2^Z) * (U1^(1-Z))
    # Beta and Sigma
    mat <- solve(t(XX)%*%(XX*UU))
    Beta <- mvrnorm(1, mat%*%t(XX)%*%((Y-RE[ID])*UU), Sig^2*mat)
    Beta.pos[k,] <- Beta
    mu <- as.vector(XX%*%Beta)
    resid <- Y - mu - RE[ID]
    Sig <- sqrt( rinvgamma(1, delta/2+N/2, delta/2+sum(resid^2*UU)/2 ) )
    Sig.pos[k] <- Sig
    
    # RE and tau
    mU <- matrix(UU, n, m)
    A <- diag( 1/Tau^2 + apply(mU, 2, sum)/Sig^2 )
    mY <- matrix((Y-as.vector(XX%*%Beta))*UU, n, m)
    B <- apply(mY, 2, sum)/Sig^2
    IA <- solve(A)
    RE <- rnorm(m, IA%*%B, IA)
    RE.pos[k,] <- RE
    Tau <- sqrt( rinvgamma(1, delta/2+m/2, delta/2+sum(RE^2)/2 ) )
    Tau.pos[k] <- Tau
    
    # Z
    mu <- as.vector(XX%*%Beta)
    logdens1 <- dnorm(Y, mu, Sig/sqrt(U1), log=T)
    logdens2 <- dnorm(Y, mu, Sig/sqrt(U2), log=T)
    pr <- 1 / ( exp( log(1-om)+logdens1-log(om)-logdens2 ) + 1 )
    Z <- rbinom(N, 1, pr)    # contamination indicator 
    
    # omega
    om <- rbeta(1, 1/2+sum(Z), 1/2+N-sum(Z))
    om.pos[k] <- om 
    
    # UU 
    resid <- Y - mu - RE[ID]
    Chi <- resid^2/Sig^2
    U2 <- Z*rgamma(N, nu/2+1/2, nu/2+Chi/2) + (1-Z)*rgamma(N, nu/2, nu/2)
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,]
  Sig.pos <- Sig.pos[-om]
  Tau.pos <- Tau.pos[-om]
  RE.pos <- RE.pos[-om,]
  Res <- list(Beta=Beta.pos, RE=RE.pos, Tau=Tau.pos, Sig=Sig.pos)
  return(Res)
}






##-----------------------------------------------------------##
##      random intercept model with normal error             ##
##-----------------------------------------------------------##
BR.RI <- function(Y, X, ID, mc=3000, burn=1000){
  ## settings
  m <- length(unique(ID))
  n <- sum(ID==1)
  N <- length(Y)
  
  ## MCMC box
  Int.pos <- c()
  Beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  tau.pos <- c()
  RE.pos <- matrix(NA, mc, m)
  
  ## initial values
  Int <- 1
  Beta <- rep(0,p)
  sig <- 1
  tau <- 0.5
  RE <- rep(0, m)
  V <- rep(1,N)  
  cc <- 1   # prior for sigma^2
  
  ## MCMC
  for(k in 1:mc){
    # Int
    resid <- as.vector(Y-X%*%Beta-RE[ID])
    Int <- rnorm(1,sum(V*resid)/sum(V),sig/sqrt(sum(V)))
    Int.pos[k] <- Int
    
    # Beta
    Yt <- Y-Int-RE[ID]
    W <- diag(V)
    invA <- ginv(t(X)%*%W%*%X/sig^2)
    Beta <- mvrnorm(1,invA%*%t(X)%*%W%*%Yt/sig^2,invA)
    Beta.pos[k,] <- Beta
    
    # sigma
    resid <- Y-Int-as.vector(X%*%Beta)-RE[ID]
    sig2 <- rinvgamma(1,N/2+cc,sum(V*resid^2)/2+cc)
    sig <- sqrt(sig2)
    sig.pos[k] <- sig
    
    # RE
    mV <- matrix(V, n, m)
    A <- diag( 1/tau^2 + apply(mV, 2, sum)/sig^2 )
    mY <- matrix((Y-Int-as.vector(X%*%Beta))*V, n, m)
    B <- apply(mY, 2, sum)/sig^2
    IA <- solve(A)
    RE <- rnorm(m, IA%*%B, IA)
    RE.pos[k,] <- RE
    
    # tau
    tau <- sqrt( rinvgamma(1, cc+m/2, cc+sum(RE^2)/2 ) )
    tau.pos[k] <- tau
  }
  
  ## Summary
  om <- 1:burn
  Beta.pos <- cbind(Int.pos, Beta.pos)[-om,]
  Sig.pos <- sig.pos[-om]
  Tau.pos <- tau.pos[-om]
  RE.pos <- RE.pos[-om,]
  Res <- list(Beta=Beta.pos, RE=RE.pos, Tau=Tau.pos, Sig=Sig.pos)
  return(Res)
}







library(MASS)
library(MCMCpack)
library(bayesm)

source("Function-RI.R")
set.seed(123)

R <- 500    # number of Monte Carlo replications 

om <- 0.05  # contamination ratio
aa <- 15    # outlier location: 5 or 10 or 15 or 20

p <- 10    # the number of predictor variables: 10 or 20
m <- 50   # the number of subjects
nn <- 10
N <- m*nn
ID <- rep(1:m, rep(nn, m))
mc <- 3000 
bn <- 1000




# settings
tau <- 0.5
sigma <- 1
Beta <- rep(0,p)  
Beta[c(1,4,7,10)] <- c(0.3, 0.3, 2, 2)   # non-zero coeffieicnts
int <- 0.5   # intercept
Para <- c(int, Beta)

quant <- function(x){ quantile(x, prob=c(0.025,0.975)) }



# covariates
rho <- 0.2
mat <- matrix(NA,p,p)
for(k in 1:p){
  for(j in 1:p){ mat[k,j] <- rho^(abs(k-j)) }
}
X <- mvrnorm(N, rep(0, p), mat) 



# result box
meth <- c("EHE", "aEHE","Cauchy", "aT", "Mix-T", "Normal")
L <- length(meth)
MSE <- matrix(NA, R, L)
MSE.RE <- matrix(NA, R, L)
CP <- matrix(NA, R, L)
AL <- matrix(NA, R, L)
dimnames(MSE)[[2]] <- dimnames(MSE.RE)[[2]] <- meth
dimnames(CP)[[2]] <- dimnames(AL)[[2]] <- meth

data <- list()

# replication 
for(r in 1:R){
  ch <- rbinom(N, 1, om)
  noise <- (1-ch)*rnorm(N, 0, 1) + ch*rnorm(N, aa, 1)
  RE <- rnorm(m, 0, tau)
  Y <- int + as.vector(X%*%Beta) + RE[ID] + sigma*noise
  plot(X[,1], Y)
  data[[r]] <- Y
  
  # EHE
  fit.EHE <- EHE.RI(Y, X, ID, gam.est=F, mc=mc, burn=bn)
  est.EHE <- apply(fit.EHE$Beta, 2, mean)
  CI.EHE <- apply(fit.EHE$Beta, 2, quant)
  RE.EHE <- apply(fit.EHE$RE, 2, mean)
  
  # aEHE
  fit.aEHE <- EHE.RI(Y, X, ID, gam.est=T, mc=mc, burn=bn)
  plot(fit.aEHE$Gam)
  est.aEHE <- apply(fit.aEHE$Beta, 2, mean)
  CI.aEHE <- apply(fit.aEHE$Beta, 2, quant)
  RE.aEHE <- apply(fit.aEHE$RE, 2, mean)
  
  # Cauchy
  fit.c <- tBR.RI(Y, X, ID, mc=mc, burn=bn, nu=1)
  est.c <- apply(fit.c$Beta, 2, mean)
  CI.c <- apply(fit.c$Beta, 2, quant)
  RE.c <- apply(fit.c$RE, 2, mean)
  
  # T (estimated)
  fit.at <- tBR.RI(Y, X, ID, mc=mc, burn=bn, nu=3, estimation=T)
  est.at <- apply(fit.at$Beta, 2, mean)
  CI.at <- apply(fit.at$Beta, 2, quant)
  RE.at <- apply(fit.at$RE, 2, mean)
  
  # mix-T 
  fit.mixt <- mix.tBR.RI(Y, X, ID, mc=mc, burn=bn, nu=1/2)
  est.mixt <- apply(fit.mixt$Beta, 2, mean)
  CI.mixt <- apply(fit.mixt$Beta, 2, quant)
  RE.mixt <- apply(fit.mixt$RE, 2, mean)
  
  # Linear regression (Benchmark)
  fit.lm <- BR.RI(Y, X, ID, mc=mc, burn=bn)
  est.lm <- apply(fit.lm$Beta, 2, mean)
  CI.lm <- apply(fit.lm$Beta, 2, quant)
  RE.lm <- apply(fit.lm$RE, 2, mean)
  
  # MSE
  Est <- cbind(est.EHE, est.aEHE, est.c, est.at, est.mixt, est.lm)
  MSE[r,] <- apply((Est-Para)^2, 2, mean) 
  
  # MSE-RE
  Est <- cbind(RE.EHE, RE.aEHE, RE.c, RE.at, RE.mixt, RE.lm)
  MSE.RE[r,] <- apply((Est-RE)^2, 2, mean) 
  
  # Coverage prob
  CP[r,1] <- mean( ifelse(CI.EHE[1,]<Para & CI.EHE[2,]>Para, 1, 0) )
  CP[r,2] <- mean( ifelse(CI.aEHE[1,]<Para & CI.aEHE[2,]>Para, 1, 0) )
  CP[r,3] <- mean( ifelse(CI.c[1,]<Para & CI.c[2,]>Para, 1, 0) )
  CP[r,4] <- mean( ifelse(CI.at[1,]<Para & CI.at[2,]>Para, 1, 0) )
  CP[r,5] <- mean( ifelse(CI.mixt[1,]<Para & CI.mixt[2,]>Para, 1, 0) )
  CP[r,6] <- mean( ifelse(CI.lm[1,]<Para & CI.lm[2,]>Para, 1, 0) )
  
  # Averae length
  AL[r,1] <- mean(CI.EHE[2,]-CI.EHE[1,])
  AL[r,2] <- mean(CI.aEHE[2,]-CI.aEHE[1,])
  AL[r,3] <- mean(CI.c[2,]-CI.c[1,])
  AL[r,4] <- mean(CI.at[2,]-CI.at[1,])
  AL[r,5] <- mean(CI.mixt[2,]-CI.mixt[1,])
  AL[r,6] <- mean(CI.lm[2,]-CI.lm[1,])
  
  if(round(r/100)==(r/100)){ print(r) }
  print(100*MSE.RE[r,])
}


save(list=ls(), file=paste0("sim-RI(a=",aa,",om=",om,").RData"))




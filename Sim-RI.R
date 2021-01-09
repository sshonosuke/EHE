rm(list=ls())
set.seed(1234)

# load R packages
library(MASS)
library(MCMCpack)
library(bayesm)

# load R function
source("EHE-RI-Function.R")

quant <- function(x){ quantile(x, prob=c(0.025, 0.975)) }


# settings
p <- 10    # the number of predictor variables
m <- 50   # the number of subjects
nn <- 10   # the number of repeated measurements
aa <- 10    # outlier location: 5 or 10 or 15 or 20
om <- 0.05  # contamination ratio:  0 or 0.05 or 0.1 ot 0.15
mc <- 3000  # length of MCMC
bn <- 1000  # length of burn-in 
N <- m*nn
ID <- rep(1:m, rep(nn, m))


# true parameter
tau <- 0.5
sigma <- 1
Beta <- rep(0,p)  
Beta[c(1,4,7,10)] <- c(0.3, 0.3, 2, 2)   # non-zero coeffieicnts
int <- 0.5   # intercept
Para <- c(int, Beta)


# covariates
rho <- 0.2
mat <- matrix(NA,p,p)
for(k in 1:p){
  for(j in 1:p){ mat[k,j] <- rho^(abs(k-j)) }
}
X <- mvrnorm(N, rep(0, p), mat) 


# data generation
ch <- rbinom(N, 1, om)
noise <- (1-ch)*rnorm(N, 0, 1) + ch*rnorm(N, aa, 1)
RE <- rnorm(m, 0, tau)
Y <- int + as.vector(X%*%Beta) + RE[ID] + sigma*noise



# proposed method (fixed gamma)
fit.EHE <- EHE.RI(Y, X, ID, gam.est=F, mc=mc, burn=bn)
est.EHE <- apply(fit.EHE$Beta, 2, mean)
CI.EHE <- apply(fit.EHE$Beta, 2, quant)
RE.EHE <- apply(fit.EHE$RE, 2, mean)

# proposed method (estimated gamma)
fit.aEHE <- EHE.RI(Y, X, ID, gam.est=T, mc=mc, burn=bn)
est.aEHE <- apply(fit.aEHE$Beta, 2, mean)
CI.aEHE <- apply(fit.aEHE$Beta, 2, quant)
RE.aEHE <- apply(fit.aEHE$RE, 2, mean)

# Cauchy error 
fit.c <- TBR.RI(Y, X, ID, mc=mc, burn=bn, nu=1)
est.c <- apply(fit.c$Beta, 2, mean)
CI.c <- apply(fit.c$Beta, 2, quant)
RE.c <- apply(fit.c$RE, 2, mean)

# t-distribution error
fit.at <- TBR.RI(Y, X, ID, mc=mc, burn=bn, nu=3, estimation=T)
est.at <- apply(fit.at$Beta, 2, mean)
CI.at <- apply(fit.at$Beta, 2, quant)
RE.at <- apply(fit.at$RE, 2, mean)

# two-component of normal and t-distribution  
fit.mixt <- MTBR.RI(Y, X, ID, mc=mc, burn=bn, nu=1/2)
est.mixt <- apply(fit.mixt$Beta, 2, mean)
CI.mixt <- apply(fit.mixt$Beta, 2, quant)
RE.mixt <- apply(fit.mixt$RE, 2, mean)

# normal error 
fit.lm <- BR.RI(Y, X, ID, mc=mc, burn=bn)
est.lm <- apply(fit.lm$Beta, 2, mean)
CI.lm <- apply(fit.lm$Beta, 2, quant)
RE.lm <- apply(fit.lm$RE, 2, mean)


# estimates of regression coeffieicnts
Est <- cbind(est.EHE, est.aEHE, est.c, est.at, est.mixt, est.lm)
100*sqrt(apply((Est-Para)^2, 2, mean))   # MSE
  
# estimates of random effects
hRE <- cbind(RE.EHE, RE.aEHE, RE.c, RE.at, RE.mixt, RE.lm)
100*sqrt(apply((hRE-RE)^2, 2, mean))    # MSE
  
# credible intervals
CI.EHE
CI.aEHE
CI.c
CI.at
CI.mixt
CI.lm

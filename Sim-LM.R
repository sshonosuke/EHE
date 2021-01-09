rm(list=ls())
set.seed(1234)


# load R packages
library(MASS)
library(MCMCpack)
library(bayesm)

# load R function
source("EHE-LM-Function.R")


# settings
p <- 20    # the number of predictor variables
n <- 300   # the number of observations
aa <- 10    # outlier location: 5 or 10 or 15 or 20
om <- 0.1  # contamination ratio:  0 or 0.05 or 0.1 ot 0.15
mc <- 3000  # length of MCMC
bn <- 1000  # length of burn-in 


# true parameter
sigma <- 0.5   # error standard deviation 
Beta <- rep(0,p)  
Beta[c(1,4,7,10)] <- c(0.3, 0.3, 2, 2)   # non-zero coeffieicnts
int <- 0.5   # intercept

# covariates
rho <- 0.2
mat <- matrix(NA,p,p)
for(k in 1:p){
  for(j in 1:p){ mat[k,j] <- rho^(abs(k-j)) }
}
X <- mvrnorm(n, rep(0, p), mat) 


# function to compute IF
IF.compute <- function(x){ as.vector(numEff(x)$f) }

# data generation
ch <- rbinom(n, 1, om)
noise <- (1-ch)*rnorm(n, 0, 1) + ch*rnorm(n, aa, 1)
Y <- int + as.vector(X%*%Beta) + noise


# fitting the proposed methods
fit.EHE <- EHE(Y, X, gam.est=F, gam=1, mc=mc, burn=bn)
fit.aEHE <- EHE(Y, X, gam.est=T, mc=mc, burn=bn, gam.hp=c(100, 100))

# fitting the alternative methods
fit.lp1 <- LPTN(Y, X, rho=0.9, mc=mc, burn=bn)
fit.lp2 <- LPTN(Y, X, rho=0.7, mc=mc, burn=bn)
fit.t <- TBR(Y, X, mc=mc, burn=bn, nu=3)
fit.c <- TBR(Y, X, mc=mc, burn=bn, nu=1)
fit.at <- TBR(Y, X, mc=mc, burn=bn, nu=3, estimation=T)
fit.mt <- MTBR(Y, X, mc=mc, burn=bn)
fit.lm <- BR(Y, X, mc=mc, burn=bn)

   
# posterior mean
est.EHE <- apply(fit.EHE$Beta, 2, mean)
est.aEHE <- apply(fit.aEHE$Beta, 2, mean)
est.lp1 <- apply(fit.lp1$Beta, 2, mean)
est.lp2 <- apply(fit.lp2$Beta, 2, mean)
est.t <- apply(fit.t$Beta, 2, mean)
est.c <- apply(fit.c$Beta, 2, mean)
est.at <- apply(fit.at$Beta, 2, mean)
est.mt <- apply(fit.mt$Beta, 2, mean)
est.lm <- apply(fit.lm$Beta, 2, mean)


# averaged squared error loss for the regression parameters
Para <- c(int, Beta)
Est <- cbind(est.EHE, est.aEHE, est.lp1, est.lp2, est.t, est.c, est.at, est.mt, est.lm)
100*sqrt( apply((Est-Para)^2, 2, mean) )  


# credible interval
CI95 <- function(x){ quantile(x, prob=c(0.025, 0.975)) }
CI.EHE <- apply(fit.EHE$Beta, 2, CI95)
CI.aEHE <- apply(fit.EHE$Beta, 2, CI95)
CI.lp1 <- apply(fit.lp1$Beta, 2, CI95)
CI.lp2 <- apply(fit.lp2$Beta, 2, CI95)
CI.t <- apply(fit.t$Beta, 2, CI95)
CI.c <- apply(fit.c$Beta, 2, CI95)
CI.at <- apply(fit.at$Beta, 2, CI95)
CI.mt <- apply(fit.mt$Beta, 2, CI95)
CI.lm <- apply(fit.lm$Beta, 2, CI95)
 
# average IF
mean( apply(fit.EHE$Beta, 2, IF.compute) )
mean( apply(fit.aEHE$Beta, 2, IF.compute) )
mean( apply(fit.lp1$Beta, 2, IF.compute) )
mean( apply(fit.lp2$Beta, 2, IF.compute) )
mean( apply(fit.c$Beta, 2, IF.compute) )
mean( apply(fit.t$Beta, 2, IF.compute) )
mean( apply(fit.at$Beta, 2, IF.compute) )
mean( apply(fit.mt$Beta, 2, IF.compute) )
mean( apply(fit.lm$Beta, 2, IF.compute) )


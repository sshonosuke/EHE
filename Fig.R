##--------------------------------------##
##    Figures for density fucntions     ##
##--------------------------------------##
EHE.dens <- function(s=0.1, gam=1, ep=seq(-4, 4, by=0.1)){
  a <- 10^8
  mc <- 10000
  W <- rgamma(mc, gam, 1)
  V <- rgamma(mc, W, 1)
  V[V<0.001] <- 0.001
  U2 <- rgamma(mc, 1, V)
  U1 <- rgamma(mc, a, a)
  Z <- rbinom(mc, 1, s)
  UU <- (1-Z)*U1 + Z*U2
  
  L <- length(ep)
  U.mat <- matrix(rep(UU, L), mc, L)
  ep.mat <- t( matrix(rep(ep, mc), L, mc) )
  dens <- apply(dnorm(ep.mat, 0, sqrt(U.mat)), 2, mean)
  return(dens)
}


ep <- seq(-10, 10, by=0.1)
s.set <- c(0.05, 0.1, 0.2)
J <- length(s.set)

dens <- matrix(NA, length(ep), J)
for(k in 1:J){
  dens[,k] <- EHE.dens(s=s.set[k], gam=1, ep=ep)
}

N.dens <- dnorm(ep)
Dens <- cbind(N.dens, dens)


# figure 1 in the main document
pdf("density.pdf", width=12, height=8, pointsize=13)
par(mfcol=c(1,2))
matplot(ep[abs(ep)<4], Dens[abs(ep)<4,], type="l", lty=1, col=1:(J+1), xlab="x", ylab="density")
legend("topright", legend=c("Normal", paste0("EHE (s=",s.set,")")), lty=1, col=1:(J+1))
matplot(ep[ep>4], Dens[ep>4,], type="l", lty=1, col=1:(J+1), xlab="x", ylab="density")
legend("topright", legend=c("Normal", paste0("EHE (s=",s.set,")")), lty=1, col=1:(J+1))
dev.off()




##--------------------------------------##
##          Figures for CDF             ##
##--------------------------------------##
H.dist <- function(gam=1, ep=seq(0, 4, by=0.1)){
  dist <- 1-( 1 + log(1+ep) )^(-gam)
}


rEHE <- function(mc=100000, s=0.1, gam=1){
  a <- 10^8
  W <- rgamma(mc, gam, 1)
  V <- rgamma(mc, W, 1)
  V[V<0.001] <- 0.001
  U2 <- rgamma(mc, 1, V)
  U1 <- rgamma(mc, a, a)
  Z <- rbinom(mc, 1, s)
  UU <- (1-Z)*U1 + Z*U2
  X <- rnorm(mc,0,sqrt(UU))
  return(X)
}


ep <- seq(0, 10, by=0.01)
s.set <- c(0.5, 1, 2)
J <- length(s.set)

dist <- matrix(NA, length(ep), J+1)
dist[,1] <- 1 - pgamma(1/ep, 0.5, rate = 0.5)
for(k in 1:J){
  dist[,k+1] <- H.dist(gam=s.set[k], ep=ep)
}


ep2 <- seq(-10, 10, by=0.1)
s.set2 <- c(0.1, 0.5, 0.8)
J2 <- length(s.set2)

DIST <- matrix(NA, length(ep2), J2+1)
DIST[,1] <-  pcauchy(ep2) #pnorm(ep2)
for(k in 1:J2){
  DIST[,k+1] <- ecdf(rEHE(s=s.set2[k], gam=1))(ep2)
}


Lab1 <- c("IG(1/2,1/2)", expression(H (gamma==0.5)), expression(H (gamma==1)), expression(H (gamma==2)))
Lab2 <- c("Cauchy", expression(EHE (s==0.1)), expression(EHE (s==0.5)), expression(EHE (s==0.8)))


# figure 2 in the main document
pdf("dist.pdf", width=12, height=8, pointsize=13)
par(mfcol=c(1,2))
matplot(ep[abs(ep)<6], dist[abs(ep)<6,], type="l", lty=1, col=1:(J+1), xlab="x", ylab="CDF")
legend("bottomright", Lab1, lty=1, col=1:(J+1))
matplot(ep2[abs(ep2)<8], DIST[abs(ep2)<8,], type="l", lty=1, col=1:(J2+1), xlab="x", ylab="CDF")
legend("bottomright", Lab2, lty=1, col=1:(J2+1))
dev.off()



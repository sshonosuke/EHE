## It will take a considerable time to run this code.

rm(list=ls())
set.seed(1)

source("EHE-sp-function.R")
MC <- 5000    # length of MCMC
Bn <- 1000    # length of burn-in 


## data
library(mlbench)
data(BostonHousing2)
data <- BostonHousing2


Lab <- c("lon","lat","crim","zn","indus","chas","nox","rm","age","dis","rad","tax","ptratio","b","lstat")
Y <- data$cmedv
X <- subset(data, select=Lab)
n <- length(Y)
pp <- dim(X)[2]

Sp <- cbind(data$lon, data$lat)   # location information 

X <- matrix(as.numeric(as.matrix(X)), n, pp)
X <- apply(X, 2, scale)
dimnames(X)[[2]] <- Lab
EX <- cbind(X, X[,-6]^2)


## (proposed) robust method 
rfit <- EHE.HS.sp(Y, EX, Sp, mc=MC, burn=Bn)

# trace plot of mixing proportion
plot(rfit$Om, type="l", ylab="mixing proportion", xlab="replication")

## non-robust method
fit <- HS.sp(Y, EX, Sp, mc=MC, burn=Bn)

## (alternative) robust method
mt.fit <- MT.HS.sp(Y, EX, Sp, mc=MC, burn=Bn)





## function for spatial plot
Plot <- function(Sp, value, ran=NULL, title=""){ 
  if(is.null(ran)){ ran <- range(value) }
  value <- (value - ran[1]) / diff(ran)
  cs <- colorRamp( c("blue", "green", "yellow", "red"), space="rgb")
  cols <- rgb( cs(value), maxColorValue=256 )
  plot(Sp, col=cols, xlim=c(-71.3, -70.6), xlab="Latitude", ylab="Longitude", main=title, pch=20, cex=1)
  ran[1] <- floor(ran[1])
  ran[2] <- ceiling(ran[2])
  cs <- colorRamp( c("blue", "green", "yellow", "red"), space="rgb")
  cols <- rgb(cs(0:1000/1000), maxColorValue=256)
  rect(-70.75, seq(42.05, 42.35, length=1001), -70.6, seq(42.05, 42.35, length=1001), col=cols, border=cols)
  tx <- seq(ran[1], ran[2], length=5)
  text(x=-70.75, y=seq(42.05, 42.35, length=5), tx, cex=0.7)
  yy <- seq(42.05, 42.35, length=5)
  for (i in 1:5){
    segments(-70.75, yy[i], -70.6, yy[i], col="white")
  }
}


## estimates of spatial effects
eta1 <- apply(rfit$Eta, 2, mean)
eta2 <- apply(fit$Eta, 2, mean)
eta3 <- apply(mt.fit$Eta, 2, mean)


ran <- c(-7, 7)
eta2[eta2<ran[1]] <- ran[1]
eta2[eta2>ran[2]] <- ran[2]
eta1[eta1<ran[1]] <- ran[1]
eta1[eta1>ran[2]] <- ran[2]


##  Figure 4 in the main document
pdf("spatial-effect.pdf", width=14, height=8, pointsize=13)
par(mfcol=c(1,2))
Plot(Sp, eta1, ran=ran, title="EH")
Plot(Sp, eta2, ran=ran, title="N")
dev.off()



## coefficients
est0 <- apply(fit$Beta, 2, mean)[-1]
est1 <- apply(rfit$Beta, 2, mean)[-1]
est2 <- apply(mt.fit$Beta, 2, mean)[-1]
Est <- cbind(est0, est1, est2)

quant <- function(x){ quantile(x, prob=c(0.025, 0.975)) }
CI0 <- apply(fit$Beta, 2, quant)[,-1]
CI1 <- apply(rfit$Beta, 2, quant)[,-1]
CI2 <- apply(mt.fit$Beta, 2, quant)[,-1]

CI <- list()
CI[[1]] <- CI0 
CI[[2]] <- CI1
CI[[3]] <- CI2
ran <- range(CI)
p <- dim(CI0)[2]

bb <- seq(-0.15, 0.15, length=3)
col <- c(1, 2, 4)


##  Figure 5 in the main document
pdf("Boston-coef.pdf", height=8, width=11, pointsize=13)
plot(NA, xlim=c(1,p), ylim=ran, xlab="covariate index", ylab="coefficient", main="")
for(j in 1:3){
  points(cbind((1:p)+bb[j], Est[,j]), ylim=ran, col=col[j], pch=20)
  for(k in 1:p){
    lines(rep(k+bb[j],2), CI[[j]][,k], col=col[j])
  }
}
abline(h=0, lwd=2)
legend("bottomright", legend=c("N", "EH", "MT"), col=col, lty=1, pch=20)
dev.off()

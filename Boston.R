# data
library(mlbench)
data(BostonHousing2)
data <- BostonHousing2


Lab <- c("lon","lat","crim","zn","indus","chas","nox","rm","age","dis","rad","tax","ptratio","b","lstat")
Y <- data$cmedv
X <- subset(data, select=Lab)
n <- length(Y)
pp <- dim(X)[2]

X <- matrix(as.numeric(as.matrix(X)), n, pp)
X <- apply(X, 2, scale)
dimnames(X)[[2]] <- Lab
EX <- cbind(X, X[,-6]^2)     # design matrix


# simple linear regression 
fit.lm <- lm(Y~EX)
hsig <- summary(fit.lm)$sigma


# method
set.seed(1)
source("EHE-function.R")

MC <- 5000
Bn <- 2000

fit <- EHE.HS(Y, EX, mc=MC, burn=Bn)

apply(fit$Beta, 2, mean)
apply(fit$Beta, 2, quantile, prob=c(0.025, 0.975))



#  Figure 3 in the main document
pdf("Boston.pdf", height=8, width=12, pointsize=15)
par(mfcol=c(1,2))
plot(fit.lm$residuals/hsig, ylab="standardized residual", xlab="sample index")
abline(h=c(-1,1)*qnorm(0.975))
abline(h=c(-1,1)*qnorm(0.995), lty=2)
legend("topleft", legend=c("95% interval", "99% interval"), lty=1:2)
plot(fit$Om, type="l",ylab="mixing proportion", xlab="replication")
dev.off()


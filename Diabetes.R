# data
library(lars)
data(diabetes)

Y <- diabetes$y
X <- diabetes$x2

# simple linear regression
fit.lm <- lm(Y~X)
hsig <- summary(fit.lm)$sigma


# methods
set.seed(1)
source("EHE-function.R")

MC <- 5000
Bn <- 2000

fit <- EHE.HS(Y, X, mc=MC, burn=Bn)

apply(fit$Beta, 2, mean)
apply(fit$Beta, 2, quantile, prob=c(0.025, 0.975))



#  Figure 4 in the main document
pdf("Diabetes.pdf", height=8, width=12, pointsize=15)
par(mfcol=c(1,2))
plot(fit.lm$residuals/hsig, ylab="standardized residual", xlab="sample index", ylim=c(-5,5))
abline(h=c(-1,1)*qnorm(0.975))
abline(h=c(-1,1)*qnorm(0.995), lty=2)
legend("topleft", legend=c("95% interval", "99% interval"), lty=1:2)
plot(fit$Om, type="l",ylab="mixing proportion", xlab="replication", ylim=c(0, 0.2))
dev.off()



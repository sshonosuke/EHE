rm(list=ls())

mse <- c()
mse.sigma <- c()
cp <- c()
al <- c()
fac <- c()

load("sim-a10.RData")
mse <- cbind(mse, apply(MSE, c(2,3), mean))
mse.sigma <- cbind(mse.sigma, apply(MSE.sigma, c(2,3), mean))
cp <- cbind(cp, apply(CP, c(2,3), mean))
al <- cbind(al, apply(AL, c(2,3), mean))
IF[is.na(IF)] <- 50
fac <- cbind(fac, apply(IF, c(2,3), mean))

load("sim-a20.RData")
mse <- cbind(mse, apply(MSE, c(2,3), mean))
mse.sigma <- cbind(mse.sigma, apply(MSE.sigma, c(2,3), mean))
cp <- cbind(cp, apply(CP, c(2,3), mean))
al <- cbind(al, apply(AL, c(2,3), mean))
IF[is.na(IF)] <- 50
fac <- cbind(fac, apply(IF, c(2,3), mean))

res <- 10*rbind(t(sqrt(mse)), t(sqrt(mse.sigma)), t(cp)*10, t(al), t(fac)/10)
write.csv(res, file="res.csv")




mse <- c()
mse.sigma <- c()
cp <- c()
al <- c()
fac <- c()

load("sim-a10-out.RData")
mse <- cbind(mse, apply(MSE, c(2,3), mean))
mse.sigma <- cbind(mse.sigma, apply(MSE.sigma, c(2,3), mean))
cp <- cbind(cp, apply(CP, c(2,3), mean))
al <- cbind(al, apply(AL, c(2,3), mean))
IF[is.na(IF)] <- 50
fac <- cbind(fac, apply(IF, c(2,3), mean))

load("sim-a20-out.RData")
mse <- cbind(mse, apply(MSE, c(2,3), mean))
mse.sigma <- cbind(mse.sigma, apply(MSE.sigma, c(2,3), mean))
cp <- cbind(cp, apply(CP, c(2,3), mean))
al <- cbind(al, apply(AL, c(2,3), mean))
IF[is.na(IF)] <- 50
fac <- cbind(fac, apply(IF, c(2,3), mean))

res <- 10*rbind(t(sqrt(mse)), t(sqrt(mse.sigma)), t(cp)*10, t(al), t(fac)/10)
write.csv(res, file="res-out.csv")




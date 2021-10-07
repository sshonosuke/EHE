mse <- c()
mse2 <- c()
cp <- c()
al <- c()

load("sim-RI(a=5,om=0.05).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))

load("sim-RI(a=5,om=0.1).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))

load("sim-RI(a=10,om=0.05).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))

load("sim-RI(a=10,om=0.1).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))

load("sim-RI(a=15,om=0.05).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))

load("sim-RI(a=15,om=0.1).RData")
mse <- cbind(mse, apply(MSE, 2, mean))
mse2 <- cbind(mse2, apply(MSE.RE, 2, mean))
cp <- cbind(cp, apply(CP, 2, mean))
al <- cbind(al, apply(AL, 2, mean))


res <- 100*rbind(t(sqrt(mse)), t(cp), t(al), t(sqrt(mse2)))


write.csv(res, file=paste0("result-RI.csv"))

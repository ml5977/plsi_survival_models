#random variable from uniform
u <- runif(n)
#extreme value function
ee <- log(-log(1-u))
#standard logistic function
#ee <- -log(1/u-1)

#generating survival time
surv.t0 <- exp(-(beta0[1]*x1 + beta0[2]*x2 + beta0[3]*x3 + beta0[4]*x4 + beta0[5]*x5 + alpha0[1]*z1 + alpha0[2]*z2) + ee)
#surv.t0 <- -(beta0[1]*x1 + beta0[2]*x2 + beta0[3]*x3 + beta0[4]*x4 + beta0[5]*x5 + alpha0[1]*z1 + alpha0[2]*z2) + ee
cen.t0 <- rexp(n, rate=1/12)
ind0 <- ifelse(surv.t0<=cen.t0, 1, 0)
obs.t0 <- ifelse(surv.t0<=cen.t0, surv.t0, cen.t0)
#summary(surv.t0)
#summary(obs.t0)
#event.rate[i] <- sum(ind0)/n

sim.data <- data.frame(id=c(1:n), time=obs.t0, status=ind0, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, z1=z1, z2=z2)
write.csv(sim.data, file="/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_transformation/gen_data_linear.csv")


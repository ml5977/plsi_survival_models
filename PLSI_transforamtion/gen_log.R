#Generating survival time
n <- 1000
alpha0 <- c(0.3, -0.3)
beta0 <- c(1, -1, 1, -1, 1); beta0 <- beta0/sqrt(sum(beta0^2))

bound <- 4
tr <- 4
x1 <- runif(n, -bound, bound)
x2 <- runif(n, -bound, bound)
x3 <- runif(n, -bound, bound)
x4 <- rnorm(n, 0, sqrt(2)); x4 <- ifelse(x4 > tr, tr, ifelse(x4 < -tr, -tr, x4))
x5 <- rnorm(n, 0, sqrt(2)); x5 <- ifelse(x5 > tr, tr, ifelse(x5 < -tr, -tr, x5))
z1 <- rnorm(n)
z2 <- rbinom(n, size=1, p=0.5)

#random variable from uniform
u <- runif(n)
#extreme value function
ee <- log(-log(1-u))
#standard logistic function
#ee <- -log(1/u-1)

#generating survival time
surv.t0 <- exp(-(log(1+(beta0[1]*x1 + beta0[2]*x2 + beta0[3]*x3 + beta0[4]*x4 + beta0[5]*x5)^2) + alpha0[1]*z1 + alpha0[2]*z2) + ee)
#surv.t0 <- -(log(1+(beta0[1]*x1 + beta0[2]*x2 + beta0[3]*x3 + beta0[4]*x4 + beta0[5]*x5)^2) + alpha0[1]*z1 + alpha0[2]*z2) + ee
cen.t0 <- rexp(n, rate=1/5)
ind0 <- ifelse(surv.t0<=cen.t0, 1, 0)
obs.t0 <- ifelse(surv.t0<=cen.t0, surv.t0, cen.t0)

sim.data <- data.frame(id=c(1:n), time=obs.t0, status=ind0, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, z1=z1, z2=z2)
write.csv(sim.data, file="/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_transformation/gen_data_log.csv")

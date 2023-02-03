###Required R packages
require(survival)
require(splines2)

###simulation data load
data <- read.csv("/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_Cox/gen_data_log.csv")
#Please look at README file for details of simulation data

###NOTE
#For time-independent covariates, we can simply use Surv(time, status) instead of Surv(start, last, status) in coxph function below.

###Time-dependent Cox regression
tdcm.fit <- coxph(Surv(start, last, status) ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z1 + z2 + z3, data=data)

###PLSI Cox regression
#initial parameters for \beta
ini.b <- tdcm.fit$coefficients[c(1:8)]; ini.b <- ini.b/sqrt(sum(ini.b^2))
tol <- 0.0001
max.iter <- 100

#iterative estimating procedure
for(ii in 1:max.iter){
  #grid of index
  data$index <- as.numeric(as.matrix(data[,c("x1","x2","x3","x4","x5","x6","x7","x8")])%*%ini.b)
  index <- c(data$index, 0)
  
  #B-spline for nonparametric link function
  knot <- 3
  sp <- ibs(index, degree=2, df=knot)
  sp <- sp - matrix(rep(sp[length(index),], length(index)), ncol=knot, byrow=TRUE)
  dsp <-  bSpline(index,degree=2,df=knot)
  
  #Step 1. estimating gamma & alpha
  data$b1 <- sp[-length(index),1]
  data$b2 <- sp[-length(index),2]
  data$b3 <- sp[-length(index),3]
  fit.gamma <- coxph(Surv(start, last, status) ~ b1 + b2 + b3 + z1 + z2 + z3, data=data)
  gamma <- fit.gamma$coefficients[1:3]
  alpha <- fit.gamma$coefficients[4:6]
  
  #Step 2. estimating beta
  data$g.curve <- dsp[-length(index),]%*%gamma
  data$gg.curve <- sp[-length(index),]%*%gamma
  data$a.curve <- as.matrix(data[,c("z1","z2","z3")])%*%alpha
  fit.beta <- coxph(Surv(start, last, status) ~ I(g.curve*x1) + I(g.curve*x2) + I(g.curve*x3) + I(g.curve*x4) + 
                      I(g.curve*x5) + I(g.curve*x6) + I(g.curve*x7) + I(g.curve*x8) + offset(gg.curve-g.curve*index+a.curve), data=data)
  beta <- fit.beta$coefficients
  
  #constraint: b1>0
  if(beta[1]<0){
    beta <- -beta
  }
  
  #iteration print out
  print(paste("PLSI Cox model - iter:",ii,"out of",max.iter))
  #print(beta/sqrt(sum(beta^2)))
  
  if(ii==max.iter){
    print(paste("Not converged within",max.iter,"iterations"))
    break
  }
  
  if(sum(abs(beta/sqrt(sum(beta^2))-ini.b)) < tol){
    print("converged!")
    break
  }
  
  ini.b <- beta/sqrt(sum(beta^2))
}

#estimations for beta
fit.beta$coefficients

#estimations for alpha and gamma
fit.gamma$coefficients

#estimated link function
plot(data$index[order(data$index)], data$gg.curve[order(data$index)], type="l") #estimated 
lines(data$index[order(data$index)], log(1+data$index[order(data$index)]^2), col="red") #true

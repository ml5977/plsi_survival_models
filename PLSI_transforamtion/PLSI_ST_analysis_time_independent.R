###packages
require(splines2)

###functions
source("/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_transformation/PLSI_ST_functions_time_independent.R")

###simulation data load
data <- read.csv("/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_transformation/gen_data_log.csv")
#Please look at README file for details of simulation data

###setting
n <- nrow(data)
rho <- 1; r<- 1 # spcify the model from ST model (rho=1:PH, rho=0:PO)
df <- 4 # number of knots

###STEP1: standard ST model for initial value
#EM algorithm
temp.mat <- matrix(data$time, nrow=n, ncol=n)
indY <- (temp.mat>=t(temp.mat))*1
m <- sum(data$status);m

#setting
nbeta <- 7; oldbeta <- rep(0, nbeta); oldq <- rep(0, n)+1/m 
epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;

#iteration between E and M steps
while(error >= epsilon && iter <= maxiter){
  #Estep
  Exi <- Estep(oldbeta=oldbeta, oldq=oldq, Y=data$time, Delta=data$status, 
               Z=data[,c("x1","x2","x3","x4","x5","z1","z2")], indY=indY, n=n, rho=rho, r=r)
  
  #Mstep
  Moutput <- Mstep(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$time, 
                   Delta=data$status, Z=data[,c("x1","x2","x3","x4","x5","z1","z2")], indY=indY, n=n)
  newbeta <- Moutput$newbeta
  newq <- Moutput$newq
  error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
  iter <- iter+1
  oldbeta <- newbeta
  oldq <- newq
}

#output
newbeta #coefficients
data$q.st <- newq

###STEP2: PLSI-ST model 
ini.parm <- newbeta[c(1:5)]
ini.parm <- ini.parm/sqrt(sum(ini.parm^2))
temp.alpha <- newbeta[c(6:7)]
tol <- 0.0001
max.iter <- 100
data <- data

j <- 1
for(j in 1:max.iter){
  #single index construction
  data$index <- as.numeric(as.matrix(data[,c("x1","x2","x3","x4","x5")])%*%ini.parm)
  index <- c(data$index, 0)
  #summary(index)
  
  #B-spline with 3 knots
  sp <- ibs(index, degree=2, df=df)
  sp <- sp - matrix(rep(sp[length(index),], length(index)), ncol=df, byrow=TRUE)
  dsp <- bSpline(index, degree=2, df=df)
  
  #estimating gamma & alpha
  data$b1 <- sp[-length(index),1]
  data$b2 <- sp[-length(index),2]
  data$b3 <- sp[-length(index),3]
  data$b4 <- sp[-length(index),4]
  
  #Transformation model for gamma and alpha
  #EM algorithm
  temp.mat <- matrix(data$time, nrow=n, ncol=n)
  indY <- (temp.mat>=t(temp.mat))*1
  m <- sum(data$status);m
  
  #setting
  nbeta <- df+length(temp.alpha); oldbeta <- c(rep(0, df),temp.alpha) ; oldq <- oldq
  epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;
  
  while(error >= epsilon && iter <= maxiter){
    #Estep
    Exi <- Estep(oldbeta=oldbeta, oldq=oldq, Y=data$time, Delta=data$status, 
                 Z=data[,c("b1","b2","b3","b4","z1","z2")], indY=indY, n=n, rho=rho, r=r)
    
    #Mstep
    Moutput <- Mstep(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$time, 
                     Delta=data$status, Z=data[,c("b1","b2","b3","b4","z1","z2")], indY=indY, n=n)
    newbeta <- Moutput$newbeta
    newq <- Moutput$newq
    error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
    iter <- iter+1
    oldbeta <- newbeta
    oldq <- newq
  }
  if(iter==maxiter){
    iter.output.plsi.1st[i] <- 1
  }
  
  temp.q <- newq
  temp.gamma <- newbeta[c(1:4)]; temp.alpha <- newbeta[c(5:6)]
  
  #estimating beta
  data$g.curve <- dsp[-length(index),]%*%temp.gamma
  data$gg.curve <- sp[-length(index),]%*%temp.gamma
  data$a.curve <- as.matrix(data[,c("z1","z2")])%*%temp.alpha
  data$offset <- data$gg.curve-data$g.curve*data$index+data$a.curve
  data$x1.new <- data$g.curve*data$x1
  data$x2.new <- data$g.curve*data$x2
  data$x3.new <- data$g.curve*data$x3
  data$x4.new <- data$g.curve*data$x4
  data$x5.new <- data$g.curve*data$x5
  
  #Transformation model for beta
  #EM algorithm
  temp.mat <- matrix(data$time, nrow=n, ncol=n)
  indY <- (temp.mat>=t(temp.mat))*1
  m <- sum(data$status);m
  
  nbeta <- length(ini.parm); oldbeta <- ini.parm; oldq <- oldq
  epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;
  
  while(error >= epsilon && iter <= maxiter){
    #Estep
    Exi <- Estep.offset(oldbeta=oldbeta, oldq=oldq, Y=data$time, Delta=data$status, 
                        Z=data[,c("x1.new","x2.new","x3.new","x4.new","x5.new")], indY=indY, n=n, rho=rho, r=r, offset=data[,c("offset")])
    
    #Mstep
    Moutput <- Mstep.offset(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$time, 
                            Delta=data$status, Z=data[,c("x1.new","x2.new","x3.new","x4.new","x5.new")], indY=indY, n=n, offset=data[,c("offset")])
    newbeta <- Moutput$newbeta
    #newq <- Moutput$newq
    error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
    iter <- iter+1
    oldbeta <- newbeta
    #oldq <- newq
  }
  if(iter==maxiter){
    iter.output.plsi.2nd[i] <- 1
  }
  temp.beta <- newbeta[c(1:5)]
  
  if(temp.beta[1]<0){
    temp.beta <- -temp.beta
  }
  
  #iteration print out
  print(paste("PLSI-trans method - iter:",j))
  #print(temp.beta/sqrt(sum(temp.beta^2)))
  if(j==max.iter){
    print(paste("Not converged within",max.iter,"iterations"))
    break
  }
  
  if(sum(abs(temp.beta/sqrt(sum(temp.beta^2))-ini.parm)) < tol){
    break
  }
  ini.parm <- temp.beta/sqrt(sum(temp.beta^2))
  
}

#estimation results
temp.beta #coefficeints for X for nonlinear effects
temp.alpha #coefficients for Z for linear effects
temp.gamma #parameter estimations for B-spline basis

#estimated link function
grid.xbeta <- seq(-2,2, by=0.01)
grid.sp <- ibs(grid.xbeta, degree=2, df=df)
grid.sp <- grid.sp - matrix(rep(grid.sp[which(grid.xbeta==0),], length(grid.xbeta)), ncol=4, byrow=TRUE)
g.function.est <- grid.sp%*%temp.gamma
#plot
#estimated link function plot
plot(grid.xbeta, g.function.est, type="l", lwd=1, col="black", ylim=c(0,2), xlab="Xbeta", ylab="Link function") #estimated
lines(grid.xbeta, log(1+grid.xbeta^2), col="red") #true



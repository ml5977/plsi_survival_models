###packages
library(MASS)
library(survival)
library(splines2)

###functions
source("/Users/leem26/Desktop/MLEE/NIEHS/Research/PhD_work/Github/PLSI_transformation/PLSI_ST_functions_time_dependent.R")

#Load data
data(VA)
head(VA)

#dummy variable
data <- VA
table(data$cell)
table(data$treat)
table(data$prior)
data$cell1 <- ifelse(data$cell==1,1,0)
data$cell2 <- ifelse(data$cell==2,1,0)
data$cell3 <- ifelse(data$cell==3,1,0)
data$treat <- as.integer(data$treat)-1
data$prior <- ifelse(data$prior==10,1,0)
data$prior <- data$prior*10
#standardized
data$age <- (data$age-mean(data$age))/sd(data$age)
data$Karn <- (data$Karn-mean(data$Karn))/sd(data$Karn)
data$diag.time <- (data$diag.time-mean(data$diag.time))/sd(data$diag.time)
head(data)
summary(data)


#cox regression
cox.fit <- coxph(Surv(stime, status) ~ age + Karn + diag.time + 
                   treat + cell1 + cell2 + cell3 + prior, data=data)
summary(cox.fit)
round(cox.fit$coefficients,3)
round(summary(cox.fit)$coefficients,3)

#previous studies
#Huang and Liu (2006): all variables were included in the link function
#Sun et al. (2008): Treatment, prior, and cell were modelled for linear, 
#while age, Karn, and diag.time were modelled by nonlinear

#setting for B-spline approximation
df <- 4
p <- 5 #linear
q <- 3 #nonlinear
#candidate rho
#cand.rho <- seq(0,2, by=0.1)
cand.rho <- c(0, 0.1, 0.5, 1.0)
MM <- length(cand.rho)

#storage
loglike.value <- rep(NA, MM)
standard.st.coef <- matrix(NA, nrow=MM, ncol=(p+q))
standard.st.q <- matrix(NA, nrow=MM, ncol=nrow(data))
standard.loglike.value <- rep(NA, MM)
standard.st.var <- list()
plsi.st.beta <- matrix(NA, nrow=MM, ncol=q)
plsi.st.alpha <- matrix(NA, nrow=MM, ncol=p)
plsi.st.gamma <- matrix(NA, nrow=MM, ncol=df)
plsi.st.q <- matrix(NA, nrow=MM, ncol=nrow(data))

ll <- 1
for(ll in 1:MM){
  ####initial value from standard ST
  #setting
  n <- nrow(data)
  rho <- cand.rho[ll]; r<- 1 # spcify the model from ST model (rho=1:PH, rho=0:PO)
  m <- sum(data$status);m
  temp.q <- rep(0, n)+1/m
  #EM algorithm
  temp.mat <- matrix(data$stime, nrow=n, ncol=n)
  indY <- (temp.mat>=t(temp.mat))*1
  m <- sum(data$status);m
  
  #setting
  nbeta <- p+q; oldbeta <- rep(0, nbeta) ; oldq <- rep(0, n)+1/m
  epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;
  
  while(error >= epsilon && iter <= maxiter){
    #Estep
    Exi <- Estep(oldbeta=oldbeta, oldq=oldq, Y=data$stime, Delta=data$status, 
                 Z=data[,c("age", "Karn", "diag.time", "treat","cell1", "cell2", "cell3", "prior")], indY=indY, n=n, rho=rho, r=r)
    
    #Mstep
    Moutput <- Mstep(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$stime, 
                     Delta=data$status, Z=data[,c("age", "Karn", "diag.time", "treat","cell1", "cell2", "cell3", "prior")], indY=indY, n=n)
    newbeta <- Moutput$newbeta
    newq <- Moutput$newq
    error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
    iter <- iter+1
    oldbeta <- newbeta
    oldq <- newq
    #print(newbeta)
    #print(iter)
  }
  
  standard.st.coef[ll,] <- newbeta
  standard.st.q[ll,] <- newq
  
  standard.loglike.value[ll] <- loglik.fun(final.gamma=newbeta[c(1:3)], final.alpha=newbeta[c(4:8)],
                                           final.q=newq, Z=data[,c("age","Karn","diag.time","treat","cell1", "cell2", "cell3", "prior")], indY=indY, Delta=data$status,
                                           n=n, rho=rho, r=r)
  standard.st.var[[ll]] <- Covest(newbeta=newbeta, newq=newq, Y=data$stime, Delta=data$status, 
                                  Z=data[,c("age", "Karn", "diag.time", "treat","cell1", "cell2", "cell3", "prior")], indY=indY, n=n, rho=rho, r=r)
  #newbeta
  #cox.fit$coefficients
  #head(newq)
  
  
  #criterion for PLSI-ST convergence
  tol <- 0.0001
  max.iter <- 1000
  
  #initial and parameter setting from standard method
  ini.parm <- newbeta[c(2,1,3)]
  ini.parm <- ini.parm/sqrt(sum(ini.parm^2)); ini.parm
  temp.alpha <- newbeta[c((q+1):(p+q))];temp.alpha
  temp.gamma <- rep(0, df)
  
  #Algorithm start!
  j <- 1
  for(j in 1:max.iter){
    #single index construction
    data$index <- as.numeric(as.matrix(data[,c("Karn", "age", "diag.time")])%*%ini.parm)
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
    temp.mat <- matrix(data$stime, nrow=n, ncol=n)
    indY <- (temp.mat>=t(temp.mat))*1
    m <- sum(data$status);m
    
    #setting
    nbeta <- df+length(temp.alpha); oldbeta <- c(temp.gamma,temp.alpha) ; oldq <- temp.q
    epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;
    
    while(error >= epsilon && iter <= maxiter){
      #Estep
      Exi <- Estep(oldbeta=oldbeta, oldq=oldq, Y=data$stime, Delta=data$status, 
                   Z=data[,c("b1","b2","b3","b4","treat","cell1", "cell2", "cell3", "prior")], indY=indY, n=n, rho=rho, r=r)
      
      #Mstep
      Moutput <- Mstep(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$stime, 
                       Delta=data$status, Z=data[,c("b1","b2","b3","b4","treat","cell1", "cell2", "cell3", "prior")], indY=indY, n=n)
      newbeta <- Moutput$newbeta
      newq <- Moutput$newq
      error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
      iter <- iter+1
      oldbeta <- newbeta
      oldq <- newq
    }
    
    iter
    temp.q <- newq
    temp.gamma <- newbeta[c(1:df)]; temp.alpha <- newbeta[c((df+1):(df+p))]
    
    #estimating beta
    data$g.curve <- dsp[-length(index),]%*%temp.gamma
    data$gg.curve <- sp[-length(index),]%*%temp.gamma
    data$a.curve <- as.matrix(data[,c("treat","cell1", "cell2", "cell3", "prior")])%*%temp.alpha
    data$offset <- data$gg.curve-data$g.curve*data$index+data$a.curve
    data$age.n <- data$g.curve*data$age
    data$Karn.n <- data$g.curve*data$Karn
    data$diag.time.n <- data$g.curve*data$diag.time
    
    
    #Transformation model for beta
    #EM algorithm
    temp.mat <- matrix(data$stime, nrow=n, ncol=n)
    indY <- (temp.mat>=t(temp.mat))*1
    m <- sum(data$status);m
    
    nbeta <- length(ini.parm); oldbeta <- ini.parm; oldq <- temp.q
    epsilon <- 0.0001; maxiter <- 500; error <- 1; iter <- 0;
    
    while(error >= epsilon && iter <= maxiter){
      #Estep
      Exi <- Estep.offset(oldbeta=oldbeta, oldq=oldq, Y=data$stime, Delta=data$status, 
                          Z=data[,c("Karn.n", "age.n", "diag.time.n")], indY=indY, n=n, rho=rho, r=r, offset=data[,c("offset")])
      
      #Mstep
      Moutput <- Mstep.offset(oldbeta=oldbeta, oldq=oldq, Exi=Exi,Y=data$stime, 
                              Delta=data$status, Z=data[,c("Karn.n", "age.n", "diag.time.n")], indY=indY, n=n, offset=data[,c("offset")])
      newbeta <- Moutput$newbeta
      #newq <- Moutput$newq
      error <- sum(abs(newbeta-oldbeta)+sum(abs(newq-oldq)))
      iter <- iter+1
      oldbeta <- newbeta
      #oldq <- newq
    }
    
    iter
    temp.beta <- newbeta[c(1:q)]
    
    if(temp.beta[1]>0){
      temp.beta <- -temp.beta
    }
    
    #iteration print out
    print(paste("PLSI-trans method - iter:",j))
    print(temp.beta/sqrt(sum(temp.beta^2)))
    print(max(abs(temp.beta/sqrt(sum(temp.beta^2))-ini.parm)))
    
    if(max(abs(temp.beta/sqrt(sum(temp.beta^2))-ini.parm)) < tol){
      break
    }
    ini.parm <- temp.beta/sqrt(sum(temp.beta^2))
    
  }
  
  #head(temp.q)
  #round(temp.beta/sqrt(sum(temp.beta^2)),3)
  #round(temp.gamma,3)
  #round(temp.alpha,3)
  #round(exp(temp.alpha),3) #treatment effect
  
  plsi.st.alpha[ll,] <- temp.alpha
  plsi.st.gamma[ll,] <- temp.gamma
  plsi.st.beta[ll,] <- temp.beta/sqrt(sum(temp.beta^2))
  plsi.st.q[ll,] <- temp.q
  
  #loglikelihood value
  
  loglike.value[ll] <- loglik.fun(final.gamma=temp.gamma, final.alpha=temp.alpha,
                                  final.q=temp.q, Z=data[,c("b1","b2","b3","b4","treat","cell1", "cell2", "cell3", "prior")], indY=indY, Delta=data$status,
                                  n=n, rho=rho, r=r)
  
  print(ll)
}


###
plot(cand.rho, standard.loglike.value, xlab=expression(rho), ylab="log-likelihood", type="l")
plot(cand.rho, loglike.value, xlab=expression(rho), ylab="log-likelihood", type="l")

#cox output
round(cox.fit$coefficients,3)
round(summary(cox.fit)$coefficients,3)

#transformation output for PH assumption (rho=1.0)
round(standard.st.coef[3,],3)

#plsi transformation output
round(plsi.st.beta[3,],3) #for "Karn", "age", "diag.time"
round(plsi.st.alpha[3,],3) #for "treat","cell1", "cell2", "cell3", "prior"

###PLOT
#optimal
ii <- 3
index <- seq(-3, 3, by=0.1)
index <- c(index, 0)
#B-spline with 3 knots
sp <- ibs(index, degree=2, df=df)
sp <- sp - matrix(rep(sp[length(index),], length(index)), ncol=df, byrow=TRUE)
dsp <- bSpline(index, degree=2, df=df)
plot.gg.curve <- sp[-length(index),]%*%plsi.st.gamma[ii,]
index <- index[-length(index)]
summary(index)
summary(plot.gg.curve)
plot(index, plot.gg.curve, 
     type="l", xlab="Single index", lwd=1.5,lty=1, ylab="Link function", 
     xlim=c(-3, 3), 
     ylim=c(-6, 10))
#title("PLSI-ST with PH")
abline(a=0, b=1, col="grey")



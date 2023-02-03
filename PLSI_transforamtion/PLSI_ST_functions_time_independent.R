simulation <- function(beta0, lambda0, n, rho, r, tau){
  tempZ1 <- rbinom(n,1,p=0.5)
  tempZ2 <- runif(n)
  Z <- cbind(tempZ1, tempZ2)
  ezbeta <- exp(Z%*%beta0)
  
  if(rho>0){
    tempT <- -log(runif(n))*rho+1
    trueT <- (tempT^(1/rho)-1)/(lambda0*ezbeta)
  }
  
  if(rho==0){
    newr <- 1/r
    tempT <- exp(-log(runif(n))/newr)-1
    trueT <- tempT/(1/newr)/(lambda0*ezbeta)
  }
  C <- runif(n)*5; C <- ifelse(C<=tau, C, tau)
  Y <- ifelse(trueT<=C,trueT, C); delta <- ifelse(trueT<=C,1,0)
  
  return(data.frame(Y=Y, delta=delta, Z1=tempZ1, Z2=tempZ2))
}

########
Gtransform <- function(temp, rho, r){
  if(rho>0){
    G <- ((1+temp)^rho-1)/rho
    dG <- 1+temp^(rho-1)
    ddG <- (rho-1)*(1+temp)^(rho-2)
    d3G <- (rho-1)*(rho-2)*(1+temp)^(rho-3)
  }
  
  if(rho==0){
    newr <- 1/r
    G <- newr*log(1+temp/newr)
    dG <- 1/(1+temp/newr)
    ddG <- -1/newr/((1+temp/newr)^2)
    d3G <- 2/(newr^2)/((1+temp/newr)^3)
  }
  return(data.frame(G=G, dG=dG, ddG=ddG, d3G=d3G))
}

###########
# Estep
# rho=0, (r1,r2) corresponds to gamma frailty (r1, r2):
#x^{r1-1}exp\{-x/r2\}
# rho>0 corresponds to 
Estep <- function(oldbeta, oldq, Y, Delta,Z, indY, n, rho, r){
  oldQ <- indY%*%oldq
  lambdabetaz <- oldQ*exp(as.matrix(Z)%*%oldbeta)
  goutput <- Gtransform(temp=lambdabetaz, rho=rho, r=r)
  Exi <- -Delta*goutput$ddG/goutput$dG+goutput$dG
  
  return(Exi)
  
}

#############
# Data Structure: Y1, D1, Z1, Y2, D2, Z2
# Mstep
#Z <- simdata[,c("Z1","Z2")]
#Delta <- simdata$delta
Mstep <- function(oldbeta, oldq, Exi,Y,Delta,Z, indY, n){
  nbeta <- length(oldbeta); EbetaZ <- exp(as.matrix(Z)%*%oldbeta); ww <- EbetaZ*Exi
  nom <- t(as.matrix(indY))%*%as.matrix((Z*matrix(ww, ncol=nbeta, nrow=length(ww))))
  denom <- t(indY)%*%ww
  g <- t(nom)/t(matrix(denom, ncol=nbeta, nrow=length(denom)))
  
  tempnom2 <- matrix(0, nrow=nbeta, ncol=nbeta)
  
  for(k in 1:nbeta){
    for(l in k:nbeta){
      tempnom2[k,l] <- sum(Delta*(t(as.matrix(indY))%*%(Z[,k]*Z[,l]*ww)/denom))
      tempnom2[l,k] <- tempnom2[k,l]
    }
  }
  
  score <- t(as.matrix(Z)) %*% as.matrix(Delta) - g %*% as.matrix(Delta)
  dscore <- -tempnom2+g%*%(matrix(Delta, ncol=nbeta, nrow=length(Delta))*t(g))
  
  
  newbeta <- oldbeta-solve(dscore)%*%score
  EnewbetaZ <- exp(as.matrix(Z)%*%newbeta)
  newq <- Delta/(t(indY)%*%as.matrix((EnewbetaZ*Exi)))
  
  return(list(newbeta=newbeta, newq=newq))
  
}

############
# covariance estimation using the Louis-formula
# rho=0, (r1,r2) corresponds to gamma frailty (r1, r2):
#x^{r1-1}exp\{-x/r2\}
# rho>0 corresponds to
#Y <- simdata$Y
#Z <- simdata[,c("Z1","Z2")]
#Delta <- simdata$delta

Covest <- function(newbeta, newq, Y, Delta, Z, indY, n, rho, r){
  nbeta <- length(newbeta)
  templambda <- newq; templambda[Delta==0] <- 1
  Lambda <- indY%*%newq; EZbeta <- exp(as.matrix(Z)%*%newbeta); LambdabetaZ <- Lambda*EZbeta
  
  goutput <- Gtransform(temp=LambdabetaZ, rho=rho, r=r)
  
  Exi1 <- -Delta*goutput$ddG/goutput$dG+goutput$dG
  Exi2 <- Delta*(goutput$d3G-3*goutput$ddG*goutput$dG+goutput$dG^3)/goutput$dG-(1-Delta)*(goutput$ddG-goutput$dG^2)
  
  S1 <- cbind(Z*matrix(Delta, ncol=nbeta, nrow=length(Delta)),
              diag(as.vector(Delta/templambda)))
  S2 <- cbind(Z*matrix(LambdabetaZ, ncol=nbeta, nrow=length(LambdabetaZ)),
              indY*matrix(EZbeta, ncol=n, nrow=length(EZbeta)))
  
  
  l1 <- S1-matrix(Exi1, ncol=(n+nbeta), nrow=length(Exi1))*S2
  l11 <- t(S1)%*%as.matrix(S1)-t(S1)%*%as.matrix((matrix(Exi1, ncol=(n+nbeta), nrow=length(Exi1))*S2)) -
    t(S2)%*%as.matrix(matrix(Exi1, ncol=(n+nbeta), nrow=length(Exi1))*S1) +
    t(S2)%*%as.matrix(matrix(Exi2, ncol=(n+nbeta), nrow=length(Exi2))*S2)
  
  l2 <- matrix(0, ncol=(n+nbeta), nrow=(n+nbeta))
  Icoef <- c(1:nbeta); I1 <- nbeta+c(1:n)
  
  l2[Icoef, Icoef] <- -t(Z)%*%as.matrix(matrix(Exi1*LambdabetaZ, ncol=nbeta, nrow=length(Exi1*LambdabetaZ))*Z)
  l2[Icoef, I1] <- -t(Z)%*%as.matrix(matrix(Exi1*EZbeta, ncol=n, nrow=length(Exi1*EZbeta))*indY)
  l2[I1, Icoef] <- t(l2[Icoef, I1])
  l2[I1, I1] <- -diag(as.vector(Delta/templambda^2))
  
  information <- -l2-(l11-t(l1)%*%as.matrix(l1))
  index <- c(1+matrix(0, nrow=nbeta), Delta)
  information <- information[index==1, index==1]
  cov <- solve(information)
  return(cov)
}


###########
# Estep
# rho=0, (r1,r2) corresponds to gamma frailty (r1, r2):
#x^{r1-1}exp\{-x/r2\}
# rho>0 corresponds to 
#Z=data[,c("x1.new","x2.new","x3.new","x4.new","x5.new")]
#offset=data[,c("offset")]
Estep.offset <- function(oldbeta, oldq, Y, Delta,Z, indY, n, rho, r, offset){
  oldQ <- indY%*%oldq
  lambdabetaz <- oldQ*exp(as.matrix(cbind(Z,offset))%*%c(oldbeta,1))
  goutput <- Gtransform(temp=lambdabetaz, rho=rho, r=r)
  Exi <- -Delta*goutput$ddG/goutput$dG+goutput$dG
  
  return(Exi)
  
}


#############
# Data Structure: Y1, D1, Z1, Y2, D2, Z2
# Mstep
#Z <- simdata[,c("Z1","Z2")]
#Delta <- simdata$delta
Mstep.offset <- function(oldbeta, oldq, Exi,Y,Delta,Z, indY, n, offset){
  nbeta <- length(oldbeta); EbetaZ <- exp(as.matrix(cbind(Z,offset))%*%c(oldbeta,1)); ww <- EbetaZ*Exi
  nom <- t(as.matrix(indY))%*%as.matrix((Z*matrix(ww, ncol=nbeta, nrow=length(ww))))
  denom <- t(indY)%*%ww
  g <- t(nom)/t(matrix(denom, ncol=nbeta, nrow=length(denom)))
  
  tempnom2 <- matrix(0, nrow=nbeta, ncol=nbeta)
  
  for(k in 1:nbeta){
    for(l in k:nbeta){
      tempnom2[k,l] <- sum(Delta*(t(as.matrix(indY))%*%(Z[,k]*Z[,l]*ww)/denom))
      tempnom2[l,k] <- tempnom2[k,l]
    }
  }
  
  score <- t(as.matrix(Z)) %*% as.matrix(Delta) - g %*% as.matrix(Delta)
  dscore <- -tempnom2+g%*%(matrix(Delta, ncol=nbeta, nrow=length(Delta))*t(g))
  
  
  newbeta <- oldbeta-solve(dscore)%*%score
  EnewbetaZ <- exp(as.matrix(Z)%*%newbeta)
  newq <- Delta/(t(indY)%*%as.matrix((EnewbetaZ*Exi)))
  
  return(list(newbeta=newbeta, newq=newq))
  
}

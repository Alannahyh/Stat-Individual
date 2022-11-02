#Alannah Hounat S2434943

library(MASS)
x<-mcycle$times
y<-mcycle$accel


pspline<-function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  
  
  lsp=seq(logsp[1],logsp[2],length=ngrid)
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
  #finding the Qr decomposition of X
  qrx<-qr(X)
  #finding the crossproduct of D for the b.hat formula 
  crossD<-crossprod(D)
  #formula for b.hat
  formula.bhat<- t(solve(qr.R(qrx)))%*%crossD%*%solve(qr.R(qrx))
  #eigen decomposition, used for finding lambda and U
  #in finding the 
  e.decomp<-eigen(formula.bhat)
  
  #lam<-e.decomp$values
  #vector containing eigenvals in the diagonal
  lambda<-diag(e.decomp$values)
  #matrix of eigen vectors
  U<-e.decomp$vectors
  
  edf <- gcv <- sig <- lsp*0 ## vectors for edf and gcv and sigma squared
  for (i in 1:length(lsp)){
    ## loop over log smoothing parameters
    fr <- get.gcv(y,X,exp(lsp[i]),D,bord,pord,qrx,lambda,U) ## fit
    gcv[i] <- fr$gcv
    edf[i] <- fr$edf
    sig[i] <- fr$sig
    
  }
  #find the optimum value for gcv
  i.opt <- max(which(gcv==min(gcv))) 
  get.gcv(y,X,exp(lsp[i.opt]),D,bord,pord,qrx,lambda,U)
  val<-get.gcv(y,X,exp(lsp[i.opt]),D,bord,pord,qrx,lambda,U)
  class(val)<-"pspline"
  return(val)
}


get.gcv <- function (y,X,sp,D,bord,pord,qrx,lambda,U){
  #number of rows and colums of matrix X
  p <- ncol(X);n <- nrow(X)
  
  #A is the matrix (I+λΛ)^-1
  A<-solve(diag(p)+sp*lambda)
  #effective degree of freedom of the model
  trA <- sum(diag(A)) 

  # computing coeff estimates,b.hat.
  b.hat<-backsolve(qr.R(qrx),U%*%solve(diag(p)+lambda*sp)%*%t(U)%*% qr.qty(qrx,y)[1:p])
  #mu.hat(fitted values)
  mu.hat <- X %*% b.hat 
  #sigma squared 
  sig<-sum((y-mu.hat)^2)/(n-trA)
  
  res<-sqrt(sig)
  V<-solve(t(X)%*%X+sp*t(D)%*%D)*sig
  #standard error 
  stan_error<-rowSums(X*(X%*%V))^0.5
  #generalized cross validation
  #gcv is used to help us find the most optimal value for λ 
  gcv <- sig/(n-trA)
  #empty vector for holding the values of y[i]-mean(y)
  sm<-c()
  #takes the elements of the y vector and subtracts mean(y)
  #this is needed for the formula of r squared
  for(elem in 1:n){
    values<-(y[elem]-mean(y))^2
    sm<-c(sm,values)
    
  }
  #summing up the values of y[i]-mean(y)
  #which needed to be unlisted prior to summation
  sum.a<-sum(unlist(sm))
  ##??
  r_sq<-1-((n-1)*sig)/sum.a
  
  list(b.hat=b.hat,mu.hat=mu.hat,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res,bord=bord,pord=pord,V=V,X=X,D=D,sp=sp)
  vals<-list(b.hat=b.hat,mu.hat=mu.hat,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res, bord=bord,pord=pord,V=V,X=X,D=D,sp=sp)
  #class(list(b.hat=b.hat,mu.hat=mu.hat,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res, bord=bord,pord=pord,V=V,X=X,D=D,sp=sp))<-"pspline"
  #return(vals)
}

#INPUT
#OUTPUT
#PURPOSE:takes in the values produced by pspline and returns 
#
print.pspline<- function(m){
  cat('Order of' ,m$bord, 'p-spline with order',m$pord, 'penalty',"\n")
  cat('Effective degrees of freedom:',m$edf,' ')
  cat('Coefficients:',m$k,"\n")
  cat('residual std dev:',m$res," ")
  cat('r-squared:',m$r_sq," ")
  cat('GCV:',m$gcv,"\n")
  lists<-list(gcv=m$gcv,edf=m$edf,r_sq=m$r_sq)
  invisible(lists)
}
#
new_x_val<-seq(min(x),max(x),by=1)

predict.pspline<- function(m,x,se=TRUE){
  

  bord<-3
  k<-20
  dk <- diff(range(new_x_val))/(k-bord) ## knot spacing
  knots <- seq(min(new_x_val)-dk*bord,by=dk,length=k+bord+1)
  Xp <- splines::splineDesign(knots,new_x_val,ord=bord+1,outer.ok=TRUE)
  
  if(se==TRUE){
    pred_y<-Xp%*%m$b.hat
    se<-rowSums(Xp*(Xp%*%m$V))^0.5
    listc<-list(pred_y=pred_y,se=se)
    return(listc)
  }
  else{
    pred_y<-Xp%*%m$b.hat
    return(pred_y)
  }
  
}


plot.pspline<-function(m){
  #plot 1
  plot(x,y,xlab='x',ylab='mu.hat')
  lines(x,m$mu.hat)
  #print(m$mu.hat)
  
  V<-solve(t(m$X)%*%m$X+m$sp*t(m$D)%*%m$D)*m$sig
  
  stan_error<-rowSums(m$X*(m$X%*%m$V))^0.5
  
  upperbound<-c(m$mu.hat+1.96*sqrt(m$sig))
  lowerbound<-c(m$mu.hat-1.96*sqrt(m$sig))
  lines(x,upperbound,lty=2)
  lines(x,lowerbound,lty=2)
  resid<-c()
  
  for(elem in 1:length(y)){
    values2<-(y[elem]-m$mu.hat[elem])
    resid<-c(resid,values2)
    
  }
  plot(m$mu.hat,resid,xlab='mu',ylab='residulas')
  
  qqnorm(resid,xlab='residuals',ylab='residuals')
  
  list2<-list(ll=lowerbound,ul=upperbound,x=x)
  invisible(list2)
  
  
}

a<-pspline(x,y)


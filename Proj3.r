#Alannah Hounat S2434943

#x and y the vectors of x, y data to smooth
#k – the number of basis functions to use
#logsp – the ends of the interval over which to search for the smoothing 
#parameter (log λ scale).
#bord – the B-spline order to use: 3 corresponds to cubic
#pord – the order of difference to use in the penalty. 2 is for the 
#penalty given above.
#ngrid – the number of smoothing parameter values to try. You should use even 
#spacing of values on the log scale.

library(MASS)
x<-mcycle$times
y<-mcycle$accel


#new_x_val<-seq(min(x),max(x),1)

pspline<-function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  
  
  lsp=seq(logsp[1],logsp[2],length=ngrid)
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
  qrx<-qr(X)
  crossD<-crossprod(D)
  formula.bhat<- t(solve(qr.R(qrx)))%*%crossD%*%solve(qr.R(qrx))
  e.decomp<-eigen(formula.bhat)
  
  lam<-e.decomp$values
  lambda<-diag(lam)
  U<-e.decomp$vectors
  
  ## function for ridge regression of y on X with generalized
  ## cross validation lsp is the set of smoothing parameters to try
  edf <- gcv <- sig <-b.hat<-fv <- lsp*0 ## vectors for edf and gcv
  for (i in 1:length(lsp)){
    ## loop over log smoothing parameters
    fr <- get.gcv(y,X,exp(lsp[i]),D,bord,pord,qrx,lambda,U) ## fit
    gcv[i] <- fr$gcv
    edf[i] <- fr$edf
    sig[i] <- fr$sig
    
  }
  #plot.gcv(edf,lsp,gcv) ## plot results
  i.opt <- max(which(gcv==min(gcv))) ## locate minimum
  
  get.gcv(y,X,exp(lsp[i.opt]),D,bord,pord,qrx,lambda,U)## return fit at minimum
  val<-get.gcv(y,X,exp(lsp[i.opt]),D,bord,pord,qrx,lambda,U)
  #print(fv[i.opt])
  #print(b.hat[i.opt])
  #print(sig[i.opt])
}
#fit.ridge
get.gcv <- function (y,X,sp,D,bord,pord,qrx,lambda,U){
  p <- ncol(X);n <- nrow(X)
  ## Compute hat matrix.
  #t(D)%*%D same as crossprod(D)
  
  A<-solve(diag(p)+sp*lambda)
  trA <- sum(diag(A)) 
  #print(length(trA))## Effective degrees of fredom
  ## compute coeff estimates.
  
  
  b.hat<-backsolve(qr.R(qrx),U%*%solve(diag(p)+lambda*sp)%*%t(U)%*% qr.qty(qrx,y)[1:p])
  fv <- X %*% b.hat ## fitted values
  
  sig<-sum((y-fv)^2)/(n-trA)
  
  res<-sqrt(sig)
  V<-solve(t(X)%*%X+sp*t(D)%*%D)*sig
  se2<-rowSums(X*(X%*%V))^0.5
  sm<-c()
  
  #print(b.hat)
  gcv <- sig/(n-trA)
  
  for(elem in 1:n){
    values<-(y[elem]-mean(y))^2
    sm<-c(sm,values)
    
  }
  sum.a<-sum(unlist(sm))
  r_sq<-1-((n-1)*sig)/sum.a
  
  list(b.hat=b.hat,fv=fv,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res, bord=bord,pord=pord,V=V,X=X,D=D,sp=sp)
  vals<-list(b.hat=b.hat,fv=fv,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res, bord=bord,pord=pord,V=V,X=X,D=D,sp=sp)
  #class(list(b.hat=b.hat,mu.hat=fv,sig=sig))<-"pspline"
}



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



predict.pspline<- function(m,x,se=TRUE){
  
  
  new_x_val<-seq(min(x),max(x),by=1)
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
  lines(x,m$fv)
  #print(m$fv)
  
  V<-solve(t(m$X)%*%m$X+m$sp*t(m$D)%*%m$D)*m$sig
  
  se2<-rowSums(m$X*(m$X%*%m$V))^0.5
  
  upperbound<-c(m$fv+1.96*sqrt(m$sig))
  lowerbound<-c(m$fv-1.96*sqrt(m$sig))
  lines(x,upperbound,lty=2)
  lines(x,lowerbound,lty=2)
  resid<-c()
  
  for(elem in 1:length(y)){
    values2<-(y[elem]-m$fv[elem])
    resid<-c(resid,values2)
    
  }
  plot(m$fv,resid,xlab='mu',ylab='residulas')
  
  qqnorm(resid,xlab='residuals',ylab='residuals')
  
  list2<-list(ll=lowerbound,ul=upperbound,x=x)
  invisible(list2)
  
  
}

a<-pspline(x,y)


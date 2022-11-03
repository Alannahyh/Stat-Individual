#Alannah Hounat S2434943

#OVERVIEW
#"Smoothing with basis expansions and penalties"
#P-Splines are regression splines used to smooth data. This is done by fitting
#them by least-squares(b.hat) and a roughness penalty.
#In this code we will be:
#-> Smoothing data x,y with generalised cross validation
# smoothing parameter selection using the pspline function
#-> Reporting details of the model fit using the print.pspline function
#-> Creating new x data and finding new predictions based off these using the 
# predict.pspline function
#-> Plotting: 1. the original x,y data with 95% confidence intervals for the data
#2. the model residuals vs the fitted values
#3. a qqplot of the residuals, all using the plot.pspline function


#INPUT:x<-x data to smooth,y<-y data to smooth,k<-number of basis functions
#logsp<-ends of interval for smoothing parameter,bord<-B-spline order,
#pord<-order of difference to use in the penalty,
#ngrid<-number of smoothing parameters
#OUTPUT:a list of class pspline
#PURPOSE:to find b.hat<-penalized least squares,
#muhat<- estimates of the residuals
#sig<-residual variance
#all of which is used to find gcv<-generalized cross validation and then
#find the optimal gcv score
pspline<-function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  
  #sequence of values for smoothing parameter
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
  
  #eigen decomposition, used for finding lambda(matrix of eigenvalues on diagonal)
  #and U (matrix of eigen vectors)
  e.decomp<-eigen(formula.bhat)
  
  #matrix containing eigenvals in the diagonal
  lambda<-diag(e.decomp$values)
  
  #matrix of eigen vectors
  U<-e.decomp$vectors
  
  edf <- gcv <- sig <- lsp*0 ## vectors for edf and gcv and sigma squared
  for (i in 1:length(lsp)){
    
    ## loop over log smoothing parameters
    gg <- get.gcv(y,X,exp(lsp[i]),D,bord,pord,qrx,lambda,U) ## fit
    gcv[i] <- gg$gcv
    edf[i] <- gg$edf
    sig[i] <- gg$sig
    
  }
  #find the optimum value for gcv
  i.opt <- max(which(gcv==min(gcv))) 

  val<-get.gcv(y,X,exp(lsp[i.opt]),D,bord,pord,qrx,lambda,U)
  class(val)<-"pspline"
  return(val)
  print(val)
}

#INPUT: y<-values, X<-matrix,sp<-smoothing parameters,bord<-,pord<-,
#qrx<- QR decomp of X matrix, lambda<-matrix with eigenvals of X in the diagonal
#U<- matrix of eigenvectors of X
#OUTPUT:GCV,b.hat,mu.hat,r_sq
#PURPOSE:takes in the values produced by pspline and returns 
get.gcv <- function (y,X,sp,D,bord,pord,qrx,lambda,U){
  #number of rows and colums of matrix X
  p <- ncol(X);n <- nrow(X)
  
  #A is the matrix (I+λΛ)^-1
  A<-solve(diag(p)+sp*lambda)
  
  #effective degree of freedom of the model
  trA <- sum(diag(A)) 

  # computing coeff estimates,b.hat
  b.hat<-backsolve(qr.R(qrx),U%*%solve(diag(p)+lambda*sp)%*%t(U)%*% qr.qty(qrx,y)[1:p])
  
  #mu.hat(fitted values)
  mu.hat <- X %*% b.hat 
  
  #sigma squared ie residual variance
  sig<-sum((y-mu.hat)^2)/(n-trA)
  
  #residual standard deviation 
  res<-sqrt(sig)
  
  #crossproduct of X
  crossX<-crossprod(X)
  
  #crossproduct of D
  crossD<-crossprod(D)

  #covariance matrix for the coefficients
  V<-solve(crossX+sp*crossD)*sig

  #standard error 
  stan_error<-rowSums(X*(X%*%V))^0.5
  
  #generalized cross validation(gcv)
  #gcv is used to help us find the most optimal value for λ 
  gcv <- sig/(n-trA)
  
  #empty vector 'sm' for holding the values of y[i]-mean(y)
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
  
  #how well data fits regression model
  r_sq<-1-((n-1)*sig)/sum.a
  
  #list of values to call in later functions using m$
  alist<-list(b.hat=b.hat,mu.hat=mu.hat,gcv=gcv,edf=trA,sig=sig,r_sq=r_sq,k=p,res=res,
       bord=bord,pord=pord,X=X,V=V,sp=sp,stan_error=stan_error,knots=knots)
  class(alist)<-"pspline"
  return(alist)
}

#INPUT:m<- model fit
#OUTPUT:reports detains of model fit m
#PURPOSE:takes in the values produced by pspline and returns 
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


#INPUT:m<- model fit ,x<- new x values,se<-standard error
#OUTPUT:if se==TRUE: list of predicted y values and corresponding standard errors
#PURPOSE:makes predictions for the smooth fit using new_x
predict.pspline<- function(m,x,se=TRUE){

  dk <- diff(range(x))/(m$k-m$bord) ## knot spacing
  knots <- seq(min(x)-dk*m$bord,by=dk,length=m$k+m$bord+1)
  
  #new model matrix for the new data supplied
  Xp <- splines::splineDesign(knots,x,ord=m$bord+1,outer.ok=TRUE)
  
  #if se=TRUE compute and return a list containing the y predictions and 
  #standard error (se) values
  if(se==TRUE){
    #y predictions found from multi0plying the model matrix by b.hat values
    pred_y<-Xp%*%m$b.hat
    
    #standard error 
    se<-rowSums(Xp*(Xp%*%m$V))^0.5
    
    #list containing predicted y values and corresponding standard errors
    listc<-list(pred_y=pred_y,se=se)
    return(listc)
  }
  #if se=FALSE compute and return y predictions
  else{
    pred_y<-Xp%*%m$b.hat
    return(pred_y)
  }
  
}

#INPUT:m<-model fit
#OUTPUT:3 plots: 
#->1.the original x,y data with 95% confidence intervals for the data
#->2.the model residuals vs the fitted values
#->3.a qqplot of the residuals
#PURPOSE:takes in the values produced by pspline and returns plots as described above
plot.pspline<-function(m){
  #Plot 1:the original x,y data with 95% confidence intervals for the data
  plot(x,y,main='Plot of original x,y data',xlab='x',ylab='mu.hat')
  lines(x,m$mu.hat,col='red')
  
  #upperbound for the confidence interval
  upperbound<-c(m$mu.hat+1.96*m$stan_error)
  
  #lowerbound values for the confidence interval
  lowerbound<-c(m$mu.hat-1.96*m$stan_error)
  lines(x,upperbound,col='#FF33CC',lty=2)
  lines(x,lowerbound,col='#6633CC',lty=2)
  
  #empty vector for the residuals
  resid<-c()
  
  #calculating the residuals of the data ie. residual=y[i]-mu.hat[i]
  for(elem in 1:length(y)){
    values2<-(y[elem]-m$mu.hat[elem])
    resid<-c(resid,values2)
    
  }
  #Plot 2: the model residuals vs the fitted values
  plot(m$mu.hat,resid,main='Model Residuals vs Fitted Values',xlab='mu',ylab='residulas')
  
  #Plot 3: qqplot of the residuals
  qqnorm(resid,main="QQPlot of Residuals")
  
  #list containing the upperbound, lowerbound and x values
  list2<-list(ll=lowerbound,ul=upperbound,x=x)
  invisible(list2)
}



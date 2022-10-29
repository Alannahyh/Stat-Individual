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


pspline<- function(x,y,k,logsp,bord,pord,ngrid){
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  #print(D)
  #print('hi')
  #print(X)
  #here Q is the orthogonal matrix we use to find the lambdas
  qrx<-qr(X)
  Q<-qr.Q(qrx)
  R<-qr.R(qrx)
  print(Q)
  print(R)
  e<- eigen(R)
  print(e)
  #p<-ncol(X)
  
  #beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p])
  #print(beta)
  #Q<-qr.Q(qrx)
  #
  #print(qrx)
  #eig<-eigen(Q)
  #print(eig)
  #plot(x,X[,1])
  #plot(x,X[,2])
  #plot(x,X[,3])
  #plot(x,X[,4])
  #plot(x,X[,5])
}

#Note that the mcycle data from the MASS library provides one set of x, y 
#(times, accel) data to try
#rdm data
pspline(x,y,6,c(-5,5),3,2,100)

#pspline(b,b,6,c(-5,5),2,2,100)
#to choose lambda (the smoothing parameter) we need to use some function
#
#sequence of log smoothing parameters - lsp
#ridge<-function(y,X,lsp){
  
#}

#x1<-rnorm(500,0,2)

#print(x1)

#plot(x,X)

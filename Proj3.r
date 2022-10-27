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
pspline<- function(x,y,k,logsp,bord,pord,ngrid){
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  print(x)
  print(X[,1])
  plot(x,X[,1])
}

print(pspline(c(1,2,3,4,5,6,7,8,9,10,11,12),c(5,6,7,3,4,5,6,7,8,9,10,11,12),4,c(-5,5),3,2,100))
#to choose lambda (the smoothing parameter) we need to use some function
#
#sequence of log smoothing parameters - lsp
#ridge<-function(y,X,lsp){
  
#}

#x1<-rnorm(500,0,2)

#print(x1)

#plot(x,X)

gg

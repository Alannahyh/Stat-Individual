#Alannah Houna S2434943



pspline<- function(x,y,k,logsp,bord,pord,ngrid){
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
}

print(pspline(c(1,2,3,4),c(5,6,7,8),5,c(-5,5),3,2,100))
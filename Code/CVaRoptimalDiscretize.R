# Set working directory as the main folder of the code
# Change to your respective directory location
setwd("~/Desktop/GITHUB/Distribution-RA-MDP")
rm(list=ls())
source("Code/UMDP.R")

X = data.frame(V = c(5,10,30),prob = c(0.1,0.6,0.3))
Y = data.frame(V = c(4,12,18),prob = c(0.1,0.6,0.3))
Z = data.frame(V = c(1,13,14),prob = c(0.1,0.6,0.3))
Q = c(iterSum(c(0.1,0.6,0.3)),1)
  
X_cvar= CVAR_multi(X$V,alphas = Q,prob = X$prob)
Y_cvar= CVAR_multi(Y$V,alphas = Q,prob = Y$prob)
Z_cvar= CVAR_multi(Z$V,alphas = Q,prob = Z$prob)

QV = list()
QV[[1]] = Q * X_cvar
QV[[2]] = Q * Y_cvar
QV[[3]] = Q * Z_cvar

plot(Q,QV[[1]],type='l')
lines(Q,QV[[2]],col="red")
lines(Q,QV[[3]],col="blue")

intersect = function(x1,x2,y1,y2,y3,y4){
  x3 = x1
  x4 = x2
  x = ( (x1*y2 - y1*x2 )*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  y = ( (x1*y2 - y1*x2 )*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  return(list(x=x,y=y))
}


ITSC = list()
for (q in 2:length(Q) ){
  ITSC[[paste0(Q[q])]] = list(x = c(),y = c())
  for (a1 in 1:2){
    for (a2 in (a1+1):3){
      list[x,y] = intersect(Q[q-1],Q[q],QV[[a1]][q-1],QV[[a1]][q],QV[[a2]][q-1],QV[[a2]][q])
      ITSC[[paste0(Q[q])]][["x"]] = c(ITSC[[paste0(Q[q])]][["x"]],x)
      ITSC[[paste0(Q[q])]][["y"]] = c(ITSC[[paste0(Q[q])]][["y"]],y)
    }
  }
}
x2 = c(sapply(ITSC,function(L) L$x))
y2 = c(sapply(ITSC,function(L) L$y))
points(x2,y2)



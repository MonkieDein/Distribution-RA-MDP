# Set working directory as the main folder of the code
# Change to your respective directory location
# setwd("~/Desktop/GITHUB/Distribution-RA-MDP")
rm(list=ls())
source("Code/UMDP.R")

X = data.frame(V = c(5,10,30),prob = c(0.1,0.6,0.3))
Y = data.frame(V = c(4,12,18),prob = c(0.1,0.6,0.3))
Z = data.frame(V = c(1,13,14),prob = c(0.1,0.6,0.3))
Q = c(0,iterSum(c(0.1,0.6,0.3)))
  
X_cvar= CVAR_multi(X$V,alphas = Q,prob = X$prob)
Y_cvar= CVAR_multi(Y$V,alphas = Q,prob = Y$prob)
Z_cvar= CVAR_multi(Z$V,alphas = Q,prob = Z$prob)

QV = list()
QV[[1]] = Q * X_cvar
QV[[2]] = Q * Y_cvar
QV[[3]] = Q * Z_cvar

plot(Q,QV[[1]],type='l',ylim = c(0,30),lwd =  1,col="darkgreen",
     main="Quantile VS QCVaR",xlab = "Q",ylab="QCVaR")
lines(Q,QV[[2]],col="red",lwd = 1)
lines(Q,QV[[3]],col="blue",lwd= 1)
for (i in 1:length(X$V)){
  lines(c(Q[i],Q[i+1]),c(X$V[i],X$V[i]),lty =2,col="darkgreen")
  if (i < length(X$V)){ points(Q[i+1],X$V[i+1],col="darkgreen")}
  points(Q[i+1],X$V[i],pch=16,col="darkgreen")
}
for (i in 1:length(Y$V)){
  lines(c(Q[i],Q[i+1]),c(Y$V[i],Y$V[i]),lty =2,col="red")
  if (i < length(Y$V)){ points(Q[i+1],Y$V[i+1],col="red")}
  points(Q[i+1],Y$V[i],pch=16,col="red")
}
for (i in 1:length(Z$V)){
  lines(c(Q[i],Q[i+1]),c(Z$V[i],Z$V[i]),lty =2,col="blue")
  if (i < length(Z$V)){ points(Q[i+1],Z$V[i+1],col="blue")}
  points(Q[i+1],Z$V[i],pch=16,col="blue")
}

intersect = function(x1,x2,y1,y2,y3,y4){
  x3 = x1
  x4 = x2
  x = ( (x1*y2 - y1*x2 )*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  y = ( (x1*y2 - y1*x2 )*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  return(list(x=x,y=y))
}

# 
# ITSC = list()
# for (q in 2:length(Q) ){
#   ITSC[[paste0(Q[q])]] = list(x = c(),y = c())
#   for (a1 in 1:2){
#     for (a2 in (a1+1):3){
#       list[x,y] = intersect(Q[q-1],Q[q],QV[[a1]][q-1],QV[[a1]][q],QV[[a2]][q-1],QV[[a2]][q])
#       ITSC[[paste0(Q[q])]][["x"]] = c(ITSC[[paste0(Q[q])]][["x"]],x)
#       ITSC[[paste0(Q[q])]][["y"]] = c(ITSC[[paste0(Q[q])]][["y"]],y)
#     }
#   }
# }
# x2 = c(sapply(ITSC,function(L) L$x))
# y2 = c(sapply(ITSC,function(L) L$y))
# points(x2,y2)
# 

ITSC2 = list()
for (q in 2:length(Q) ){
  a1 = which.max(c(QV[[1]][q-1],QV[[2]][q-1],QV[[3]][q-1]))
  minq = Inf
  minq_i = a1
  for (a2 in order(sapply(1:3,function(i) QV[[i]][q]),decreasing = TRUE) ){
    if (QV[[a2]][q] > QV[[a1]][q]){
      list[x,y] = intersect(Q[q-1],Q[q],QV[[a1]][q-1],QV[[a1]][q],QV[[a2]][q-1],QV[[a2]][q])
      if (x < minq){
        ITSC2[[paste0(Q[q])]][[paste(a1,a2)]] = list(x=x,y=y)
        ITSC2[[paste0(Q[q])]][[paste(a1,minq_i)]] = NULL
        ITSC2[[paste0(Q[q])]][[paste(minq_i,a2)]] = intersect(Q[q-1],Q[q],QV[[minq_i]][q-1],QV[[minq_i]][q],QV[[a2]][q-1],QV[[a2]][q])
        minq = x
        minq_i = a2
      }
    }
  }
}
x3 = c(sapply(ITSC2,function(L) sapply(L,function(L2) L2$x)))
y3 = c(sapply(ITSC2,function(L) sapply(L,function(L2) L2$y)))
points(x3,y3)

Qcvar= sort(c(Q,x3-1e-15))

QV2 = list()
QV2[[1]] = Qcvar * CVAR_multi(X$V,alphas = Qcvar,prob = X$prob)
QV2[[2]] = Qcvar * CVAR_multi(Y$V,alphas = Qcvar,prob = Y$prob)
QV2[[3]] = Qcvar * CVAR_multi(Z$V,alphas = Qcvar,prob = Z$prob)
QVmax = Reduce(pmax,QV2)
QVwhichmax = apply(sapply(QV2,function(x) x),1,which.max)

Vvar = list()
Vvar[[1]] = VAR_multi(X$V,alphas = Qcvar,prob = X$prob)
Vvar[[2]] = VAR_multi(Y$V,alphas = Qcvar,prob = Y$prob)
Vvar[[3]] = VAR_multi(Z$V,alphas = Qcvar,prob = Z$prob)
for (i in 2:(length(QVwhichmax))){
  lines(c(Qcvar[i-1],Qcvar[i]),c(Vvar[[QVwhichmax[i]]][i],Vvar[[QVwhichmax[i]]][i]),lwd=2,lty=2 )
  points(Qcvar[i-1],Vvar[[QVwhichmax[i]]][i])
  points(Qcvar[i],Vvar[[QVwhichmax[i]]][i],pch=16)
}

Qcvarmax = sapply(2:(length(QVwhichmax)),function(i) Vvar[[QVwhichmax[i]]][i])
probcvar = sapply(2:(length(QVwhichmax)),function(i) Qcvar[i]-Qcvar[i-1])
cvar2 = Qcvar * CVAR_multi(Qcvarmax,alphas = Qcvar,prob = probcvar)
lines(Qcvar,cvar2,lwd=2)


# Set working directory as the main folder of the code
# Change to your respective directory location
# setwd("~/Desktop/GITHUB/Distribution-RA-MDP")
rm(list=ls())
source("Code/UMDP.R")

X = list()
X[[1]] = data.frame(V = c(10,30),prob = c(0.8,0.2))
X[[2]] = data.frame(V = c(4,12,18),prob = c(0.1,0.6,0.3))
X[[3]] = data.frame(V = c(3,15),prob = c(0.2,0.8))

Q = lapply(X,function(x_i) c(0,iterSum(x_i$prob)))
UQ = sort( unique(unlist(Q)) )
ymin = min(unlist(lapply(X,function(x) x$V)))
ymax = max(unlist(lapply(X,function(x) x$V)))
cvarX = lapply(X,function(x_i) CVAR_multi(x_i$V,alphas = UQ,prob = x_i$prob))
varX = sapply(X,function(x_i) VAR_multi(x_i$V,UQ,x_i$prob))
maxvar = apply(varX,1,max)
  
Qcvar = lapply(cvarX,function(x) UQ*x)

color = list("darkgreen","red","blue")

plot(c(),c(),ylim = c(0,ymax),
     xlim=c(0,1),main="Quantile VS QCVaR",xlab = "Q",ylab="QCVaR")
for (a in 1:length(Qcvar)){
  lines(UQ,Qcvar[[a]],col=color[[a]],lwd =  1)
}

for (a in 1:length(Qcvar)){
  for (q in 1:(length(Q[[a]])-1)){
    lines(c(Q[[a]][q],Q[[a]][q+1]),c(VAR_multi(X[[a]]$V,Q[[a]][q+1],X[[a]]$prob),
                             VAR_multi(X[[a]]$V,Q[[a]][q+1],X[[a]]$prob)),
          lty =2,col=color[[a]])
    points(Q[[a]][q],VAR_multi(X[[a]]$V,Q[[a]][q+1],X[[a]]$prob),col=color[[a]])
    points(Q[[a]][q+1],VAR_multi(X[[a]]$V,Q[[a]][q+1],X[[a]]$prob),pch=16,col=color[[a]])
  }
}
# optimize var
for (q in 1:(length(maxvar)-1)){
  lines(c(UQ[q],UQ[q+1]),c(maxvar[q+1],maxvar[q+1]),lty =2,lwd=2)
  points(UQ[q],maxvar[q+1])
  points(UQ[q+1],maxvar[q+1],pch=16)
}


ITSC2 = list()
for (q in 2:length(UQ) ){
  a1 = which.max(c(Qcvar[[1]][q-1],Qcvar[[2]][q-1],Qcvar[[3]][q-1]))
  minq = Inf
  minq_i = a1
  for (a2 in order(sapply(1:3,function(i) Qcvar[[i]][q]),decreasing = TRUE) ){
    if (Qcvar[[a2]][q] > Qcvar[[a1]][q]){
      list[x,y] = intersect(UQ[q-1],UQ[q],Qcvar[[a1]][q-1],Qcvar[[a1]][q],
                            Qcvar[[a2]][q-1],Qcvar[[a2]][q])
      if (x < minq){
        ITSC2[[paste0(UQ[q])]][[paste(a1,a2)]] = list(x=x,y=y)
        ITSC2[[paste0(UQ[q])]][[paste(a1,minq_i)]] = NULL
        ITSC2[[paste0(UQ[q])]][[paste(minq_i,a2)]] = intersect(UQ[q-1],UQ[q],
                                                               Qcvar[[minq_i]][q-1],Qcvar[[minq_i]][q],
                                                               Qcvar[[a2]][q-1],Qcvar[[a2]][q])
        minq = x
        minq_i = a2
      }
    }
  }
}
x3 = c(sapply(ITSC2,function(L) sapply(L,function(L2) L2$x)))
y3 = c(sapply(ITSC2,function(L) sapply(L,function(L2) L2$y)))
points(x3,y3)

UQ2 = sort(unique(c(UQ,x3-1e-15)))
Qcvar2 = lapply(X,function(x_i)UQ2 * CVAR_multi(x_i$V,alphas = UQ2,prob = x_i$prob))

Qcvarmax = Reduce(pmax,Qcvar2)
Qcvarwhichmax = apply(sapply(Qcvar2,function(x) x),1,which.max)
Vvar = lapply(X, function(x_i)VAR_multi(x_i$V,alphas = UQ2,prob = x_i$prob))

for (q in 2:(length(Qcvarwhichmax))){
  lines(c(UQ2[q-1],UQ2[q]),c(Vvar[[Qcvarwhichmax[q]]][q],Vvar[[Qcvarwhichmax[q]]][q]),lwd=2,lty=2 )
  points(UQ2[q-1],Vvar[[Qcvarwhichmax[q]]][q])
  points(UQ2[q],Vvar[[Qcvarwhichmax[q]]][q],pch=16)
}
# Use our new distribution to construct the optimal CVaR curve.
Qvaluemax = sapply(2:(length(Qcvarwhichmax)),function(i) Vvar[[Qcvarwhichmax[i]]][i])
probcvar = sapply(2:(length(Qcvarwhichmax)),function(i) UQ2[i]-UQ2[i-1])
cvar2 = UQ2 * CVAR_multi(Qvaluemax,alphas = UQ2,prob = probcvar)
lines(UQ2,cvar2,lwd=2)


# Set working directory as the main folder of the code
# Change to your respective directory location
# setwd("~/Desktop/GITHUB/Distribution-RA-MDP")
rm(list=ls())
source("Code/UMDP.R")

plotVAR = function(X,Q,lty=1,lwd=1,col="black"){
  points(Q[1],X[1],pch=16,col=col)
  for (q in 1:(length(Q)-1)){
    lines(c(Q[q],Q[q+1]),c(X[q+1],X[q+1]),lty =lty,lwd=lwd,col=col)
    points(Q[q],X[q+1],col=col)
    points(Q[q+1],X[q+1],pch=16,col=col)
  }
}

plotoptVAR = function(X,Q,optcol,lty=1,lwd=1){
  points(Q[1],X[1],pch=16,col=optcol[1])
  for (q in 1:(length(Q)-1)){
    lines(c(Q[q],Q[q+1]),c(X[q+1],X[q+1]),lty =lty,lwd=lwd,col=optcol[q+1])
    points(Q[q],X[q+1],col=optcol[q+1])
    points(Q[q+1],X[q+1],pch=16,col=optcol[q+1])
  }
}

S = 1:3
A = 1:2

Vt2 = list(data.frame(V = c(8,10),prob = c(0.5,0.5)),
           data.frame(V = c(8,12),prob = c(0.4,0.6)),
           data.frame(V = c(9,13),prob = c(0.7,0.3)))

for (s in S){
  plot(c(),c(),ylim = c(0,max(unlist(sapply(S, function(s) Vt2[[s]]$V)))),
     xlim=c(0,1),main=paste0("V2-s",s),xlab = "Q",ylab="V")

  plotVAR(Vt2[[s]]$V[c(1,1:length(Vt2[[s]]$V))],c(0,iterSum(Vt2[[s]]$prob)),lwd = 2)
}

R = list()
for (s in S){
  R[[s]] = list()
  for (a in A){
    R[[s]][[a]] = c(s+a,s-a,s)
  }
}

P = list()
for (s in S){
  P[[s]] = list()
  P[[s]][[1]] = c(0.2,0.5,0.3)
  P[[s]][[2]] = c((s+1)/10,(s+2)/10, (7-2*s)/10)
}


Vsa = lapply(S , function(s) lapply(A, function(a) NULL ) )
for (s in S){
  for (a in A){ 
    V = c(sapply(S, function(s2) Vt2[[s2]]$V + R[[s]][[a]][s2] ))
    prob = c(sapply(S, function(s2) P[[s]][[a]][s2]*Vt2[[s2]]$prob) )
    o = order(V)
    tmp = data.frame(V = V[o],prob = prob[o])
    Vsa[[s]][[a]] = aggregate(x = tmp["prob"], by = tmp[c("V")], FUN = sum)
  }
}

Qt1 = lapply(S,function(s) lapply(Vsa[[s]],function(x_i) c(0,iterSum(x_i$prob))))
UQ1 = lapply(S,function(s) sort( unique(unlist(Qt1[[s]])) ) )
ymin1 = lapply(S,function(s) min(unlist(lapply(Vsa[[s]],function(x) x$V))))
ymax1 = lapply(S,function(s) max(unlist(lapply(Vsa[[s]],function(x) x$V))))
varX1 = lapply(S,function(s) sapply(Vsa[[s]],function(x_i) VAR_multi(x_i$V,UQ1[[s]],x_i$prob)))
maxvar1 = lapply(S,function(s) apply(varX1[[s]],1,max))
maxwhich.var1 = lapply(S,function(s) apply(varX1[[s]],1,which.max))

Vt1 = list()
for (s in S){
  tmp = data.frame(V = maxvar1[[s]][-1],
                   prob = UQ1[[s]][-1]-UQ1[[s]][-length(UQ1[[s]])]) 
  Vt1[[s]] = aggregate(x = tmp["prob"], by = tmp[c("V")], FUN = sum)
}

color = list("darkgreen","blue")

for (s in S){
  plot(c(),c(),ylim = c(0,ymax1[[s]]),
       xlim=c(0,1),main=paste0("V1-s",s),xlab = "Q",ylab="V")
  
  for (a in A){
    plotVAR(VAR_multi(Vsa[[s]][[a]]$V,Qt1[[s]][[a]],Vsa[[s]][[a]]$prob),
            Qt1[[s]][[a]],col=color[[a]],lwd=3,lty=2)
  }
  plotVAR(maxvar1[[s]],UQ1[[s]],lwd = 1,lty=1)
  legend(0.8, 10, legend=c("a1", "a2","a*"),
         col=c("darkgreen","blue","black"), lty=c(2,2,1),cex=0.5)
}

optcol = lapply(S,function(s) sapply(maxwhich.var1[[s]],function(i) color[[i]]))

for (s in S){
  plot(c(),c(),ylim = c(0,ymax1[[s]]),
       xlim=c(0,1),main=paste0("V1-s",s),xlab = "Q",ylab="V")
  plotoptVAR(maxvar1[[s]],UQ1[[s]],optcol[[s]],lwd = 1,lty=1)
  legend(0.8, 10, legend=c("a1", "a2"),
         col=c("darkgreen","blue"), lty=c(1,1),cex=0.5)
}


# The zeroth iteration
Vsa = lapply(S , function(s) lapply(A, function(a) NULL ) )
for (s in S){
  for (a in A){ 
    V = unlist(sapply(S, function(s2) Vt1[[s2]]$V + R[[s]][[a]][s2] ))
    prob = unlist(sapply(S, function(s2) P[[s]][[a]][s2]*Vt1[[s2]]$prob) )
    o = order(V)
    tmp = data.frame(V = V[o],prob = prob[o])
    Vsa[[s]][[a]] = aggregate(x = tmp["prob"], by = tmp[c("V")], FUN = sum)
  }
}
Qt0 = lapply(S,function(s) lapply(Vsa[[s]],function(x_i) c(0,iterSum(x_i$prob))))
UQ0 = lapply(S,function(s) sort( unique(unlist(Qt0[[s]])) ) )
ymin0 = lapply(S,function(s) min(unlist(lapply(Vsa[[s]],function(x) x$V))))
ymax0 = lapply(S,function(s) max(unlist(lapply(Vsa[[s]],function(x) x$V))))
varX0 = lapply(S,function(s) sapply(Vsa[[s]],function(x_i) VAR_multi(x_i$V,UQ0[[s]],x_i$prob)))
maxvar0 = lapply(S,function(s) apply(varX0[[s]],1,max))
maxwhich.var0 = lapply(S,function(s) apply(varX0[[s]],1,which.max))

color = list("darkgreen","blue")

for (s in S){
  plot(c(),c(),ylim = c(0,ymax0[[s]]),
       xlim=c(0,1),main=paste0("V0-s",s),xlab = "Q",ylab="V")
  
  for (a in A){
    plotVAR(VAR_multi(Vsa[[s]][[a]]$V,Qt0[[s]][[a]],Vsa[[s]][[a]]$prob),
            Qt0[[s]][[a]],col=color[[a]],lwd=3,lty=2)
  }
  plotVAR(maxvar0[[s]],UQ0[[s]],lwd = 1)
  legend(0.8, 11, legend=c("a1", "a2","a*"),
         col=c("darkgreen","blue","black"), lty=c(2,2,1),cex=0.5)
}

optcol = lapply(S,function(s) sapply(maxwhich.var0[[s]],function(i) color[[i]]))

for (s in S){
  plot(c(),c(),ylim = c(0,ymax0[[s]]),
       xlim=c(0,1),main=paste0("V0-s",s),xlab = "Q",ylab="V")
  plotoptVAR(maxvar0[[s]],UQ0[[s]],optcol[[s]],lwd = 1,lty=1)
  legend(0.8, 10, legend=c("a1", "a2"),
         col=c("darkgreen","blue"), lty=c(1,1),cex=0.5)
}


# for (a in A){
#   for (q in 1:( length(Qt1[[s]][[a]])-1) ){
#     lines(c(Qt1[[s]][[a]][q],Qt1[[s]][[a]][q+1]),c(VAR_multi(Vsa[[s]][[a]]$V,Qt1[[s]][[a]][q+1],Vsa[[s]][[a]]$prob),
#                                      VAR_multi(Vsa[[s]][[a]]$V,Qt1[[s]][[a]][q+1],Vsa[[s]][[a]]$prob)),
#           lty =2,col=color[[a]])
#     points(Qt1[[s]][[a]][q],VAR_multi(Vsa[[s]][[a]]$V,Qt1[[s]][[a]][q+1],Vsa[[s]][[a]]$prob),col=color[[a]])
#     points(Qt1[[s]][[a]][q+1],VAR_multi(Vsa[[s]][[a]]$V,Qt1[[s]][[a]][q+1],Vsa[[s]][[a]]$prob),pch=16,col=color[[a]])
#   }
# }
# optimize var
# for (q in 1:(length(maxvar1[[s]])-1)){
#   lines(c(UQ1[[s]][q],UQ1[[s]][q+1]),c(maxvar1[[s]][q+1],maxvar1[[s]][q+1]),lty =2,lwd=2)
#   points(UQ1[[s]][q],maxvar1[[s]][q+1])
#   points(UQ1[[s]][q+1],maxvar1[[s]][q+1],pch=16)
# }

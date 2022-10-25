
solveQMDPvar = function(MDP, decimal = 0, horizon = 1000){
  V = list()
  V[[t+1]] = lapply(MDP$S,function(s) list(qBegin = c(0) ,prob = c(1) ,v = c(0) ) )
  for (t in horizon:1){
    V[[t]] = lapply(MDP$S,function(s) list(qBegin = NULL ,prob = NULL ,v = NULL ) )
    for (s in MDP$S){
      V_ts = lapply(MDP$A, function(a) Distribution(MDP,V[[t+1]],s,a))
      Q = sort(unique(unlist(sapply(MDP$A, function(a) V_ts[[a]][["qBegin"]]))))
      VAR_tsa = lapply(MDP$A, function(a) VAR_multi( V_ts[[a]][["v"]] , Q , V_ts[[a]][["prob"]] ) )
      VAR_max = Reduce(pmax,VAR_tsa)
    }
  }
}

Distribution = function(MDP, v, s, a){
  # Compute reward and probability
  X = data.frame(
    v = unlist(sapply(MDP$S,function(s2) MDP$R[s,a,s2] + MDP$gamma * v[[s2]][["v"]])),
    prob = unlist(sapply(MDP$S,function(s2) MDP$P[s,a,s2] * v[[s2]][["prob"]]))
  )
  
  # remove zeros probability
  X = X[X$prob>0,]
  # sum probability that has the same value function
  X = aggregate(x = X["prob"], by = X["v"], FUN = sum)
  
  # sum_{i=1:(n-1)} prob[i] for quantile begin at [n] 
  X[["qBegin"]] = rep(0,nrow(X))
  for (i in 2:nrow(X)){
    X[["qBegin"]][i] = X[["qBegin"]][i-1] + X$prob[i-1]
  }

  return(X)
}







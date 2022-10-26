
# compute sum_{i=1:(n-1)} prob[i] for quantile begin at [n] 
iterSum = function(prob){
  sumprob = prob*0
  if (length(prob) > 1){
    for (i in 2:length(prob)){
      sumprob[i] = sumprob[i-1] + prob[i-1]
    }
  }
  return(sumprob)
}


# Given value function distribution of previous time and current state action pair
# compute value function distribution of next time step.
Distribution = function(MDP, v, s, a, precision = 5){
  # Compute reward and probability distribution
  X = data.frame(
    v = unlist(sapply(MDP$S,function(s2) round(MDP$R[s,a,s2] + MDP$gamma * v[[s2]][["v"]],precision ) )),
    prob = unlist(sapply(MDP$S,function(s2) MDP$P[s,a,s2] * v[[s2]][["prob"]]))
  )
  
  # remove zeros probability
  X = X[X$prob>0,]
  # sum probability that has the same value function
  X = aggregate(x = X["prob"], by = X["v"], FUN = sum)
  
  # sum_{i=1:(n-1)} prob[i] for quantile begin at [n] 
  X[["qBegin"]] = iterSum(X$prob)
  
  return(X)
}



solveQMDPvar = function(MDP, decimal = 5, horizon = 1000){
  V = list()
  PI = list()
  V[[horizon+1]] = lapply(MDP$S,function(s) data.frame(qBegin = c(0) ,prob = c(1) ,v = c(0) ) )
  for (t in horizon:1){
    V[[t]] = list()
    PI[[t]] = list()
    for (s in MDP$S){
      # Retrieve distribution of value function for each s,a
      V_ts = lapply(MDP$A, function(a) Distribution(MDP,V[[t+1]],s,a,decimal))
      # Extract unique quantile and compute value
      Q = sort(unique(unlist(sapply(MDP$A, function(a) V_ts[[a]][["qBegin"]]))))
      if (length(Q) == 1){ prob = 1 } else { prob = c(Q[2:length(Q)],1)-Q }
      VAR_tsa = matrix(sapply(MDP$A, function(a) VAR_multi( V_ts[[a]][["v"]] , Q+1e-8 , V_ts[[a]][["prob"]] ) )
                       ,ncol = MDP$lAl)
      # Optimal policy PI and optimal value at risk
      tmp = data.frame(prob = prob, v =apply(VAR_tsa,1,max),
                       pi =  apply(VAR_tsa,1,which.max))
      tmp = aggregate(x = tmp["prob"], by = tmp[c("pi","v")], FUN = sum)
      # store optimal value function
      PI[[t]][[s]] = tmp$pi
      V[[t]][[s]] = data.frame(qBegin = iterSum(tmp$prob) ,prob = tmp["prob"] ,
                               v = tmp["v"])
    }
  }
  return(list(V=V,PI=PI))
}

solveQMDPcvar = function(MDP, decimal = 5, horizon = 1000){
  V = list()
  PI = list()
  V[[horizon+1]] = lapply(MDP$S,function(s) data.frame(qBegin = c(0) ,prob = c(1) ,v = c(0) ) )
  for (t in horizon:1){
    V[[t]] = list()
    PI[[t]] = list()
    for (s in MDP$S){
      # Retrieve distribution of value function for each s,a
      V_ts = lapply(MDP$A, function(a) Distribution(MDP,V[[t+1]],s,a,decimal))
      # Extract unique quantile and compute value
      Q = sort(unique(unlist(sapply(MDP$A, function(a) V_ts[[a]][["qBegin"]]))))
      if (length(Q) == 1){ prob = 1 } else { prob = c(Q[2:length(Q)],1)-Q }
      VAR_tsa = matrix(sapply(MDP$A, function(a) VAR_multi( V_ts[[a]][["v"]] , Q+1e-8 , V_ts[[a]][["prob"]] ) )
                       ,ncol = MDP$lAl)
      CVAR_tsa = matrix(sapply(MDP$A, function(a) CVAR_multi( V_ts[[a]][["v"]] , Q+1e-8 , V_ts[[a]][["prob"]] ) )
                        ,ncol = MDP$lAl)
      PI[[t]][[s]] = apply(CVAR_tsa,1,which.max)
      # Optimal policy PI and optimal value at risk
      tmp = data.frame(prob = prob, v = sapply(1:length(PI[[t]][[s]]), function(i) VAR_tsa[ i,PI[[t]][[s]][i] ]),
                       pi = PI[[t]][[s]])
      tmp = aggregate(x = tmp["prob"], by = tmp[c("pi","v")], FUN = sum)
      # store optimal value function
      PI[[t]][[s]] = tmp$pi
      V[[t]][[s]] = data.frame(qBegin = iterSum(tmp$prob) ,prob = tmp$prob ,
                               v = tmp$v)
    }
  }
  return(list(V=V,PI=PI))
}

# start from s'(s_0,q_0), v_0 = v(s_0,q_0) then pick a_1 = pi(s_0,q_0), 
# get transition to s_1 obtain some r_0 = r(s_0,a_1,s_1). 
# v_1 = (v_0 - r_0) / \gamma , use v_1 as sufficient statistics to map 
# to our s'.








solveE = function(MDP, horizon = 1000){
  V = list()
  PI = list()
  V[[horizon+1]] = sapply(MDP$S,function(s) 0 ) 
  for (t in horizon:1){
    V[[t]] = rep(-Inf,MDP$lSl)
    PI[[t]] = rep(0,MDP$lSl)
    for (s in MDP$S){
      # Retrieve distribution of value function for each s,a
      V_ts = sapply(MDP$A, function(a) 
        sum(MDP$P[s,a,] * (MDP$R[s,a,] + MDP$gamma * V[[t+1]]) )
        )
      V[[t]][s] = max(V_ts)
      PI[[t]][s] = which.max(V_ts)
    }
  }
  return(list(V=V,PI=PI))
}

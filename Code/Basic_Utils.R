library(gsubfn)  # need 0.7-0 or later
library(stringr)
library(zeallot)

# TURN single outcome mdp.data.frame into it's respective reward,transition (S,A,S') matrix.
prep_MDP = function(mdp.df,MDP){
  P <- array(0, c(MDP$lSl,MDP$lAl,MDP$lSl)) # transition prob dim
  R = array(0, c(MDP$lSl,MDP$lAl,MDP$lSl)) # reward dim 
  dimnames(P) = list(MDP$S,MDP$A,MDP$S)
  
  n_rows = nrow(mdp.df)
  # Assign transition prob and reward matrix w.r.t(S,A,S')
  for (r in 1:n_rows){
    R[MDP$S == mdp.df$idstatefrom[r],MDP$A == mdp.df$idaction[r],MDP$S == mdp.df$idstateto[r]] = mdp.df$reward[r]
    P[MDP$S == mdp.df$idstatefrom[r],MDP$A == mdp.df$idaction[r],MDP$S == mdp.df$idstateto[r]] = mdp.df$probability[r]
  }
  return(list(P=P,R=R))
}

csvToMDP = function(file_name){
  mdp.df = read.csv( paste0(file_name)  ,header = TRUE)
  if ((nrow(unique(mdp.df[,-4])) - nrow(unique(mdp.df[,-c(4,5)])))==0){
    mdp.df = aggregate(x = mdp.df["probability"], by = mdp.df[c("idstatefrom","idaction","idstateto","reward" )], FUN = sum)
  } else {
    warning("error: mdp is not unique, exist reward uncertainty")
  }
  MDP = list()
  # Parse in Basic parameter States, Actions, Outcomes space and their length.
  MDP$S = sort(unique(unique(c(mdp.df$idstatefrom,mdp.df$idstateto)))) 
  MDP$A = sort(unique(mdp.df$idaction))
  MDP$lSl = length(MDP$S) 
  MDP$lAl = length(MDP$A) 
  # Extract Reward and Transition Matrix.
  list[MDP$P,MDP$R] = prep_MDP(mdp.df,MDP)
  return(MDP)
}


# compute intersection for two lines L1 = ((x1,y1),(x2,y2)), L2 = ((x3,y3),(x4,y4))
intersect = function(x1,x2,y1,y2,y3,y4,x3 = x1,x4 = x2){
  x = ( (x1*y2 - y1*x2 )*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  y = ( (x1*y2 - y1*x2 )*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4) )/( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
  return(list(x=x,y=y))
}

# Expectation function which calculate the mean for vector-reward X and probability of occurrence.
E = function(X, prob = NULL){
  if (!is.null(prob)){
    if (length(prob) != length(X)){
      stop("Mismatch Dimensions of prob and value")
    }
    return(sum(X*prob))
  }
  return(mean(X))
}

# Entropic Risk Measure
# Given a set of values, the function evaluate the appropriate entropic risk value
# Alpha refer to the risk aversion parameter where 
# ERM(beta = lim -> 0) is the expectation and ERM(beta = lim -> infty) is worst case
# Entropic Risk Measure with Log-Sum-Exp Trick with Prob
ERM = function(X, beta=0.9, prob = NULL){
  if (beta == 0){
    return(E(X,prob))
  } else if (beta == Inf){
    return(ifelse(is.null(prob),min(X),min(X[prob>0])))
  }
  Y = -beta*X
  C = max(Y)+1
  if (!is.null(prob)){
    if (length(prob) != length(X)){
      stop("Mismatch Dimensions of prob and value")
    }
    # The hat is avoid underflow, we need to shift the value to max of non-zero prob
    Yhat = Y[prob!=0]
    Chat = max(Yhat)
    probhat = prob[prob!=0]
    return(-(Chat+log(sum(exp(Yhat-Chat)*probhat)))/beta)
  }
  return(-(C+log(mean(exp(Y-C))))/beta)
}

# (alpha = lim -> 0) is the worst case, (alpha = lim -> 1) is the mean  
EVAR = function(X,betas,alpha=0.05, prob=NULL, ERMs = NULL){
  if (abs(alpha)<1e-10){
    return(list(value = ifelse(is.null(prob),min(X),min(X[prob>0])), b_index = which.max(betas), b = max(betas)))
  } else{
    if (is.null(ERMs)){
      ERMs = sapply(betas, function(b) ERM(X,beta = b,prob = prob))
    }
    values = ERMs + log(alpha)/betas
    return(list(value = max(values), b_index = which.max(values), b = betas[which.max(values)]  ))
  }
}

# CVAR single alpha, uniform distribute
CVAR = function(X,alpha=0.05){
  X = sort(X)
  n = length(X)
  portion = alpha*n
  x_a = ceiling(portion)
  return((mean(X[1:x_a]) + X[x_a])*(1-(x_a/portion)) )
}

# remove values with 0 probability, sort value and its probability
removeNsort = function(X,prob){
  # remove values with 0 probability
  remain = prob>0
  X = X[remain]
  prob = prob[remain]
  # sort value and its probability
  ord = order(X)
  X = X[ord]
  prob = prob[ord]
  return(list(X,prob))
}
# CVAR function for multiple alphas
# alpha is used as significant level, alpha = 0 equivalent to minimum, alpha = 1 equivalent to average
CVAR_multi = function(X,alphas,prob = NULL){
  n = length(X)
  if (is.null(prob)){
    prob = rep(1/n,n)
  } else if (length(prob) != length(X)){
    stop("CVAR: Mismatch Dimensions of prob and value")
  } else if (abs(sum(prob) - 1) > 1e-15){
    stop("CVAR: Distribution probability does not sum to one (1)")
  }

  # remove values with 0 probability, sort value and its probability
  list[X,prob] = removeNsort(X,prob)

  # sort alphas
  alphas = sort(alphas)
  if (max(alphas)>1 || min(alphas)<0){
    stop("CVAR: Undefined significant level. (alphas) should be an array between 0 and 1.")
  }
  lQl = length(alphas)
  v = alphas*0
  names(v) = alphas

  # Initialize parameter for loop
  Psum = 0
  Vsum = 0
  index = 1
  # sequentially solves for all alpha in alphas
  for (l in 1:lQl){
    k = alphas[l]
    if (k == 0){
      v[l] = X[1]
    } else {
      while (k > (Psum + prob[index] +1e-15) ){
        Psum = sum(prob[1:index])# (Psum + prob[index])
        Vsum = Vsum + prob[index]*X[index]
        index = index + 1
      }
      v[l] = (Vsum + (k-Psum)*X[index])/k  # This is Piecewise Linear
    }
  }
  return(v)
}

# Total Discounted Return (Cost)
# TDR take in vector of returns V[t=0,t=1,...] and discount factor 
# Calculated the total discounted return
TDR = function(V, discount=0.9){
  D = sapply(1:length(V),function(t) discount^(t-1))
  return(sum(D * V))
}

# check directory if DNE, create folder.
wdir = function(directory_name){
  if (!dir.exists(directory_name)){
    cat("Directory",directory_name,"Not Exist. Creating Directory...\n")
    dir.create(directory_name)
  }
  return(directory_name)
}

# Generate Monte Carlo sampling instances for every time step. 
MonteCarloSamplingS_ = function(n,t,folder_name,seed = 0){
  set.seed(seed)
  for (i in 1:n){
    write.csv(data.frame(S_ = runif(t)),paste0(folder_name,"/instance_",i,".csv"))
  }
}

# Draw Next State S'
drawS_ = function(weights,choice){
  choiceIndex = 1
  for (w in weights){
    choice = choice - w
    if (choice <= 1e-15){
      return(choiceIndex)
    }
    choiceIndex = choiceIndex + 1
  }
}

# Value at Risk
VAR = function(X,alpha = 0.05,prob = NULL){
  if (is.null(prob)){
    return(quantile(X,alpha,type = 1,names = FALSE))
  }
  if (abs(sum(prob) - 1) > 1e-15){
    stop("VAR: Distribution probability does not sum to one (1)")
  }
  # sort value and its probability
  ord = order(X)
  X = X[ord]
  prob = prob[ord]
  return(X[drawS_(prob,alpha)])
}

# VAR function for multiple alphas
VAR_multi = function(X,alphas,prob = NULL){
  if (is.null(prob)){
    v = quantile(X,alphas,type = 1,names = FALSE)
    names(v) = alphas
    return(v)
  } else if (length(prob) != length(X)){
    stop("VAR: Mismatch Dimensions of prob and value")
  } else if (abs(sum(prob) - 1) > 1e-15){
    stop("VAR: Distribution probability does not sum to one (1)")
  }
  
  # remove values with 0 probability, sort value and its probability
  list[X,prob] = removeNsort(X,prob)
  
  # sort alphas
  alphas = sort(alphas)
  if (max(alphas)>1 || min(alphas)<0){
    stop("VAR: Undefined significant level. (alphas) should be an array between 0 and 1.")
  }
  lQl = length(alphas)
  v = alphas*0
  names(v) = alphas
  
  index = 1
  Psum = 0
  for (l in 1:lQl){
    k = alphas[l]
    while (k > (Psum + prob[index] + 1e-15) ){
      Psum = sum(prob[1:index]) # Psum + prob[index]
      index = index + 1
    }
    v[l] = X[index]  # This is Piecewise Linear
  }
  return(v)
}

code = function(codetext){
  return( eval(parse(text=codetext)) )
}

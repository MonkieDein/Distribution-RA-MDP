# Set working directory as the main folder of the code
# Change to your respective directory location
setwd("~/Desktop/GITHUB/Distribution-RA-MDP")
rm(list=ls())
source("Code/UMDP.R")

# Domain put to test
domains = c("inventory_generic","ruin","riverswim","inventory","machine","population" ,"cancer") #  : file too big cannot push to git.
# Folders
TrainFold = wdir(paste0(wdir("Eval/"),"train/"))
TestFold = wdir(paste0(wdir("Eval/"),"test/"))
folder_name = "domains_output/"

Time = list()
TrainOut = list()

domain = "riverswim"
# for (domain in domains){
  cat("domain",domain)
  
  # Define parameter specified for the domain
  Time[[domain]]      <- list()
  TrainOut[[domain]]  <- list()
  saveFold            <- wdir(paste0(TrainFold,domain,"/"))
  
  # If trained before then use cache
  if (file.exists(paste0(saveFold,"Train.RData"))){
    cat(paste0("Use Cache Trained Output \n",domain))
    next
  }
  
  # Parse MDP information given csv folder
  MDP = csvToMDP(paste0(folder_name,"output_",domain,".csv"))
  # Parse in other MDP parameters.
  param = read.csv(paste0(folder_name,"/parameters.csv"),header = TRUE,stringsAsFactors=F)
  list[,MDP$risk_evar,MDP$risk_erm,MDP$gamma,MDP$tolerance,MDP$vspan,MDP$horizon, MDP$S_0]=param[param$domain == domain,] 
  
  # Value function stores V[[t]][[s]][c("qBegin","prob","V")] as distribution of the value function
  # pi on another hand stores pi[[t]][[s]]["qBegin"] as policy for the quantile level.
  TrainOutVAR <- solveQMDPvar(MDP,decimal=1,horizon = 30)
  # CVaR is discretize via optimal var discretization therefore is not accurate.
  TrainOutCVAR <- solveQMDPcvar(MDP,decimal=1,horizon = 50)
  # For expectation V[[t]][s] and pi[[t]][s] stores the optimal value function and policy.
  TrainOutE <- solveE(MDP,horizon = 50)
  # solve MDP with each algorithms
  # for (algo in algorithms){
  #   Time[[domain]][[algo]]            <- list()
  #   Time[[domain]][[algo]][["Start"]] <- Sys.time()
  #   TrainOut[[domain]][[algo]]        <- code( paste0("solve",algo,"(MDP, algoParam[['",algo,"']])") )
  #   Time[[domain]][[algo]][["End"]]   <- Sys.time()
  # }

  # rm(list = ls.str(mode = 'numeric'))
  # rm(list = ls.str(mode = 'function'))
  # save.image(file = paste0(saveFold,"Train.RData"))
  # 
  # rm(list = ls.str(mode = 'list'))
# }






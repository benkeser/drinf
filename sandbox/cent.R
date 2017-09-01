#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

ns <- c(500, 1000, 5000)
bigB <- 2000
parm <- expand.grid(seed=1:bigB,
                    n=ns)

# source in simulation Functions
source("~/drinf/makeData.R")
# load drinf
library(drinf, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2/")

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  for(i in 1:nrow(parm)){
     set.seed(parm$seed[i])
     dat <- makeData(n=parm$n[i])
     save(dat, file=paste0("~/drinf/scratch/dataList_n=",parm$n[i],
                           "_seed=",parm$seed[i],".RData"))
   }
   print(paste0('initial datasets saved to: ~/drinf/scratch/dataList ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # load data
    load(paste0("~/drinf/scratch/dataList_n=",parm$n[i],
                "_seed=",parm$seed[i], ".RData"))
    
    # set seed
    set.seed(parm$seed[i])
    
    # faster to call mean.tmle
    object <- drinf.tmle(
    L0 = dat$L0, L1 = dat$L1, L2 = dat$L2, A0 = dat$A0, A1 = dat$A1, 
    abar = c(1,1), 
    SL.Q = "SL.hal9001",
    SL.g = "SL.glm", 
    SL.Qr = "SL.hal9001",
    SL.gr = "SL.hal9001",
    flucOrd = c("targetg0","targetg1","redReg",
               "targetQ2","targetQ1","redReg"),
    universal = FALSE,
    universalStepSize = 1e-4,  
    return.models = FALSE,
    verbose = TRUE,
    maxIter = 1,
    return.ltmle = TRUE,
    allatonce = FALSE,
    tolg = 1e-2,
    tolQ = 1e-2
    )
    drtmle_ci <- rep(object$est,2) + c(-1.96, 1.96) * rep(object$se,2)
    ltmle_ci <- rep(object$est.ltmle,2) + c(-1.96, 1.96) * rep(object$se.ltmle,2)

    # computed locally
    truth <- 1.16667    
    # output should look like 
    # seed, n, truth
    # drtmle est, ci, coverage
    # ltmle est, ci, coverage 
    out <- c(parm$seed[i], parm$n[i], truth, 
             object$est, drtmle_ci,
             as.numeric(drtmle_ci[1] < truth & drtmle_ci[2] > truth),
             object$est.ltmle, ltmle_ci,
             as.numeric(ltmle_ci[1] < truth & ltmle_ci[2] > truth))

    # save output 
    save(out, file = paste0("~/drinf/scratch/out_n=",
                            parm$n[i],"_seed=",parm$seed[i],".RData.tmp"))
    file.rename(paste0("~/drinf/scratch/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],".RData.tmp"),
                paste0("~/drinf/scratch/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
    ns <- c(500, 1000, 5000)
    bigB <- 2000
    parm <- expand.grid(seed=1:bigB,
                        n=ns)


    rslt <- NULL
    for(i in 1:nrow(parm)){
        tmp <- tryCatch({
            load(paste0("~/drinf/scratch/out_n=",
                        parm$n[i],"_seed=",parm$seed[i],".RData"))
            out
        }, error=function(e){
          c(parm$seed[i], parm$n[i], rep(NA,11))
        })
        rslt <- rbind(rslt, tmp)
    }
    # format
    out <- data.frame(rslt)
    colnames(out) <- c("seed","n","truth","drtmle","drtmle_cil","drtmle_ciu","drtmle_cov",
                       "ltmle","ltmle_cil","ltmle_ciu","ltmle_cov")
    save(out, file=paste0('~/drinf/out/allOut.RData'))
}
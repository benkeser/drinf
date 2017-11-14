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

ns <- c(1000)
bigB <- 1000
g <- c("SL.hal9001","SL.glm")
Q <- c("SL.hal9001","SL.glm")
cv <- c(1,5)
# g <- c("SL.glm.interaction")
# Q <- c("SL.glm.interaction")
parm <- expand.grid(seed=1:bigB,
                    n=ns, g = g, Q = Q, cv = cv, 
                    stringsAsFactors = FALSE)
# n = 500 => up to 49 seconds => get about 20 done per run
# n = 1000 => up to 147 seconds = 2.5 minutes => get about 10 done per run
# n = 5000 => just do one per run
parm$g[(parm$g == "SL.glm" & parm$Q == "SL.glm")] <- "SL.glm.interaction"
parm$Q[(parm$g == "SL.glm.interaction" & parm$Q == "SL.glm")] <- "SL.glm.interaction"

# parm <- parm[1,,drop=FALSE]
# source in simulation Functions
source("~/drinf/makeData.R")
# load drinf
library(drinf, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2/")
library(gam, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2/")
library(hal9001, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2/")
library(SuperLearner)
library(methods)


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
    
    system.time(
    # faster to call mean.tmle
    object <- drinf.tmle(
    L0 = dat$L0, L1 = dat$L1, L2 = dat$L2, 
    A0 = dat$A0, A1 = dat$A1, 
    abar = c(1,1), 
    SL.Q = parm$Q[i],
    SL.g = parm$g[i], 
    cvFolds = parm$cv[i],
    SL.Qr = "SL.gam",
    SL.gr = "SL.gam",
    flucOrd = c("targetg0","targetg1","redReg",
               "targetQ2","targetQ1","redReg"),
    return.models = FALSE,
    verbose = TRUE,
    maxIter = 25,
    return.ltmle = TRUE,
    allatonce = FALSE,
    tolg = 1e-2,
    tolQ = 1e-2, stratify = TRUE
    )
    )
    drtmle_ci <- rep(object$est,2) + c(-1.96, 1.96) * rep(object$se,2)
    ltmle_ci <- rep(object$est.ltmle,2) + c(-1.96, 1.96) * rep(object$se.ltmle,2)

    # computed locally
    truth <- 2.700148    
    # output should look like 
    # seed, n, truth
    # drtmle est, ci, coverage
    # ltmle est, ci, coverage 
    # drtmle iterations
    # drtmle meanIC
    out <- c(parm$seed[i], parm$n[i], truth, 
             parm$Q[i], parm$g[i],
             object$est, drtmle_ci,
             object$est_trace, # tmles with maxIter 1:10
             object$se_trace,
             as.numeric(drtmle_ci[1] < truth & drtmle_ci[2] > truth),
             object$est.ltmle, ltmle_ci,
             as.numeric(ltmle_ci[1] < truth & ltmle_ci[2] > truth),
             object$iter, unlist(object$ic))

    # save output 
    save(out, file = paste0("~/drinf/out/out_n=",
                            parm$n[i],"_seed=",parm$seed[i],
                            "_Q=",parm$Q[i],"_g=",parm$g[i],
                            "_cvFolds=",parm$cv[i],".RData.tmp"))
    file.rename(paste0("~/drinf/out/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],
                       "_Q=",parm$Q[i],"_g=",parm$g[i],"_cvFolds=",parm$cv[i],
                       ".RData.tmp"),
                paste0("~/drinf/out/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],
                       "_Q=",parm$Q[i],"_g=",parm$g[i],"_cvFolds=",parm$cv[i],
                       ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
    # ns <- c(500, 1000, 5000)
    # bigB <- 1000
    # g <- c("SL.hal9001","SL.glm")
    # Q <- c("SL.hal9001","SL.glm")
    # # g <- c("SL.glm.interaction")
    # # Q <- c("SL.glm.interaction")
    # parm <- expand.grid(seed=1:bigB,
    #                     n=ns, g = g, Q = Q,
    #                     stringsAsFactors = FALSE)
    # rslt <- NULL
    # for(i in 1:nrow(parm)){
    #     tmp_1 <- tryCatch({
    #         load(paste0("~/drinf/out/out_n=",
    #                     parm$n[i],"_seed=",parm$seed[i],
    #                    "_Q=",parm$Q[i],"_g=",parm$g[i],
    #                    "_cvFolds=1.RData"))
    #         out
    #     }, error=function(e){
    #       rep(NA,17 + 25*2)
    #     })
    #     tmp_5 <- tryCatch({
    #         load(paste0("~/drinf/out/out_n=",
    #                     parm$n[i],"_seed=",parm$seed[i],
    #                    "_Q=",parm$Q[i],"_g=",parm$g[i],
    #                    "_cvFolds=5.RData"))
    #         out[-(1:5)]
    #     }, error=function(e){
    #       rep(NA, 17 + 25*2 - 5)
    #     })
    #     tmp <- c(tmp_1, tmp_5)
    #     rslt <- rbind(rslt, tmp)
    # }
    # # format
    # out <- data.frame(rslt)
    # sim_names <- c("drtmle", "drtmle_cil","drtmle_ciu",
    #                    paste0("drtmle_maxIter",1:10),
    #                    "drtmle_cov",
    #                    "ltmle","ltmle_cil","ltmle_ciu","ltmle_cov",
    #                    "drtmle_iter", "origIC","missQIC","missgIC")
    # colnames(out) <- c("seed","n","truth","Q","g", sim_names,
    #                    paste0("cv_", sim_names))

    # out[,(1:ncol(out))[c(-4,-5)]] <- apply(out[,(1:ncol(out))[c(-4,-5)]], 2, as.numeric)
    # save(out, file=paste0('~/drinf/out/allOut_pluscv.RData'))

    # # post processing
    # getBias <- function(out, n, Q, g){
    #   rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
    #   bias <- by(rslt, rslt$n, function(x){
    #     bias_drtmle <- mean(x$drtmle - x$truth, na.rm = TRUE)
    #     bias_drtmle_1 <- mean(x$drtmle_maxIter1 - x$truth, na.rm = TRUE)
    #     bias_ltmle <- mean(x$ltmle - x$truth, na.rm = TRUE)
    #     c(nrow(x), bias_drtmle, bias_drtmle_1, bias_ltmle)
    #   })
    #   bias
    # }
    # getBias(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.glm")
    # getBias(out, n = c(500,1000,5000), g = "SL.hal9001", Q = "SL.glm")
    # getBias(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.hal9001")
    # getBias(out, n = c(500,1000,5000), Q = "SL.glm.interaction", g = "SL.glm.interaction")

    # getRootNBias <- function(out, n, Q, g){
    #   rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
    #   rootn_bias <- by(rslt, rslt$n, function(x){
    #     bias_drtmle <- sqrt(x$n[1])*mean(x$drtmle - x$truth, na.rm = TRUE)
    #     bias_drtmle_1 <- mean(x$drtmle_maxIter1 - x$truth, na.rm = TRUE)        
    #     bias_ltmle <- sqrt(x$n[1])*mean(x$ltmle - x$truth, na.rm = TRUE)
    #     c(nrow(x), bias_drtmle, bias_drtmle_1, bias_ltmle)
    #   })
    #   rootn_bias
    # }
    # getRootNBias(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.glm")
    # getRootNBias(out, n = c(500,1000,5000), g = "SL.hal9001", Q = "SL.glm")
    # getRootNBias(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.hal9001")
    # getRootNBias(out, n = c(500,1000,5000), Q = "SL.glm.interaction", g = "SL.glm.interaction")

    # getCov <- function(out, n, Q, g){
    #   rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
    #   cov <-  by(rslt, rslt$n, function(x){
    #     cov_drtmle <- mean(x$drtmle_cov, na.rm = TRUE)
    #     cov_ltmle <- mean(x$ltmle_cov, na.rm = TRUE)
    #     c(cov_drtmle, cov_ltmle)
    #   })
    #   cov
    # }
    # getCov(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.glm")
    # getCov(out, n = c(500,1000,5000), g = "SL.hal9001", Q = "SL.glm")
    # getCov(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.hal9001")
    # getCov(out, n = c(500,1000,5000), Q = "SL.glm.interaction", g = "SL.glm.interaction")

    # getIC <- function(out, n, Q, g){
    #   rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
    #   ic <-  by(rslt, rslt$n, function(x){
    #     colMeans(x[ , grepl("IC", colnames(x))])
    #   })
    #   ic
    # }
    # getIC(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.glm")
    # getIC(out, n = c(500,1000,5000), g = "SL.hal9001", Q = "SL.glm")
    # getIC(out, n = c(500,1000,5000), Q = "SL.hal9001", g = "SL.hal9001")
    # getIC(out, n = c(500,1000,5000), Q = "SL.glm.interaction", g = "SL.glm.interaction")


}
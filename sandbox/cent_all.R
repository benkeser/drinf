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

ns <- c(500,1000,5000,7000)
bigB <- 500
g <- c("SL.hal9001","SL.glm")
Q <- c("SL.hal9001","SL.glm")
cv <- 1
# g <- c("SL.glm.interaction")
# Q <- c("SL.glm.interaction")
parm <- expand.grid(seed=1:bigB,
                    n=ns, g = g, Q = Q, cv = cv, 
                    stringsAsFactors = FALSE)
# # n = 500 => up to 49 seconds => get about 20 done per run
# # n = 1000 => up to 147 seconds = 2.5 minutes => get about 10 done per run
# # n = 5000 => just do one per run
# parm$g[(parm$g == "SL.glm" & parm$Q == "SL.glm")] <- "SL.glm.interaction"
# parm$Q[(parm$g == "SL.glm.interaction" & parm$Q == "SL.glm")] <- "SL.glm.interaction"

# load('~/drinf/out/noboot_allOut_nocv_newest.RData')
# redo_parm <- out[is.na(out$truth),c(1,2,4,5)]
# redo_parm$cv <- 1
# redo_parm$Q <- as.character(redo_parm$Q)
# redo_parm$g <- as.character(redo_parm$g)
# # only do seed 1:100 
# parm <- redo_parm[redo_parm$seed <= 500, ]

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
  # for(i in 1:nrow(parm)){
  #    set.seed(parm$seed[i])
  #    dat <- makeData(n=parm$n[i])
  #    save(dat, file=paste0("~/drinf/scratch/dataList_n=",parm$n[i],
  #                          "_seed=",parm$seed[i],".RData"))
  #  }
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
    maxIter <- 5

    tm <- system.time(
      # faster to call mean.tmle
      object <- drinf.tmle(
      L0 = dat$L0, L1 = dat$L1, L2 = dat$L2, 
      A0 = dat$A0, A1 = dat$A1, 
      abar = c(1,1), 
      SL.Q = parm$Q[i],
      SL.g = parm$g[i], 
      cvFolds = parm$cv[i],
      SL.Qr = c("SL.gam","SL.glm","SL.mean","SL.earth","SL.hal9001"),
      SL.gr = c("SL.gam","SL.glm","SL.mean","SL.earth","SL.hal9001"),
      flucOrd = c("targetg0","targetg1",
                  "targetQ2","targetQ1"),
      # flucOrd = c("targetQg", "redReg"),
      return.models = FALSE,
      verbose = FALSE,
      maxIter = maxIter,
      return.ltmle = TRUE,
      allatonce = FALSE,
      tolg = 1e-2,
      tolQ = 1e-2, stratify = TRUE
      )
    )

    # bootstrap
    # nBoot <- 200
    # estMatrix <- rep(NA, nBoot)
    # for(j in 1:nBoot){
    #     idx <- sample(1:parm$n[i], replace = TRUE)
    #     # faster to call mean.tmle
    #     boot <- drinf.tmle(
    #         L0 = dat$L0[idx,,drop=FALSE], 
    #         L1 = dat$L1[idx,,drop=FALSE], 
    #         L2 = dat$L2[idx], 
    #         A0 = dat$A0[idx], 
    #         A1 = dat$A1[idx], 
    #         abar = c(1,1), 
    #         SL.Q = parm$Q[i],
    #         SL.g = parm$g[i], 
    #         cvFolds = parm$cv[i],
    #         SL.Qr = c("SL.gam","SL.glm","SL.mean","SL.hal9001"),
    #         SL.gr = c("SL.gam","SL.glm","SL.mean","SL.hal9001"),
    #         flucOrd = c("redReg","targetg0","targetg1",
    #                     "redReg","targetQ2","targetQ1"),
    #         # flucOrd = c("targetQg", "redReg"),
    #         return.models = FALSE,
    #         verbose = FALSE,
    #         maxIter = maxIter,
    #         return.ltmle = TRUE,
    #         only.ltmle = TRUE,
    #         allatonce = FALSE,
    #         tolg = 1e-2,
    #         tolQ = 1e-2, stratify = TRUE
    #     )
    #     estMatrix[j] <- boot$est.ltmle
    # }
    # ltmle_boot_ci <- quantile(estMatrix, c(0.025, 0.975))

    # drtmle with max(if) < 1/n stopping criteria
    # drtmle_ci <- rep(object$est,2) + c(-1.96, 1.96) * rep(object$se,2)
    # drtmle with norm(if) < 1/n stopping criteria
    drtmle_norm_n <- object$est_trace[object$n_norm_iter]
    se_drtmle_norm_n <- object$se_trace[object$n_norm_iter]
    drtmle_norm_n_ci <- rep(drtmle_norm_n,2) + c(-1.96, 1.96) * rep(se_drtmle_norm_n,2)
    # drtmle with norm(if) < 1/sqrt(n) stopping criteria
    drtmle_norm_sqrt_n <- object$est_trace[object$sqrt_n_norm_iter]
    se_drtmle_norm_sqrt_n <- object$se_trace[object$sqrt_n_norm_iter]
    drtmle_norm_sqrt_n_ci <- rep(drtmle_norm_sqrt_n,2) + c(-1.96, 1.96) * rep(se_drtmle_norm_sqrt_n,2)
    
    # drtmle with max(if) < 1/sqrt(n) stopping criteria
    drtmle_max_sqrt_n <- object$est_trace[object$sqrt_n_max_iter]
    se_drtmle_max_sqrt_n <- object$se_trace[object$sqrt_n_max_iter]
    drtmle_max_sqrt_n_ci <- rep(drtmle_max_sqrt_n,2) + c(-1.96, 1.96) * rep(se_drtmle_max_sqrt_n,2)
    
    # drtmle with max(if) < 1/n stopping criteria
    drtmle_max_n <- object$est_trace[object$n_max_iter]
    se_drtmle_max_n <- object$se_trace[object$n_max_iter]
    drtmle_max_n_ci <- rep(drtmle_max_n,2) + c(-1.96, 1.96) * rep(se_drtmle_max_n,2)
    
    drtmle_ci_trace <- rep(object$est_trace,2) + c(rep(-1.96, maxIter), rep(1.96, maxIter)) * rep(object$se_trace, 2)
    ltmle_ci <- rep(object$est.ltmle,2) + c(-1.96, 1.96) * rep(object$se.ltmle,2)

    # computed locally
    truth <- 1.3    
    # output should look like 
    # seed, n, truth
    # drtmle est, ci, coverage
    # ltmle est, ci, coverage 
    # drtmle iterations
    # drtmle meanIC
    out <- c(parm$seed[i], parm$n[i], truth, 
             parm$Q[i], parm$g[i],
             drtmle_max_n, drtmle_max_n_ci,
             as.numeric(drtmle_max_n_ci[1] < truth & drtmle_max_n_ci[2] > truth),
             drtmle_max_sqrt_n, drtmle_max_sqrt_n_ci,
             as.numeric(drtmle_max_sqrt_n_ci[1] < truth & drtmle_max_sqrt_n_ci[2] > truth),
             drtmle_norm_n, drtmle_norm_n_ci,
             as.numeric(drtmle_norm_n_ci[1] < truth & drtmle_norm_n_ci[2] > truth),
             drtmle_norm_sqrt_n, drtmle_norm_sqrt_n_ci,
             as.numeric(drtmle_norm_sqrt_n_ci[1] < truth & drtmle_norm_sqrt_n_ci[2] > truth),
             object$est_trace, # tmles with maxIter 1:25
             object$se_trace,
             # add confidence intervals for other estimators 
             object$est.ltmle, ltmle_ci, # ltmle_boot_ci,
             as.numeric(ltmle_ci[1] < truth & ltmle_ci[2] > truth),
             # as.numeric(ltmle_boot_ci[1] < truth & ltmle_boot_ci[2] > truth),
             object$iter, object$sqrt_n_max_iter, object$n_max_iter,
             object$n_norm_iter, object$sqrt_n_norm_iter,
             unlist(object$ic), tm)

    # save output 
    save(out, file = paste0("~/drinf/out/outincomp_n=",
                            parm$n[i],"_seed=",parm$seed[i],
                            "_Q=",parm$Q[i],"_g=",parm$g[i],
                            "_cvFolds=",parm$cv[i],".RData.tmp"))
    file.rename(paste0("~/drinf/out/outincomp_n=",
                       parm$n[i],"_seed=",parm$seed[i],
                       "_Q=",parm$Q[i],"_g=",parm$g[i],"_cvFolds=",parm$cv[i],
                       ".RData.tmp"),
                paste0("~/drinf/out/outincomp_n=",
                       parm$n[i],"_seed=",parm$seed[i],
                       "_Q=",parm$Q[i],"_g=",parm$g[i],"_cvFolds=",parm$cv[i],
                       ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
  ns <- c(500,1000,5000,7000)
  bigB <- 500
  g <- c("SL.hal9001","SL.glm")
  Q <- c("SL.hal9001","SL.glm")
  cv <- c(1)
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

    rslt <- matrix(NA, nrow = nrow(parm), ncol = 43)
    for(i in 1:nrow(parm)){
        tmp_1 <- tryCatch({
            load(paste0("~/drinf/out/outincomp_n=",
                        parm$n[i],"_seed=",parm$seed[i],
                       "_Q=",parm$Q[i],"_g=",parm$g[i],
                       "_cvFolds=1.RData"))
            out
        }, error=function(e){
          c(parm$seed[i], parm$n[i], NA, parm$Q[i], parm$g[i], rep(NA, 43 - 5))
        })
        # tmp_5 <- tryCatch({
        #     load(paste0("~/drinf/out/out_n=",
        #                 parm$n[i],"_seed=",parm$seed[i],
        #                "_Q=",parm$Q[i],"_g=",parm$g[i],
        #                "_cvFolds=5.RData"))
        #     out[-(1:5)]
        # }, error=function(e){
        #   rep(NA, 17 + 25*2 - 5)
        # })
        # tmp <- c(tmp_1, tmp_5)
        rslt[i,] <- tmp_1
    }
    # # format
    out <- data.frame(rslt)
    sim_names <- c("seed","n","truth","Q","g",
                   paste0("max_n_", c("est","cil","ciu","cov")),
                   paste0("max_sqrt_n_", c("est","cil","ciu","cov")),
                   paste0("norm_n_", c("est","cil","ciu","cov")),
                   paste0("norm_sqrt_n_", c("est","cil","ciu","cov")),
                   paste0("drtmle_maxIter",1:5),
                   paste0("se_drtmle_maxIter",1:5),
                   "ltmle","ltmle_cil","ltmle_ciu",
                   # "ltmleboot_cil", 
                   # "ltmleboot_ciu", 
                   "ltmle_cov", 
                   # "ltmleboot_cov",
                   paste0(c("total_","sqrt_n_max_","n_max_","n_norm_","sqrt_n_norm_"),"iter"),
                   "origIC","missQIC","missgIC")
    colnames(out) <- sim_names

    out[,(1:ncol(out))[c(-4,-5)]] <- apply(out[,(1:ncol(out))[c(-4,-5)]], 2, as.numeric)
    save(out, file=paste0('~/drinf/out/allOut_incomp.RData'))

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
    getRootNBias <- function(out, n, Q, g, est = c("max_sqrt_n_est",
                                                   "norm_sqrt_n_est",
                                                   paste0("drtmle_maxIter",1:5),
                                                   "ltmle")){
      rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
      rootn_bias <- by(rslt, rslt$n, function(x){
        o <- matrix(c(sum(!is.na(x$truth)), rep(NA, length(est))), nrow = 1)
        ct <- 1
        for(e in est){
          # browser()
          ct <- ct + 1
          o[ct] <- sqrt(x$n[1])*mean(x[,e] - x$truth, na.rm = TRUE)
        }
        colnames(o) <- c("nsim", est)
        o
      })
      ou <- Reduce(rbind, rootn_bias)
      ou <- cbind(unique(rslt$n), ou)
      ou
    }
    getRootNBias(out, n = c(500,1000,5000,7000), Q = "SL.hal9001", g = "SL.glm")
                 # est = paste0("cv_drtmle_maxIter", 1:25))
    getRootNBias(out, n = c(500,1000,5000,7000), g = "SL.hal9001", Q = "SL.glm")
                 # est = paste0("cv_drtmle_maxIter", 1:25))
    getRootNBias(out, n = c(500,1000,5000,7000), Q = "SL.hal9001", g = "SL.hal9001")
                 # est = paste0("cv_drtmle_maxIter", 1:25))
    getRootNBias(out, n = c(500,1000,5000,7000), Q = "SL.glm.interaction", g = "SL.glm.interaction")
                 # est = paste0("cv_drtmle_maxIter", 1:25))

    getCov <- function(out, n, Q, g,est = c("max_sqrt_n",
                                             "norm_sqrt_n",
                                             paste0("drtmle_maxIter",1:5),
                                             "ltmle")){
      rslt <- out[out$n %in% n & out$Q %in% Q & out$g %in% g, ]
      cov <- by(rslt, rslt$n, function(x){
        o <- matrix(c(sum(!is.na(x$truth)), rep(NA, length(est) + 1)), nrow = 1)
        ct <- 1
        for(e in est){
          # browser()
          ct <- ct + 1
          cov_avail <- any(grepl(paste0(e,"_cov"), colnames(rslt)))
          if(cov_avail){
            o[,ct] <- mean(x[,paste0(e,"_cov")], na.rm = TRUE)
          }else{
            this_est <- x[,paste0(e)]
            this_se <- x[,paste0("se_",e)]
            cil <- this_est - 1.96*this_se; ciu <- this_est + 1.96*this_se
            o[,ct] <- mean(cil < x$truth[1] & ciu > x$truth[1], na.rm = TRUE)
          }
        }
        # add in ltmle with mc standard deviation interval
        sd_ltmle <- sd(x$ltmle, na.rm = TRUE)
        cil <- x$ltmle - 1.96 * sd_ltmle
        ciu <- x$ltmle + 1.96 * sd_ltmle
        o[,ct + 1] <- mean(cil < x$truth[1] & ciu > x$truth[1], na.rm = TRUE)

        colnames(o) <- c("nsim", est, "ltmle_mc")
        o
      })
      ou <- Reduce(rbind, cov)
      ou <- cbind(unique(rslt$n), ou)
      ou
    }
    getCov(out, n = c(500,1000,5000,7000), Q = "SL.hal9001", g = "SL.glm")
    getCov(out, n = c(500,1000,5000,7000), g = "SL.hal9001", Q = "SL.glm")
    getCov(out, n = c(500,1000,5000,7000), Q = "SL.hal9001", g = "SL.hal9001")
    getCov(out, n = c(500,1000,5000,7000), Q = "SL.glm.interaction", g = "SL.glm.interaction")

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
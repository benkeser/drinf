library(devtools)
install_github("benkeser/modifySL")


# load packages
library(modifySL)
library(haltmle.sim)
library(drtmle)

# parameters
ns <- c(100,500,1000,2000)
bigB <- 2000

# directories to save in 
saveDir <- "~/haltmle_sim/out"
scratchDir <- "~/haltmle_sim/scratch"

# # simulation parameters
parm <- expand.grid(seed=1:bigB,
                    n=ns)

i <- 1

set.seed(parm$seed[i])
    dat <- haltmle.sim:::makeRandomData(n=parm$n[i], maxD = 8)
 algo <- c("SL.glm","SL.bayesglm", 
              "SL.earth")

# fit super learner with all algorithms
set.seed(parm$seed[i])
X <- data.frame(dat$A, dat$W)
fullSL <- SuperLearner(Y=dat$Y$Y,X=X,family=gaussian(), SL.library=algo,
               verbose=TRUE,method="method.CC_LS")
# modify SL to drop HAL
dropSL <- modifySL::modifySL(fit = fullSL, Y = dat$Y$Y, library = fullSL$libraryNames[-length(algo)])
                 
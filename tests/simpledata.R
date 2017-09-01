library(roxygen2)
library(devtools)
setwd("~/Dropbox/R")
# document("drinf")
# build("drinf")
install.packages("~/Dropbox/R/drinf",repos = NULL, type = "source")
devtools::install_github("nhejazi/lassi")
library(hal9001)
library(drinf)
library(SuperLearner)
# debug(targetQg)
#debug(drinf.tmle)
#debug(universalStep)
#; debug(targetg1)
set.seed(12346272)
makeData <- function(n = n){
    L0 <- data.frame(x.0 = runif(n,-1,1))
    A0 <- rbinom(n, 1, plogis(L0$x.0^2))
    L1 <- data.frame(x.1 = L0$x.0^2*A0 + runif(n))
    A1 <- rbinom(n, 1, plogis(L0$x.0*L1$x.1))
    L2 <- rnorm(n, L0$x.0^2*A0*A1 + L1$x.1)
    return(list(L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1))
}

getTruth <- function(n = 1e6, abar = c(1,1)){
    L0 <- data.frame(x.0 = runif(n,-1,1))
    L1 <- data.frame(x.1 = L0$x.0^2*abar[1] + runif(n))
    L2 <- rnorm(n, L0$x.0^2*abar[1]*abar[2] + L1$x.1)
    return(mean(L2))
}

dat <- makeData(n = 500)
truth <- getTruth()

test <- drinf.tmle(
    L0 = dat$L0, L1 = dat$L1, L2 = dat$L2, A0 = dat$A0, A1 = dat$A1, 
    abar = c(1,1), 
    SL.Q = "SL.hal9001",
    SL.g = "SL.glm", 
    SL.Qr = "SL.hal9001",
    SL.gr = "SL.hal9001",
#    flucOrd = c("targetg0","targetg1","targetQ2","targetQ1"),
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


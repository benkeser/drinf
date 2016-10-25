L0 <- data.frame(x.0 = rnorm(100))
A0 <- rbinom(100, 1, plogis(L0$x.0))
L1 <- data.frame(x.1 = rnorm(100, L0$x.0 + A0))
A1 <- rbinom(100, 1, plogis(L0$x.0 + L1$x.1))
L2 <- rnorm(100, L0$x.0 + A0 + A1 + L1$x.1)

library(roxygen2)
library(devtools)
setwd("~/Dropbox/")
document("drinf")
build("drinf")
library(drinf)
# debug(targetQg)
# debug(drinf.tmle)
# debug(targetg0); debug(targetg1)
test <- drinf.tmle(
    L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1, 
    abar = c(1,1), SL.Q = c("SL.glm","SL.mean","SL.gam"),
    SL.g = c("SL.glm","SL.mean","SL.gam"), 
    SL.Qr = c("SL.glm","SL.mean","SL.gam"),
    SL.gr = c("SL.glm","SL.mean","SL.gam"),
    #flucOrd = c("targetg0","targetg1","targetQ2","targetQ1"),
    flucOrd = c("targetg0","redReg","targetg1","redReg",
                "targetQ2","redReg","targetQ1"),
    return.models = TRUE, 
    verbose = FALSE,
    maxIter = 100,
    return.ltmle = FALSE,
    allatonce = TRUE,
    tolg = 1e-4,
    tolQ = 1e-4
)




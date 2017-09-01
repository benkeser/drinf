makeData <- function(n, b = 0.25, setA = NULL){
    L0.1 <- runif(n, -4, 4)
    L0.2 <- rbinom(n, 1, 0.5)
    
    if(all(is.null(setA))){
        A0 <- rbinom(n, 1, plogis(1 + b*L0.1 - 2*b*L0.2*L0.1))
    }else{
        A0 <- rep(setA[1],n)
    }
    
    L1.1 <- runif(n, -4, 4)
    L1.2 <- rbinom(n, 1, 0.5)
    
    if(all(is.null(setA))){
        A1 <- rbinom(n, 1, plogis(2 + b*L0.1 - 2*b*L0.2*L0.1+ 
                                      b/2*L1.1 - 2*b/2*L1.2*L1.1))
    }else{
        A1 <- rep(setA[2],n)
    }
    
    L2 <- rbinom(n, 1, plogis(-0.25*A0 - 0.25*A1 + 
                                  b*L0.1 - 2*b*L0.2*L0.1 + 
                                  b*L1.1 - 2*b*L1.2*L1.1))
    
    return(list(
      L0 = data.frame(L0.1 = L0.1, L0.2 = L0.2),
      A0 = A0,
      L1 = data.frame(L1.1 = L1.1, L1.2 = L1.2),
      A1 = A1, 
      L2 = L2
    ))
}

bigData <- makeData(n = 1e6, setA = c(1,1))
mean(bigData$L2)

## function to simulate from truth
truthQ <- function(
    Y, X, newX, family, obsWeights, b=0.25, ...
){
    # for Q1n the number of columns will only be 2
    if(ncol(X) == 2){
        n <- length(Y)
        err <- runif(nrow(newX), min = -n^(-1/4), max = n^(-1/4))
        pred <- 1/(8*b) * log( 
            (1+ exp(-0.5 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + 4*b))/
                (1+ exp(-0.5 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 - 4*b))
        )  
    }else{
        n <- length(Y)
        pred <- plogis(-0.5 + 
                           b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + 
                           b*newX$L1.1 - 2*b*newX$L1.2*newX$L1.1+ 
                           + runif(nrow(newX), min = -n^(-1/4), max = n^(-1/4)))
    }
    return(list(fit = n, pred = pred))
}



## function to simulate from truth
truthG <- function(
    Y, X, newX, family, obsWeights, b=1, ...
){
    n <- length(Y)
    err <- runif(nrow(newX), min = -n^(-1/4), max = n^(-1/4))
    
    # for g0n the number of columns will only be 2
    if(ncol(X) == 2){
        pred <- plogis(2 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + err)
    }else{
        n <- length(Y)
        pred <- plogis(2 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1+ 
                   b/2*newX$L1.1 - 2*b/2*newX$L1.2*newX$L1.1 + err)
    }
    return(list(fit = n, pred = pred))
}

## try to confirm my truth is correct
# makeData2 <- function(n, b = 0.25, setA = NULL){
#     L0.1 <- runif(n, -4, 4)
#     L0.2 <- rbinom(n, 1, 0.5)
#     
#     Q1 <- 1/(8*b) * log( 
#         (1+ exp(-0.5 + b*L0.1 - 2*b*L0.2*L0.1 + 4*b))/
#             (1+ exp(-0.5 + b*L0.1 - 2*b*L0.2*L0.1 - 4*b))
#         )
#     return(Q1)
# }
# bigData2 <- makeData2(n = 1e6, setA = c(1,1))
# mean(bigData2)
# 


library(drinf)
# g wrong, Q right
set.seed(1234)
grbg <- makeData(5000)
test <- drinf.tmle(
    L0 = grbg$L0, L1 = grbg$L1, L2 = grbg$L2, A0 = grbg$A0, A1 = grbg$A1, 
    abar = c(1,1), 
    #SL.Q = c("SL.glm","SL.mean","SL.gam","SL.step.interaction"),
    SL.Q = "truthQ",
    SL.g = c("SL.glm","SL.mean","SL.gam"), 
    SL.Qr = getFromNamespace("SL.earth","SuperLearner"),
    SL.gr = getFromNamespace("SL.earth","SuperLearner"),
    #    flucOrd = c("targetg0","targetg1","targetQ2","targetQ1"),
    flucOrd = c("targetg0","redReg","targetg1","redReg",
                "targetQ2","redReg","targetQ1"),
    return.models = FALSE,
    verbose = FALSE,
    maxIter = 20,
    return.ltmle = TRUE,
    tolg = 1e-4,
    tolQ = 1e-4,
    SL.Q.options = list(family = binomial())
)

# g right, Q wrong
set.seed(1234)
grbg <- makeData(1000)
test <- drinf.tmle(
    L0 = grbg$L0, L1 = grbg$L1, L2 = grbg$L2, A0 = grbg$A0, A1 = grbg$A1, 
    abar = c(1,1), 
    #SL.Q = c("SL.glm","SL.mean","SL.gam","SL.step.interaction"),
    SL.g = "truthG",
    SL.Q = c("SL.glm","SL.mean","SL.gam"), 
    SL.Qr = getFromNamespace("SL.earth","SuperLearner"),
    SL.gr = getFromNamespace("SL.earth","SuperLearner"),
    #    flucOrd = c("targetg0","targetg1","targetQ2","targetQ1"),
    flucOrd = c("targetg0","redReg","targetg1","redReg",
                "targetQ2","redReg","targetQ1"),
    return.models = FALSE,
    verbose = FALSE,
    maxIter = 100,
    return.ltmle = TRUE,
    tolg = 1e-1,
    tolQ = 1e-4,
    tolIF = 1/1000,
    SL.Q.options = list(family = binomial())
)


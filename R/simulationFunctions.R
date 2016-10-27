#' makeData
#' 
#' A function to simulate data for use in the simulation study. The function
#' can be used to generate observed data or data from an intervened SCM; the latter
#' is achieved via the \code{setA} option. The data are simulated as follows: 
#' \itemize{
#' \item L0.1 ~ Uniform(-4,4)
#' \item L0.2 ~ Bernoulli(1/2)
#' \item A0 ~ Bernoulli(expit(b0 + b*L0.1 - 2b*L0.2*L0.1))
#' \item L1.1 ~ Uniform(-4,4)
#' \item L1.2 ~ Bernoulli(1/2)
#' \item A1 ~ Bernoulli(expit(2*b0 + b*L0.1 - 2*b*L0.2*L0.1 + b/2*L1.1 - b*L1.2*L1.1))
#' \item L2 ~ Bernoulli(expit(ba*A0 + ba*A1 + b*L0.1 - 2*b*L0.2*L0.1 + b*L1.1 - 2*b*L1.2*L1.1))
#' }
#' 
#' @param Y The outcome
#' @param X The predictors
#' @param newX The prediction frame
#' @param family Not used, just for compatibility
#' @param obsWeights Not used, just for compatibility
#' @param b Interaction coefficient for L01 and L11
#' @param ba The effect of each treatment on the outcome. See \code{makeData}
#' to see how it is used. 
#' @param r The exponent of the rate at which error goes to zero as function of n.
#' 
#' @return A \code{list} with L0, A0, L1, A1, L2 to input into \code{drinf.tmle}.

makeData <- function(n, b = 0.25, b0 = 1, ba = -0.25, setA = NULL){
    L0.1 <- runif(n, -4, 4)
    L0.2 <- rbinom(n, 1, 0.5)
    
    if(all(is.null(setA))){
        A0 <- rbinom(n, 1, plogis(b0 + b*L0.1 - 2*b*L0.2*L0.1))
    }else{
        A0 <- rep(setA[1],n)
    }
    
    L1.1 <- runif(n, -4, 4)
    L1.2 <- rbinom(n, 1, 0.5)
    
    if(all(is.null(setA))){
        A1 <- rbinom(n, 1, plogis(2*b0 + b*L0.1 - 2*b*L0.2*L0.1+ 
                                      b/2*L1.1 - 2*b/2*L1.2*L1.1))
    }else{
        A1 <- rep(setA[2],n)
    }
    
    L2 <- rbinom(n, 1, plogis(ba*A0 + ba*A1 + 
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

#' truthQ
#' 
#' A function to generate an estimate of g0n and g1n that is equal to 
#' Q10 and Q20 plus an error term. The function is written for compatibility
#' with the drinf.tmle function, i.e., as though it were a SuperLearner wrapper.
#' If \code{newX} is two columns, returns Q1n, if four columns, returns Q2n. 
#' 
#' @param Y The outcome
#' @param X The predictors
#' @param newX The prediction frame
#' @param family Not used, just for compatibility
#' @param obsWeights Not used, just for compatibility
#' @param b Interaction coefficient for L01 and L11
#' @param ba The effect of each treatment on the outcome. See \code{makeData}
#' to see how it is used. 
#' @param r The exponent of the rate at which error goes to zero as function of n.
#' 
#' @return A \code{list} with predictions equal to truth + error

truthQ <- function(
    Y, X, newX, family, obsWeights, ba = -0.25, b=0.25, r=-1/4, ...
){
    # for Q1n the number of columns will only be 2
    if(ncol(X) == 2){
        n <- length(Y)
        err <- runif(nrow(newX), min = -n^(r), max = n^(r))
        pred <- 1/(8*b) * log( 
            (1+ exp(2*ba + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + 4*b))/
                (1+ exp(2*ba + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 - 4*b))
        )  
    }else{
        n <- length(Y)
        pred <- plogis(2*ba + 
                           b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + 
                           b*newX$L1.1 - 2*b*newX$L1.2*newX$L1.1+ 
                           + runif(nrow(newX), min = -n^(r), max = n^(r)))
    }
    return(list(fit = n, pred = pred))
}

#' truthG
#' 
#' A function to generate an estimate of g0n and g1n that is equal to 
#' g00 and g10 plus an error term. The function is written for compatibility
#' with the drinf.tmle function, i.e., as though it were a SuperLearner wrapper.
#' If \code{newX} is two columns, returns g0n, if four columns, returns g1n. 
#' 
#' @param Y The outcome
#' @param X The predictors
#' @param newX The prediction frame
#' @param family Not used, just for compatibility
#' @param obsWeights Not used, just for compatibility
#' @param b Interaction coefficient for L01 and L11
#' @param b0 The intercept
#' @param r The exponent of the rate at which error goes to zero as function of n.
#' 
#' @return A \code{list} with predictions equal to truth + error

truthG <- function(
    Y, X, newX, family, obsWeights, b=1,  b0 = 1, r = -1/4, ...
){
    n <- length(Y)
    err <- runif(nrow(newX), min = -n^(r), max = n^(r))
    
    # for g0n the number of columns will only be 2
    if(ncol(X) == 2){
        pred <- plogis(b0 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1 + err)
    }else{
        n <- length(Y)
        pred <- plogis(2*b0 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1+ 
                           b/2*newX$L1.1 - 2*b/2*newX$L1.2*newX$L1.1 + err)
    }
    return(list(fit = n, pred = pred))
}

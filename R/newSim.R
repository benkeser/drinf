#' makeData
#' 
#' A function to simulate data for use in the simulation study. The function
#' can be used to generate observed data or data from an intervened SCM; the latter
#' is achieved via the \code{setA} option. 
#' 
#' @export
#' 
#' @return A \code{list} with L0, A0, L1, A1, L2 to input into \code{drinf.tmle}.

makeData <- function(n, b = 0.5, b0 = 1, ba = -0.25, rho=0.05, setA = NULL){
    L0.1 <- runif(n, -4, 4)
    L0.2 <- rbinom(n, 1, 0.5)
    
    c00 <- 2*b0 + b*L0.1 - 2*b*L0.2*L0.1

    if(all(is.null(setA))){
        A0 <- rbinom(n, 1, plogis(c00))
    }else{
        A0 <- rep(setA[1],n)
    }
    
    f0 <- L0.1/2 + 0.25*A0
    L1.1 <-  f0 + runif(n,-rho,rho)

    p0 <- 0.1 + 0.8*L0.2 + 0.05*A0
    L1.2 <- rbinom(n, 1, p0)
    
    c1 <- b*L1.1 - 2*b*L1.2*L1.1
    if(all(is.null(setA))){
        A1 <- rbinom(n, 1, plogis(c00 + c1 + 0.05*A0))
    }else{
        A1 <- rep(setA[2],n)
    }
    
    c0 <- ba*A0 + ba*A1 + c00
    L2 <- c0 + c1 + runif(n,-1,1)
   
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
#' @export 
#' 
#' @return A \code{list} with predictions equal to truth + error

truthQ <- function(
    Y, X, newX, family, obsWeights, b0 = 1, ba = -0.25, b=0.5, r=-1/4, cons = 2, ...
){  
    c00 <- 2*b0 + b*newX$L0.1 - 2*b*newX$L0.2*newX$L0.1
    c0 <- ba*1 + ba*1 + c00
    f0 <- newX$L0.1/2 + 0.25
    p0 <- 0.1 + 0.8*newX$L0.2 + 0.05*1

    # for Q1n the number of columns will only be 2
    if(ncol(X) == 2){
        n <- length(Y)
        err <- cons*runif(nrow(newX), min = -n^(r), max = n^(r))
        pred <- (1-p0)*(c0 + b*f0) + p0*(c0 - b*f0) + err
    }else{
        n <- length(Y)
        c1 <- b*newX$L1.1 - 2*b*newX$L1.2*newX$L1.1
        err <- cons*runif(nrow(newX), min = -n^(r), max = n^(r))
        pred <- c0 + c1 + err
    }
    return(list(fit = n, pred = pred))
}




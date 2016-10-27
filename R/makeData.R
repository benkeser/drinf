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
#' @export
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

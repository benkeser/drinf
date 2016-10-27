
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

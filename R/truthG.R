
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
#' @export 
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

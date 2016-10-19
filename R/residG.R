#' residG
#' 
#' This function computes the residuals that are used as outcomes in the 
#' reduced dimension regressions for g. 
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' @param ... Other arguments (not currently used)
#' 
#' @return A list with elements rg1 and rg2 that are the outcomes in the
#' reduced dimension regressions on Q1n and Q2n, respectively. 

residG <- function(
    A0, A1, g0n, g1n, abar, ...
){
    rg0 <- (as.numeric(A0==abar[1])-g0n)/g0n
    rg1 <- as.numeric(A0==abar[1])/g0n * (as.numeric(A1==abar[2]) - g1n)/g1n
    return(list(
        rg0 = rg0, rg1 = rg1
    ))
}
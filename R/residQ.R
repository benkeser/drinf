#' residQ
#' 
#' This function computes the residuals that are used as outcomes in the 
#' reduced dimension regressions for Q. 
#' 
#' @param L2 A \code{vector} outcome of interest
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' 
#' @return A list with elements rQ1 and rQ2 that are the outcomes in the
#' reduced dimension regressions on g0n and g1n, respectively. 

residQ <- function(
    L2, A0, A1, Q2n, Q1n, g0n, g1n, abar, ...
){
    rQ1 <- as.numeric(A0==abar[1] & A1 == abar[2]) / (g0n*g1n) * (L2 - Q2n) + 
        as.numeric(A0==abar[1])/g0n * (Q2n - Q1n)
    rQ2 <- as.numeric(A0==abar[1] & A1 == abar[2])  / (g0n*g1n) * (L2 - Q2n)
    return(list(
        rQ1 = rQ1, rQ2 = rQ2
    ))
}
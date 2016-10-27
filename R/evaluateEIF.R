
#' evaluateEIF
#' 
#' Function that returns estimated efficient influence function contributions.
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest. 
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' @param ... Additional arguments (not currently used)
#' 
#' @return A list with named entries corresponding to the different pieces of the 
#' EIF at the nuisance parameter estimates evaluated at the observations. 

evaluateEIF <- function(
    A0, A1, L2, Q2n, Q1n, g1n, g0n, abar, ...
){
    # usual pieces of the EIF
    Dstar2 <- as.numeric(A0==abar[1] & A1==abar[2]) / (g0n * g1n) * (L2 - Q2n)
    Dstar1 <- as.numeric(A0==abar[1]) / g0n * (Q2n - Q1n)
    Dstar0 <- Q1n - mean(Q1n)
    
    #--------------
    # return 
    #--------------
    return(list(
        Dstar0 = Dstar0, Dstar1 = Dstar1, Dstar2 = Dstar2
    ))
}

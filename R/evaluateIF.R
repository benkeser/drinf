#' evaluateIF
#' 
#' Function that returns estimated influence function contributions.
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest. 
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param Q2nr.obsa A \code{vector} of estimates of Q_{2,0}^r at the observed values of A0.
#' @param Q1nr1 A \code{vector} of estimates of Q_{1,0,1}^r
#' @param Q1nr2 A \code{vector} of estimates of Q_{1,0,2}^r
#' @param g0nr A \code{vector} of estimates of g_{0,0}^r
#' @param g1nr A \code{vector} of estimates of g_{0,0}^r
#' @param h0nr A \code{vector} of estimates of h_{0,0}^r
#' @param h1nr A \code{vector} of estimates of h_{1,0}^r
#' @param hbarnr A \code{vector} of estimates of h^r, the iterated reduced
#' dimension regression.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.  
#' @return A list with named entries corresponding to the different pieces of the 
#' EIF at the nuisance parameter estimates evaluated at the observations. 

evaluateIF <- function(
    A0, A1, L2, Q2n, Q1n, g1n, g0n, Q2nr.obsa, Q1nr1, Q1nr2,
    g0nr, g1nr, h0nr, h1nr, hbarnr, abar
){
    # usual pieces of the EIF
    Dstar2 <- as.numeric(A0==abar[1] & A1==abar[2]) / (g0n * g1n) * (L2 - Q2n)
    Dstar1 <- as.numeric(A0==abar[1]) / g0n * (Q2n - Q1n)
    Dstar0 <- Q1n - mean(Q1n)
    
    #--------------------------------------------------
    # pieces resulting from misspecification of Q
    # targeting of g should knock these equations out
    #--------------------------------------------------
    # targeting of g1 to get rid of this one
    Dg1.Q2 <- Q2nr.obsa / g1n^2 * (as.numeric(A1 == abar[2]) - g1n)
    
    # targeting of g0 to get rid of this one
    Dg0.Q1 <- (Q1nr1 + Q1nr2) / g0n^2 * (as.numeric(A0 == abar[1]) - g0n)
    
    #--------------------------------------------------
    # pieces resulting from misspecification of g
    # targeting of Q should knock these equations out
    #--------------------------------------------------
    # targeting of Q2 to get rid of these two
    # this one from misspecification of g1
    DQ2.g1 <- as.numeric(A0==abar[1] & A1==abar[2]) / g1nr * h1nr * (L2 - Q2n)
    # this one from misspecification of g0 -- involves iterated reduction
    DQ2.g0 <- as.numeric(A0==abar[1] & A1==abar[2]) / g1nr * hbarnr * (L2 - Q2n)
    
    # targeting of Q1 to get rid of this one
    DQ1.g0 <- as.numeric(A0==abar[1])/g0nr * h0nr * (Q2n - Q1n)
    
    #--------------
    # return 
    #--------------
    return(list(
        Dstar0 = Dstar0, Dstar1 = Dstar1, Dstar2 = Dstar2, 
        Dg1.Q2 = Dg1.Q2, Dg0.Q1 = Dg0.Q1, DQ2.g1 = DQ2.g1, DQ2.g0 = DQ2.g0, 
        DQ1.g0 = DQ1.g0
    ))
}

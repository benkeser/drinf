#' evaluateRisk
#' 
#' Function that returns empirical risk
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
#' @return A \code{numeric} value for the empirical risk

evaluateRisk <- function(
    A0, A1, L2, Q2n, Q1n, g1n, g0n, abar, tolg
){
    # scale L2, Q2n, Q1n to be in (0,1)
    L2.min <- min(L2); L2.max <- max(L2)
    # scale Q2n,Q1n
    Q2ns <- (Q2n - L2.min)/(L2.max - L2.min)
    Q1ns <- (Q1n - L2.min)/(L2.max - L2.min)
    L2s <- (L2 - L2.min)/(L2.max - L2.min)
    # scale g0n g1n
    A0s <- (as.numeric(A0 == abar[1]) - tolg)/(1 - 2*tolg)
    A1s <- (as.numeric(A1 == abar[2]) - tolg)/(1 - 2*tolg)
    g1ns <- (g1n - tolg)/(1 - 2*tolg)
    g0ns <- (g0n - tolg)/(1 - 2*tolg)
    g1ns[g1ns == 0] <- .Machine$double.neg.eps
    g1ns[g1ns == 1] <- 1 - .Machine$double.neg.eps
    g0ns[g0ns == 0] <- .Machine$double.neg.eps
    g0ns[g0ns == 1] <- 1 - .Machine$double.neg.eps

    # Loss for Q2
    LQ2 <- mean(as.numeric(A0==abar[1] & A1==abar[2])*
        (-L2s * log(Q2ns) - (1-L2s) * log(1-Q2ns)))
    # Loss for Q1
    LQ1 <- mean(as.numeric(A0==abar[1])*
        (-Q2ns * log(Q1ns) - (1-Q2ns) * log(1-Q1ns)))
    # Loss for g1
    Lg1 <- mean(as.numeric(A0==abar[1])*
        (-as.numeric(A1s) * log(g1ns) - 
        (1 - as.numeric(A1s)) * log(1-g1ns)))
    # Loss for g0
    Lg0 <- mean(as.numeric(
        -as.numeric(A0s) * log(g0ns) - 
        (1 - as.numeric(A0s)) * log(1-g0ns)))
    
    # return risk
    return(list(sum = LQ2 + LQ1 + Lg1 + Lg0,
                Q2 = LQ2, Q1 = LQ1, Lg1 = Lg1, Lg0 = Lg0))
}

#' redReg
#' 
#' Function that computes the reduced dimension regressions needed for the 
#' extra targeting steps of the tmle.
#' 
#' @return A list with Qnr and gnr objects. See \code{estimateQr} and \code{estimategr}
#' for details. 


redReg <- function(
    L2, A0, A1, Q2n, Q1n, g0n, g1n, abar, verbose, return.models,
    SL.Qr, SL.gr, ...
){
    #---------------------------------------
    # residual outcomes for Qr regressions
    #---------------------------------------
    rQ <- residQ(
        L2 = L2, A0 = A0, A1 = A1, Q2n = Q2n, Q1n = Q1n, 
        g0n = g0n, g1n = g1n, abar = abar, ...
    )
    
    #---------------------------
    # estimate Qr regressions
    #---------------------------
    Qnr <- estimateQr(
        rQ1 = rQ$rQ1, rQ2 = rQ$rQ2, g0n = g0n, g1n = g1n, 
        A0 = A0, A1 = A1, SL.Qr = SL.Qr, abar = abar, 
        return.models = return.models, verbose = verbose, ...
    )
    
    #---------------------------------------
    # residual outcomes for gr regressions
    #---------------------------------------
    rg <- residG(
        A0 = A0, A1 = A1, g0n = g0n, g1n = g1n, abar = abar, ...
    )
    
    #---------------------------
    # estimate gr regressions
    #---------------------------
    gnr <- estimategr(
        rg0 = rg$rg0, rg1 = rg$rg1, g0n = g0n, 
        g1n = g1n, A0 = A0, A1 = A1, Q2n = Q2n, Q1n = Q1n, 
        SL.gr = SL.gr, abar = abar, return.models = return.models, 
        tolg = tolg, verbose = verbose, ...
    )
    
    return(list(
        Qnr = Qnr, gnr = gnr
    ))
}
#' redReg
#' 
#' Function that computes the reduced dimension regressions needed for the 
#' extra targeting steps of the tmle.
#' 
#' @return A list with Qnr and gnr objects. See \code{estimateQr} and \code{estimategr}
#' for details. 


redReg <- function(
    L2, A0, A1, Qn, gn, folds, validFold, abar, verbose, return.models,
    SL.Qr, SL.gr, tolg, ...
){  
    # combine full Qn and gn across folds
    full_Qn <- do.call(Map, c(c, Qn))
    full_gn <- do.call(Map, c(c, gn))

    if(all(folds == 1)){
        train_L2 <- valid_L2 <- L2
        train_A0 <- valid_A0 <- A0
        train_A1 <- valid_A1 <- A1
        train_Q1n <- valid_Q1n <- full_Qn$Q1n
        train_Q2n <- valid_Q2n <- full_Qn$Q2n
        train_g0n <- valid_g0n <- full_gn$g0n
        train_g1n <- valid_g1n <- full_gn$g1n
    }else{
        # training data
        train_L2 <- L2[folds != validFold]
        train_A0 <- A0[folds != validFold]
        train_A1 <- A1[folds != validFold]
        train_Q1n <- full_Qn$Q1n[folds != validFold]
        train_Q2n <- full_Qn$Q2n[folds != validFold]
        train_g0n <- full_gn$g0n[folds != validFold]
        train_g1n <- full_gn$g1n[folds != validFold]
    
        # validation data
        valid_L2 <- L2[folds == validFold]
        valid_A0 <- A0[folds == validFold]
        valid_A1 <- A1[folds == validFold]
        valid_Q1n <- full_Qn$Q1n[folds == validFold]
        valid_Q2n <- full_Qn$Q2n[folds == validFold]
        valid_g0n <- full_gn$g0n[folds == validFold]
        valid_g1n <- full_gn$g1n[folds == validFold]
    }
    #---------------------------------------
    # residual outcomes for Qr regressions
    #---------------------------------------
    rQ <- residQ(
        L2 = train_L2, A0 = train_A0, A1 = train_A1, Q2n = train_Q2n, 
        Q1n = train_Q1n, g0n = train_g0n, g1n = train_g1n, abar = abar, ...
    )
    
    #---------------------------
    # estimate Qr regressions
    #---------------------------
    Qnr <- estimateQr(
        rQ1_1 = rQ$rQ1_1, rQ1_2 = rQ$rQ1_2, rQ2 = rQ$rQ2, 
        g0n = full_gn$g0n, g1n = full_gn$g1n, folds = folds,
        validFold = validFold, A0 = A0, A1 = A1, SL.Qr = SL.Qr, abar = abar, 
        return.models = return.models, verbose = verbose, ...
    )
    
    #---------------------------------------
    # residual outcomes for gr regressions
    #---------------------------------------
    rg <- residG(
        A0 = train_A0, A1 = train_A1, g0n = train_g0n, g1n = train_g1n, 
        abar = abar, ...
    )
    
    #---------------------------
    # estimate gr regressions
    #---------------------------
    gnr <- estimategr(
        rg0 = rg$rg0, rg1 = rg$rg1, g0n = full_gn$g0n, 
        g1n = full_gn$g1n, folds = folds, validFold = validFold,
        A0 = A0, A1 = A1, Q2n = full_Qn$Q2n, 
        Q1n = full_Qn$Q1n, SL.gr = SL.gr, abar = abar, 
        return.models = return.models, tolg = tolg, verbose = verbose, ...
    )
    
    return(list(
        Qnr = Qnr, gnr = gnr
    ))
}
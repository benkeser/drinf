#' estimateQr
#' 
#' A function used to estimate the reduced dimension regressions for Q
#' 
#' @param rQ1 The "residual" for the first reduced dimension regression (on g0n), 
#' equal to sum of the EIF at time 1 and 2. 
#' @param rQ2 The "residual" for the second reduced dimension regression (on g1n),
#' equal to the EIF at time 2. 
#' @param g1n A \code{vector} of estimates of g_{1,0}.
#' @param g0n A \code{vector} of estimates of g_{0,0}.
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param SL.Qr A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the reduced-dimension regression to protect against misspecification of the
#' outcome regressions.  See \code{SuperLearner} package for details.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' @param returnModels  A \code{boolean} indicating whether the models for Qr0 should be 
#' returned with the output. 
#'  
#' @return A list with elements Q2nr, Q1nr, Q2mod, and Q1mod, corresponding to the
#' predicted values of the reduced dimension regressions and the models used to fit
#' them, respectively. 


estimateQr <- function(
    rQ1, rQ2, g0n, g1n, A0, A1, SL.Qr, abar, returnModels, ...
){
    multiAlgos <- (length(SL.Qr) > 1 | is.list(SL.Qr))
    
    #--------
    # Q2nr 
    #--------
    Q2rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.Qr),args=list(
        Y=rQ2[A0==abar[1]], # fit only using A0==abar[1] obs. 
        X=data.frame(g1n = g1n[A0==abar[1]]),
        obsWeights = rep(1, length(g1n)),
        family = gaussian(),
        SL.library=SL.Qr,
        verbose=verbose))
    # empty vector of predictions
    Q2nr <- rep(NA, length(A0))
    # zeros for these obs. 
    Q2nr[A0 != abar[1]] <- 0
    if(multiAlgos){
        # Super Learner predictions for A0==abar[1] obs. 
        Q2nr[A0 == abar[1]] <- Q2mod$SL.predict
    }else{
        Qrn[A0 == abar[1]] <- Q2mod$pred
    }
    #------
    # Q1nr 
    #------
    Q1rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.Qr),args=list(
        Y=rQ1, # fit using all obs. 
        X=data.frame(g0n = g0n),
        SL.library=SL.Qr,
        obsWeights = rep(1, length(g0n)),
        family = gaussian(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions for A0==abar[1] obs. 
        Q1nr[A0 == abar[1]] <- Q2mod$SL.predict
    }else{
        Q1nr[A0 == abar[1]] <- Q2mod$pred
    }

    #--------
    # return
    #--------
    out <- list(Q2nr = Q2nr, Q1nr = Q1nr,
                Q2rmod = NULL, Q1rmod = NULL)
    
    if(returnModels){
        Q2rmod$call <- Q1rmod$call <- NULL
        out$Q2rmod <- Q2rmod
        out$Q1rmod <- Q1rmod
    }
    
    return(out)
}
#' estimategr
#' 
#' A function used to estimate the reduced dimension regressions for g
#' 
#' @param rg0 The "residual" for the first reduced dimension regression (on Q1n).
#' @param rg1 The "residual" for the second reduced dimension regression (on Q2n).
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param SL.gr A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the reduced-dimension regression to protect against misspecification of the
#' treatment regressions.  See \code{SuperLearner} package for details.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' @param returnModels  A \code{boolean} indicating whether the models for Qr0 should be 
#' returned with the output. 
#'  
#' @return A list with elements g0nr, g1nr, h0nr, h1nr, and hbarnr, corresponding to the
#' predicted values of the reduced dimension regressions. Also included in output are the
#' models used to obtain these predicted values (set to \code{NULL} if \code{returnModels = FALSE})


# TO DO : All this code for one dimensional regressions
# could be streamlined into a single function that would take as input
# the outcome variable and the regressor and return SL predictions. 

estimategr <- function(
    rg0, rg1, g0n, g1n, A0, A1, Q2n, Q1n, SL.gr, abar, returnModels, ...
){
    multiAlgos <- (length(SL.gr) > 1 | is.list(SL.gr))
    #-------------
    # g0nr 
    # A0 ~ Q1n 
    #-------------
    g0rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.gr),args=list(
        Y=as.numeric(A0==abar[1]),  
        X=data.frame(Q1n = Q1n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(Q1n)),
        family = binomial(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions 
        g0nr <- g0rmod$SL.predict
    }else{
        g0nr <- g0rmod$pred
    }

    #----------------
    # g1nr
    # A0*A1 ~ Q2n
    #----------------
    g1rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.gr),args=list(
        Y=as.numeric(A0==abar[1] & A1==abar[2]), 
        X=data.frame(Q2n = Q2n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(Q2n)),
        family = binomial(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions 
        g1nr <- g1rmod$SL.predict
    }else{
        g1nr <- g1rmod$pred
    }
    
    #-------------------------------
    # h0nr
    # rg0 [= (A0 - g0n)/g0n] ~ Q1n
    #-------------------------------
    h0rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.gr),args=list(
        Y=rg0,
        X=data.frame(Q1n = Q1n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(Q1n)),
        family = gaussian(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions 
        h0nr<- h0rmod$SL.predict
    }else{
        h0nr <- h0rmod$pred
    }
    
    #---------------------------------------
    # h1rn
    # rg1 [= A0/g0n * (A1 - g1n)/g1n] ~ Q2n
    #---------------------------------------
    h1rmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.gr),args=list(
        Y=rg1,
        X=data.frame(Q2n = Q2n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(Q2n)),
        family = gaussian(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions 
        h1nr <- h1rmod$SL.predict
    }else{
        h1nr <- h1rmod$pred
    }
    
    #--------------------------------------
    # hbarnr
    # A0/g0nr * h0nr ~ Q2n 
    #--------------------------------------
    hbarrmod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.gr),args=list(
        Y=A0/g0nr * h0nr,
        X=data.frame(Q2n = Q2n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(Q2n)),
        family = gaussian(),
        verbose=verbose))
    if(multiAlgos){
        # Super Learner predictions 
        hbarnr <- hbarrmod$SL.predict
    }else{
        hbarnr <- hbarrmod$pred
    }
    
    #--------
    # return
    #--------
    out <- list(g0nr = g0nr, g1nr = g1nr, h0nr = h0nr, h1nr = h1nr,
                hbarnr = hbarnr, g0rmod = NULL, g1rmod = NULL, h0rmod = NULL, 
                h1rmod = NULL, hbarrmod = NULL) 
    
    if(returnModels){
        out$g0rmod <- g0rmod
        out$g1rmod <- g1rmod
        out$h0rmod <- h0rmod
        out$h1rmod <- h1rmod
        hbarrmod <- hbarrmod
    }
    
    return(out)
}
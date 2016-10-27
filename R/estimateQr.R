#' estimateQr
#' 
#' A function used to estimate the reduced dimension regressions for Q. The regression 
#' can be computed using a user specified function, passed through \code{SL.Qr} or using
#' \code{SuperLearner} when \code{length(SL.Qr) == 1} or \code{is.list(SL.Qr)}. There is 
#' an error proofing of the \code{SuperLearner} implementation that deals with situations where
#' the \code{NNLS} procedure in the Super Learner ensemble fails and so the function returns 
#' zero weights for every coefficient. In this case, the code will default to using the discrete
#' Super Learner; that is, the learner with lowest CV-risk. 
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
#' @param return.models  A \code{boolean} indicating whether the models for Qr0 should be 
#' returned with the output. 
#' @param verbose A \code{boolean} indicating whether messages should be printed to indicate progress.
#'  
#' @importFrom SuperLearner SuperLearner
#' 
#' @return A list with elements Q2nr.obsa, Q2rn.seta, Q1nr, Q2mod, 
#' and Q1mod. Q2nr.obsa corresponds to the predicted value of the reduced dimension
#' regression where A0 is its observed value, while Q2nr.seta is the reduced dimension 
#' regression where A0 is set to abar[1]. 


estimateQr <- function(
    rQ1, rQ2, g0n, g1n, A0, A1, SL.Qr, abar, return.models, verbose, ...
){
    multiAlgos <- (length(SL.Qr) > 1 | is.list(SL.Qr))
    #--------
    # Q2nr 
    #--------
    Q2rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Qr),args=list(
        Y=rQ2[A0==abar[1]], # fit only using A0==abar[1] obs. 
        X=data.frame(g1n = g1n[A0==abar[1]]),
        newX = data.frame(g1n = g1n[A0==abar[1]]),
        obsWeights = rep(1, length(g1n[A0==abar[1]])),
        family = gaussian(),
        SL.library=SL.Qr,
        verbose=verbose))
    # empty vector of predictions
    Q2nr.obsa <- rep(NA, length(A0))
    # zeros for these obs. 
    Q2nr.obsa[A0 != abar[1]] <- 0
    if(multiAlgos){
        weightfail <- all(Q2rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q2nr.obsa[A0 == abar[1]] <- Q2rmod$SL.predict
            # get prediction also setting A0 = abar[1] that is needed later
            # for targeting step
            Q2nr.seta <- predict(Q2rmod, newdata = data.frame(g1n = g1n))[[1]]
        }else{
            # find dsl
            dslcol <- which(Q2rmod$cvRisk == min(Q2rmod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q2nr.obsa[A0 == abar[1]] <- Q2rmod$library.predict[,dslcol]
            Q2nr.seta <- predict(Q2rmod, newdata = data.frame(g1n = g1n))[[2]][,dslcol]
        }
    }else{
        Q2nr.obsa[A0 == abar[1]] <- Q2rmod$pred
        Q2nr.seta <- predict(Q2rmod$fit, newdata = data.frame(g1n = g1n))
    }
    #------
    # Q1nr 
    #------
    Q1rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Qr),args=list(
        Y=rQ1, # fit using all obs. 
        X=data.frame(g0n = g0n),
        newX=data.frame(g0n = g0n),
        SL.library=SL.Qr,
        obsWeights = rep(1, length(g0n)),
        family = gaussian(),
        verbose=verbose))
    if(multiAlgos){
        weightfail <- all(Q1rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q1nr <- Q1rmod$SL.predict
        }else{
            # find dsl
            dslcol <- which(Q1rmod$cvRisk == min(Q1rmod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q1nr <- Q1rmod$library.predict[,dslcol]
        }
    }else{
        Q1nr <- Q1rmod$pred
    }

    #--------
    # return
    #--------
    out <- list(Q2nr.obsa = Q2nr.obsa, Q2nr.seta = Q2nr.seta, Q1nr = Q1nr,
                Q2rmod = NULL, Q1rmod = NULL)
    
    if(return.models){
        Q2rmod$call <- Q1rmod$call <- NULL
        out$Q2rmod <- Q2rmod
        out$Q1rmod <- Q1rmod
    }
    
    return(out)
}
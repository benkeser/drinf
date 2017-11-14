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
#' @param rQ1_1 The "residual" for the first of two reduced-dimension regressions (on g0n). 
#' @param rQ1_2 The "residual" for the second of two reduced-dimension regressions (on g0n). 
#' @param rQ2 The "residual" for the second reduced dimension regression (on g1n),
#' equal to the EIF at time 2. 
#' @param g1n A \code{vector} of estimates of g_{1,0}.
#' @param g0n A \code{vector} of estimates of g_{0,0}.
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param folds Vector of cross-validation folds
#' @param validFold Which fold is the validation fold
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
    rQ1_1, rQ1_2, rQ2, g0n, g1n, A0, A1, SL.Qr, folds, validFold, abar, return.models, verbose, ...
){
    if(all(folds == 1)){
        train_g0n <- valid_g0n <- g0n
        train_g1n <- valid_g1n <- g1n
        train_A0 <- valid_A0 <- A0
        train_A1 <- valid_A1 <- A1
    }else{
        # training data
        train_g0n <- g0n[folds != validFold]
        train_g1n <- g1n[folds != validFold]
        train_A0 <- A0[folds != validFold]
        train_A1 <- A1[folds != validFold]

        # validation data
        valid_g0n <- g0n[folds == validFold]
        valid_g1n <- g1n[folds == validFold]
        valid_A0 <- A0[folds == validFold]
        valid_A1 <- A1[folds == validFold]
    }
    multiAlgos <- (length(SL.Qr) > 1 | is.list(SL.Qr))

    #--------
    # Q2nr(A0, L0, L1)
    # = A0*A1 / g0n * (L2 - Q2n) ~ g1n | A0 = 1
    # = 0 | A0 = 0
    #--------
    Q2rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Qr),args=list(
        Y=rQ2[train_A0==abar[1]], # fit only using A0==abar[1] obs. 
        X=data.frame(g1n = train_g1n[train_A0==abar[1]]),
        newX = data.frame(g1n = valid_g1n[valid_A0==abar[1]]),
        obsWeights = rep(1, length(train_g1n[train_A0==abar[1]])),
        family = gaussian(),
        SL.library=SL.Qr,
        method = "method.CC_LS",
        verbose=verbose))
    # empty vector of predictions
    Q2nr.obsa <- rep(NA, length(valid_A0))
    # zeros for these obs. 
    Q2nr.obsa[valid_A0 != abar[1]] <- 0
    if(multiAlgos){
        weightfail <- all(Q2rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q2nr.obsa[valid_A0 == abar[1]] <- Q2rmod$SL.predict
            # get prediction also setting A0 = abar[1] that is needed later
            # for targeting step
            Q2nr.seta <- predict(Q2rmod, newdata = data.frame(g1n = valid_g1n))[[1]]
        }else{
            # find dsl
            dslcol <- which(Q2rmod$cvRisk == min(Q2rmod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q2nr.obsa[valid_A0 == abar[1]] <- Q2rmod$library.predict[,dslcol]
            Q2nr.seta <- predict(Q2rmod, newdata = data.frame(g1n = valid_g1n))[[2]][,dslcol]
        }
    }else{
        Q2nr.obsa[valid_A0 == abar[1]] <- Q2rmod$pred
        Q2nr.seta <- predict(Q2rmod$fit, newdata = data.frame(g1n = valid_g1n))
    }
    #------
    # Q1nr1
    # = A0 * (Q2n - Q1n) ~ g0n 
    #------
    Q1r1mod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Qr),args=list(
        Y=rQ1_1, # fit using all obs. 
        X=data.frame(g0n = train_g0n),
        newX=data.frame(g0n = valid_g0n),
        SL.library=SL.Qr,
        obsWeights = rep(1, length(train_g0n)),
        family = gaussian(),
        method = "method.CC_LS",
        verbose = verbose))
    if(multiAlgos){
        weightfail <- all(Q1r1mod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q1nr1 <- Q1r1mod$SL.predict
        }else{
            # find dsl
            dslcol <- which(Q1r1mod$cvRisk == min(Q1r1mod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q1nr1 <- Q1r1mod$library.predict[,dslcol]
        }
    }else{
        Q1nr1 <- Q1r1mod$pred
    }

    #------
    # Q1nr2
    # = A0*A1 / g1n * (L2 - Q2n) ~ g0n 
    #------
    Q1r2mod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Qr),args=list(
        Y = rQ1_2, # fit using all obs. 
        X = data.frame(g0n = train_g0n),
        newX = data.frame(g0n = valid_g0n),
        SL.library = SL.Qr,
        obsWeights = rep(1, length(train_g0n)),
        family = gaussian(),
        method = "method.CC_LS",
        verbose=verbose))
    if(multiAlgos){
        weightfail <- all(Q1r2mod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q1nr2 <- Q1r2mod$SL.predict
        }else{
            # find dsl
            dslcol <- which(Q1r2mod$cvRisk == min(Q1r2mod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q1nr2 <- Q1r2mod$library.predict[,dslcol]
        }
    }else{
        Q1nr2 <- Q1r2mod$pred
    }
    #--------
    # return
    #--------
    out <- list(Q2nr.obsa = Q2nr.obsa, Q2nr.seta = Q2nr.seta, Q1nr1 = Q1nr1,
                Q1nr2 = Q1nr2, Q2rmod = NULL, Q1r1mod = NULL, Q1r2mod = NULL)
    
    if(return.models){
        Q2rmod$call <- Q1rmod$call <- NULL
        out$Q2rmod <- Q2rmod
        out$Q1r1mod <- Q1r1mod
        out$Q1r2mod <- Q1r2mod
    }
    
    return(out)
}
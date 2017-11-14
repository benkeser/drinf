#' estimategr
#' 
#' A function used to estimate the reduced dimension regressions for g. The regression 
#' can be computed using a user specified function, passed through \code{SL.gr} or using
#' \code{SuperLearner} when \code{length(SL.gr) == 1} or \code{is.list(SL.gr)}. There is 
#' an error proofing of the \code{SuperLearner} implementation that deals with situations where
#' the \code{NNLS} procedure in the Super Learner ensemble fails and so the function returns 
#' zero weights for every coefficient. In this case, the code will default to using the discrete
#' Super Learner; that is, the learner with lowest CV-risk. 
#' 
#' @param rg0 The "residual" for the first reduced dimension regression (on Q1n).
#' @param rg1 The "residual" for the second reduced dimension regression (on Q2n).
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param folds Vector of cross-validation folds
#' @param validFold Which fold is the validation fold
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param SL.gr A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the reduced-dimension regression to protect against misspecification of the
#' treatment regressions.  See \code{SuperLearner} package for details.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest.
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param return.models  A \code{boolean} indicating whether the models for Qr0 should be 
#' returned with the output. 
#'  
#' @importFrom SuperLearner SuperLearner

#' @return A list with elements g0nr, g1nr, h0nr, h1nr, and hbarnr, corresponding to the
#' predicted values of the reduced dimension regressions. Also included in output are the
#' models used to obtain these predicted values (set to \code{NULL} if \code{return.models = FALSE})


# TO DO : All this code for one dimensional regressions
# could be streamlined into a single function that would take as input
# the outcome variable and the regressor and return SL predictions. 

estimategr <- function(
    rg0, rg1, g0n, g1n, A0, A1, folds, validFold, Q2n, Q1n, 
    SL.gr, abar, return.models, tolg, verbose, ...
){
    all1 <- all(folds ==1)
    if(all1){
        train_g0n <- valid_g0n <- g0n
        train_g1n <- valid_g1n <- g1n
        train_Q1n <- valid_Q1n <- Q1n
        train_Q2n <- valid_Q2n <- Q2n
        train_A0 <- valid_A0 <- A0
        train_A1 <- valid_A1 <- A1
    }else{
        # training data
        train_g0n <- g0n[folds != validFold]
        train_g1n <- g1n[folds != validFold]    
        train_Q1n <- Q1n[folds != validFold]
        train_Q2n <- Q2n[folds != validFold]
        train_A0 <- A0[folds != validFold]
        train_A1 <- A1[folds != validFold]

        # validation data
        valid_g0n <- g0n[folds == validFold]
        valid_g1n <- g1n[folds == validFold]
        valid_A0 <- A0[folds == validFold]
        valid_A1 <- A1[folds == validFold]
        valid_Q1n <- Q1n[folds == validFold]
        valid_Q2n <- Q2n[folds == validFold]
    }
    multiAlgos <- (length(SL.gr) > 1 | is.list(SL.gr))
    #-------------
    # g0nr 
    # A0 ~ Q1n 
    #-------------
    g0rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.gr),args=list(
        Y=as.numeric(train_A0==abar[1]),  
        X=data.frame(Q1n = train_Q1n),
        newX=data.frame(Q1n = Q1n), # need predictions on full data 
        SL.library=SL.gr,
        obsWeights = rep(1, length(train_Q1n)),
        family = binomial(),
        method = "method.CC_nloglik",
        verbose=verbose))
    if(multiAlgos){
        weightfail <- all(g0rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions 
            g0nr <- g0rmod$SL.predict[folds == validFold]
            if(!all1){
                train_g0nr <- g0rmod$SL.predict[folds != validFold]                
            }else{
                train_g0nr <- g0nr
            }
        }else{
            dslcol <- which(g0rmod$cvRisk == min(g0rmod$cvRisk, na.rm = TRUE))
            g0nr <- g0rmod$library.predict[folds == validFold, dslcol]
            if(!all1){
                train_g0nr <- g0rmod$library.predict[folds != validFold, dslcol]
            }else{
                train_g0nr <- g0nr
            }
        }
    }else{
        g0nr <- g0rmod$pred[folds == validFold]
        if(!all1){
            train_g0nr <- g0rmod$pred[folds != validFold]
        }else{
            train_g0nr <- g0nr
        }
    }
    g0nr[g0nr < tolg] <- tolg
    train_g0nr[train_g0nr < tolg] <- tolg

    #----------------
    # g1nr
    # A0*A1 ~ Q2n
    #----------------
    g1rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.gr),args=list(
        Y=as.numeric(train_A0==abar[1] & train_A1==abar[2]), 
        X=data.frame(Q2n = train_Q2n),
        newX = data.frame(Q2n = valid_Q2n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(train_Q2n)),
        family = binomial(),
        method = "method.CC_nloglik", 
        verbose=verbose))
    if(multiAlgos){
        weightfail <- all(g1rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions 
            g1nr <- g1rmod$SL.predict
        }else{
            dslcol <- which(g1rmod$cvRisk == min(g1rmod$cvRisk, na.rm = TRUE))
            g1nr <- g1rmod$library.predict[,dslcol]
        }
    }else{
        g1nr <- g1rmod$pred
    }
    # trim small values
    g1nr[g1nr < tolg] <- tolg
    
    #-------------------------------
    # h0nr
    # rg0 [= (A0 - g0n)/g0n] ~ Q1n
    #-------------------------------
    h0rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.gr),args=list(
        Y = rg0,
        X = data.frame(Q1n = train_Q1n),
        newX = data.frame(Q1n = Q1n), # need predictions on full data here
        SL.library = SL.gr,
        obsWeights = rep(1, length(train_Q1n)),
        family = gaussian(),
        method = "method.CC_LS",
        verbose = verbose))
    if(multiAlgos){
        weightfail <- all(h0rmod$coef == 0)
        if(!weightfail){
            # Super Learner predictions 
            h0nr<- h0rmod$SL.predict[folds == validFold]
            if(!all1){
                train_h0nr<- h0rmod$SL.predict[folds != validFold]
            }else{
                train_h0nr <- h0nr
            }
        }else{
            dslcol <- which(h0rmod$cvRisk == min(h0rmod$cvRisk, na.rm = TRUE))
            h0nr <- h0rmod$library.predict[folds == validFold, dslcol]
            if(!all1){
                train_h0nr <- h0rmod$library.predict[folds != validFold, dslcol]                
            }else{
                train_h0nr <- h0nr
            }
        }
    }else{
        h0nr <- h0rmod$pred[folds == validFold]
        if(!all1){
            train_h0nr <- h0rmod$pred[folds != validFold]            
        }else{
            train_h0nr <- h0nr
        }
    }
    
    #---------------------------------------
    # h1rn
    # rg1 [= A0/g0n * (A1 - g1n)/g1n] ~ Q2n
    #---------------------------------------
    h1rmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.gr),args=list(
        Y=rg1,
        X=data.frame(Q2n = train_Q2n),
        newX = data.frame(Q2n = valid_Q2n),
        SL.library=SL.gr,
        obsWeights = rep(1, length(train_Q2n)),
        family = gaussian(),
        method = "method.CC_LS",
        verbose=verbose))
    if(multiAlgos){
        weightfail <- all(h1rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions 
            h1nr <- h1rmod$SL.predict
        }else{
            dslcol <- which(h1rmod$cvRisk == min(h1rmod$cvRisk, na.rm = TRUE))
            h1nr <- h1rmod$library.predict[,dslcol]
        }
    }else{
        h1nr <- h1rmod$pred
    }
    
    #--------------------------------------
    # hbarnr
    # A0/g0nr * h0nr ~ Q2n 
    #--------------------------------------
    hbarrmod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.gr),args=list(
        Y = train_A0 / train_g0nr * train_h0nr,
        X = data.frame(Q2n = train_Q2n),
        newX = data.frame(Q2n = valid_Q2n),
        SL.library = SL.gr,
        obsWeights = rep(1, length(train_Q2n)),
        family = gaussian(),
        method = "method.CC_LS", 
        verbose = verbose))
    if(multiAlgos){
        weightfail <- all(hbarrmod$coef==0)
        if(!weightfail){
            # Super Learner predictions 
            hbarnr <- hbarrmod$SL.predict
        }else{
            dslcol <- which(hbarrmod$cvRisk == min(hbarrmod$cvRisk, na.rm = TRUE))
            hbarnr <- hbarrmod$library.predict[,dslcol]
        }
    }else{
        hbarnr <- hbarrmod$pred
    }
    
    #--------
    # return
    #--------
    out <- list(g0nr = g0nr, g1nr = g1nr, h0nr = h0nr, h1nr = h1nr,
                hbarnr = hbarnr, g0rmod = NULL, g1rmod = NULL, h0rmod = NULL, 
                h1rmod = NULL, hbarrmod = NULL) 
    
    if(return.models){
        out$g0rmod <- g0rmod
        out$g1rmod <- g1rmod
        out$h0rmod <- h0rmod
        out$h1rmod <- h1rmod
        out$hbarrmod <- hbarrmod
    }
    
    return(out)
}
#' estimateQ
#' 
#' This function computes the conditional treatment probabilities at both timepoints.
#' 
#' @param L0 A \code{data.frame} featuring covariates measured at baseline.
#' @param L1 A \code{data.frame} featuring time-varying covariates measured at 
#' the first timepoint.
#' @param L2 A \code{vector} outcome of interest
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest. 
#' @param stratify A \code{boolean} indicating whether to pool across treatment
#' nodes or to estimate outcome regression separately in each category.
#' @param SL.Q A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the conditional probability of treatment at each time point.  
#' See \code{SuperLearner} package for details.
#' @param glm.Q A \code{character} specifying the right-hand side of the \code{glm} 
#' formula used to estimate the conditional probability of treatment at each time point. 
#' Only used if \code{SL.Q = NULL}.
#' @param SL.Q.options A \code{list} of additional arguments passed to \code{SuperLearner} 
#' for outcome regression fits.
#' @param glm.Q.options A \code{list} of additional arguments passed to \code{glm} for the
#' outcome regression fits. Typically, the \code{family} argument.
#' @param return.models  A \code{boolean} indicating whether the models for g00 should be 
#' returned with the output. 
#' @param ... Other arguments (not currently used).
#' 
#' @return Returns a list with \code{Q2n}, \code{Q1n}, and the estimated model objects if
#' \code{return.models = TRUE}
#' 
#' @importFrom SuperLearner SuperLearner
#' 
#' @export
#' 
#' @examples TO DO : add examples


estimateQ <- function(
    validFold, folds, L0, L1, L2, A0, A1, abar, SL.Q, SL.Q.options, glm.Q, 
    glm.Q.options, return.models, verbose, stratify, 
    ...
){  
    all1 <- all(folds == 1)
    if(all1){
      train_L0 <- valid_L0 <- L0
      train_L1 <- valid_L1 <- L1
      train_A0 <- valid_A0 <- A0
      train_A1 <- valid_A1 <- A1
      train_L2 <- valid_L2 <- L2
    }else{
      # training data
      train_L0 <- L0[folds != validFold, , drop = FALSE]
      train_L1 <- L1[folds != validFold, , drop = FALSE]
      train_A0 <- A0[folds != validFold]
      train_A1 <- A1[folds != validFold]
      train_L2 <- L2[folds != validFold]

      # validation data
      valid_L0 <- L0[folds == validFold, , drop = FALSE]
      valid_L1 <- L1[folds == validFold, , drop = FALSE]
      valid_A0 <- A0[folds == validFold]
      valid_A1 <- A1[folds == validFold]
      valid_L2 <- L2[folds == validFold]
    }
    if(is.null(SL.Q) & is.null(glm.Q)){ 
        stop("Specify Super Learner library or GLM formula for g")
    }
    if(!is.null(SL.Q) & !is.null(glm.Q)){
        warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
        glm.Q <- NULL
    }
    # Super Learner
    if(!is.null(SL.Q)){
        # if multiple learners specified, call SuperLearner
        multiAlgos <- (length(SL.Q) > 1 | is.list(SL.Q))
        #--------
        # Q2n
        #--------
        Q2mod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Q),args=c(list(
            Y=eval(parse(text=paste0(
                ifelse(stratify,
                       "train_L2[train_A0==abar[1] & train_A1==abar[2]]",
                       "train_L2"
                       )
            ))), 
            X=eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(train_L0,train_L1)[train_A0==abar[1] & train_A1==abar[2],]",
                       "cbind(train_L0,train_L1,train_A0,train_A1)"
                       )
            ))), 
            newX = eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(L0,L1)",
                       "cbind(L0,L1,A0=abar[1],A1=abar[2])"
                )
            ))),
            SL.library=SL.Q,
            obsWeights = rep(1, ifelse(stratify, length(train_L2[train_A0==abar[1] & train_A1==abar[2]]),length(train_L2))),
            verbose=verbose), SL.Q.options
        ))
        
        if(multiAlgos){
            # Super Learner predictions
            Q2n_all <- as.numeric(Q2mod$SL.predict)
        }else{
            Q2n_all <- as.numeric(Q2mod$pred)
        }
        # replace extrapolated predictions with 
        # smallest/largest value
        Q2n_all[Q2n_all < min(L2)] <- min(L2)
        Q2n_all[Q2n_all > max(L2)] <- max(L2)
        
        # on validation
        Q2n <- Q2n_all[folds == validFold]
        if(!all1){
          train_Q2n <- Q2n_all[folds != validFold]
        }else{
          train_Q2n <- Q2n
        }

        #-----------
        # Q1n
        #-----------
        Q1mod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.Q),args = c(list(
            Y=eval(parse(text=paste0(
                ifelse(stratify,
                       "train_Q2n[train_A0==abar[1]]",
                       "train_Q2n"
                )
            ))), 
            X=eval(parse(text=paste0(
                ifelse(stratify,
                       "train_L0[train_A0==abar[1],,drop=FALSE]",
                       "cbind(train_L0,train_A0)"
                )
            ))), 
            newX = eval(parse(text=paste0(
                ifelse(stratify,
                       "valid_L0",
                       "cbind(valid_L0,A0=abar[1])"
                )
            ))),
            SL.library = SL.Q,
            obsWeights = rep(1, ifelse(stratify, length(train_Q2n[train_A0==abar[1]]),length(train_Q2n))),
            verbose = verbose), SL.Q.options
        ))
        if(multiAlgos){
            # Super Learner predictions
            Q1n <- as.numeric(Q1mod$SL.predict)
        }else{
            Q1n <- as.numeric(Q1mod$pred)
        }
        # replace extrapolated predictions with 
        # smallest/largest value
        Q1n[Q1n < min(Q2n)] <- min(Q2n)
        Q1n[Q1n > max(Q2n)] <- max(Q2n)
    } # end Super Learner g call
    
    # if SL.Q == NULL then call glm
    if(!is.null(glm.Q)){
        #-------
        # Q2n 
        #-------
        Q2mod <- do.call("glm", c(list(formula = as.formula(paste0(
            ifelse(stratify,
                   "train_L2[train_A0==abar[1] & train_A1==abar[2]] ~",
                   "train_L2 ~"
                   ),glm.Q)),
            data=eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(train_L0,train_L1)[train_A0==abar[1] & train_A1==abar[2],]",
                       "cbind(train_L0,train_L1,train_A0)"
                       ))))), glm.Q.options))
        
        Q2n_all <- predict(Q2mod, type="response", 
                       newdata = eval(parse(text=paste0(
                           ifelse(stratify,
                                  "data.frame(L0, L1)",
                                  "data.frame(L0, L1, A0 = abar[1], A1 = abar[2])"
                                  ))))
                       )
        
        # replace extrapolated predictions with 
        # smallest/largest value
        Q2n_all[Q2n_all < min(L2)] <- min(L2)
        Q2n_all[Q2n_all > max(L2)] <- max(L2)
        # on validation
        Q2n <- Q2n_all[folds == validFold]
        if(!all1){
          train_Q2n <- Q2n_all[folds != validFold]
        }else{
          
          train_Q2n <- Q2n 
        }
        
        #-------
        # Q1n
        #-------
        Q1mod <- do.call("glm", c(list(formula = as.formula(paste0(
            ifelse(stratify,
                   "train_Q2n[train_A0==abar[1]] ~",
                   "train_Q2n ~"
            ),glm.Q)),
            data = eval(parse(text=paste0(
                ifelse(stratify,
                       "train_L0[train_A0==abar[1],,drop=FALSE]",
                       "cbind(train_L0,train_A0)"
                ))))), glm.Q.options))
        
        Q1n <- predict(Q1mod, type="response", 
                       newdata = eval(parse(text=paste0(
                           ifelse(stratify,
                                  "valid_L0",
                                  "data.frame(valid_L0, A0 = abar[1])"
                           ))))
        )
        
        # replace extrapolated predictions with 
        # smallest/largest value
        Q1n[Q1n < min(Q2n)] <- min(Q2n)
        Q1n[Q1n > max(Q2n)] <- max(Q2n)
    } # end glm call
    
    #-----------
    # returning
    #-----------
    # don't like the output when do.call is used in call
    Q2mod$call <- Q1mod$call <- NULL
    
    out <- list(Q2n = Q2n, Q1n = Q1n, 
                Q2mod = NULL, Q1mod = NULL)
    if(return.models){
        out$Q2mod <- Q2mod
        out$Q1mod <- Q1mod
    }
    return(out)
} # end function
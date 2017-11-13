#' estimateG
#' 
#' This function computes the conditional treatment probabilities at both timepoints.
#' 
#' @param L0 A \code{data.frame} featuring covariates measured at baseline.
#' @param L1 A \code{data.frame} featuring time-varying covariates measured at 
#' the first timepoint.
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest. 
#' @param stratify A \code{boolean} indicating whether to pool across treatment
#' nodes or to estimate outcome regression separately in each category.
#' @param SL.g A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the conditional probability of treatment at each time point.  See \code{SuperLearner}
#' package for details.
#' @param glm.g A \code{character} specifying the right-hand side of the \code{glm} 
#' formula used to estimate the conditional probability of treatment at each time point. 
#' Only used if \code{SL.g = NULL}.
#' @param SL.g.options A \code{list} of additional arguments passed to \code{SuperLearner} for condtional treatment 
#' probability fits.
#' @param return.models  A \code{boolean} indicating whether the models for g00 should be 
#' returned with the output. 
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param ... Other arguments (currently passed to \code{SuperLearner}).
#' 
#' @return Returns a list with \code{g0n}, \code{g1n}, and the estimated model objects if
#' \code{return.models = TRUE}
#' 
#' @importFrom SuperLearner SuperLearner
#' @export
#' 
#' @examples TO DO : add examples


estimateG <- function(
    validFold, folds, L0, L1, A0, A1, abar, SL.g, glm.g, stratify, 
    return.models, SL.g.options, verbose, tolg, ...
){  
    if(all(folds == 1)){
        train_L0 <- valid_L0 <- L0
        train_L1 <- valid_L1 <- L1
        train_A0 <- valid_A0 <- A0
        train_A1 <- valid_A1 <- A1
    }else{
        # training data
        train_L0 <- L0[folds != validFold, , drop = FALSE]
        train_L1 <- L1[folds != validFold, , drop = FALSE]
        train_A0 <- A0[folds != validFold]
        train_A1 <- A1[folds != validFold]

        # validation data
        valid_L0 <- L0[folds == validFold, , drop = FALSE]
        valid_L1 <- L1[folds == validFold, , drop = FALSE]
        valid_A0 <- A0[folds == validFold]
        valid_A1 <- A1[folds == validFold]
    }

    if(is.null(SL.g) & is.null(glm.g)){ 
        stop("Specify Super Learner library or GLM formula for g")
    }
    if(!is.null(SL.g) & !is.null(glm.g)){
        warning("Super Learner library and GLM formula specified. Proceeding with Super Learner only.")
        glm.g <- NULL
    }
    # Super Learner
    if(!is.null(SL.g)){
        # if multiple learners specified, call SuperLearner
        multiAlgos <- (length(SL.g) > 1 | is.list(SL.g))
        #--------
        # g0n
        #--------
        g0mod <- do.call(
            ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.g),args = c(list(
            Y = as.numeric(train_A0==abar[1]), X = train_L0, newX = valid_L0, 
            SL.library = SL.g, obsWeights = rep(1, length(train_A0)),
            verbose = verbose), SL.g.options
        ))
        if(multiAlgos){       
            # Super Learner predictions
            g0n <- g0mod$SL.predict
        }else{
            g0n <- g0mod$pred
        }
        # replace small predictions with tolg
        g0n[g0n < tolg] <- tolg
        g0n[g0n > 1 - tolg] <- 1 - tolg
        
        #-----------
        # g1n            
        #-----------
        g1mod <- do.call(ifelse(multiAlgos,getFromNamespace("SuperLearner","SuperLearner"),SL.g),args = c(list(
            Y=eval(parse(text=paste0(
                ifelse(stratify,
                   "as.numeric(train_A1[train_A0==abar[1]] == abar[2])",
                       "as.numeric(train_A1==abar[2])"
                       )
            ))),
            X=eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(train_L0,train_L1)[train_A0==abar[1],]",
                       "cbind(train_L0,train_L1,train_A0)"
            )))),
            newX=eval(parse(text=paste0(
                ifelse(stratify,
                       "data.frame(valid_L0, valid_L1)",
                       "data.frame(valid_L0, valid_L1, A0 = abar[1])"
                )))),
            SL.library=SL.g,
            obsWeights=rep(1,ifelse(stratify, sum(train_A0==abar[1]), length(train_A0))),
            verbose=verbose), SL.g.options)
        )
        if(multiAlgos){
            # Super Learner predictions
            g1n <- g1mod$SL.predict
        }else{
            g1n <- g1mod$pred
        }
        # replace small predictions with tolg
        g1n[g1n < tolg] <- tolg
        g1n[g1n > 1 - tolg] <- 1 - tolg
        # if only one learner specified, call it directly to avoid CV
    }
    # if SL.g == NULL then call glm
    if(!is.null(glm.g)){
        #-------
        # g0n
        #-------
        g0mod <- glm(as.formula(paste0("as.numeric(train_A0==abar[1])~",glm.g)), 
                  data=train_L0, family=binomial())
        g0n <- predict(g0mod, type="response", newdata = valid_L0)
        g0n[g0n < tolg] <- tolg
        g0n[g0n > 1 - tolg] <- 1 - tolg
        
        #-------
        # g1n 
        #-------
        g1mod <- glm(as.formula(paste0(
            ifelse(stratify,
                   "as.numeric(train_A1[train_A0==abar[1]]==abar[2])~",
                   "as.numeric(train_A1==abar[2])~"
            ),glm.g)), 
                     data=eval(parse(text=paste0(
                         ifelse(stratify,
                                "cbind(train_L0,train_L1)[train_A0==abar[1],]",
                                "cbind(train_L0,train_L1,train_A0)"
                         )))), family=binomial())
        g1n <- predict(g1mod, type="response", newdata = eval(parse(text=paste0(
            ifelse(stratify,
                   "data.frame(valid_L0, valid_L1)",
                   "data.frame(valid_L0, valid_L1, A0 = abar[1])"
            )))))
        g1n[g1n < tolg] <- tolg
        g1n[g1n > 1 - tolg] <- 1 - tolg
    } # end glm call
    
    #-----------
    # returning
    #-----------
    # don't like the output when do.call is used in call
    g0mod$call <- g1mod$call <- NULL
    
    out <- list(g0n = g0n, g1n = g1n, 
                g0mod = NULL, g1mod = NULL)
    if(return.models){
        out$g0mod <- g0mod
        out$g1mod <- g1mod
    }
    return(out)
} # end function
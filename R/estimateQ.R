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
#' @param ... Other arguments (not currently used).
#' 
#' @return Returns a list with \code{Q2n}, \code{Q1n}, and the estimated model objects if
#' \code{returnModels = TRUE}
#' 
#' @export
#' 
#' @examples TO DO : add examples


estimateQ <- function(
    L0, L1, L2, A0, A1, abar, SL.Q, glm.Q, glm.Q.options, ...
){
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
        Q2mod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.Q),args=c(list(
            Y=eval(parse(text=paste0(
                ifelse(stratify,
                       "L2[A0==abar[1] & A1==abar[2]]",
                       "L2"
                       )
            ))), 
            X=eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(L0,L1)[A0==abar[1] & A1==abar[2],]",
                       "cbind(L0,L1,A0,A1)"
                       )
            ))), 
            newX = eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(L0,L1)",
                       "cbind(L0,L1,A0=abar[1],A1=abar[2])"
                )
            ))),
            SL.library=SL.Q,
            verbose=verbose), SL.Q.options
        ))
        
        if(multiAlgos){
            # Super Learner predictions
            Q2n <- Q2mod$SL.predict
        }else{
            Q2n <- Q2mod$pred
        }
        # replace extrapolated predictions with 
        # smallest/largest value
        Q2n[Q2n < min(L2)] <- min(L2)
        Q2n[Q2n > max(L2)] <- max(L2)
        
        #-----------
        # Q1n
        #-----------
        Q1mod <- do.call(ifelse(multiAlgos,"SuperLearner",SL.Q),args = c(list(
            Y=eval(parse(text=paste0(
                ifelse(stratify,
                       "Q2n[A0==abar[1],]",
                       "Q2n"
                )
            ))), 
            X=eval(parse(text=paste0(
                ifelse(stratify,
                       "L0[A0==abar[1],,drop=FALSE]",
                       "cbind(L0,A0)"
                )
            ))), 
            newX = eval(parse(text=paste0(
                ifelse(stratify,
                       "L0",
                       "cbind(L0,A0=abar[1])"
                )
            ))),
            SL.library=SL.Q,
            verbose=verbose), SL.Q.options
        ))
        if(multiAlgos){
            # Super Learner predictions
            Q1n <- Q1mod$SL.predict
        }else{
            Q1n <- Q1mod$pred
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
                   "L2[A0==abar[1] & A1==abar[2]] ~",
                   "L2 ~"
                   ),glm.Q)),
            data=eval(parse(text=paste0(
                ifelse(stratify,
                       "cbind(L0,L1)[A0==abar[1] & A1==abar[2],]",
                       "cbind(L0,L1,A0)"
                       ))))), glm.Q.options))
        
        Q2n <- predict(Q2mod, type="response", 
                       newdata = eval(parse(text=paste0(
                           ifelse(stratify,
                                  "data.frame(L0, L1)",
                                  "data.frame(L0, L1, A0 = abar[1], A1 = abar[2])"
                                  ))))
                       )
        
        # replace extrapolated predictions with 
        # smallest/largest value
        Q2n[Q2n < min(L2)] <- min(L2)
        Q2n[Q2n > max(L2)] <- max(L2)
        
        #-------
        # Q1n
        #-------
        Q1mod <- do.call("glm", c(list(formula = as.formula(paste0(
            ifelse(stratify,
                   "Q2n[A0==abar[1]] ~",
                   "Q2n ~"
            ),glm.Q)),
            data = eval(parse(text=paste0(
                ifelse(stratify,
                       "L0[A0==abar[1],,drop=FALSE]",
                       "cbind(L0,A0)"
                ))))), glm.Q.options))
        
        Q1n <- predict(Q1mod, type="response", 
                       newdata = eval(parse(text=paste0(
                           ifelse(stratify,
                                  "L0",
                                  "data.frame(L0, A0 = abar[1])"
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
    
    out <- list(Q2n = Q2n, Q1n = g1n, 
                Q2mod = NULL, Q1mod = NULL)
    if(returnModels){
        out$Q2mod <- Q2mod
        out$Q1mod <- Q1mod
    }
    return(out)
} # end function
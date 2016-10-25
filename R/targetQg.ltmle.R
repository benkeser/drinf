#' targetQg.ltmle
#' 
#' Function that targets estimates of Q to solve the EIF. Fits a single logistic 
#' fluctuation regression of c(Q2n, L2) on relevant offsets and so-called
#' clever covariates. 
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest. 
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest. 
#' @param tolQ A \code{numeric}
#' @param return.models  A \code{boolean} indicating whether the fluctuation model should be 
#' returned with the output.  
#' @importFrom SuperLearner trimLogit
#' 
#' @return A list with named entries corresponding to the estimators of the 
#' fluctuated nuisance parameters evaluated at the observed data values. If 
#' \code{return.models = TRUE} output also includes the fitted fluctuation model. 

targetQg.ltmle <- function(
    A0, A1, L2, Q2n, Q1n, g1n, g0n, abar, tolQ, return.models, ...
){
    #-------------------------------------------
    # making outcomes for logistic fluctuation
    #-------------------------------------------
    # scale L2, Q2n, Q1n to be in (0,1)
    L2.min <- min(L2); L2.max <- max(L2)
    Q2n.min <- min(Q2n); Q2n.max <- max(Q2n)
    Q1n.min <- min(Q1n); Q1n.max <- max(Q1n)
    
    # scale L2
    L2s <- (L2 - L2.min)/(L2.max - L2.min)
    # scale Q2n for when it's an offset with L2 as outcome
    Q2ns.off <- (Q2n - L2.min)/(L2.max - L2.min)
    # scale Q2n for when it's the outcome with Q1n as offset
    Q2ns.out <- (Q2n - Q2n.min)/(Q2n.max - Q2n.min)
    # scall Q1n for when it's an offset with Q2n as outcome
    Q1ns.off <- (Q1n - Q2n.min)/(Q2n.max - Q2n.min)
    
    flucOut <- c(Q2ns.out, L2s)
    
    #-------------------------------------------
    # making offsets for logistic fluctuation
    #-------------------------------------------
    flucOff <- c(
        SuperLearner::trimLogit(Q1ns.off, trim = tolQ),
        SuperLearner::trimLogit(Q2ns.off, trim = tolQ)
    )
    
    #-------------------------------------------
    # making covariates for fluctuation
    #-------------------------------------------
    flucCov <- c(
        as.numeric(A0 == abar[1]) * (1/g0n), # for Q2n as outcome, logit(Q1n) offset
        as.numeric(A0==abar[1] & A1==abar[2]) / (g0n * g1n)  # for L2 as outcome, logit(Q2n) as covariate
    )
    
    #-------------------------------------------
    # making covariates for prediction
    #-------------------------------------------
    # getting the values of the clever covariates evaluated at 
    # \bar{A} = abar
    predCov <- c(
        1/g0n, # all A0 == abar[1]
        1/(g0n * g1n)  # all c(A0,A1) = abar
    )
    
    #-------------------------------------------
    # fitting fluctuation submodel
    #-------------------------------------------
    flucmod <- glm(
        formula = "out ~ -1 + offset(fo) + fc",
        data = data.frame(out = flucOut, fo = flucOff, fc = flucCov),
        family = binomial()
    )
    
    #--------------------------------------------
    # get predictions back using predCov
    #--------------------------------------------
    # vector of nuisance parameters
    # just set out = 0 in newdata so predict.glm does not complain
    etastar <- predict(
        flucmod, newdata = data.frame(out = 0, fo = flucOff, fc = predCov),
        type = "response"
    )
    
    #------------------------------------------------------------
    # assign etastar values to corresponding nuisance parameters 
    #------------------------------------------------------------
    # length of output
    n <- length(etastar)
    # first n entries are Q1nstar and transform back to original scale
    Q1nstar <- etastar[(1):(n)]*(Q2n.max - Q2n.min) + Q2n.min
    # next n entries are Q2nstar and transform back to original scale
    Q2nstar <- etastar[(n+1):(2*n)]*(L2.max - L2.min) + L2.min
    
    #----------------
    # return
    #----------------
    out <- list(
        Q1nstar = Q1nstar, Q2nstar = Q2nstar,
        flucmod = NULL
    )
    if(return.models){
        out$flucmod = flucmod
    }
    return(out)
}
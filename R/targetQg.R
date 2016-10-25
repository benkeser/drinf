#' targetQg
#' 
#' Function that targets estimates of Q and g to solve the EIF + remainder 
#' approximations that result from misspecification. Fits a single logistic 
#' fluctuation regression of c(A0, A1, Q2n, L2) on relevant offsets and so-called
#' clever covariates. 
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest. 
#' @param Q2n A \code{vector} of estimates of Q_{2,0}
#' @param Q1n A \code{vector} of estimates of Q_{1,0}
#' @param g1n A \code{vector} of estimates of g_{1,0}
#' @param g0n A \code{vector} of estimates of g_{0,0}
#' @param Q2nr.obsa A \code{vector} of estimates of Q_{2,0}^r evaluated at observed A0 value
#' @param Q2nr.seta A \code{vector} of estimates of Q_{2,0}^r evaluated at A0 = abar[1]
#' @param Q1nr A \code{vector} of estimates of Q_{1,0}^r
#' @param g0nr A \code{vector} of estimates of g_{0,0}^r
#' @param g1nr A \code{vector} of estimates of g_{0,0}^r
#' @param h0nr A \code{vector} of estimates of h_{0,0}^r
#' @param h1nr A \code{vector} of estimates of h_{1,0}^r
#' @param hbarnr A \code{vector} of estimates of h^r, the iterated reduced
#' dimension regression.
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest. 
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param tolQ A \code{numeric}
#' @param return.models  A \code{boolean} indicating whether the fluctuation model should be 
#' returned with the output.  
#' @param offset.model A \code{boolean} indicating whether to fit a fluctuation model
#' with an offset term (if \code{TRUE}) or intercept-only with weights (if \code{FALSE})
#' @param allatonce Do all the targeting with a single fit; if \code{FALSE}, then target g, followed by Q
#' 
#' @importFrom SuperLearner trimLogit
#' 
#' @return A list with named entries corresponding to the estimators of the 
#' fluctuated nuisance parameters evaluated at the observed data values. If 
#' \code{return.models = TRUE} output also includes the fitted fluctuation model. 
 
targetQg <- function(
    A0, A1, L2, Q2n, Q1n, g1n, g0n, Q2nr.obsa, Q2nr.seta, Q1nr, g0nr, g1nr, h0nr, h1nr, 
    hbarnr, abar, tolg, tolQ, return.models,offset.model = TRUE, allatonce = FALSE, 
    ...
){
    #-------------------------------------------
    # making outcomes for logistic fluctuation
    #-------------------------------------------
    # length of output
    n <- length(g0n)
    
    # scale L2, Q2n, Q1n to be in (0,1)
    L2.min <- min(L2); L2.max <- max(L2)

    # scale L2
    L2s <- (L2 - L2.min)/(L2.max - L2.min)
    # scale Q2n,Q1n
    Q2ns <- (Q2n - L2.min)/(L2.max - L2.min)
    Q1ns <- (Q1n - L2.min)/(L2.max - L2.min)
    
    flucOut <- c(A0, A1, Q2ns, L2s)
    
    #-------------------------------------------
    # making offsets for logistic fluctuation
    #-------------------------------------------
    flucOff <- c(
        SuperLearner::trimLogit(g0n, trim = tolg),
        SuperLearner::trimLogit(g1n, trim = tolg),
        SuperLearner::trimLogit(Q1ns, trim = tolQ),
        SuperLearner::trimLogit(Q2ns, trim = tolQ)
    )
    
    #-------------------------------------------
    # making covariates for fluctuation
    #-------------------------------------------
    # the original "clever covariates"
    flucCov1 <- c(
        rep(0,n), 
        rep(0,n),
        (L2.max - L2.min) * as.numeric(A0 == abar[1])/g0n, # the usual guy
        (L2.max - L2.min) * as.numeric(A0==abar[1] & A1==abar[2])/(g0n * g1n) # the usual guy
    )
    # the new "clever covariates" for Q
    flucCov2 <- c(
        rep(0, n), # not needed for A0
        rep(0, n), # not needed for A1    
        # the extra term for targeting Q1n
        (L2.max - L2.min) * as.numeric(A0 == abar[1]) * (h0nr/g0nr),
        # the sum of the extra two terms for targeting Q2n
        (L2.max - L2.min) * as.numeric(A0==abar[1] & A1==abar[2]) *  
            ((hbarnr + h1nr)/g1nr)
    )
    # the new "clever covariates" for g
    flucCov3 <- c(
        # the term for targeting g0n
        Q1nr / g0n,
        # the term for targeting g1n 
        Q2nr.obsa / g1n, 
        rep(0,n),
        rep(0,n)
    )
    
    
    #-------------------------------------------
    # making covariates for prediction
    #-------------------------------------------
    # getting the values of the clever covariates evaluated at 
    # \bar{A} = abar
    predCov1 <- c(
        rep(0,n), rep(0,n),
        (L2.max - L2.min)/g0n, # all A0 == abar[1]
        (L2.max - L2.min)/(g0n * g1n)  # all c(A0,A1) = abar
    )
    
    predCov2 <- c(
        rep(0,n), rep(0,n),
        (L2.max - L2.min) * h0nr/g0nr, # all A0 == abar[1]
        (L2.max - L2.min) * (hbarnr + h1nr)/g1nr  # all c(A0,A1) = abar
    )
    
    predCov3 <- c(
        Q1nr / g0n, # no need to change this one
        Q2nr.seta / g1n, # this one is set to Q2nr(abar[1], \bar{L}_1)
        rep(0,n), rep(0,n)
    )
    #-------------------------------------------
    # fitting fluctuation submodel
    #-------------------------------------------
    if(offset.model){
        if(allatonce){
            flucmod <- suppressWarnings(glm(
                formula = "out ~ -1 + offset(fo) + fc1 + fc2 + fc3",
                data = data.frame(out = flucOut, fo = flucOff, 
                                  fc1 = flucCov1, fc2 = flucCov2, fc3 = flucCov3),
                family = binomial(), start = rep(0,3)
            ))
            #--------------------------------------------
            # get predictions back using predCov
            #--------------------------------------------
            # vector of nuisance parameters
            # just set out = 0 in newdata so predict.glm does not complain
            etastar <- predict(
                flucmod, 
                newdata = data.frame(out = 0, fo = flucOff,
                                     fc1 = predCov1, fc2 = predCov2, 
                                     fc3 = predCov3),
                type = "response"
            )
        }else{
            # first do g's
            flucmodg <- suppressWarnings(glm(
                formula = "out ~ -1 + offset(fo) + fc",
                data = data.frame(out = flucOut, fo = flucOff, fc = flucCov)[1:(2*n),],
                family = binomial(), start = 0
            ))
            flucmodQ <- suppressWarnings(glm(
                formula = "out ~ -1 + offset(fo) + fc",
                data = data.frame(out = flucOut, fo = flucOff, fc = flucCov)[(2*n+1):(4*n),],
                family = binomial(), start = 0
            ))
            flucmod <- list(Q = flucmodQ, g = flucmodg)
            
            #--------------------------------------------
            # get predictions back using predCov
            #--------------------------------------------
            # vector of nuisance parameters
            # just set out = 0 in newdata so predict.glm does not complain
            etastar <- c(
                predict(flucmodg, newdata = data.frame(out = 0, fo = flucOff[1:(2*n)], fc = predCov[1:(2*n)]),
                        type = "response"),
                predict(flucmodg, newdata = data.frame(out = 0, fo = flucOff[(2*n+1):(4*n)], fc = predCov[(2*n+1):(4*n)]),
                        type = "response")
            )
        }
    }else{
        stop("optim code not written yet")
        # flucmod <- glm(
        #     formula = "out ~ offset(fo)", 
        #     weights = fc, 
        #     data = data.frame(out = flucOut, fo = flucOff, fc = flucCov),
        #     family = binomial()
        # )
    }

    #------------------------------------------------------------
    # assign etastar values to corresponding nuisance parameters 
    #------------------------------------------------------------
    # first n entries are g0nstar
    g0nstar <- etastar[1:n]
    # next n entries are g1nstar
    g1nstar <- etastar[(n+1):(2*n)]
    # next n entries are Q1nstar, transform back to original scale
    Q1nstar <- etastar[(2*n+1):(3*n)]*(L2.max - L2.min) + L2.min
    # next n entries are Q2nstar, transform back to original scale
    Q2nstar <- etastar[(3*n+1):(4*n)]*(L2.max - L2.min) + L2.min
    
    #----------------
    # return
    #----------------
    out <- list(
        g0nstar = g0nstar, g1nstar = g1nstar, Q1nstar = Q1nstar, Q2nstar = Q2nstar,
        flucmod = NULL
    )
    if(return.models){
        out$flucmod = flucmod
    }
    return(out)
}
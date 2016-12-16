#' targetQ1.ltmle
#' 
#' Function that targets Q2n for ltmle.
#' 
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest. 
#' @param Qn A \code{list} of current estimates of Q2n and Q1n
#' @param gn A \code{list} of current estimates of g1n and g0n
#' @param abar A \code{vector} of length 2 indicating the treatment assignment 
#' that is of interest. 
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param tolQ A \code{numeric}
#' @param return.models  A \code{boolean} indicating whether the fluctuation model should be 
#' returned with the output.  
#' @param tol.coef A \code{numeric} indicating the coefficient above which the minimization along the
#' submodel using \code{glm} is deemed to be unreasonable. In these cases \code{optim} is used 
#' instead to perform the fluctuation along the same submodel. 
#' @importFrom SuperLearner trimLogit
#' 
#' @return A list with named entries corresponding to the estimators of the 
#' fluctuated nuisance parameters evaluated at the observed data values. If 
#' \code{return.models = TRUE} output also includes the fitted fluctuation model. 

targetQ1.ltmle <- function(
    A0, A1, L2, Qn, gn, 
    abar, tolg, tolQ, return.models, tol.coef = 1e1, ...
){
    #-------------------------------------------
    # making outcomes for logistic fluctuation
    #-------------------------------------------
    # length of output
    n <- length(gn$g0n)
    
    # scale L2, Q2n, Q1n to be in (0,1)
    L2.min <- min(L2); L2.max <- max(L2)
    # scale L2
    L2s <- (L2 - L2.min)/(L2.max - L2.min)
    # scale Q2n,Q1n
    Q2ns <- (Qn$Q2n - L2.min)/(L2.max - L2.min)
    Q1ns <- (Qn$Q1n - L2.min)/(L2.max - L2.min)
    
    #-------------------------------------------
    # making offsets for logistic fluctuation
    #-------------------------------------------
    flucOff <- c(
        SuperLearner::trimLogit(Q1ns, trim = tolQ)
    )
    
    #-------------------------------------------
    # making covariates for fluctuation
    #-------------------------------------------
    # the original "clever covariates"
    flucCov1 <- c(
        (L2.max - L2.min) * as.numeric(A0==abar[1])/(gn$g0n) # the usual guy
    )
    
    #-------------------------------------------
    # making covariates for prediction
    #-------------------------------------------
    # getting the values of the clever covariates evaluated at 
    # \bar{A} = abar
    predCov1 <- c(
        (L2.max - L2.min)/(gn$g0n)  # all c(A0,A1) = abar
    )
    
    #-------------------------------------------
    # fitting fluctuation submodel
    #-------------------------------------------
    # first fluctuation submodel to solve original equations
    flucmod1 <- suppressWarnings(glm(
        formula = "out ~ -1 + offset(fo) + fc1",
        data = data.frame(out = Q2ns, fo = flucOff, 
                          fc1 = flucCov1),
        family = binomial(), start = 0
    ))
    # see if the fluctuation coefficient is reasonable
    if(abs(flucmod1$coefficients) < tol.coef){
        # get predictions 
        Q1nstar <- predict(
            flucmod1, 
            newdata = data.frame(out = 0, fo = flucOff,
                                 fc1 = predCov1),
            type = "response"
        )*(L2.max - L2.min) + L2.min
    }else{
        # use optim to try the minimization along submodel if glm 
        # looks wonky
        flucmod1 <- optim(
            par = 0, fn = offnegloglik, gr = gradient.offnegloglik,
            method = "L-BFGS-B", lower = -100, upper = 100,
            control = list(maxit = 10000),
            Y = Q2ns, offset = flucOff, weight = flucCov1
        )
        epsilon <- flucmod1$par
        Q1nstar <- plogis(flucOff + predCov1 * epsilon)*(L2.max - L2.min) + L2.min
    }
    
    
    #--------------
    # output 
    #-------------
    out <- list(
        Q1nstar = Q1nstar,
        flucmod = NULL
    )
    if(return.models){
        out$flucmod = list(flucmod1)
    }
    return(out)
}
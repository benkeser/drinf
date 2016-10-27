#' drinf.tmle
#' 
#' This function computes the time-varying covariate-adjusted mean of 
#' an outcome under a specified treatment assignment using targeted 
#' minimum loss-based estimation.  
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
#' nodes or to estimate outcome regression separately in each category. Should be 
#' kept \code{TRUE} until I have more time to think about how to pool across 
#' treatment arms?
#' @param SL.Q A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the outcome regressions at each time point. See \code{SuperLearner}
#' package for details.
#' @param SL.g A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the conditional probability of treatment at each time point.  See \code{SuperLearner}
#' package for details.
#' @param SL.Qr A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the reduced-dimension regression to protect against misspecification of the
#' outcome regressions.  See \code{SuperLearner} package for details.
#' @param SL.gr A \code{vector} or \code{list} specifying the SuperLearner library
#' to be used to estimate the reduced-dimension regression to protect against misspecification of the
#' conditional treatment probabilities. See \code{SuperLearner} package for details.
#' @param glm.Q A \code{character} specifying the right-hand side of the \code{glm} 
#' formula used to estimate the outcome regressions at each time point. Only used if \code{SL.Q = NULL}.
#' @param glm.g A \code{character} specifying the right-hand side of the \code{glm} 
#' formula used to estimate the conditional probability of treatment at each time point. 
#' Only used if \code{SL.g = NULL}.
#' @param guard A \code{vector} of \code{characters}, either \code{"Q"}, \code{"g"}, both, or neither (\code{NULL}).
#' Indicates whether to guard against misspecification of outcome or treatment regressions or both. Currently only works
#' with \code{c("Q","g")}. 
#' @param return.models A \code{boolean} indicating whether the models for Q, g, Qr, and gr should be 
#' returned with the output. 
#' @param maxIter A \code{numeric} indicating the maximum number of TMLE iterations before stopping. 
#' @param tolIF A \code{numeric} stopping criteria for the TMLE updates based on the empirical average of the 
#' estimated influence curve. 
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param verbose A \code{boolean} indicating whether messages should be printed to indicate progress.
#' @param SL.Q.options A \code{list} of additional arguments passed to \code{SuperLearner} for outcome
#' regression fits.
#' @param SL.g.options A \code{list} of additional arguments passed to \code{SuperLearner} for condtional treatment 
#' probability fits.
#' @param glm.Q.options A \code{list} of additional arguments passed to \code{glm} for outcome
#' regression fits.
#' @param return.ltmle A \code{boolean} indicating whether to compute the LTMLE estimate using a similar
#' iterative updating scheme.  
#' @param return.naive A \code{boolean} indicating whether to return the naive plug-in estimate. 
#' @param ... Other arguments (not currently used)
#' @return TO DO: Add return values
#' 
#' @export 
#' 
#' @examples 
#' TO DO : Add Examples

#-------------------------------------#
# TO DO 
#-------------------------------------#
# - Allow separate specification of SL or glm for each timepoint
# - Allow to adjust for misspecification separately at each timepoint
# - Optional guarding against Q, g, or both. 
# - Add separate verbose option for SL and for main function 

drinf.tmle <- function(L0, L1, L2, 
                       A0, A1,
                       abar=c(1,1),
                       stratify=TRUE,
                       SL.Q=NULL,
                       SL.g=NULL,
                       SL.Qr=NULL,
                       SL.gr=NULL,
                       glm.Q=NULL,
                       glm.g=NULL,
                       guard=c("Q","g"),
                       flucOrd = c("targetQ2","targetQ1","targetg1","targetg0"), ############### new option testing for flexible fluctuations
                       return.models=FALSE,
                       maxIter=20,
                       tolIF=1/(10*sqrt(length(L2))), 
                       tolg=1e-8,
                       tolQ=1e-8,
                       verbose=TRUE,
                       SL.Q.options = list(family = gaussian()),
                       SL.g.options = list(family = binomial()),
                       glm.Q.options = list(family = gaussian()),
                       return.ltmle = TRUE,
                       return.naive = TRUE,
                      ...){
    
    #---------------------
    # estimate g
    #---------------------
    gn <- estimateG(
        L0 = L0, L1 = L1, A0 = A0, A1 = A1, abar = abar, 
        SL.g = SL.g, glm.g = glm.g, stratify = stratify, 
        return.models = return.models, SL.g.options = SL.g.options, 
        verbose = verbose, tolg = tolg, ...
    )
    
    #---------------------
    # estimate Q
    #---------------------
    Qn <- estimateQ(
        L0 = L0, L1 = L1, L2 = L2, A0 = A0, A1 = A1, abar = abar,  
        SL.Q = SL.Q, SL.Q.options = SL.Q.options, return.models = return.models, 
        glm.Q = glm.Q, glm.Q.options = glm.Q.options, verbose = verbose,
        stratify = stratify, ...
    )
    
    #-----------------------------------------
    # compute reduced dimension regressions
    #----------------------------------------
    Qnr.gnr <- redReg(A0 = A0, A1 = A1, L2 = L2, abar = abar, 
                      gn = gn, Qn = Qn, verbose = verbose, tolg = tolg, 
                      SL.Qr = SL.Qr, SL.gr = SL.gr, return.models = return.models)
    
    #----------------------------------------
    # target Q and g according to flucOrd
    #----------------------------------------
    # initialize vectors in Qnstar and gnstar that will hold
    # targeted nuisance parameters by setting equal to the initial values
    Qnstar <- vector(mode = "list")
    gnstar <- vector(mode = "list")
    
    Qnstar$Q2n <- Qn$Q2n
    Qnstar$Q1n <- Qn$Q1n
    gnstar$g1n <- gn$g1n
    gnstar$g0n <- gn$g0n
    
    # flucOrd is a vector of functions to call in sequence to 
    # perform the targeting.
    for(ff in flucOrd){
        if(verbose){
            cat("Calling ", ff, " for targeting step. \n")
        }
        flucOut <- do.call(ff, args = list(
            A0 = A0, A1 = A1, L2 = L2, Qn = Qnstar, gn = gnstar, Qnr.gnr = Qnr.gnr, 
            tolg = tolg, tolQ = tolQ, abar = abar, return.models = return.models,
            SL.Qr = SL.Qr, SL.gr = SL.gr, verbose = verbose, ...
        ))
        # look for what names are in function output and assign values 
        # accordingly
        flucOutNames <- names(flucOut)
        if("Q2nstar" %in% flucOutNames){
            if(verbose) cat("Q2n was targeted by ", ff,". \n")
            Qnstar$Q2n <- flucOut$Q2nstar
        }
        if("Q1nstar" %in% flucOutNames){
            if(verbose) cat("Q1n was targeted by ", ff,". \n")
            Qnstar$Q1n <- flucOut$Q1nstar
        }
        if("g1nstar" %in% flucOutNames){
            if(verbose) cat("g1n was targeted by ", ff,". \n")
            gnstar$g1n <- flucOut$g1nstar
        }
        if("g0nstar" %in% flucOutNames){
            if(verbose) cat("g0n was targeted by ", ff,". \n")
            gnstar$g0n <- flucOut$g0nstar
        }
        # one of these functions could be redReg, to update the reduced
        # dimension regressions in between targeting steps -- if so, 
        # replace Qnr.gnr
        if("Qnr" %in% flucOutNames){
            Qnr.gnr <- flucOut
        }
    } 
    
    #-------------------------
    # evaluate IF
    #-------------------------
    if.dr <- evaluateIF(
        A0 = A0, A1 = A1, L2 = L2, 
        Q2n = Qnstar$Q2n, Q1n = Qnstar$Q1n, 
        g1n = gnstar$g1n, g0n = gnstar$g0n, 
        Q2nr.obsa = Qnr.gnr$Qnr$Q2nr.obsa, Q1nr = Qnr.gnr$Qnr$Q1nr, 
        g0nr = Qnr.gnr$gnr$g0nr, g1nr = Qnr.gnr$gnr$g1nr, 
        h0nr = Qnr.gnr$gnr$h0nr, h1nr = Qnr.gnr$gnr$h1nr, hbarnr = Qnr.gnr$gnr$hbarnr,
        abar = abar
    )
    # mean of IF -- first three terms are added, last 5 are subtracted
    meanif.dr <- c(
        # original terms
        t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[1:3])),
        # extra terms tageting g's
        t(matrix(c(1,1)))%*%colMeans(Reduce("cbind",if.dr[4:5])),
        # extra terms targeting Q's
        t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[6:8]))
    )
    
    #-----------------------------------
    # targeting loop for dr inference
    #-----------------------------------
    iter <- 1
    while(max(abs(meanif.dr)) > tolIF & iter < maxIter){
        #-----------------------------------------
        # compute reduced dimension regressions
        #----------------------------------------
        Qnr.gnr <- redReg(A0 = A0, A1 = A1, L2 = L2, abar = abar, 
                          gn = gnstar, Qn = Qnstar, verbose = verbose, tolg = tolg, 
                          SL.Qr = SL.Qr, SL.gr = SL.gr, return.models = return.models)
        
        #--------------------
        # target Q and g
        #--------------------
        for(ff in flucOrd){
            if(verbose){
                cat("Calling ", ff, " for targeting step. \n")
            }
            flucOut <- do.call(ff, args = list(
                A0 = A0, A1 = A1, L2 = L2, Qn = Qnstar, gn = gnstar, Qnr.gnr = Qnr.gnr, 
                tolg = tolg, tolQ = tolQ, abar = abar, return.models = return.models,
                SL.Qr = SL.Qr, SL.gr = SL.gr,verbose = verbose, ...
            ))
            # look for what names are in function output and assign values 
            # accordingly
            flucOutNames <- names(flucOut)
            if("Q2nstar" %in% flucOutNames){
                if(verbose) cat("Q2n was targeted by ", ff,". \n")
                Qnstar$Q2n <- flucOut$Q2nstar
            }
            if("Q1nstar" %in% flucOutNames){
                if(verbose) cat("Q1n was targeted by ", ff,". \n")
                Qnstar$Q1n <- flucOut$Q1nstar
            }
            if("g1nstar" %in% flucOutNames){
                if(verbose) cat("g1n was targeted by ", ff,". \n")
                gnstar$g1n <- flucOut$g1nstar
            }
            if("g0nstar" %in% flucOutNames){
                if(verbose) cat("g0n was targeted by ", ff,". \n")
                gnstar$g0n <- flucOut$g0nstar
            }
            if("Qnr" %in% flucOutNames){
                Qnr.gnr <- flucOut
            }
        } 
        
        #-------------------------
        # evaluate IF
        #-------------------------
        if.dr <- evaluateIF(
            A0 = A0, A1 = A1, L2 = L2, 
            Q2n = Qnstar$Q2n, Q1n = Qnstar$Q1n, 
            g1n = gnstar$g1n, g0n = gnstar$g0n, 
            Q2nr.obsa = Qnr.gnr$Qnr$Q2nr.obsa, Q1nr = Qnr.gnr$Qnr$Q1nr, 
            g0nr = Qnr.gnr$gnr$g0nr, g1nr = Qnr.gnr$gnr$g1nr, 
            h0nr = Qnr.gnr$gnr$h0nr, h1nr = Qnr.gnr$gnr$h1nr, hbarnr = Qnr.gnr$gnr$hbarnr,
            abar = abar
        )
        # mean of IF -- first three terms are added, last 5 are subtracted
        meanif.dr <- c(
            # original terms
            t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[1:3])),
            # extra terms tageting g's
            t(matrix(c(1,1)))%*%colMeans(Reduce("cbind",if.dr[4:5])),
            # extra terms targeting Q's
            t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[6:8]))
        )
        
        # update iteration
        iter <- iter + 1
        cat(#"\n epsilon = ", etastar$flucmod$coefficients[1],
            "\n mean ic = ", meanif.dr, 
             "\n")
    }
    
    # cat("out of dr tmle loop")
    #------------------------------------------
    # evaluate parameter and compute std err
    #------------------------------------------
    # evaluate parameter
    psin <- mean(Qnstar$Q1n)
    
    # compute standard errors
    se <- sqrt(mean(
        (if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2 -
            if.dr$Dg1.Q2 - if.dr$Dg0.Q1 - if.dr$DQ2.g1 - 
            if.dr$DQ2.g0 - if.dr$DQ1.g0)^2
    )/length(A0))

    #-----------------------------------
    # Computing the regular LTMLE
    #-----------------------------------
    if(return.ltmle){
        Q2nstar.ltmle <- targetQ2.ltmle(
            A0 = A0, A1 = A1, L2 = L2, Qn = Qn, gn = gn, 
            abar = abar, tolg = tolg, tolQ = tolQ, return.models = return.models
        )$Q2nstar
        
        Q1n <- estimateQ1.ltmle(
            L2 = L2, L0 = L0, Q2n = Q2nstar.ltmle, A0 = A0, A1 = A1, 
            abar = abar, SL.Q = SL.Q, SL.Q.options = SL.Q.options, 
            glm.Q = glm.Q, glm.Q.options = glm.Q.options, 
            return.models = return.models, verbose = verbose, stratify = stratify
        )$Q1n
        
        Q1nstar.ltmle <- targetQ1.ltmle(
            A0 = A0, A1 = A1, L2 = L2, Qn = list(Q2n = Q2nstar.ltmle, Q1n = Qn$Q1n), gn = gn, 
            abar = abar, tolg = tolg, tolQ = tolQ, return.models = return.models
        )$Q1nstar
        
        # point estimate
        psin.ltmle <- mean(Q1nstar.ltmle)
        
        # influence functions
        if.ltmle <- evaluateEIF(
            A0 = A0, A1 = A1, L2 = L2, Q2n = Q2nstar.ltmle, 
            Q1n = Q1nstar.ltmle, g1n = gn$g1n, g0n = gn$g0n, abar = abar
        )
        
        # se
        se.ltmle <- sqrt(mean(
            (if.ltmle$Dstar2 + if.ltmle$Dstar1 + if.ltmle$Dstar0)^2
        )/length(A0))
        
    } # end if(return.ltmle)
    
    #----------------------------
    # formatting output
    #----------------------------
    out <- vector(mode = "list")
    
    # dr tmle output
    out$est <- psin
    out$se <- se
    
    # ltmle output
    out$est.ltmle <- NULL
    out$se.ltmle <- NULL
    if(return.ltmle){
        out$est.ltmle <- psin.ltmle
        out$se.ltmle <- se.ltmle
    }
    
    # naive output
    out$est.naive <- NULL
    if(return.ltmle){
        out$est.naive <- mean(Qn$Q1n)
    }
    
    # number of iterations
    out$iter <- iter
    
    # model output
    out$Qmod <- vector(mode = "list")
    out$gmod <- vector(mode = "list")
    out$Qrmod <- vector(mode = "list")
    out$grmod <- vector(mode = "list")
    
    # setting up models for return
    if(return.models){
        out$Qmod$Q2n <- Qn$Q2nmod
        out$Qmod$Q1n <- Qn$Q1nmod
        out$gmod$g0n <- gn$g0nmod
        out$gmod$g1n <- gn$g1nmod
        out$Qrmod$Q2nr <- Qnr.gnr$Q2nrmod
        out$Qrmod$Q1nr <- Qnr.gnr$Q1nrmod
        out$grmod$g1nr <- Qnr.gnr$g1nrmod
        out$grmod$g0nr <- Qnr.gnr$g0nrmod
        out$grmod$h0nr <- Qnr.gnr$h0nrmod
        out$grmod$h1nr <- Qnr.gnr$h1nrmod
        out$grmod$hbarnr <- Qnr.gnr$hbarnrmod
        # out$flucmod <- etastar$flucmod
        if(return.ltmle){
            out$flucmod.ltmle <- etastar.ltmle$flucmod
        }
    }
    
    return(out)
}

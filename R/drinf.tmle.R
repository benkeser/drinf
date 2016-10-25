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
    
    #-------------------------------------#
    # Workflow 
    #-------------------------------------#
    # ++ Estimate g
    # ++ Estimate Q
    # ++ Compute outcomes of Qr regressions
    # ++ Estimate Qr regressions
    # ++ Estimate gr regression
    # ++ - internally, compute outcome(s) of iterated gr regressions
    # ++ Evaluate IC
    # while Pn D* < tolIF {
    #    ++ Target Q, g
    #    ++ Compute outcomes of Qr regressions
    #    ++ Estimate Qr
    #    ++ Estimate gr regression
    #    ++    - internally, compute outcomes of iterated gr regressions
    #    ++ Evaluate IC
    # }
    # ++ Evaluate parameter
    # ++ Evaluate standard error
    
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
    #-----------------------------------------
    Qrn.grn <- redReg(A0 = A0, A1 = A1, L2 = L2, abar = abar, 
                      SL.Qr = SL.Qr, SL.gr = SL.gr, return.models = return.models)
    
    #--------------------
    # target Q and g
    #--------------------
    # initialize vectors in Qnstar and gnstar that will hold
    # targeted nuisance parameters by setting equal to the initial values
    Qnstar <- vector(mode = "list")
    gnstar <- vector(mode = "list")

    Qnstar$Q2nstar <- Qn$Q2n
    Qnstar$Q1nstar <- Qn$Q1n
    gnstar$g1nstar <- gn$g1n
    gnstar$g0nstar <- gn$g0n
    
    # run the targeting function
    etastar <- targetQg(
        A0 = A0, A1 = A1, L2 = L2, Q2n = Qnstar$Q2nstar, Q1n = Qnstar$Q1nstar, 
        g1n = gnstar$g1nstar, g0n = gnstar$g0nstar, 
        Q2nr.seta = Qnr$Q2nr.seta, Q2nr.obsa = Qnr$Q2nr.obsa, 
        Q1nr = Qnr$Q1nr, g0nr = gnr$g0nr, g1nr = gnr$g1nr, 
        h0nr = gnr$h0nr, h1nr = gnr$h1nr, 
        hbarnr = gnr$hbarnr, abar = abar, tolg = tolg, tolQ = tolQ, 
        return.models = return.models, ...
    )
    #--------------------------------------------------------
    # assign new values of nuisance parameters to Qn and gn
    #--------------------------------------------------------
    Qnstar$Q2nstar <- etastar$Q2nstar
    Qnstar$Q1nstar <- etastar$Q1nstar
    gnstar$g1nstar <- etastar$g1nstar
    gnstar$g0nstar <- etastar$g0nstar
    
    #-------------------------
    # evaluate IF
    #-------------------------
    if.dr <- evaluateIF(
        A0 = A0, A1 = A1, L2 = L2, 
        Q2n = Qnstar$Q2nstar, Q1n = Qnstar$Q1nstar, 
        g1n = gnstar$g1nstar, g0n = gnstar$g0nstar, 
        Q2nr.obsa = Qnr$Q2nr.obsa, Q1nr = Qnr$Q1nr, 
        g0nr = gnr$g0nr, g1nr = gnr$g1nr, 
        h0nr = gnr$h0nr, h1nr = gnr$h1nr, hbarnr = gnr$hbarnr,
        abar = abar
    )
    # mean of IF -- first three terms are added, last 5 are subtracted
    meanif.dr <- t(matrix(c(1,1,1,-1,-1,-1,-1,-1)))%*%
        colMeans(Reduce("cbind",if.dr))
    
    #-----------------------------------
    # targeting loop for dr inference
    #-----------------------------------
    iter <- 1
    while(abs(meanif.dr) > tolIF & iter < maxIter){
        #---------------------------------------
        # residual outcomes for Qr regressions
        #---------------------------------------
        # write over original values
        # same call as earlier, but now with Qnstar and gnstar
        rQ <- residQ(
            L2 = L2, A0 = A0, A1 = A1, 
            Q2n = Qnstar$Q2nstar, Q1n = Qnstar$Q1nstar, 
            g0n = gnstar$g0nstar, g1n = gnstar$g1nstar, abar = abar, ...
        )
        
        #---------------------------
        # estimate Qr regressions
        #---------------------------
        # write over original values
        # same call as earlier, but now with Qnstar and gnstar
        Qnr <- estimateQr(
            rQ1 = rQ$rQ1, rQ2 = rQ$rQ2, 
            g0n = gnstar$g0nstar, g1n = gnstar$g1nstar, 
            A0 = A0, A1 = A1, SL.Qr = SL.Qr, abar = abar, 
            verbose = verbose,return.models = return.models, ...
        )
        
        #---------------------------------------
        # residual outcomes for gr regressions
        #---------------------------------------
        # write over original values
        # same call as earlier, but now with Qnstar and gnstar
        rg <- residG(
            A0 = A0, A1 = A1, g0n = gnstar$g0nstar, g1n = gnstar$g1nstar, abar = abar, ...
        )
        
        #---------------------------
        # estimate gr regressions
        #---------------------------
        # write over original values
        # same call as earlier, but now with Qnstar and gnstar
        gnr <- estimategr(
            rg0 = rg$rg0, rg1 = rg$rg1, g0n = gnstar$g0nstar, 
            g1n = gnstar$g1nstar, A0 = A0, A1 = A1, 
            Q2n = Qnstar$Q2nstar,Q1n = Qnstar$Q1nstar, verbose = verbose, 
            SL.gr = SL.gr, abar = abar, return.models = return.models, 
            tolg = tolg, ...
        )
        
        #--------------------
        # target Q and g
        #--------------------
        # run the targeting function
        etastar <- targetQg(
            A0 = A0, A1 = A1, L2 = L2, 
            Q2n = Qnstar$Q2nstar, Q1n = Qnstar$Q1nstar, 
            g1n = gnstar$g1nstar, g0n = gnstar$g0nstar, 
            Q2nr.seta = Qnr$Q2nr.seta, Q2nr.obsa = Qnr$Q2nr.obsa, 
            Q1nr = Qnr$Q1nr, g0nr = gnr$g0nr, g1nr = gnr$g1nr, 
            h0nr = gnr$h0nr, h1nr = gnr$h1nr, 
            hbarnr = gnr$hbarnr, abar = abar, tolg = tolg, tolQ = tolQ, 
            return.models = return.models, ...
        )
        #--------------------------------------------------------
        # assign new values of nuisance parameters to Qn and gn
        #--------------------------------------------------------
        Qnstar$Q2nstar <- etastar$Q2nstar
        Qnstar$Q1nstar <- etastar$Q1nstar
        gnstar$g1nstar <- etastar$g1nstar
        gnstar$g0nstar <- etastar$g0nstar
        
        #-------------------------
        # evaluate IC
        #-------------------------
        # write over original values
        # same call as earlier, but now with Qnstar, gnstar
        if.dr <- evaluateIF(
            A0 = A0, A1 = A1, L2 = L2, 
            Q2n = Qnstar$Q2nstar, Q1n = Qnstar$Q1nstar, 
            g1n = gnstar$g1nstar, g0n = gnstar$g0nstar, 
            Q2nr.obsa = Qnr$Q2nr.obsa, Q1nr = Qnr$Q1nr, 
            g0nr = gnr$g0nr, g1nr = gnr$g1nr, 
            h0nr = gnr$g0nr, h1nr = gnr$h1nr, hbarnr = gnr$hbarnr, abar = abar
        )
        # mean of IC
        meanif.dr <- t(matrix(c(1,1,1,-1,-1,-1,-1,-1)))%*%
            colMeans(Reduce("cbind",if.dr))
        # update iteration
        iter <- iter + 1
        cat("\n epsilon = ", etastar$flucmod$coefficients[1],
            "\n mean ic = ", meanif.dr, 
            "\n")
    }
    
    cat("out of dr tmle loop")
    #------------------------------------------
    # evaluate parameter and compute std err
    #------------------------------------------
    # evaluate parameter
    psin <- mean(Qnstar$Q1nstar)
    
    # compute standard errors
    se <- sqrt(mean(Reduce("c",if.dr)^2)/length(A0))
    
    #-----------------------------------
    # Computing the regular LTMLE
    #-----------------------------------
    if(return.ltmle){
        #-------------------------
        # evaluate IF
        #-------------------------
        # TO DO : This is redundant -- should it be renamed above to 
        # avoid this redundancy?
        if.ltmle <- evaluateEIF(
            A0 = A0, A1 = A1, L2 = L2, Q2n = Qn$Q2n, Q1n = Qn$Q1n, 
            g1n = gn$g1n, g0n = gn$g0n, abar
        )
        # mean of IC
        meanif.ltmle <- sum(colMeans(Reduce("cbind",if.ltmle)))
        
        # initialize vectors in Qnstar and gnstar that will hold
        # targeted nuisance parameters by setting equal to the initial values
        Qnstar.ltmle <- vector(mode = "list")
        gnstar.ltmle <- vector(mode = "list")
        
        Qnstar.ltmle$Q2nstar <- Qn$Q2n
        Qnstar.ltmle$Q1nstar <- Qn$Q1n
        gnstar.ltmle$g1nstar <- gn$g1n
        gnstar.ltmle$g0nstar <- gn$g0n
        
        iter.ltmle <- 1
        while((abs(meanif.ltmle) > tolIF & iter.ltmle < maxIter) | iter.ltmle == 1){
            #--------------------
            # fluctuate Q and g
            #--------------------
            etastar.ltmle <- targetQg.ltmle(
                A0 = A0, A1 = A1, L2 = L2, 
                Q2n = Qnstar.ltmle$Q2nstar, Q1n = Qnstar.ltmle$Q1nstar, 
                g1n = gnstar.ltmle$g1nstar, g0n = gnstar.ltmle$g0nstar, abar = abar, tolQ = tolQ, 
                return.models = return.models, ...
            )
            #--------------------------------------------------------
            # assign new values of nuisance parameters to Qn and gn
            #--------------------------------------------------------
            Qnstar.ltmle$Q2star <- etastar.ltmle$Q1nstar
            Qnstar.ltmle$Q1star <- etastar.ltmle$Q0nstar
            gnstar.ltmle$g1star <- etastar.ltmle$g1nstar
            gnstar.ltmle$g0star <- etastar.ltmle$g0nstar
            #-------------------------
            # evaluate EIF
            #-------------------------
            # write over original values
            # same call as earlier, but now with Qnstar, gnstar
            if.ltmle <- evaluateEIF(
                A0 = A0, A1 = A1, L2 = L2, abar = abar, 
                Q2n = Qnstar.ltmle$Q2nstar, Q1n = Qnstar.ltmle$Q1nstar, 
                g1n = gnstar.ltmle$g1nstar, g0n = gnstar.ltmle$g0nstar
            )
            # mean of EIF
            meanif.ltmle <- sum(colMeans(Reduce("cbind",if.ltmle)))
        }
        
        #------------------------------------------
        # evaluate parameter and compute std err
        #------------------------------------------
        # evaluate parameter
        psin.ltmle <- mean(Qnstar.ltmle$Q1nstar)
        
        # compute standard errors
        se.ltmle <- sqrt(mean(Reduce("c",if.ltmle)^2)/length(A0))
        # update iter
        iter <- iter + 1
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
        out$Qrmod$Q2nr <- Qnr$Q2nrmod
        out$Qrmod$Q1nr <- Qnr$Q1nrmod
        out$grmod$g1nr <- gnr$g1nrmod
        out$grmod$g0nr <- gnr$g0nrmod
        out$grmod$h0nr <- gnr$h0nrmod
        out$grmod$h1nr <- gnr$h1nrmod
        out$grmod$hbarnr <- gnr$hbarnrmod
        out$flucmod <- etastar$flucmod
        if(return.ltmle){
            out$flucmod.ltmle <- etastar.ltmle$flucmod
        }
    }
    
    return(out)
}

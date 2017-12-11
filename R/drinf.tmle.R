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
#' @param tolQ A \code{numeric} indicating the truncation level for transformed outcome regressions.
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
#' @param universal A \code{boolean} indicating whether to perform TMLE step using locally least favorable
#' parametric submodels (if \code{FALSE}) or universally least favorable submodels (if \code{TRUE})
#' @param universalStepSize A \code{numeric} indicating the step size for the recursive calculation of 
#' universally least favorable submodel. Default is \code{0.005}. 
#' @param only.ltmle Only return ltmle (for bootstrapping)
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
                       universal=FALSE,
                       universalStepSize=1e-4,
                       printFreq = 50, 
                       flucOrd = c("targetQ2","targetQ1","targetg1","targetg0"),
                       return.models=FALSE,
                       maxIter=20,
                       tolIF=1/(length(L2)), 
                       tolg=1e-8,
                       tolQ=1e-8,
                       verbose=TRUE,
                       SL.Q.options = list(family = gaussian()),
                       SL.g.options = list(family = binomial()),
                       glm.Q.options = list(family = gaussian()),
                       return.ltmle = TRUE, 
                       only.ltmle = FALSE, 
                       return.naive = TRUE,
                       cvFolds = 1, 
                      ...){
    
    #---------------------
    # estimate g
    #---------------------
    n <- length(L2)
    folds <- rep(1:cvFolds, each = ceiling(n/cvFolds), length.out=n)

    gn_list <- sapply(1:cvFolds, estimateG, 
                        L0 = L0, A0 = A0, L1 = L1, A1 = A1,
                        abar = abar, folds = folds, 
                        SL.g = SL.g, glm.g = glm.g, stratify = stratify, 
                        return.models = return.models, SL.g.options = SL.g.options, 
                        verbose = verbose, tolg = tolg, simplify = FALSE)
    
    #---------------------
    # estimate Q
    #---------------------
    Qn_list <- sapply(1:cvFolds, estimateQ, 
                      L0 = L0, A0 = A0, 
                      L1 = L1, A1 = A1, 
                      L2 = L2, FUN = estimateQ, 
                      abar = abar, folds = folds, 
                      SL.Q = SL.Q, SL.Q.options = SL.Q.options, 
                      return.models = return.models, 
                      glm.Q = glm.Q, glm.Q.options = glm.Q.options, 
                      verbose = verbose, stratify = stratify,
                      simplify = FALSE)
    
    #-----------------------------------------
    # compute reduced dimension regressions
    #----------------------------------------
    if(!only.ltmle){
    Qnr.gnr_list <- sapply(1:cvFolds, redReg, 
                      gn = gn_list, Qn = Qn_list, 
                      A0 = A0, A1 = A1, L2 = L2, 
                      folds = folds, abar = abar, 
                      verbose = verbose, 
                      tolg = tolg, SL.Qr = SL.Qr, SL.gr = SL.gr,
                      return.models = return.models,
                      simplify = FALSE)
    
    #----------------------------------------
    # targeting portion of method
    #----------------------------------------
    if(!universal){
        #----------------------------------------
        # targeting according to locally least 
        # favorable models
        #----------------------------------------
        
        # set up an empty vector of estimates
        psin_trace <- rep(NA, maxIter)
        se_trace <- rep(NA, maxIter)
        #----------------------------------------
        # target Q and g according to flucOrd
        #----------------------------------------
        # initialize vectors in Qnstar and gnstar that will hold
        # targeted nuisance parameters by setting equal to the initial values
        Qnstar_list <- Qn_list
        gnstar_list <- gn_list
        full_Qnstar <- do.call(Map, c(c, Qnstar_list))
        full_gnstar <- do.call(Map, c(c, gnstar_list))
        tmp <- do.call(Map, c(c, Qnr.gnr_list))
        full_Qnr.gnr <- vector(mode = "list")
        full_Qnr.gnr$Qnr <- sapply(unique(names(tmp$Qnr)), function(x) unname(unlist(tmp$Qnr[names(tmp$Qnr)==x])), simplify=FALSE)
        full_Qnr.gnr$gnr <- sapply(unique(names(tmp$gnr)), function(x) unname(unlist(tmp$gnr[names(tmp$gnr)==x])), simplify=FALSE)

        est.naive <- mean(full_Qnstar$Q1n)
        # flucOrd is a vector of functions to call in sequence to 
        # perform the targeting.
        # for(ff in flucOrd){
        #     if(verbose){
        #         cat("Calling ", ff, " for targeting step. \n")
        #     }
        #     flucOut <- do.call(ff, args = list(
        #         A0 = A0, A1 = A1, L2 = L2, Qn = Qnstar, gn = gnstar, Qnr.gnr = Qnr.gnr, 
        #         tolg = tolg, tolQ = tolQ, abar = abar, return.models = return.models,
        #         SL.Qr = SL.Qr, SL.gr = SL.gr, verbose = verbose, ...
        #     ))
        #     # look for what names are in function output and assign values 
        #     # accordingly
        #     flucOutNames <- names(flucOut)
        #     if("Q2nstar" %in% flucOutNames){
        #         if(verbose) cat("Q2n was targeted by ", ff,". \n")
        #         Qnstar$Q2n <- flucOut$Q2nstar
        #     }
        #     if("Q1nstar" %in% flucOutNames){
        #         if(verbose) cat("Q1n was targeted by ", ff,". \n")
        #         Qnstar$Q1n <- flucOut$Q1nstar
        #     }
        #     if("g1nstar" %in% flucOutNames){
        #         if(verbose) cat("g1n was targeted by ", ff,". \n")
        #         gnstar$g1n <- flucOut$g1nstar
        #     }
        #     if("g0nstar" %in% flucOutNames){
        #         if(verbose) cat("g0n was targeted by ", ff,". \n")
        #         gnstar$g0n <- flucOut$g0nstar
        #     }
        #     # one of these functions could be redReg, to update the reduced
        #     # dimension regressions in between targeting steps -- if so, 
        #     # replace Qnr.gnr
        #     if("Qnr" %in% flucOutNames){
        #         Qnr.gnr <- flucOut
        #     }
        # } 
        
        # #-------------------------
        # # evaluate IF
        # #-------------------------
        # if.dr <- evaluateIF(
        #     A0 = A0, A1 = A1, L2 = L2, 
        #     Q2n = Qnstar$Q2n, Q1n = Qnstar$Q1n, 
        #     g1n = gnstar$g1n, g0n = gnstar$g0n, 
        #     Q2nr.obsa = Qnr.gnr$Qnr$Q2nr.obsa, 
        #     Q1nr1 = Qnr.gnr$Qnr$Q1nr1, Q1nr2 = Qnr.gnr$Qnr$Q1nr2, 
        #     g0nr = Qnr.gnr$gnr$g0nr, g1nr = Qnr.gnr$gnr$g1nr, 
        #     h0nr = Qnr.gnr$gnr$h0nr, h1nr = Qnr.gnr$gnr$h1nr, hbarnr = Qnr.gnr$gnr$hbarnr,
        #     abar = abar
        # )
        # # mean of IF -- first three terms are added, last 5 are subtracted
        # meanif.dr <- c(
        #     # original terms
        #     t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[1:3])),
        #     # extra terms tageting g's
        #     t(matrix(c(1,1)))%*%colMeans(Reduce("cbind",if.dr[4:5])),
        #     # extra terms targeting Q's
        #     t(matrix(c(1,1,1)))%*%colMeans(Reduce("cbind",if.dr[6:8]))
        # )
        
        #-----------------------------------
        # targeting loop for dr inference
        #-----------------------------------
        iter <- 0
        meanif.dr <- Inf
        n_minus_12 <- length(L2)^(-1/2)
        sqrt_n_max_block <- n_max_block <- FALSE
        sqrt_n_norm_block <- n_norm_block <- FALSE
        sqrt_n_max_iter <- sqrt_n_norm_iter <- n_norm_iter <- n_max_iter <- maxIter
        while((max(abs(meanif.dr)) > tolIF | sqrt(sum(meanif.dr^2)) > tolIF) & 
                iter < maxIter){
            iter <- iter + 1

            
            #-----------------------------------------
            # compute reduced dimension regressions
            #----------------------------------------
            # Qnr.gnr <- redReg(A0 = A0, A1 = A1, L2 = L2, abar = abar, 
            #                   gn = gnstar, Qn = Qnstar, verbose = verbose, 
            #                   tolg = tolg, SL.Qr = SL.Qr, SL.gr = SL.gr, return.models = return.models)
            
            #--------------------
            # target Q and g
            #--------------------
            in_ct <- 0
            for(ff in flucOrd){
                in_ct <- in_ct + 1
                if(verbose){
                    cat("Calling ", ff, " for targeting step. \n")
                }
                if(ff != "redReg"){
                    flucOut <- do.call(ff, args = list(
                        A0 = A0, A1 = A1, L2 = L2, Qn = full_Qnstar, gn = full_gnstar, 
                        Qnr.gnr = full_Qnr.gnr, tolg = tolg, tolQ = tolQ, abar = abar, 
                        return.models = return.models, SL.Qr = SL.Qr, SL.gr = SL.gr, 
                        verbose = verbose, ...
                    ))
                    flucOutNames <- names(flucOut)
                }else if(ff == "redReg" & !(iter == 1 & in_ct == 1)){
                    # fit reduced regression
                    Qnr.gnr_list <- sapply(1:cvFolds, redReg, 
                      gn = gnstar_list, Qn = Qnstar_list, 
                      A0 = A0, A1 = A1, L2 = L2, 
                      folds = folds, abar = abar, 
                      verbose = verbose, 
                      tolg = tolg, SL.Qr = SL.Qr, SL.gr = SL.gr,
                      return.models = return.models,
                      simplify = FALSE)
                    tmp <- do.call(Map, c(c, Qnr.gnr_list))
                    full_Qnr.gnr <- vector(mode = "list")
                    full_Qnr.gnr$Qnr <- sapply(unique(names(tmp$Qnr)), function(x) unname(unlist(tmp$Qnr[names(tmp$Qnr)==x])), simplify=FALSE)
                    full_Qnr.gnr$gnr <- sapply(unique(names(tmp$gnr)), function(x) unname(unlist(tmp$gnr[names(tmp$gnr)==x])), simplify=FALSE)
                    if(verbose) cat("Reduced-dimension regressions updated. \n")
                    flucOutNames <- "Nothing to see here."
                }else if(ff == "redReg" & iter == 1 & in_ct ==1){
                    flucOutNames <- "Nothing to see here."
                }
                # look for what names are in function output and assign values 
                # accordingly
                if("Q2nstar" %in% flucOutNames){
                    if(verbose) cat("Q2n was targeted by ", ff,". \n")
                    split_Q2nstar <- split(flucOut$Q2nstar, folds)
                    Qnstar_list <- mapply(x = Qnstar_list, values = split_Q2nstar, 
                                     FUN = list_replace, MoreArgs = list(list = "Q2n"),
                                     SIMPLIFY = FALSE)
                    full_Qnstar <- do.call(Map, c(c, Qnstar_list))
                }
                if("Q1nstar" %in% flucOutNames){
                    if(verbose) cat("Q1n was targeted by ", ff,". \n")
                    split_Q1nstar <- split(flucOut$Q1nstar, folds)
                    Qnstar_list <- mapply(x = Qnstar_list, values = split_Q1nstar, 
                                     FUN = list_replace, MoreArgs = list(list = "Q1n"),
                                     SIMPLIFY = FALSE)
                    full_Qnstar <- do.call(Map, c(c, Qnstar_list))
                }
                if("g1nstar" %in% flucOutNames){
                    if(verbose) cat("g1n was targeted by ", ff,". \n")
                    split_g1nstar <- split(flucOut$g1nstar, folds)
                    gnstar_list <- mapply(x = gnstar_list, values = split_g1nstar, 
                                     FUN = list_replace, MoreArgs = list(list = "g1n"),
                                     SIMPLIFY = FALSE)
                    full_gnstar <- do.call(Map, c(c, gnstar_list))        
                }
                if("g0nstar" %in% flucOutNames){
                    if(verbose) cat("g0n was targeted by ", ff,". \n")
                    split_g0nstar <- split(flucOut$g0nstar, folds)
                    gnstar_list <- mapply(x = gnstar_list, values = split_g0nstar, 
                                     FUN = list_replace, MoreArgs = list(list = "g0n"),
                                     SIMPLIFY = FALSE)
                    full_gnstar <- do.call(Map, c(c, gnstar_list))             
                }
                # if("Qnr" %in% flucOutNames){
                #     Qnr.gnr <- flucOut
                # }
            } 
            
            #-------------------------
            # evaluate IF
            #-------------------------
            # full_Qnstar <- do.call(Map, c(c, Qnstar_list))
            # full_gnstar <- do.call(Map, c(c, gnstar_list))
            # tmp <- do.call(Map, c(c, Qnr.gnr_list))
            # full_Qnr.gnr <- vector(mode = "list")
            # full_Qnr.gnr$Qnr <- sapply(unique(names(tmp$Qnr)), function(x) unname(unlist(tmp$Qnr[names(tmp$Qnr)==x])), simplify=FALSE)
            # full_Qnr.gnr$gnr <- sapply(unique(names(tmp$gnr)), function(x) unname(unlist(tmp$gnr[names(tmp$gnr)==x])), simplify=FALSE)

            if.dr <- evaluateIF(
                A0 = A0, A1 = A1, L2 = L2, 
                Q2n = full_Qnstar$Q2n, Q1n = full_Qnstar$Q1n, 
                g1n = full_gnstar$g1n, g0n = full_gnstar$g0n, 
                Q2nr.obsa = full_Qnr.gnr$Qnr$Q2nr.obsa, 
                Q1nr1 = full_Qnr.gnr$Qnr$Q1nr1, Q1nr2 = full_Qnr.gnr$Qnr$Q1nr2, 
                g0nr = full_Qnr.gnr$gnr$g0nr, g1nr = full_Qnr.gnr$gnr$g1nr, 
                h0nr = full_Qnr.gnr$gnr$h0nr, h1nr = full_Qnr.gnr$gnr$h1nr, 
                hbarnr = full_Qnr.gnr$gnr$hbarnr,
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

            if(max(abs(meanif.dr)) < n_minus_12 & !sqrt_n_max_block){
                sqrt_n_max_iter <- iter
                sqrt_n_max_block <- TRUE # so it doesn't reset 
            }

            if(sqrt(sum(meanif.dr^2)) < n_minus_12 & !sqrt_n_norm_block){
                sqrt_n_norm_iter <- iter
                sqrt_n_norm_block <- TRUE
            }

            if(sqrt(sum(meanif.dr^2)) < tolIF & !n_norm_block){
                n_norm_iter <- iter
                n_norm_block <- TRUE
            }
            
            if(max(abs(meanif.dr)) < tolIF & !n_max_block){
                n_max_iter <- iter
                n_max_block <- TRUE
            }

            # update iteration
            psin_trace[iter] <- mean(full_Qnstar$Q1n)
            se_trace[iter] <- sqrt(var(
                (if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2 -
                if.dr$Dg1.Q2 - if.dr$Dg0.Q1 - if.dr$DQ2.g1 - 
                if.dr$DQ2.g0 - if.dr$DQ1.g0)
            )/length(A0))

            cat(#"\n epsilon = ", etastar$flucmod$coefficients[1],
                "\n mean ic = ", meanif.dr, 
                 "\n")
        }
    }else{
        #----------------------------------------
        # targeting according to uniformly least
        # favorable submodels
        #----------------------------------------

        #------------------------
        # Work flow
        #------------------------
        # set epsilon = 0
        # while the || P_n D*(Q,g,Qr,gr) || > 1/n
        #    {
        #    set epsilon = epsilon + dx               ====> dx should be option for function
        #    set \bar{Q}_\epsilon = expit(logit(\bar{Q}_{\epsilon-dx}) - dx H^T P_n D*(Q_{\epsilon-dx}, G_{\epsilon -dx})
        #    estimate gr and Qr
        #    stop if || P_n D*(Q_{\epsilon}, G_{\epsilon}) || < 1/n
        #    }

        # what functions will be used?
        # universalStep 
        #    takes as input current Q, g, Qr, gr, dx, P_n D* (vector)
        #    returns Q_\epsilon 
        # redReg (calls estimateGr estimateQr)
        # evaluateIF
        # compute norm of P_n D*

        # initialize list that will hold the current step of the fluctuated nuisance parameters 
        # by setting equal to the initial values
        QnEps <- vector(mode = "list")
        gnEps <- vector(mode = "list")
        QnEps$Q2n <- Qn$Q2n
        QnEps$Q1n <- Qn$Q1n
        gnEps$g1n <- gn$g1n
        gnEps$g0n <- gn$g0n

        #-------------------------
        # evaluate IF
        #-------------------------
        if.dr <- evaluateIF(
            A0 = A0, A1 = A1, L2 = L2, Q2n = QnEps$Q2n, Q1n = QnEps$Q1n, 
            g1n = gnEps$g1n, g0n = gnEps$g0n, Q2nr.obsa = Qnr.gnr$Qnr$Q2nr.obsa, 
            Q1nr1 = Qnr.gnr$Qnr$Q1nr1, Q1nr2 = Qnr.gnr$Qnr$Q1nr2,
            g0nr = Qnr.gnr$gnr$g0nr, g1nr = Qnr.gnr$gnr$g1nr, 
            h0nr = Qnr.gnr$gnr$h0nr, h1nr = Qnr.gnr$gnr$h1nr, 
            hbarnr = Qnr.gnr$gnr$hbarnr, abar = abar
        )
        # P_n D*, P_n D_g, P_n D_Q
        PnDFull <- colMeans(Reduce("cbind",if.dr))
        PnDQ2 <- matrix(c(PnDFull[3], PnDFull[6] + PnDFull[7]), ncol = 1)
        PnDQ1 <- matrix(c(PnDFull[2],PnDFull[8]), ncol = 1)
        PnDg0 <- matrix(PnDFull[5], ncol = 1)
        PnDg1 <- matrix(PnDFull[4], ncol = 1)
        # compute norm
        normPnD <- as.numeric(sqrt(sum(PnDQ2^2) + sum(PnDQ1^2) + PnDg0^2 + PnDg1^2))
        iter <- 0 
        all.risks <- evaluateRisk(L2 = L2, A0 = A0, A1 = A1, Q2n = QnEps$Q2n,
                     Q1n = QnEps$Q1n, g1n = gnEps$g1n, g0n = gnEps$g0n,
                     abar = abar, tolg = tolg)
        risk <- all.risks$sum
        badSteps <- 0
        # start taking steps
        while(normPnD > tolIF){
            iter <- iter + 1
            #--------------------
            # take small step 
            #--------------------
            stepOut <- universalStep(
                L2 = L2, Qn = QnEps, gn = gnEps, Qnr.gnr = Qnr.gnr, 
                PnDQ2 = PnDQ2, PnDQ1 = PnDQ1, PnDg0 = PnDg0, PnDg1 = PnDg1,
                normPnD = normPnD, dx = universalStepSize, tolg = tolg, tolQ = tolQ 
            )
            QnEps <- list(Q2n = stepOut$Q2n, Q1n = stepOut$Q1n)
            gnEps <- list(g1n = stepOut$g1n, g0n = stepOut$g0n)

            #----------------------------
            # update reduced regressions
            #----------------------------
            Qnr.gnr <- redReg(A0 = A0, A1 = A1, L2 = L2, abar = abar, 
                              gn = gnEps, Qn = QnEps, verbose = verbose, tolg = tolg, 
                              SL.Qr = SL.Qr, SL.gr = SL.gr, return.models = return.models)

            #-------------------------
            # evaluate IF
            #-------------------------
            if.dr <- evaluateIF(
                A0 = A0, A1 = A1, L2 = L2, Q2n = QnEps$Q2n, Q1n = QnEps$Q1n, 
                g1n = gnEps$g1n, g0n = gnEps$g0n, Q2nr.obsa = Qnr.gnr$Qnr$Q2nr.obsa, 
                Q1nr1 = Qnr.gnr$Qnr$Q1nr1, Q1nr2 = Qnr.gnr$Qnr$Q1nr2, 
                g0nr = Qnr.gnr$gnr$g0nr, g1nr = Qnr.gnr$gnr$g1nr, 
                h0nr = Qnr.gnr$gnr$h0nr, h1nr = Qnr.gnr$gnr$h1nr, 
                hbarnr = Qnr.gnr$gnr$hbarnr, abar = abar
            )
            PnDFull <- colMeans(Reduce("cbind",if.dr))
            PnDQ2 <- matrix(c(PnDFull[3], PnDFull[6]+PnDFull[7]), ncol = 1)
            PnDQ1 <- matrix(c(PnDFull[2],PnDFull[8]), ncol = 1)
            PnDg0 <- matrix(PnDFull[5], ncol = 1)
            PnDg1 <- matrix(PnDFull[4], ncol = 1)

            # compute norm
            normPnD <- as.numeric(sqrt(sum(PnDQ2^2) + sum(PnDQ1^2) + PnDg0^2 + PnDg1^2))
        
            # compute empirical risk
            tmp.risk <- evaluateRisk(L2 = L2, A0 = A0, A1 = A1, Q2n = QnEps$Q2n,
                                     Q1n = QnEps$Q1n, g1n = gnEps$g1n, g0n = gnEps$g0n,
                                     abar = abar, tolg = tolg)
            warn <- tmp.risk$sum > risk
            if(warn){
                diff <- mapply(t = tmp.risk, a = all.risks, function(t,a){
                    a - t
                })
            }
            risk <- tmp.risk$sum
            all.risks <- tmp.risk 
            if(iter%%printFreq == 0){
                cat("\n epsilon = ", iter*universalStepSize, "\n norm IC = ", normPnD, 
                    "\n risk = ", risk,"\n estimate = ", mean(QnEps$Q1n),"\n")
                if(warn){
                    cat("!!! WARNING - INCREASING RISK !!!")
                    badSteps <- badSteps + 1
                }
            }
        } # end while loop
        # assign values at approximate MLE to Qnstar and gnstar
        Qnstar <- QnEps; gnstar <- gnEps 
    } # end if universal

    #------------------------------------------
    # evaluate parameter and compute std err
    #------------------------------------------
    # evaluate parameter
    psin <- mean(full_Qnstar$Q1n)
       
    # compute standard errors
    se <- sqrt(var(
            (if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2 -
            if.dr$Dg1.Q2 - if.dr$Dg0.Q1 - if.dr$DQ2.g1 - 
            if.dr$DQ2.g0 - if.dr$DQ1.g0)
    )/length(A0))
    }
    #-----------------------------------
    # Computing the regular LTMLE
    #-----------------------------------
    if(return.ltmle){
        Qnstar_list.ltmle <- Qn_list
        full_Qn.ltmle <- do.call(Map, c(c, Qn_list))
        full_gn <- do.call(Map, c(c, gn_list))

        Q2nstar.ltmle <- targetQ2.ltmle(
            A0 = A0, A1 = A1, L2 = L2, Qn = full_Qn.ltmle, gn = full_gn, 
            abar = abar, tolg = tolg, tolQ = tolQ, return.models = return.models
        )$Q2nstar

        # replace values
        split_Q2nstar.ltmle <- split(Q2nstar.ltmle, folds)
        Qnstar_list.ltmle <- mapply(x = Qnstar_list.ltmle, values = split_Q2nstar.ltmle, 
                         FUN = list_replace, MoreArgs = list(list = "Q2n"),
                         SIMPLIFY = FALSE)
        full_Qnstar.ltmle <- do.call(Map, c(c, Qnstar_list.ltmle))

        Q1n.ltmle <- sapply(1:cvFolds, estimateQ1.ltmle,
                        L2 = L2, L0 = L0, Q2n = full_Qnstar.ltmle$Q2n, 
                        A0 = A0, A1 = A1, folds = folds, abar = abar, 
                        SL.Q = SL.Q, SL.Q.options = SL.Q.options, 
                        glm.Q = glm.Q, glm.Q.options = glm.Q.options, 
                        return.models = return.models, verbose = verbose, 
                        stratify = stratify, simplify = FALSE)
        # replace values
        Qnstar_list.ltmle <- mapply(x = Qnstar_list.ltmle, 
                                    values = lapply(Q1n.ltmle, "[[", "Q1n"), 
                 FUN = list_replace, MoreArgs = list(list = "Q1n"),
                 SIMPLIFY = FALSE)
        # collapse list 
        full_Qnstar.ltmle <- do.call(Map, c(c, Qnstar_list.ltmle))

        Q1nstar.ltmle <- targetQ1.ltmle(
            A0 = A0, A1 = A1, L2 = L2, 
            Qn = full_Qnstar.ltmle, 
            gn = full_gn, abar = abar, tolg = tolg, tolQ = tolQ, 
            return.models = return.models
        )$Q1nstar
        
        # point estimate
        psin.ltmle <- mean(Q1nstar.ltmle)
        
        # influence functions
        if.ltmle <- evaluateEIF(
            A0 = A0, A1 = A1, L2 = L2, Q2n = full_Qnstar.ltmle$Q2n, 
            Q1n = full_Qnstar.ltmle$Q1n, g1n = full_gn$g1n, g0n = full_gn$g0n, 
            abar = abar
        )
        
        # se
        se.ltmle <- sqrt(var(
            (if.ltmle$Dstar2 + if.ltmle$Dstar1 + if.ltmle$Dstar0)
        )/length(A0))
        
    } # end if(return.ltmle)
    
    #----------------------------
    # formatting output
    #----------------------------
    out <- vector(mode = "list")
    
    # ltmle output
    out$est.ltmle <- NULL
    out$se.ltmle <- NULL
    if(return.ltmle){
        out$est.ltmle <- psin.ltmle
        out$se.ltmle <- se.ltmle
    }
    # dr tmle output    
    if(!only.ltmle){
    out$est <- psin
    out$se <- as.numeric(se)

    # replace NA's in psin_trace
    n_notNA_psin_trace <- sum(!is.na(psin_trace))
    if(n_notNA_psin_trace < length(psin_trace)){
        psin_trace[(n_notNA_psin_trace+1):length(psin_trace)] <- psin_trace[n_notNA_psin_trace]
        se_trace[(n_notNA_psin_trace+1):length(psin_trace)] <- se_trace[n_notNA_psin_trace]
    }

    out$est_trace <- psin_trace
    out$se_trace <- se_trace

    # naive output
    out$est.naive <- NULL
    if(return.ltmle){
        out$est.naive <- est.naive
    }
    
    # number of iterations
    out$iter <- iter
    out$sqrt_n_max_iter <- sqrt_n_max_iter
    out$n_max_iter <- n_max_iter
    out$sqrt_n_norm_iter <- sqrt_n_norm_iter
    out$n_norm_iter <- n_norm_iter

    if(universal) out$badSteps <- badSteps
    # model output
    out$Qmod <- vector(mode = "list")
    out$gmod <- vector(mode = "list")
    out$Qrmod <- vector(mode = "list")
    out$grmod <- vector(mode = "list")
    
    # add on influence functions
    out$ic <- list(orig = mean(if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2),
                   missQ = mean(if.dr$Dg1.Q2 + if.dr$Dg0.Q1),
                   missg = mean(if.dr$DQ2.g1 + if.dr$DQ2.g0 + if.dr$DQ1.g0))

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
    }
    return(out)
}

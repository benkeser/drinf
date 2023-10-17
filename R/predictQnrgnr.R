predict_Qnrgnr <- function(Q2n, Q1n, g0n, g1n, 
                           A0, A1, L2, tolg, 
                           SL.Qr, SL.gr,
                           abar, verbose,
                           Qnr.gnr, ...){
    multiAlgos <- (length(SL.Qr) > 1 | is.list(SL.Qr))
    
    # empty vector of predictions
    Q2nr.obsa <- rep(NA, length(A0))
    # zeros for these obs. 
    Q2nr.obsa[A0 != abar[1]] <- 0

    if(multiAlgos){
        weightfail <- all(Q2rmod$coef==0)
        if(!weightfail){
            # Super Learner predictions for A0==abar[1] obs. 
            Q2nr.obsa[A0 == abar[1]] <- predict(Qnr.gnr$Qnr$Q2rmod, newdata = data.frame(g1n = g1n[A0 == abar[1]]))[[1]]
            # get prediction also setting A0 = abar[1] that is needed later
            # for targeting step
            Q2nr.seta <- predict(Qnr.gnr$Qnr$Q2rmod, newdata = data.frame(g1n = g1n))[[1]]
        }else{
            # find dsl
            dslcol <- which(Qnr.gnr$Qnr$Q2rmod$cvRisk == min(Q2rmod$cvRisk, na.rm = TRUE))
            # use discrete Super Learner predictions
            Q2nr.obsa[A0 == abar[1]] <- predict(Qnr.gnr$Qnr$Q2rmod, newdata = data.frame(g1n = g1n[A0 == abar[1]]))[[2]][,dslcol]
            Q2nr.seta <- predict(Qnr.gnr$Qnr$Q2rmod, newdata = data.frame(g1n = g1n))[[2]][,dslcol]
        }
    }else{
        Q2nr.obsa[A0 == abar[1]] <- predict(Qnr.gnr$Qnr$Q2rmod$fit, newdata = data.frame(g1n = g1n[A0 == abar[1]]))
        Q2nr.seta <- predict(Qnr.gnr$Qnr$Q2rmod$fit, newdata = data.frame(g1n = g1n))
    }

    # Q1r1
    Q1nr1 <- predict_from_fit(fit = Qnr.gnr$Qnr$Q1r1mod,
                              covariate = g0n, covariate_name = "g0n",
                              multiAlgos = multiAlgos)
    Q1nr2 <- predict_from_fit(fit = Qnr.gnr$Qnr$Q1r2mod,
                              covariate = g0n, covariate_name = "g0n",
                              multiAlgos = multiAlgos)
    g0nr <- predict_from_fit(fit = Qnr.gnr$gnr$g0rmod, 
                             covariate = Q1n, covariate_name = "Q1n",
                             multiAlgos = multiAlgos, truncate = tolg)
    g1nr <- predict_from_fit(fit = Qnr.gnr$gnr$g1rmod, 
                             covariate = Q2n, covariate_name = "Q2n",
                             multiAlgos = multiAlgos, truncate = tolg)
    h0nr <- predict_from_fit(fit = Qnr.gnr$gnr$h0rmod, 
                             covariate = Q1n, covariate_name = "Q1n",
                             multiAlgos = multiAlgos)
    h1nr <- predict_from_fit(fit = Qnr.gnr$gnr$h1rmod, 
                             covariate = Q2n, covariate_name = "Q2n",
                             multiAlgos = multiAlgos)
    hbarnr <- predict_from_fit(fit = Qnr.gnr$gnr$hbarrmod, 
                             covariate = Q2n, covariate_name = "Q2n",
                             multiAlgos = multiAlgos)

    return(list(
        Qnr = list(Q2nr.obsa = Q2nr.obsa, Q2nr.seta = Q2nr.seta, Q1nr1 = Q1nr1,
                   Q1nr2 = Q1nr2, Q2rmod = Qnr.gnr$Qnr$Q2rmod, Q1r1mod = Qnr.gnr$Qnr$Q1r1mod,
                   Q1r2mod = Qnr.gnr$Qnr$Q1r2mod),
        gnr = list(g0nr = g0nr, g1nr = g1nr, h0nr = h0nr, h1nr = h1nr,
                hbarnr = hbarnr, g0rmod = Qnr.gnr$gnr$g0rmod, g1rmod = Qnr.gnr$gnr$g1rmod,
                h0rmod = Qnr.gnr$gnr$h0rmod, h1rmod = Qnr.gnr$gnr$h1rmod,
                hbarrmod = Qnr.gnr$gnr$hbarrmod)
    ))
}



predict_from_fit <- function(fit, covariate, 
                                 covariate_name, multiAlgos, truncate = NULL){
    newdata <- data.frame(covariate)
    colnames(newdata) <- covariate_name
    if(multiAlgos){
        weightfail <- all(fit==0)
        if(!weightfail){
            # Super Learner predictions 
            pred <- predict(fit, newdata = newdata)[[1]]
        }else{
            dslcol <- which.min(fit$cvRisk, na.rm = TRUE)
            pred <- predict(fit, newdata = newdata)[[2]][, dslcol]
        }
        }else{
            pred <- predict(fit$fit, newdata = newdata)
    }
    if(!is.null(truncate)){
        pred[pred < truncate] <- truncate
    }
    return(pred)
}

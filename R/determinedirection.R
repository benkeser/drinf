determine_direction <- function(L2, A0, A1, QnEps, gnEps, Qnr.gnr, 
                PnDQ2, PnDQ1, PnDg0, PnDg1,
                normPnD, tolg, tolQ, universalStepSize, abar){
	# go negative
    stepOut <- universalStep(
        L2 = L2, Qn = QnEps, gn = gnEps, Qnr.gnr = Qnr.gnr, 
        PnDQ2 = PnDQ2, PnDQ1 = PnDQ1, PnDg0 = PnDg0, PnDg1 = PnDg1,
        normPnD = normPnD, dx = -universalStepSize, tolg = tolg, tolQ = tolQ 
    )
    QnEps <- list(Q2n = stepOut$Q2n, Q1n = stepOut$Q1n)
    gnEps <- list(g1n = stepOut$g1n, g0n = stepOut$g0n)

    # risk
    risk.neg <- evaluateRisk(L2 = L2, A0 = A0, A1 = A1, Q2n = QnEps$Q2n,
                             Q1n = QnEps$Q1n, g1n = gnEps$g1n, g0n = gnEps$g0n,
                             abar = abar, tolg = tolg)
    # go postive
    stepOut <- universalStep(
        L2 = L2, Qn = QnEps, gn = gnEps, Qnr.gnr = Qnr.gnr, 
        PnDQ2 = PnDQ2, PnDQ1 = PnDQ1, PnDg0 = PnDg0, PnDg1 = PnDg1,
        normPnD = normPnD, dx = universalStepSize, tolg = tolg, tolQ = tolQ 
    )
    QnEps <- list(Q2n = stepOut$Q2n, Q1n = stepOut$Q1n)
    gnEps <- list(g1n = stepOut$g1n, g0n = stepOut$g0n)

    # risk
    risk.pos <- evaluateRisk(L2 = L2, A0 = A0, A1 = A1, Q2n = QnEps$Q2n,
                             Q1n = QnEps$Q1n, g1n = gnEps$g1n, g0n = gnEps$g0n,
                             abar = abar, tolg = tolg)
    # if negative better risk, return -1,
    # if positive better risk, return 1
    return(ifelse(risk.neg$sum < risk.pos$sum, -1, 1))
}

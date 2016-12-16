#' universalStep
#' 
#' Take one step in the practical construction of the universally 
#' least favorable submodel. 
#' 
#' @param Qn A \code{list} of current estimates of Q2n and Q1n
#' @param gn A \code{list} of current estimates of g0n and g1n
#' @param Qnr.gnr A \code{list} of current estimates of reduced dim. regressions
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities.
#' @param tolQ A \code{numeric} indicating the truncation level for transformed outcome regressions.
#' @param PnD A \code{matrix} of dimension 3x1 that is the empirical mean of the relevant
#' portions of the influence function (one for original EIF, one for extra G terms, 
#' one for extra Q terms). 
#' @param normPnD A \code{numeric} giving the norm of PnD
#' @param dx A \code{numeric} step size


universalStep <- function(
	L2, Qn, gn, Qnr.gnr, 
	PnDQ2, PnDQ1, PnDg1, PnDg0, normPnD, dx = 0.005, tolg, tolQ, ...
){

	# compute min and max of L2 for scaling
    L2.min <- min(L2); L2.max <- max(L2)
    # scale Q2n and Q1n
    Q2ns <- (Qn$Q2n - L2.min)/(L2.max - L2.min)
	Q1ns <- (Qn$Q1n - L2.min)/(L2.max - L2.min)
    
	# compute vector of clever covariates
	HQ2 <- cbind(
		(L2.max - L2.min) / (gn$g0n * gn$g1n),
		(L2.max - L2.min) * (Qnr.gnr$gnr$hbarnr + Qnr.gnr$gnr$h1nr)/Qnr.gnr$gnr$g1nr
	)
	HQ1 <- cbind(
		(L2.max - L2.min) / gn$g0n,
		(L2.max - L2.min) * Qnr.gnr$gnr$h0nr/Qnr.gnr$gnr$g0nr
	)
	Hg1 <- cbind(
		Qnr.gnr$Qnr$Q2nr.seta / gn$g1n
	)
	Hg0 <- cbind(
		Qnr.gnr$Qnr$Q1nr / gn$g0n
	)

	# take a step for Q2n 
	Q2nEps <- (L2.max - L2.min)*plogis(SuperLearner::trimLogit(Q2ns,trim=tolQ) + dx*(HQ2%*%PnDQ2)/normPnD) + L2.min
	Q1nEps <- (L2.max - L2.min)*plogis(SuperLearner::trimLogit(Q1ns,trim=tolQ) + dx*(HQ1%*%PnDQ1)/normPnD) + L2.min
	g0nEps <- plogis(SuperLearner::trimLogit(gn$g0n,trim=tolg) + dx*(Hg0%*%PnDg0)/normPnD)
	g1nEps <- plogis(SuperLearner::trimLogit(gn$g1n,trim=tolg) + dx*(Hg1%*%PnDg1)/normPnD)
	g0nEps[g0nEps < tolg] <- tolg
	g0nEps[g0nEps > 1-tolg] <- 1-tolg
	g1nEps[g1nEps < tolg] <- tolg
	g1nEps[g1nEps > 1-tolg] <- 1-tolg

	return(list(
		Q2n = Q2nEps, Q1n = Q1nEps, g0n = g0nEps, g1n = g1nEps
	))
}
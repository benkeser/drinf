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
#' Indicates whether to guard against misspecification of outcome or treatment regressions or both. 
#' @param returnModels A \code{boolean} indicating whether the models for Q, g, Qr, and gr should be 
#' returned with the output. 
#' @param maxIter A \code{numeric} indicating the maximum number of TMLE iterations before stopping. 
#' @param tolEps A \code{numeric} stopping criteria for the TMLE updates based on the value of the fluctuation 
#' parameter of the parametric submodels. 
#' @param tolIC A \code{numeric} stopping criteria for the TMLE updates based on the empirical average of the 
#' estimated influence curve. 
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities. 
#' @param verbose A \code{boolean} indicating whether messages should be printed to indicate progress.
#' @param Qsteps A \code{numeric} equal to 1 or 2, indicating whether the fluctuation to solve score equations 
#' for submodels with iterated means as outcomes should be performed in a single or multiple steps. 
#' @param SL.Q.options A \code{list} of additional arguments passed to \code{SuperLearner} for outcome
#' regression fits.
#' @param SL.g.options A \code{list} of additional arguments passed to \code{SuperLearner} for condtional treatment 
#' probability fits.
#' @param glm.Q.options A \code{list} of additional arguments passed to \code{glm} for outcome
#' regression fits.
#' 
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
                       reduction="univariate",
                       returnModels=FALSE,
                       maxIter=100,
                       tolEps=1e-4, 
                       tolIC=1/(10*sqrt(length(L2))), 
                       tolg=1e-8,
                       verbose=TRUE,
                       Qsteps=2,
                       SL.Q.options = list(family = gaussian()),
                       SL.g.options = list(family = binomial()),
                       glm.Q.options = list(family = gaussian()),
                      ...){
    
    #-------------------------------------#
    # Workflow 
    #-------------------------------------#
    # ++ Estimate g
    # ++ Estimate Q
    # ++ Compute outcomes of Qr regressions
    # ++ Estimate Qr regressions
    # ++ Estimate gr regression
    # ++ - internally, compute outcomes of iterated gr regressions
    # Evaluate IC
    # while Pn D* < tolIC {
    #    Target Q
    #    Compute outcomes of Qr regressions
    #    Estimate Qr
    #    Target g
    # Estimate gr regression
    #  - internally, compute outcomes of iterated gr regressions
    #    Evaluate IC
    # }
    # Evaluate parameter
  
}

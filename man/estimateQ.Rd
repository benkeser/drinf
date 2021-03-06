% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateQ.R
\name{estimateQ}
\alias{estimateQ}
\title{estimateQ}
\usage{
estimateQ(validFold, folds, L0, L1, L2, A0, A1, abar, SL.Q, SL.Q.options, glm.Q,
  glm.Q.options, return.models, verbose, stratify, ...)
}
\arguments{
\item{L0}{A \code{data.frame} featuring covariates measured at baseline.}

\item{L1}{A \code{data.frame} featuring time-varying covariates measured at 
the first timepoint.}

\item{L2}{A \code{vector} outcome of interest}

\item{A0}{A \code{vector} treatment delivered at baseline.}

\item{A1}{A \code{vector} treatment deliver after \code{L1} is measured.}

\item{abar}{A \code{vector} of length 2 indicating the treatment assignment 
that is of interest.}

\item{SL.Q}{A \code{vector} or \code{list} specifying the SuperLearner library
to be used to estimate the conditional probability of treatment at each time point.  
See \code{SuperLearner} package for details.}

\item{SL.Q.options}{A \code{list} of additional arguments passed to \code{SuperLearner} 
for outcome regression fits.}

\item{glm.Q}{A \code{character} specifying the right-hand side of the \code{glm} 
formula used to estimate the conditional probability of treatment at each time point. 
Only used if \code{SL.Q = NULL}.}

\item{glm.Q.options}{A \code{list} of additional arguments passed to \code{glm} for the
outcome regression fits. Typically, the \code{family} argument.}

\item{return.models}{A \code{boolean} indicating whether the models for g00 should be 
returned with the output.}

\item{stratify}{A \code{boolean} indicating whether to pool across treatment
nodes or to estimate outcome regression separately in each category.}

\item{...}{Other arguments (not currently used).}
}
\value{
Returns a list with \code{Q2n}, \code{Q1n}, and the estimated model objects if
\code{return.models = TRUE}
}
\description{
This function computes the conditional treatment probabilities at both timepoints.
}
\examples{
TO DO : add examples
}

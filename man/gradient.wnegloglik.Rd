% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wnegloglik.R
\name{gradient.wnegloglik}
\alias{gradient.wnegloglik}
\title{gradient.wnegloglik}
\usage{
gradient.wnegloglik(epsilon, weight, Y, offset)
}
\arguments{
\item{epsilon}{The scalar parameter of the fluctuation submodel.}

\item{weight}{The \code{vector} of weights, i.e., the clever covariates.}

\item{Y}{The \code{vector} of regression outcomes.}

\item{offset}{The \code{vector} of offsets.}
}
\value{
A \code{vector} of the gradient of the loss evaluated at the 
data inputs.
}
\description{
gradient.wnegloglik
}

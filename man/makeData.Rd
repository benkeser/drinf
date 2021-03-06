% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newSim.R
\name{makeData}
\alias{makeData}
\title{makeData}
\usage{
makeData(n, b = 0.5, b0 = 1, ba = -0.25, rho = 0.05, setA = NULL)
}
\value{
A \code{list} with L0, A0, L1, A1, L2 to input into \code{drinf.tmle}.
}
\description{
A function to simulate data for use in the simulation study. The function
can be used to generate observed data or data from an intervened SCM; the latter
is achieved via the \code{setA} option.
}

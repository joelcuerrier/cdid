% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMM_convert_result.R
\name{gmm_convert_result}
\alias{gmm_convert_result}
\title{GMM_convert_result}
\usage{
gmm_convert_result(dp, type)
}
\arguments{
\item{dp}{a dp object}

\item{type}{1 for 2step weighting, 2 for identity weighting}
}
\value{
a \code{\link{DIDparams}} object
}
\description{
Function to convert results so they can be used by the did package developed by Brantly Callaway. For more details on the methodology, see:
Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences",
Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
}

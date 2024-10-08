% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fonction_simu_attrition.R
\docType{data}
\name{fonction_simu_attrition}
\alias{fonction_simu_attrition}
\title{Cross-sectional data simulation}
\format{
data: a data frame
}
\usage{
data=data_sim=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
}
\description{
This function generates a simulated dataset with attrition. The dataset is used to estimate the treatment effect using the chained method. See the simulation in "The Chained Difference-in-Differences" paper.
}
\examples{
data=data_sim=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
}
\keyword{datasets}

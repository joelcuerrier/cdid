% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DIDparams.R
\name{DIDparams}
\alias{DIDparams}
\title{DIDparams}
\usage{
DIDparams(
  yname,
  tname,
  idname = NULL,
  gname,
  xformla = NULL,
  data,
  control_group,
  anticipation = 0,
  weightsname = NULL,
  alp = 0.05,
  bstrap = TRUE,
  biters = 1000,
  clustervars = NULL,
  cband = TRUE,
  print_details = TRUE,
  pl = FALSE,
  cores = 1,
  est_method = "chained",
  base_period = "varying",
  panel = TRUE,
  true_repeated_cross_sections,
  n = NULL,
  nG = NULL,
  nT = NULL,
  tlist = NULL,
  glist = NULL,
  call = NULL
)
}
\arguments{
\item{true_repeated_cross_sections}{Whether or not the data really
is repeated cross sections.  (We include this because unbalanced
panel code runs through the repeated cross sections code)}

\item{n}{The number of observations.  This is equal to the
number of units (which may be different from the number
of rows in a panel dataset).}

\item{nG}{The number of groups}

\item{nT}{The number of time periods}

\item{tlist}{a vector containing each time period}

\item{glist}{a vector containing each group}

\item{call}{Function call to att_gt}
}
\description{
Object to hold did parameters that are passed across functions
}

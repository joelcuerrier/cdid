# Adapted from Callaway, B. (2024). The did Library. Department of Economics, University of Georgia. Available at: https://github.com/bcallaway11/did

#' @title MP
#' @importFrom methods is
#' @description Multi-period objects that hold results for group-time average treatment effects
#'
#' @param group which group (defined by period first treated) an group-time average treatment effect is for
#' @param t which time period a group-time average treatment effect is for
#' @param att the group-average treatment effect for group \code{group} and time
#'  period \code{t}
#' @param c simultaneous critical value if one is obtaining simultaneous confidence
#'  bands. Otherwise it reports the critical value based on pointwise normal
#'  approximation.
#' @param V_analytical Analytical estimator for the asymptotic variance-covariance matrix for group-time average treatment effects
#' @param se standard errors for group-time average treatment effects. If bootstrap is set to TRUE, this provides bootstrap-based se.
#' @param inffunc the influence function for estimating group-time average treatment effects
#' @param n the number of unique cross-sectional units (unique values of idname)
#' @param W the Wald statistic for pre-testing the common trends assumption
#' @param Wpval the p-value of the Wald statistic for pre-testing the
#'  common trends assumption
#' @param aggte an aggregate treatment effects object
#' @param alp the significance level, default is 0.05
#' @param debT first time period
#' @param DIDparams a [`DIDparams`] object.
#'
#' @return MP object
#' @export
MP <- function(group, t, att, V_analytical, se, c, inffunc, n=NULL, W=NULL, Wpval=NULL, aggte=NULL, alp = 0.05, DIDparams=NULL,debT) {

  out <- list(group=group, t=t, att=att, V_analytical=V_analytical, se=se, c=c,
  inffunc=inffunc, n=n, W=W, Wpval=Wpval, aggte=aggte, alp = alp,
  DIDparams=DIDparams, call=DIDparams$call)

  class(out) <- "MP"
  return(out)
}

#' @title summary.MP
#'
#' @description Prints a detailed summary of an \code{MP} object. The function outputs key details of
#' the group-time average treatment effects, such as estimation method, control group,
#' and pre-test results for parallel trends.
#'
#' @param object An \code{MP} object, representing the results of a multi-period analysis.
#' @param ... Additional arguments passed to the function.
#'
#' @return No return value. This function is called for its side effects of
#' printing a summary of the \code{MP} object to the console, including:
#' \itemize{
#'   \item Call: The call used to create the \code{MP} object.
#'   \item Group-Time Average Treatment Effects: A table of estimates with confidence bands.
#'   \item Control Group: Information about the chosen control group (e.g., "Never Treated").
#'   \item Anticipation Periods: Number of periods used to account for anticipation effects.
#'   \item Estimation Method: Method used for treatment effect estimation.
#'   \item Pre-Test Results: p-values for the test of parallel trends assumption, if available.
#' }
#'
#' @seealso \code{\link{MP}}, \code{\link{print.MP}}
#'
#' @export
summary.MP <- function(object, ...) {
  mpobj <- object

  # call
  cat("\n")
  cat("Call:\n")
  print(mpobj$DIDparams$call)
  cat("\n")

  # citation
 cat("Bellego C., Benatia D., and V. Dortet-Bernadet (2024). \"The Chained Difference-in-Differences.\" Journal of Econometrics. doi:10.1016/j.jeconom.2023.11.002.")
  cat("\n")

  # group time average treatment effects
  cat("Group-Time Average Treatment Effects:\n")

  cband_text1a <- paste0(100*(1-mpobj$alp),"% ")
  cband_text1b <- ifelse(mpobj$DIDparams$bstrap,
                         ifelse(mpobj$DIDparams$cband, "Simult. ", "Pointwise "),
                         "Pointwise ")
  cband_text1 <- paste0("[", cband_text1a, cband_text1b)

  cband_lower <- mpobj$att - mpobj$c*mpobj$se
  cband_upper <- mpobj$att + mpobj$c*mpobj$se

  sig <- (cband_upper < 0) | (cband_lower > 0)
  sig[is.na(sig)] <- FALSE
  sig_text <- ifelse(sig, "*", "")

  out <- cbind.data.frame(mpobj$group, mpobj$t, mpobj$att, mpobj$se, cband_lower, cband_upper)
  out <- round(out,4)
  out <- cbind.data.frame(out, sig_text)


  colnames(out) <- c("Group", "Time", "ATT(g,t)","Std. Error", cband_text1, "Conf. Band]", "")
  print(out, row.names=FALSE)
  cat("---\n")
  cat("Signif. codes: `*' confidence band does not cover 0")
  cat("\n\n")

  # report pre-test
  if (!is.null(mpobj$Wpval)) {
    cat("P-value for pre-test of parallel trends assumption:  ")
    cat(as.character(mpobj$Wpval))
    cat("\n")
  }

  # set control group text
  control_group <- mpobj$DIDparams$control_group
  control_group_text <- NULL
  if (control_group == "nevertreated") {
    control_group_text <- "Never Treated"
  } else if (control_group == "notyettreated") {
    control_group_text <- "Not Yet Treated"
  }

  if (!is.null(control_group)) {
    cat("Control Group:  ")
    cat(control_group_text)
    cat(",  ")
  }

  # anticipation periods
  cat("Anticipation Periods:  ")
  cat(mpobj$DIDparams$anticipation)
  cat("\n")

  # estimation method text
  est_method <- mpobj$DIDparams$est_method
  if ( is(est_method,"character") ) {
    est_method_text <- est_method
    if (est_method == "dr") {
      est_method_text <- "Doubly Robust"
    } else if (est_method == "ipw") {
      est_method_text <- "Inverse Probability Weighting"
    } else if (est_method == "reg") {
      est_method_text <- "Outcome Regression"
    }

    cat("Estimation Method:  ")
    cat(est_method_text)
    cat("\n")
  }
}

#' @title print.MP
#'
#' @description Prints a summary of the results contained in an \code{MP} object.
#' This function calls \code{summary.MP} to display the details of the multi-period
#' analysis results in a user-friendly format.
#'
#' @param x An \code{MP} object, representing the results of multi-period analysis.
#' @param ... Additional arguments passed to \code{summary.MP}.
#'
#' @return No return value. This function is called for its side effects of
#' printing the summary of the \code{MP} object to the console.
#'
#' @seealso \code{\link{summary.MP}}
#'
#' @export
print.MP <- function(x, ...) {
  summary.MP(x, ...)
}


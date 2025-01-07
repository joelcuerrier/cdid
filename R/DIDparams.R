# Callaway, B. (2024). The did Library. Department of Economics, University of Georgia. Available at: https://github.com/bcallaway11/did

#' @title DIDparams
#'
#' @description Creates a \code{DIDparams} object to hold parameters for difference-in-differences analysis,
#' including data structure details and user-specified options. This object is designed to streamline
#' parameter passing across functions in the `cdid` package.#'
#'
#' @inheritParams pre_process_cdid
#' @param n The number of observations.  This is equal to the
#'  number of units (which may be different from the number
#'  of rows in a panel dataset).
#' @param nG The number of groups
#' @param nT The number of time periods
#' @param tlist a vector containing each time period
#' @param glist a vector containing each group
#' @param true_repeated_cross_sections Whether or not the data really
#'  is repeated cross sections.  (We include this because unbalanced
#'  panel code runs through the repeated cross sections code)
#'
#' @return A \code{DIDparams} object, which is a list containing the following elements:
#' \itemize{
#'   \item \code{yname}: The name of the outcome variable.
#'   \item \code{tname}: The name of the time variable.
#'   \item \code{idname}: The name of the unit identifier variable (if applicable).
#'   \item \code{gname}: The name of the group variable (e.g., treatment group).
#'   \item \code{xformla}: A formula specifying covariates for the model.
#'   \item \code{data}: The dataset used for analysis.
#'   \item \code{control_group}: The type of control group (e.g., "never treated" or "not yet treated").
#'   \item \code{anticipation}: The number of periods of anticipation before treatment.
#'   \item \code{weightsname}: The name of the variable containing sampling weights (if applicable).
#'   \item \code{alp}: The significance level (default is 0.05).
#'   \item \code{bstrap}: Logical. Indicates whether bootstrap is used for standard errors.
#'   \item \code{biters}: The number of bootstrap iterations (if bootstrap is enabled).
#'   \item \code{clustervars}: Variables used for clustering standard errors.
#'   \item \code{cband}: Logical. Indicates whether simultaneous confidence bands are computed.
#'   \item \code{print_details}: Logical. Indicates whether detailed results should be printed.
#'   \item \code{pl}: Logical. Parallelization flag for computations.
#'   \item \code{cores}: The number of cores to use for parallelization (if enabled).
#'   \item \code{est_method}: The estimation method (e.g., "chained").
#'   \item \code{base_period}: The base period used for comparison (e.g., "varying").
#'   \item \code{panel}: Logical. Indicates whether the data is a panel dataset.
#'   \item \code{true_repeated_cross_sections}: Logical. Indicates whether the data is truly repeated cross-sections.
#'   \item \code{n}: The number of observations (units).
#'   \item \code{nG}: The number of groups.
#'   \item \code{nT}: The number of time periods.
#'   \item \code{tlist}: A vector containing all time periods.
#'   \item \code{glist}: A vector containing all groups.
#'   \item \code{call}: The call that generated the \code{DIDparams} object.
#' }
#'
#' @seealso \code{\link{pre_process_cdid}}
#'
#' @export
DIDparams <- function(yname,
                   tname,
                   idname=NULL,
                   gname,
                   xformla=NULL,
                   data,
                   control_group,
                   anticipation=0,
                   weightsname=NULL,
                   alp=0.05,
                   bstrap=TRUE,
                   biters=1000,
                   clustervars=NULL,
                   cband=TRUE,
                   print_details=TRUE,
                   pl=FALSE,
                   cores=1,
                   est_method="chained",
                   base_period="varying",
                   panel=TRUE,
                   true_repeated_cross_sections,
                   n=NULL,
                   nG=NULL,
                   nT=NULL,
                   tlist=NULL,
                   glist=NULL,
                   call=NULL) {

  out <- list(yname=yname,
              tname=tname,
              idname=idname,
              gname=gname,
              xformla=xformla,
              data=data,
              control_group=control_group,
              anticipation=anticipation,
              weightsname=weightsname,
              alp=alp,
              bstrap=bstrap,
              biters=biters,
              clustervars=clustervars,
              cband=cband,
              print_details=print_details,
              pl=pl,
              cores=cores,
              est_method=est_method,
              base_period=base_period,
              panel=panel,
              true_repeated_cross_sections=true_repeated_cross_sections,
              n=n,
              nG=nG,
              nT=nT,
              tlist=tlist,
              glist=glist,
              call=call)
  class(out) <- "DIDparams"
  out
}

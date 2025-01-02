#' @title att_gt_cdid
#' @description `att_gt_cdid` computes average treatment effects.
#' Our estimator accommodates (1) multiple time
#' periods, (2) variation in treatment timing, (3) treatment effect heterogeneity,
#' and (4) general missing data patterns. For more details on the methodology, see:
#' Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences",
#' Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
#' @param yname The name of the outcome variable
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param idname The individual (cross-sectional unit) id name
#' @param gname The name of the variable in `data` that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
#' @param weightsname The name of the column containing weights.
#'  If not set, all observations have same weight.
#' @param alp the significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using
#'  the multiplier bootstrap.  If standard errors are clustered, then one
#'  must set `bstrap=TRUE`. Default is `TRUE` (in addition, cband
#'  is also by default `TRUE` indicating that uniform confidence bands
#'  will be returned.  If bstrap is `FALSE`, then analytical
#'  standard errors are reported.
#' @param biters The number of bootstrap iterations to use.  The default is 1000,
#'  and this is only applicable if `bstrap=TRUE`.
#' @param clustervars A vector of variables names to cluster on.  At most, there
#'  can be two variables (otherwise will throw an error) and one of these
#'  must be the same as idname which allows for clustering at the individual
#'  level. By default, we cluster at individual level (when `bstrap=TRUE`).
#' @param cband Boolean for whether or not to compute a uniform confidence
#'  band that covers all of the group-time average treatment effects
#'  with fixed probability `1-alp`.  In order to compute uniform confidence
#'  bands, `bstrap` must also be set to `TRUE`.  The default is
#' `TRUE`.
#' @param print_details Whether or not to show details/progress of computations.
#'   Default is `FALSE`.
#' @param pl Whether or not to use parallel processing
#' @param cores The number of cores to use for parallel processing
#' @param est_method the method to compute group-time average treatment effects.  At the moment, one can only use the IPW estimator
#' with either "2-step" or "Identity" weighting matrix to aggregate Delta ATT into ATT.
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.
#' @param xformla A formula for the covariates to include in the
#'  model.  It should be of the form `~ X1 + X2`.  Default
#'  is NULL which is equivalent to `xformla=~1`.  This is
#'  used to create a matrix of covariates which is then passed
#'  to the 2x2 DID estimator chosen in `est_method`.
#'  X's are assumed fixed across the time dimension in this version.
#'  Use different columns Xt, Xt+1 if time-varying covariates are needed.
#' @param panel (Not used) This is not used as balanced and unbalanced panel data is treated similarly.
#' @param allow_unbalanced_panel (Not used) This is not used as balanced and unbalanced panel data is treated similarly.
#' @param control_group Which units to use the control group.
#'  The default is "nevertreated" which sets the control group
#'  to be the group of units that never participate in the
#'  treatment.  This group does not change across groups or
#'  time periods.  The other option is to set
#'  `group="notyettreated"`.  In this case, the control group
#'  is set to the group of units that have not yet participated
#'  in the treatment in that time period.  This includes all
#'  never treated units, but it includes additional units that
#'  eventually participate in the treatment, but have not
#'  participated yet.
#' @param anticipation (Not used) The number of time periods before participating
#'  in the treatment where units can anticipate participating in the
#'  treatment and therefore it can affect their untreated potential outcomes
#' @param base_period (Not used) The cdid package only uses the g-1 base period for the moment. Whether to use a "varying" base period or a
#'  "universal" base period.  Either choice results in the same
#'  post-treatment estimates of ATT(g,t)'s.  In pre-treatment
#'  periods, using a varying base period amounts to computing a
#'  pseudo-ATT in each treatment period by comparing the change
#'  in outcomes for a particular group relative to its comparison
#'  group in the pre-treatment periods (i.e., in pre-treatment
#'  periods this setting computes changes from period t-1 to period
#'  t, but repeatedly changes the value of t)
#'
#'  A universal base period fixes the base period to always be
#'  (g-anticipation-1).  This does not compute
#'  pseudo-ATT(g,t)'s in pre-treatment periods, but rather
#'  reports average changes in outcomes from period t to
#'  (g-anticipation-1) for a particular group relative to its comparison
#'  group.  This is analogous to what is often reported in event
#'  study regressions.
#'
#'  Using a varying base period results in an estimate of
#'  ATT(g,t) being reported in the period immediately before
#'  treatment.  Using a universal base period normalizes the
#'  estimate in the period right before treatment (or earlier when
#'  the user allows for anticipation) to be equal to 0, but one
#'  extra estimate in an earlier period.
#'
#' @references Bellego, Benatia, and Dortet-Bernadet (2024) \"The Chained Difference-in-Differences",
#' Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
#'
#' @return an [`MP`] object containing all the results for group-time average
#'  treatment effects
#' @export

att_gt_cdid <- function(yname,
                    tname,
                    idname=NULL,
                    gname,
                    xformla=NULL,
                    data,
                    panel=TRUE,
                    allow_unbalanced_panel=TRUE,
                    control_group,
                    anticipation=0,
                    weightsname=NULL,
                    alp=0.05,
                    bstrap=TRUE,
                    cband=TRUE,
                    biters=1000,
                    clustervars=NULL,
                    est_method="2-step",
                    base_period="varying",
                    print_details=FALSE,
                    pl=FALSE,
                    cores=1){

  #Part1. Pre-process step can be useful to carry over the parameters of the functions
  dp=pre_process_cdid(yname,
                      tname,
                      idname,
                      gname,
                      xformla,
                      data,
                      panel=TRUE,
                      allow_unbalanced_panel=TRUE,
                      control_group,
                      anticipation=0,
                      weightsname,
                      alp=0.05,
                      bstrap=TRUE,
                      biters=1000,
                      clustervars,
                      cband=TRUE,
                      est_method,
                      base_period="varying",
                      print_details=FALSE,
                      pl=FALSE,
                      cores=1,
                      call=match.call())

  ## Part 2. Prelim checks and compute Delta ATT
  result <- gmm_compute_delta_att(dp)

  # Part 3. Post-estimation aggregation step. Converts delta ATT into aggregated ATT.
  result <- gmm_convert_delta_to_att(result)

  # # Part 4. Result must be converted to be used in the aggte function.
  if   (dp$est_method == "2-step") {
  if (!exists("result") || is.null(result)) stop("Error: 'result' is NULL or does not exist before gmm_convert_result for 2-step.")
  result <- gmm_convert_result(result,1)}
  if (!exists("result") || is.null(result)) stop("Error: 'result' is NULL or does not exist before gmm_convert_result for identity.")
  if (dp$est_method == "Identity")  { #if not 2-step, then identity
  result <- gmm_convert_result(result,2)}

  return(result)
}



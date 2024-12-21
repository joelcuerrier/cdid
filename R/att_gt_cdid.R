

#' att_gt_cdid Function:
#' @description `att_gt` computes average treatment effects.
#' Our estimator accommodates (1) multiple time
#' periods, (2) variation in treatment timing, (3) treatment effect heterogeneity,
#' and (4) general missing data patterns
#' Bellego, Benatia,Dortet-Bernardet (2024) for a detailed description.
#' @param yname The name of the outcome variable
#' @param data The name of the data.frame that contains the data
#' @param tname The name of the column containing the time periods
#' @param idname The individual (cross-sectional unit) id name
#' @param gname The name of the variable in `data` that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "group" a unit belongs to.  It should be 0 for units
#'  in the untreated group.
#' @param weightsname The name of the column containing the sampling weights.
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
#' @param est_method the method to compute group-time average treatment effects.  The default is "dr" which uses the doubly robust
#' approach in the `DRDID` package.  Other built-in methods
#' include "ipw" for inverse probability weighting and "reg" for
#' first step regression estimators.  The user can also pass their
#' own function for estimating group time average treatment
#' effects.  This should be a function
#' `f(Y1,Y0,treat,covariates)` where `Y1` is an
#' `n` x `1` vector of outcomes in the post-treatment
#' outcomes, `Y0` is an `n` x `1` vector of
#' pre-treatment outcomes, `treat` is a vector indicating
#' whether or not an individual participates in the treatment,
#' and `covariates` is an `n` x `k` matrix of
#' covariates.  The function should return a list that includes
#' `ATT` (an estimated average treatment effect), and
#' `inf.func` (an `n` x `1` influence function).
#' The function can return other things as well, but these are
#' the only two that are required. `est_method` is only used
#' if covariates are included.
#' @param xformla A formula for the covariates to include in the
#'  model.  It should be of the form `~ X1 + X2`.  Default
#'  is NULL which is equivalent to `xformla=~1`.  This is
#'  used to create a matrix of covariates which is then passed
#'  to the 2x2 DID estimator chosen in `est_method`.
#' @param panel Whether or not the data is a panel dataset.
#'  The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a
#'  particular point in time.  The default is TRUE.  When
#'  is using a panel dataset, the variable `idname` must
#'  be set.  When `panel=FALSE`, the data is treated
#'  as repeated cross sections.
#' @param allow_unbalanced_panel Whether or not function should
#'  "balance" the panel with respect to time and id.  The default
#'  values if `FALSE` which means that [att_gt()] will drop
#'  all units where data is not observed in all periods.
#'  The advantage of this is that the computations are faster
#'  (sometimes substantially).
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
#' @param anticipation The number of time periods before participating
#'  in the treatment where units can anticipate participating in the
#'  treatment and therefore it can affect their untreated potential outcomes
#' @param faster_mode This option enables a faster version of `did`, optimizing
#' computation time for large datasets by improving data management within the package.
#' The default is set to `FALSE`. While the difference is minimal for small datasets,
#' it is recommended for use with large datasets.
#' @param base_period Whether to use a "varying" base period or a
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
#' @references Bellego, Benatia,Dortet-Bernardet (2022)  \"The Chained Difference-in-Differences\" , \url{https://www.davidbenatia.com/publication/chaineddid/}
#'
#' @return an [`MP`] object containing all the results for group-time average
#'  treatment effects
#' @export

att_gt_cdid <- function(yname,
                    tname,
                    idname=NULL,
                    gname,
                    xformla=NULL,
                    # propensityformla,
                    data,
                    chained = FALSE, #True to use chained.
                    panel=TRUE, 
                    allow_unbalanced_panel=TRUE, 
                    control_group, 
                    anticipation=0, 
                    weightsname,
                    weight_assumption=NULL,
                    alp=0.05,
                    bstrap=TRUE,
                    cband=TRUE, 
                    biters=1000,
                    clustervars=NULL, 
                    est_method, 
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
                      clustervars=NULL,
                      cband=TRUE,
                      est_method, 
                      base_period="varying",
                      print_details=FALSE,
                      pl=FALSE,
                      cores=1,
                      call=match.call())
  
  #gmm.R calls GMM_estimPeriod_Boot dans fonctions_estimation_Boot.R
  ########################################
  ## Part 2. Prelim checks (previously done with mp.spatt.GMM)
  ## And Compute Delta ATT (previously done with compute.mp.spatt.GMM)
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


 
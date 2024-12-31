#' @name pre_process_cdid
#' @title Process `cdid` Function Arguments
#' @importFrom stats complete.cases
#' @importFrom utils globalVariables
utils::globalVariables(c("..idname", "x"))
#' @description Function to process arguments passed to the main methods in the
#'  `cdid` package as well as conducting some tests to ensure
#'  data is in proper format and provides helpful error messages.
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
#' @param call (Not used) a call control var
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
#' @return a [`DIDparams`] object
#' @references Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences",
#' Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
#' @export
#'
pre_process_cdid <- function(yname,
                            tname,
                            idname,
                            gname,
                            xformla = NULL,
                            data,
                            panel = TRUE, #
                            allow_unbalanced_panel,
                            control_group = c("nevertreated","notyettreated"),
                            anticipation = 0,
                            weightsname = NULL,
                            alp = 0.05,
                            bstrap = FALSE,
                            cband=FALSE, #
                            biters = 1000,
                            clustervars = NULL,
                            est_method = "dr", #
                            base_period = "varying",
                            print_details = FALSE,
                            pl = FALSE,
                            cores = 1,
                            call = NULL) {
  #-----------------------------------------------------------------------------
  # Data pre-processing and error checking
  #-----------------------------------------------------------------------------

  # set control group
  #First the column nevertreated is created in the data set

  control_group <- control_group[1]
  data$control_group <- ifelse(data[[gname]] > 0, 0, 1)
  if(!(control_group %in% c("nevertreated","notyettreated"))){
    stop("control_group must be either 'nevertreated' or 'notyettreated'")
  }
  # make sure dataset is a data.frame
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(data) == "data.frame")) {
    data <- as.data.frame(data)
  }

  # make sure time periods are numeric
  if (! (is.numeric(data[, tname])) ) stop("data[, tname] must be numeric")

  #  make sure gname is numeric
  if (! (is.numeric(data[, gname])) ) stop("data[, gname] must be numeric")

  # put in blank xformla if no covariates
  if (is.null(xformla)) {
    xformla <- ~1
  }

  # drop irrelevant columns from data
  # data <- cbind.data.frame(data[,c(idname, tname, yname, gname, weightsname)], model.frame(xformla, data=data, na.action=na.pass))

  # check if any covariates were missing
  n_orig <- nrow(data)
  data <- data[complete.cases(data),]
  n_diff <- n_orig - nrow(data)
  if (n_diff != 0) {
    warning(paste0("dropped ", n_diff, " rows from original data due to missing data"))
  }

  # weights if null
  ifelse(is.null(weightsname), w <- rep(1, nrow(data)), w <- data[,weightsname])

  if (".w" %in% colnames(data)) stop("`did` tried to use column named \".w\" internally, but there was already a column with this name")
  data$.w <- w #sampling weights.

  # Outcome variable will be denoted by y
  # data$.y <- data[, yname]

  # figure out the dates
  # list of dates from smallest to largest
  tlist <- unique(data[,tname])[order(unique(data[,tname]))]

  # Groups with treatment time bigger than max time period are considered to be never treated
  asif_never_treated <- (data[,gname] > max(tlist, na.rm = TRUE))
  asif_never_treated[is.na(asif_never_treated)] <- FALSE
  data[asif_never_treated, gname] <- 0

  # list of treated groups (by time) from smallest to largest
  glist <- unique(data[,gname], )[order(unique(data[,gname]))]


  # Check if there is a never treated group
  if ( length(glist[glist==0]) == 0) {
    if(control_group=="nevertreated"){

      stop("There is no available never-treated group")
    } else {
      # Drop all time periods with time periods >= latest treated
      data <- subset(data,(data[,tname] < (max(glist)-anticipation)))
      # Replace last treated time with zero
      # lines.gmax <- data[,gname]==max(glist, na.rm = TRUE)
      # data[lines.gmax,gname] <- 0

      tlist <- sort(unique(data[,tname]))
      glist <- sort(unique(data[,gname]))

      # don't comput ATT(g,t) for groups that are only treated at end
      # and only play a role as a comparison group
      glist <- glist[ glist < max(glist)]
    }
  }

  # Only the treated groups
  glist <- glist[glist>0]

  # drop groups treated in the first period or before
  first.period <- tlist[1]
  glist <- glist[glist > first.period + anticipation]

  # check for groups treated in the first period and drop these
  # nfirstperiod <- length(unique(data[ !((data[,gname] > first.period) | (data[,gname]==0)), ] )[,idname])
  treated_first_period <- ( data[,gname] <= first.period ) & ( !(data[,gname]==0) )
  treated_first_period[is.na(treated_first_period)] <- FALSE
  nfirstperiod <- ifelse(panel, length(unique(data[treated_first_period,][,idname])), nrow(data[treated_first_period,]))
  if ( nfirstperiod > 0 ) {
    warning(paste0("Dropped ", nfirstperiod, " units that were already treated in the first period."))
    data <- data[ data[,gname] %in% c(0,glist), ]
    # update tlist and glist
    tlist <- unique(data[,tname])[order(unique(data[,tname]))]
    glist <- unique(data[,gname], )[order(unique(data[,gname]))]
    glist <- glist[glist>0]

    # drop groups treated in the first period or before
    first.period <- tlist[1]
    glist <- glist[glist > first.period + anticipation]

  }

  #  make sure id is numeric
  if (! is.null(idname)){
    #  make sure id is numeric
    if (! (is.numeric(data[, idname])) ) stop("data[, idname] must be numeric")

    ## # checks below are useful, but removing due to being slow
    ## # might ought to figure out a way to do these faster later
    ## # these checks are also closely related to making sure
    ## # that we have a well-balanced panel, so it might make
    ## # sense to move them over to the BMisc package

    ## # Check if idname is unique by tname
    ## n_id_year = all( table(data[, idname], data[, tname]) <= 1)
    ## if (! n_id_year) stop("The value of idname must be the unique (by tname)")

    ## # make sure gname doesn't change across periods for particular individuals
    ## if (!all(sapply( split(data, data[,idname]), function(df) {
    ##   length(unique(df[,gname]))==1
    ## }))) {
    ##   stop("The value of gname must be the same across all periods for each particular individual.")
    ## }
  }



  # if user specifies repeated cross sections,
  # set that it really is repeated cross sections
  true_repeated_cross_sections <- FALSE
  if (!panel) {
    true_repeated_cross_sections <- TRUE
  }

  #-----------------------------------------------------------------------------
  # setup data in panel case
  #-----------------------------------------------------------------------------
  if (panel) {

    # check for unbalanced panel
    if (allow_unbalanced_panel) {

      # code will run through repeated cross sections, so set panel to be FALSE
      panel <- FALSE
      true_repeated_cross_sections <- FALSE

      if (!is.numeric(data[,idname])) {
        stop("Must provide a numeric id")
      }

    } else {

      # this is the case where we coerce balanced panel

      # check for complete cases
      keepers <- complete.cases(data)
      n <- length(unique(data[,idname]))
      n.keep <- length(unique(data[keepers,idname]))
      if (nrow(data[keepers,]) < nrow(data)) {
        warning(paste0("Dropped ", (n-n.keep), " observations that had missing data."))
        data <- data[keepers,]
      }

      # make it a balanced data set
      n.old <- length(unique(data[,idname]))
      data <- BMisc::makeBalancedPanel(data, idname, tname)
      n <- length(unique(data[,..idname]))
      if (n < n.old) {
        warning(paste0("Dropped ", n.old-n, " observations while converting to balanced panel."))
      }

      # If drop all data, you do not have a panel.
      if (nrow(data)==0) {
        stop("All observations dropped to converted data to balanced panel. Consider setting `panel = FALSE' and/or revisit 'idname'.")
      }

      n <- nrow(data[ data[,tname]==tlist[1], ])

      # slow, repeated check here...
      ## # check that first.treat doesn't change across periods for particular individuals
      ## if (!all(sapply( split(data, data[,idname]), function(df) {
      ##   length(unique(df[,gname]))==1
      ## }))) {
      ##   stop("The value of gname must be the same across all periods for each particular individual.")
      ## }

    }
  }

  #-----------------------------------------------------------------------------
  # code for setting up repeated cross sections (and unbalanced panel)
  #-----------------------------------------------------------------------------

  if (!panel) {
    # check for complete cases


    keepers <- complete.cases(data)
    # keepers <- na.omit(data)
    if (nrow(data[keepers,]) < nrow(data)) {
      warning(paste0("Dropped ", nrow(data) - nrow(data[keepers,]), " observations that had missing data."))
      data <- data[keepers,]
    }

    # If drop all data, you do not have a panel.
    if (nrow(data)==0) {
      stop("All observations dropped due to missing data problems.")
    }

    # n-row data.frame to hold the influence function
    if (true_repeated_cross_sections) {
      data$.rowid <- seq(1:nrow(data))
      idname <- ".rowid"
    } else {
      # set rowid to idname for repeated cross section/unbalanced
      data$.rowid <- data[, idname]
    }

    # n is unique number of cross section observations
    # this is different for repeated cross sections and unbalanced panel
    n <- length(unique(data[,idname]))
  }

  ## # Update tlist and glist because of data handling
  ## # figure out the dates
  ## # list of dates from smallest to largest
  ## tlist <- unique(data[,tname])[order(unique(data[,tname]))]
  ## # list of treated groups (by time) from smallest to largest
  ## glist <- unique(data[,gname])[order(unique(data[,gname]))]

  ## # Only the treated groups
  ## glist <- glist[glist>0]

  ## # drop groups treated in the first period or before
  ## first.period <- tlist[1]
  ## glist <- glist[glist > first.period + anticipation]

  # Check if groups is empty (usually a problem with the way people defined groups)
  if(length(glist)==0){
    stop("No valid groups. The variable in 'gname' should be expressed as the time a unit is first treated (0 if never-treated).")
  }

  # if there are only two time periods, then uniform confidence
  # bands are the same as pointwise confidence intervals
  if (length(tlist)==2) {
    cband <- FALSE
  }

  #-----------------------------------------------------------------------------
  # more error handling after we have balanced the panel

  # check against very small groups
  gsize <- aggregate(data[,gname], by=list(data[,gname]), function(x) length(x)/length(tlist))

  # how many in each group before give warning
  # 5 is just a buffer, could pick something else, but seems to work fine
  reqsize <- length(BMisc::rhs.vars(xformla)) + 5

  # which groups to warn about
  gsize <- subset(gsize, x < reqsize) # x is name of column from aggregate

  # warn if some groups are small
  if (nrow(gsize) > 0) {
    gpaste <-  paste(gsize[,1], collapse=",")
    warning(paste0("Be aware that there are some small groups in your dataset.\n  Check groups: ", gpaste, "."))

    if ( (0 %in% gsize[,1]) & (control_group == "nevertreated") ) {
      stop("never treated group is too small, try setting control_group=\"notyettreated\"")
    }
  }
  #----------------------------------------------------------------------------

  # How many time periods
  nT <- length(tlist)
  # How many treated groups
  nG <- length(glist)

  # order dataset wrt idname and tname
  data <- data[order(data[,idname], data[,tname]),]

  # store parameters for passing around later
  dp <- DIDparams(yname=yname,
                  tname=tname,
                  idname=idname,
                  gname=gname,
                  xformla=xformla,
                  data=as.data.frame(data),
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
}

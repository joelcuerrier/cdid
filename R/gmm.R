# The vignettes included in this library are based on those authored by [Brantly Callaway and Pedro H.C. Sant'Anna] for the 'did' package. We have adapted the content for use in this library to provide users with accessible and informative documentation. The original vignettes authored by [Brantly Callaway and Pedro H.C. Sant'Anna] can be found in the 'did' package documentation [https://bcallaway11.github.io/did/articles/did-basics.html].
#' @title chained GMM 
#' @description The chained difference-in-differences, leverages the overlapping structure of many unbalanced panel data set.
#' @param tname The name of the column containing the time periods
#' @param yname The name of the outcome variable
#' @param idname The individual (cross-sectional unit) id name
#' @param gname The name of the variable in data that contains the first period when a particular observation is treated. This should be a positive number for all observations in treated groups. It defines which "group" a unit belongs to. It should be 0 for units in the untreated group.
#' @param xformla A formula for the covariates to include in the model. It should be of the form ~ X1 + X2. Default is NULL which is equivalent to xformla=~1. This is used to create a matrix of covariates which is then passed to the 2x2 DID estimator chosen in est_method.
#' @param propensityformla A formula for the covariates to include in the calculation of the propensity scores for the weights. This should be used only if weight_assumption is not NULL. Must be panel data.
#' @param data The name of the data.frame that contains the data
#' @param anticipation The number of time periods before participating in the treatment where units can anticipate participating in the treatment and therefore it can affect their untreated potential outcomes
#' @param weightsname The name of the column containing the sampling weights. If not set, all observations have same weight.
#' @param weight_assumption This is used to specify the assumption about the weights. The default is NULL which means that the weights are assumed to be constant over time. The other options are "missing_trends" and "sequential_missing". The first option assumes that the weights are constant over time, but that there are missing trends in the data. The second option assumes that the weights are constant over time, but that there are missing observations in the data. The last option is "varying" which assumes that the weights are varying over time. This is not yet implemented.
#' @param alp The significance level, default is 0.05
#' @param bstrap Boolean for whether or not to compute standard errors using the multiplier bootstrap. If standard errors are clustered, then one must set bstrap=TRUE. Default is TRUE (in addition, cband is also by default TRUE indicating that uniform confidence bands will be returned. If bstrap is FALSE, then analytical standard errors are reported.
#' @param biters The number of bootstrap iterations to use. The default is 1000, and this is only applicable if bstrap=TRUE.
#' @param base_period Whether to use a "varying" base period or a "universal" base period. Either choice results in the same post-treatment estimates of ATT(g,t)'s. In pre-treatment periods, using a varying base period amounts to computing a pseudo-ATT in each treatment period by comparing the change in outcomes for a particular group relative to its comparison group in the pre-treatment periods (i.e., in pre-treatment periods this setting computes changes from period t-1 to period t, but repeatedly changes the value of t). A universal base period fixes the base period to always be (g-anticipation-1). This does not compute pseudo-ATT(g,t)'s in pre-treatment periods, but rather reports average changes in outcomes from period t to (g-anticipation-1) for a particular group relative to its comparison group. This is analogous to what is often reported in event study regressions. Using a varying base period results in an estimate of ATT(g,t) being reported in the period immediately before treatment. Using a universal base period normalizes the estimate in the period right before treatment (or earlier when the user allows for anticipation) to be equal to 0, but one extra estimate in an earlier period.
#' @param print_details Whether or not to show details/progress of computations. Default is FALSE.
#' @param pl Whether or not to use parallel processing. Not yet functionning.
#' @param cores The number of cores to use for parallel processing. Not yet in function.
#' @param debT The first period of treatment.
#' @param finT The last period of treatment.
#' @param deb The first period of the data.
#' @param fin The last period of the data.
#' @param select A vector of observations to select.
#' @param treated The name of the column containing the treatment indicator. 
#' @param clustervars The name of the column containing the clustering variable. Default is NULL.
#' @param cband Boolean for whether or not to compute a uniform confidence band that covers all of the group-time average treatment effects with fixed probability `1-alp`.  In order to compute uniform confidence bands, `bstrap` must also be set to `TRUE`.  The default is `TRUE`.
#' @name gmm
#' @export
#' @examples
#' data=data_sim = fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
#' results=gmm(
#'                 yname="Y1_chaine",
#'                 tname="annee",
#'                 idname="id",
#'                 gname="annee_G",
#'                 xformla=~X, 
#'                 propensityformla=c("X"), 
#'                 data=data,      
#'                 anticipation=0,      
#'                 weightsname=c("P_Y1_chaine"), 
#'                 # weight_assumption="missing_trends", #Not yet implemented.
#'                 weight_assumption=NULL,
#'                 bstrap=FALSE, 
#'                 biters=1000,
#'                 debT=3,
#'                 finT=8,
#'                 deb=1,
#'                 fin=8,
#'                 select='select',
#'                 treated='traite_G',
#'                 pl=FALSE,
#'                 cores=1,
#'                 cband=TRUE,
#'                 clustervars=NULL)
#' results
#' 
#' 

#' @return  \item{att}{The average treatment effect on the treated}



  library(devtools)
  library(openxlsx)
  library(data.table) 
  library(BMisc)
  library(did)
  library(jsonlite)
  library(devtools)
  library(tidyr)
  set.seed(123)

#Compare the results gmm vs chained 
  # source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition.R")
  # source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition_nofe.R")
  # source("C:/Users/cuerr/Documents/cdid/R/fonctions_estimation_Boot.R")
  # source("C:/Users/cuerr/Documents/cdid/R/mp_spatt_Boot.R")
  # source("C:/Users/cuerr/Documents/cdid/R/compute_mp_spatt_Boot_alt.R")
  # source("C:/Users/cuerr/Documents/cdid/R/panelDiffV.R")
  # source("C:/Users/cuerr/Documents/cdid/R/gg.R")
  # source("C:/Users/cuerr/Documents/cdid/R/agregat.R")
  # source("C:/Users/cuerr/Documents/cdid/R/process_attgt.R")
  # source("C:/Users/cuerr/Documents/cdid/R/pre_process_did.R")
  # source("C:/Users/cuerr/Documents/cdid/R/DIDparams.R")
  # source("C:/Users/cuerr/Documents/cdid/R/mboot.R")
  # source("C:/Users/cuerr/Documents/cdid/R/MP.R")
  # source("C:/Users/cuerr/Documents/cdid/R/chained.R")
  # source("C:/Users/cuerr/Documents/cdid/R/compute.aggte.R")
  # source("C:/Users/cuerr/Documents/cdid/R/aggte.R")
  # source("C:/Users/cuerr/Documents/cdid/R/process_attgt_gmm.R")
  # source("C:/Users/cuerr/Documents/cdid/R/MP.R")
  
 
  gmm <-function(yname,
                    tname,
                    idname=NULL,
                    gname,
                    xformla=NULL,
                    propensityformla,
                    data,
                    panel=TRUE,  
                    allow_unbalanced_panel=TRUE, 
                    control_group=c("nevertreated","notyettreated"), 
                    anticipation=0, 
                    weightsname=NULL,
                    weight_assumption=NULL,
                    alp=0.05,
                    bstrap=TRUE,
                    cband=TRUE, 
                    biters=1000,
                    clustervars=NULL, 
                    est_method="chained", 
                    base_period="varying",
                    print_details=FALSE, 
                    pl=FALSE, 
                    cores=1, 
                    debT,
                    finT,
                    deb,
                    fin,
                    select,
                    treated){

  dp=pre_process_did(yname=yname,
                    tname=tname,
                    idname=idname,
                    gname=gname,
                    xformla=xformla,
                    data=data,
                    panel=panel,
                    allow_unbalanced_panel=allow_unbalanced_panel,
                    control_group=control_group,
                    anticipation=anticipation,
                    weightsname=weightsname,
                    alp=alp,
                    bstrap=bstrap,
                    biters=biters,
                    clustervars=clustervars,
                    cband=cband,
                    est_method=est_method,
                    base_period=base_period,
                    print_details=print_details,
                    pl=pl,
                    cores=cores,
                    call=match.call()
                    ,treated=treated
                    )

    #-----------------------------------------------------------------------------
    # Compute all ATT(g,t)
    #-----------------------------------------------------------------------------

  results=GMM_estimPeriod_Boot(yname=yname,
                            tname=tname,
                            idname=idname,
                            gname=gname,
                            xformla=xformla,
                            propensityformla=propensityformla,  
                            data=data,
                            debT=debT,
                            finT=finT,
                            deb=deb,
                            fin=fin,
                            select=select,
                            weightsname=weightsname,
                            weight_assumption=weight_assumption,
                            cband=cband,
                            alp=0.05,
                            bstrap=bstrap,
                            biters=biters,
                            treated=treated)    



    attgt_list=results[[1]][1]
    attgt_list=as.data.frame(attgt_list)

    attgt.list <- lapply(1:nrow(attgt_list), function(i) { #168,6
    list(
    att = as.numeric(attgt_list[i, yname]),#42
    group = as.numeric(attgt_list[i, gname]), #42
    
    year = as.numeric(attgt_list[i, tname]),
    post = ifelse(attgt_list[i, gname] > attgt_list[i, tname], 1, 0))})
    
    
    #Il y a un ajustement des résultats GMM pour les rendre compatibles avec les résultats de cDiD.
    #L'ajustement consiste à ajouter la dimension tbis. On a donc un ATT associé à chaque g,t et tbis. Sans quoi on a plusieurs ATT pour un seul g,t. 168 au lieu de 42.

    inffunc=results[[2]] #6597x42
    inffunc=as.data.frame(inffunc)
    inffunc=as.matrix(inffunc)
    attgt.results=process_attgt_gmm(attgt.list)
    
    #This patch is used to aggregate the results by g,t. It allows to have 8 groups instead of 2 JC 30-06-2024
    attgt.results<- aggregate(att ~ group + tt, data = as.data.frame(attgt.results), FUN = mean)
      
    group=attgt.results$group
    att=attgt.results$att
    tt=attgt.results$tt

    #Cette partie du code provient de la librarie cdid. L'output de cdid est adapté. 
    # analytical standard errors
    # estimate variance
    # this is analogous to cluster robust standard errors that
    # are clustered at the unit level

    # note to self: this def. won't work with unbalanced panel,
    # same with clustered standard errors
    # but it is always ignored b/c bstrap has to be true in that case

    n <- length(unique(data[,idname]))
    V <- Matrix::t(inffunc)%*%inffunc/n
    se <- sqrt(Matrix::diag(V)/n)  
    se[se <= sqrt(.Machine$double.eps)*10] <- NA

    # #muted for now because no clustvars argument in the current function  
    # # # if clustering along another dimension...we require using the
    # # # bootstrap (in principle, could come up with analytical standard
    # # # errors here though)
    # if ( (length(clustervars) > 0) & !bstrap) {
    #   warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
    # }

    # Identify entries of main diagonal V that are zero or NA
    zero_na_sd_entry <- unique(which(is.na(se))) 

    # bootstrap variance matrix
    if (bstrap) {
    bout <- mboot(inffunc, DIDparams=dp, pl=pl, cores=cores)
    bres <- bout$bres
    if(length(zero_na_sd_entry)>0) {
        se[-zero_na_sd_entry] <- bout$se[-zero_na_sd_entry]
    } else {
        se <- bout$se
    }
    }
    # Zero standard error replaced by NA
    se[se <= sqrt(.Machine$double.eps)*10] <- NA

    #-----------------------------------------------------------------------------
    # compute Wald pre-test
    #-----------------------------------------------------------------------------

    # select which periods are pre-treatment
    pre <- which(group > tt)

    #En utilisant process_attgt_gmm il y a 168 ATT. Il y a 42 groupes et 4 périodes. J'ajoute une patch pour retirer les valeurs supérieures à 42.
    pre <- pre[pre<= length(unique(group))*length(unique(tt))]

    # Drop group-periods that have variance equal to zero (singularity problems)
    if(length(zero_na_sd_entry)>0){
    pre <- pre[!(pre %in% zero_na_sd_entry)]
    }
    # pseudo-atts in pre-treatment periods
    preatt <- as.matrix(att[pre])

    # covariance matrix of pre-treatment atts
    preV <- as.matrix(V[pre,pre])

    # check if there are actually any pre-treatment periods
    if (length(preV) == 0) {
    message("No pre-treatment periods to test")
    W  <- NULL
    Wpval <- NULL
    } else if(sum(is.na(preV))) {
    warning("Not returning pre-test Wald statistic due to NA pre-treatment values")
    W <- NULL
    Wpval <- NULL
    } else if (rcond(preV) <= .Machine$double.eps) {
    # singluar covariance matrix for pre-treatment periods
    warning("Not returning pre-test Wald statistic due to singular covariance matrix")
    W <- NULL
    Wpval <- NULL
    } else {
    
    W <- n*t(preatt)%*%solve(preV)%*%preatt
    q <- length(pre) # number of restrictions
    Wpval <- round(1-pchisq(W,q),10)
    }


    #-----------------------------------------------------------------------------
    # compute confidence intervals / bands
    #-----------------------------------------------------------------------------

    # critical value from N(0,1), for pointwise
    cval <- qnorm(1-alp/2)

    # in order to get uniform confidence bands
    # HAVE to use the bootstrap
    if (bstrap){
    if (cband) {
        # for uniform confidence band
        # compute new critical value
        # see paper for details
        bSigma <- apply(bres, 2,
                        function(b) (quantile(b, .75, type=1, na.rm = T) -
                                        quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

        bSigma[bSigma <= sqrt(.Machine$double.eps)*10] <- NA

        # sup-t confidence band
        bT <- apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = TRUE))
        cval <- quantile(bT, 1-alp, type=1, na.rm = T)
        if(cval >= 7){
        warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
        }
    }
    }

    results=MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp,debT=debT)
      
    unique_id <- unique(results$DIDparams$data$id)
    unique_group <- unique(results$group)
    unique_t <- unique(results$t)

    inffunc_data <- expand.grid(id = unique_id,
                                group = unique_group,
                                t = unique_t)

    inffunc_data$gt <- as.integer(factor(paste0(inffunc_data$group, inffunc_data$t), levels = unique(paste0(inffunc_data$group, inffunc_data$t))))
    
    # inffunc_data$x=results$inffunc
    inffunc_data$x <- if (length(results$inffunc) >= nrow(inffunc_data)) results$inffunc[1:nrow(inffunc_data)] else c(results$inffunc, rep(NA, nrow(inffunc_data) - length(results$inffunc)))
    
    sparse_matrix <- sparseMatrix(
    i = inffunc_data$id,
    j = inffunc_data$gt,
    x = inffunc_data$x,
    dimnames = list(NULL, NULL)
    )

    # results$inffunc=sparse_matrix
    results$inffunc=inffunc
    return(results)
    }  
            
        

#' @export

GMM <- function(    yname,
                    tname,
                    idname=NULL,
                    gname,
                    xformla=NULL,
                    propensityformla,
                    data,
                    panel=TRUE, 
                    allow_unbalanced_panel=TRUE, 
                    control_group=c("nevertreated"), 
                    anticipation=0, 
                    weightsname,
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

  #Part1. Pre-process step can be useful to carry over the parameters of the functions
  dp=pre_process_cdid(yname="Y",
                     tname="date",
                     idname="id",
                     gname="date_G",
                     xformla=~X,
                     data=data0,
                     panel=TRUE,
                     allow_unbalanced_panel=TRUE,
                     control_group=c("notyettreated"), #either "nevertreated" or "notyettreated"
                     anticipation=0,
                     weightsname="S",
                     alp=0.05,
                     bstrap=TRUE,
                     biters=1000,
                     clustervars=NULL,
                     cband=TRUE,
                     est_method="chained",
                     base_period="varying",
                     print_details=FALSE,
                     pl=FALSE,
                     cores=1,
                     call=match.call()
                     ,treated="treated")
  

  #gmm.R calls GMM_estimPeriod_Boot dans fonctions_estimation_Boot.R
  ########################################
  ## Part 2. Prelim checks (previously done with mp.spatt.GMM)
  ## And Compute Delta ATT (previously done with compute.mp.spatt.GMM)
  
  result = gmm_compute_delta_att(dp)

  # Part 3. Post-estimation aggregation step. Converts delta ATT into aggregated ATT. 
  # (to be cleaned)
  result = gmm_convert_delta_to_att(result) 


  # # Part 4. Result must be converted to be used in the aggte function.  
  # # (to be cleaned)
  result = gmm_convert_result(result)
  
  return(result)
}


 
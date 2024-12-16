
#' @export

# # 0. We generate data in several steps.
# # Simulate data
# data=fonction_simu_attrition(theta2_alpha_Gg=0, lambda1_alpha_St=0, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=1)
# # data <- subset(data, select = -c(Y1_CS, Y1_longDID))  # delte this. and check i think in the actual data process there is observations that are dropped.

# # Filter indivs never observed across all years
# data0 <- data[data$P_Y1_longDID == 1,]
# data0 <- data0 %>%
#   group_by(id) %>%
#   filter(!(all(P_Y1_chaine == 0 & annee %in% 1:8))) %>%
#   ungroup()

# data0 <- data0[,c(1,2,3,6,7,8)]
# data0 <- rename(data0,Y = Y1_chaine)
# data0 <- rename(data0,date = annee)
# data0 <- rename(data0,date_G = annee_G)
# data0$P_Y1_chaine <- NULL

# # Sort data0 by id and date
# data0 <- data0 %>%
#   arrange(id, date)

# ### Correct id's
# list_id <-as.data.frame(unique(data0[,"id"]))
# list_id$iden_num<-1:dim(list_id)[1] 
# data0 <- merge(data0, list_id, by.x = "id", by.y = "id", all.x = TRUE)
# data0$id <- data0$iden_num
# data0$iden_num <- NULL

# # Add a treatment var
# data0$treated <- 0
# data0$treated[data0$date_G!=0] <- 1

# # Add a new variable S
# data0 <- data0 %>%
#   group_by(id) %>%                                # Group data by id
#   mutate(S = ifelse(date == min(date), 1, 0)) %>% # S == 1 for lowest year, 0 otherwise
#   ungroup()                                       # Ungroup data

# #### CUE: data0 is the typical dataset. The script starts here. The above script is just to make sure we have the expected data form
# # Data must have:
# # binary treatment var (treated)
# # sampling indicator for being observed in date t and t+1 (S)
# # outcome variable (Y)
# # predictor variable (X): several predictors should work
# # id for individuals (id)
# # dates (date)
# # treatment dates / cohorts (date_G)

chained <- function(yname,
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

  # Part 1. Pre-process step can be useful to carry over the parameters of the functions
  dp=pre_process_cdid(yname=yname,
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
                    est_method="chained",
                    base_period="varying",
                    print_details=print_details,
                    pl=pl,
                    cores=cores,
                    call=match.call(),
                    treated=treated)

  ## Part 2. Compute Delta ATT (previously done with chained.compute.mp.spatt.Boot)
  result = compute_delta_att(dp) 
  
  # Part 3. Post-estimation aggregation step. Converts delta ATT into aggregated ATT. 
  # (to be cleaned)
  result = convert_delta_to_att(result) 

  # Part 4. Result must be converted to be used in the aggte function.  
  # (to be cleaned)
  result = convert_result(result)
  
  return(result)
}

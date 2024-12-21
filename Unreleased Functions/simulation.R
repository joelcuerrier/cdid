
source("R/fonction_simu_attrition_params.R")
library(BMisc)
library(dplyr)
library(Matrix)
#########################################
# Code pour la Simulation
nsims=1
beta_hat_chaine  = matrix(NA,nsims,6)
IC_inf_chaine  = matrix(NA,nsims,6)
IC_sup_chaine  = matrix(NA,nsims,6)

set.seed(123)
for (simu_i in 1:nsims){
    print(simu_i)
    
    data=fonction_simu_attrition_params(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=0.6)
    
    # Filter indivs never observed across all years
    data0 <- data[data$P_Y1_longDID==1,]
    data0 <- data0 %>%
      group_by(id) %>%
      filter(!(all(P_Y1_chaine == 0 & annee %in% 1:8))) %>%
      ungroup()

    data0 <- data0[,c(1,2,3,6,7,8)]
    data0 <- rename(data0,Y = Y1_chaine)
    data0 <- rename(data0,date = annee)
    data0 <- rename(data0,date_G = annee_G)
    data0$P_Y1_chaine <- NULL

    # Sort data0 by id and date
    data0 <- data0 %>%
      arrange(id, date)

    ### Correct id's
    list_id <-as.data.frame(unique(data0[,"id"]))
    list_id$iden_num<-1:dim(list_id)[1] 
    data0 <- merge(data0, list_id, by.x = "id", by.y = "id", all.x = TRUE)
    data0$id <- data0$iden_num
    data0$iden_num <- NULL

    # Add a treatment var
    data0$treated <- 0
    data0$treated[data0$date_G!=0] <- 1

    # Add a new variable S
    data0 <- data0 %>%
      group_by(id) %>%                                # Group data by id
      # mutate(S = ifelse(date == min(date), 1, 0)) %>% # S == 1 for lowest year, 0 otherwise
      mutate(S = ifelse(date == min(date), 1, 1)) %>% # S == 1 for lowest year, 0 otherwise
      ungroup()   
    
    chained.results = GMM(
                    yname="Y",
                    tname="date",
                    idname="id",
                    gname="date_G",
                    xformla=~X, 
                    # propensityformla=c("X"), 
                    data=data0,      
                    weightsname="S",  
                    bstrap=FALSE, 
                    biters=1000,
                    # debT=3,
                    # finT=8,
                    # deb=1,
                    # fin=8,
                    # select='select',
                    treated='treated',
                    cband=TRUE)
      
      agg.es.chained <- aggte(MP = chained.results, type = 'dynamic')


      betas = tail(agg.es.chained$att.egt, 6)
      beta_hat_chaine[simu_i,1:length(betas)] = betas

      cband_lower <- agg.es.chained$att.egt - agg.es.chained$crit.val.egt*agg.es.chained$se.egt
      cband_lower = tail(cband_lower, 6)
      IC_inf_chaine[simu_i, 1:length(cband_lower)] = cband_lower
      
      cband_upper <- agg.es.chained$att.egt + agg.es.chained$crit.val.egt*agg.es.chained$se.egt
      cband_upper = tail(cband_upper, 6)
      IC_sup_chaine[simu_i, 1:length(cband_upper)] = cband_upper

      simu_attri=list(beta_hat_chaine,IC_inf_chaine,IC_sup_chaine)
      simu_attri
}

result_sim_attri = data.frame()

nb_estimateur=1
for (j in 1:nb_estimateur){
  for (i in 1:6){    
    result_sim_attri[(i*2-1),j]=round(mean(round(mean(beta_hat_chaine[,i]),digits=3)),digits=3)
    result_sim_attri[(i*2),1]=paste0(paste0("(",round(sd(simu_attri[[1]][,i]),digits=3)),")")
  }
}


simu_attri
result_sim_attri

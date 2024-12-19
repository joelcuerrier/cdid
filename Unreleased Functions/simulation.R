
  set.seed(123)

  source("R/fonction_simu_attrition.R")
  source("R/fonction_simu_attrition.R")
  source("R/fonction_simu_attrition_nofe.R")
  source("R/fonctions_estimation_Boot.R")
  source("R/mp_spatt_Boot.R")
  source("R/compute_mp_spatt_Boot_alt.R")
  source("R/panelDiffV.R")
  source("R/gg.R")
  source("R/agregat.R")
  source("R/process_attgt.R")
  source("R/pre_process_did.R")
  source("R/DIDparams.R")
  source("R/mboot.R")
  source("R/MP.R")
  source("R/chained.R")
  source("R/compute.aggte.R")
  source("R/aggte.R")
  source("R/compute.aggte.R")
  source("R/gmm.R")


install.packages("devtools")
library(devtools)
remotes::install_github("joelcuerrier/cdid", ref = "main", force = TRUE, dependencies=TRUE)
library(cdid)

# ===============================
#  Chained
# ===============================
# Nouveauté dans la simulation et notes importantes:
# Les influ obtenu avec WildBoostrap sont utilisés pour calculer les ICs. Pour cette raison
# j'utilise la fonction chained_estimPeriod_Boot, pour accéder aux résultats att et influence 
# aggrégés. La fonction chained n'utilise pas les aggrégations. Les codes de la vignettes eux, 
# utilisent les fonctions d'agrégat de Callaway qui utilise un simple Bootstrap. Il y a possibilité d'intégrer le wild boostrap dans cette fonction.

nsims=10

beta_hat_chaine  = matrix(NA,nsims,6)
IC_inf_chaine  = matrix(NA,nsims,6)
IC_sup_chaine  = matrix(NA,nsims,6)

set.seed(123)
for (simu_i in 1:nsims){
    print(simu_i)
    
    data=fonction_simu_attrition(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=0.6)
    data <- subset(data, select = -c(Y1_CS, Y1_longDID)) 
    chained.results=chained_estimPeriod_Boot(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                # anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption=NULL,
                bstrap=FALSE, 
                biters=1000,
                debT=3,
                finT=8,
                deb=1,
                fin=8,
                select='select',
                treated='traite_G'
                # pl=FALSE,
                # cores=1,
                # cband=TRUE,
                # clustervars=NULL
                )

    yname = "Y1_chaine"
    # beta_hat_chaine[simu_i, 1:length(chained.results[[1]][chained.results[[1]][,1]>0, 2])] = chained.results[[1]][chained.results[[1]][,1]>0,2]
    beta_hat_chaine[simu_i, 1:length(chained.results[[4]][chained.results[[4]][,1]>0, 2])] = chained.results[[4]][chained.results[[4]][,1]>0,2]
    # IC<-wild_bootstrap(chained.results[[1]],chained.results[[2]],chained.results[[3]],nom_outcome=yname,biters=1000,nivtest=0.05,seedvec=NULL)
    
    
    
    IC<-wild_bootstrap(chained.results[[4]],chained.results[[5]],chained.results[[3]],nom_outcome=yname,biters=1000,nivtest=0.05,seedvec=NULL)
    IC_inf_chaine[simu_i, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
    IC_sup_chaine[simu_i, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
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

# export to csv
library(readxl)
write.csv(simu_attri[[1]], "simu_attri1.csv") #ATT g,t
write.csv(simu_attri[[2]], "simu_attri2.csv") #IC inf
write.csv(simu_attri[[3]], "simu_attri3.csv") #IC sup
write.csv(result_sim_attri, "simu_attri4.csv") #result


  source("R/fonction_simu_attrition.R")
  source("R/fonction_simu_attrition.R")
  source("R/fonction_simu_attrition_nofe.R")
  source("R/fonctions_estimation_Boot.R")
  source("R/mp_spatt_Boot.R")
  source("R/compute_mp_spatt_Boot_alt.R")
  source("R/panelDiffV.R")
  source("R/gg.R")
  source("R/agregat.R")
  source("R/process_attgt.R")
  source("R/pre_process_did.R")
  source("R/DIDparams.R")
  source("R/mboot.R")
  source("R/MP.R")
  source("R/chained.R")
  source("R/compute.aggte.R")
  source("R/aggte.R")
  source("R/compute.aggte.R")
  source("R/gmm.R")

# ===============================
#  Gmm
# ===============================
nsims=100

set.seed(123)
beta_hat_chaine  = matrix(NA,nsims,6)
IC_inf_chaine  = matrix(NA,nsims,6)
IC_sup_chaine  = matrix(NA,nsims,6)

for (simu_i in 1:nsims){
    print(simu_i)
    
    data=fonction_simu_attrition(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=0.6)
    
    
    chained.results=GMM_estimPeriod_Boot(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                # anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption=NULL,
                bstrap=FALSE, 
                biters=1000,
                debT=3,
                finT=8,
                deb=1,
                fin=8,
                select='select',
                treated='traite_G'
                # pl=FALSE,
                # cores=1,
                # cband=TRUE,
                # clustervars=NULL
            )

    yname = "Y1_chaine" #I kepts the same variable name but they are defined as GMM outputs.
    # beta_hat_chaine[simu_i, 1:length(chained.results[[1]][chained.results[[1]][,1]>0, 2])] = chained.results[[1]][chained.results[[1]][,1]>0,2]
    beta_hat_chaine[simu_i, 1:length(chained.results[[4]][chained.results[[4]][,1]>0, 2])] = chained.results[[4]][chained.results[[4]][,1]>0,2]

    # IC<-wild_bootstrap(chained.results[[1]],chained.results[[2]],chained.results[[3]],nom_outcome=yname,biters=1000,nivtest=0.05,seedvec=NULL)
    IC<-wild_bootstrap(chained.results[[4]],chained.results[[5]],chained.results[[3]],nom_outcome=yname,biters=1000,nivtest=0.05,seedvec=NULL)
    IC_inf_chaine[simu_i, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
    IC_sup_chaine[simu_i, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
    simu_attri=list(beta_hat_chaine,IC_inf_chaine,IC_sup_chaine)

}

result_sim_attri = data.frame()
nb_estimateur=1
for (j in 1:nb_estimateur){
  for (i in 1:6){

    result_sim_attri[(i*2-1),j]=round(mean(round(mean(beta_hat_chaine[,i]),digits=3)),digits=3)
    result_sim_attri[(i*2),1]=paste0(paste0("(",round(sd(simu_attri[[1]][,i]),digits=3)),")")
  }
}

# export to csv
library(readxl)
write.csv(simu_attri[[1]], "simu_attri1gmm.csv") #ATT g,t
write.csv(simu_attri[[2]], "simu_attri2gmm.csv") #IC inf
write.csv(simu_attri[[3]], "simu_attri3gmm.csv") #IC sup
write.csv(result_sim_attri, "simu_attri4gmm.csv") #result

result_sim_attri
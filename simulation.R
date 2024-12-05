nsims=1000
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


beta_hat_chaine  = matrix(NA,nsims,6)
IC_inf_chaine  = matrix(NA,nsims,6)
IC_sup_chaine  = matrix(NA,nsims,6)
# mat.att <- matrix(nrow = 0, ncol = 0)
# mat.se <- matrix(nrow = 0, ncol = 0)

for (simu_i in 1:nsims){
    print(simu_i)
    data=fonction_simu_attrition(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=0.9)
    
    chained.results=chained(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption=NULL,
                bstrap=FALSE, 
                biters=1000,
                debT=3,
                finT=8,
                deb=1,
                fin=8,
                select='select',
                treated='traite_G',
                pl=FALSE,
                cores=1,
                cband=TRUE,
                clustervars=NULL)

    yname = "Y1_chaine"
    beta_hat_chaine[simu_i, 1:length(chained.results[[1]][chained.results[[1]][,1]>0, 2])] = chained.results[[1]][chained.results[[1]][,1]>0,2]
    IC<-wild_bootstrap(chained.results[[1]],chained.results[[2]],chained.results[[3]],nom_outcome=yname,biters=1000,nivtest=0.05,seedvec=NULL)
    IC_inf_chaine[simu_i, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
    IC_sup_chaine[simu_i, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
    simu_attri=list(beta_hat_chaine,IC_inf_chaine,IC_sup_chaine)
    # agg <- aggte(MP = chained.results, type = 'dynamic')
    # if (nrow(mat.att) == 0) {
    # mat.att <- matrix(nrow = length(agg$att.egt), ncol = 0)} #we adjust the size of the matrix for cbind.
    # if (nrow(mat.se) == 0) {
    # mat.se <- matrix(nrow = length(agg$se.egt), ncol = 0)} #we adjust the size of the matrix for cbind.
    
    
    # mat.att <- cbind(mat.att, agg$att.egt)
    # mat.se <- cbind(mat.se, unlist(agg$se.egt))
}

# row_averages = rowMeans(mat.att, na.rm = TRUE)
# results.att <- matrix(row_averages, ncol = 1)
# row_averages = rowMeans(mat.se, na.rm = TRUE)
# results.se <- matrix(row_averages, ncol = 1)
# results.final <- cbind(results.att, results.se)

# num_rows <- nrow(results.final)
# tt <- c(-(num_rows %/% 2):((num_rows - 1) %/% 2))
# results.final <- cbind(tt, results.final)
# results.final


result_sim_attri = data.frame()

nb_estimateur=1
for (j in 1:nb_estimateur){
  for (i in 1:6){
    
    result_sim_attri[(i*2-1),j]=round(mean(round(mean(beta_hat_chaine[,i]),digits=3)),digits=3)
    result_sim_attri[(i*2),1]=paste0(paste0("(",round(sd(simu_attri[[1]][,i]),digits=3)),")")
    
  }
}


result_sim_attri


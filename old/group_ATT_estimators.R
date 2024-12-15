
# install.packages("openxlsx")
# install.packages("data.table")
library(openxlsx)
library(data.table)
library(BMisc)
set.seed(123)

source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/base_de_donnees.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/fonctions_estimation_Boot.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/mp_spatt_Boot.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/compute_mp_spatt_Boot_alt.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/panelDiffV.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/gg.R")
source("C:/Users/cuerr/OneDrive - HEC Montréal/Desktop/research_assistant/économétrie/did_cdd/did_R_4.3.1/R/chained/agregat.R")
data=data_sim=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)

# View(data)
group_att_estimators <-function(yname,tname,treated,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp,bstrap,biters=1000,STRATE,pond_RD=NULL,bootStrap=FALSE){
#' @param yname : nom des variables dependantes
#' @param tname : nom de la variable temporelle
#' @param treated : nom de la variable traitee
#' @param idname : nom de la variable id
#' @param gname : nom de la variable de groupe
#' @param xformla : formule des variables explicatives
#' @param data : base de donnees
#' @param debT : debut de la periode d'estimation
#' @param finT : fin de la periode d'estimation
#' @param deb : debut de la periode d'analyse
#' @param fin : fin de la periode d'analyse
#' @param select : type de selection
#' @param weightsname : nom des poids
#' @param alp : niveau de test
#' @param bstrap : boolean pour le bootstrap
#' @param biters : nombre d'iterations du bootstrap
#' @param STRATE : type de stratification
#' @param pond_RD : ponderation RD
#' @param bootStrap : boolean pour le bootstrap

IC_inf_chaine  = matrix(NA,1,6)
IC_inf_CS      = matrix(NA,1,6)
IC_inf_GMM     = matrix(NA,1,6)
IC_inf_longDID = matrix(NA,1,6)
IC_sup_chaine  = matrix(NA,1,6)
IC_sup_CS      = matrix(NA,1,6)
IC_sup_longDID = matrix(NA,1,6)
IC_sup_GMM     = matrix(NA,1,6)

# print("chaine")
# chaine=estimPeriod_Boot(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp,bstrap,biters,STRATE,pond_RD,treated)    
# beta_hat_chaine  = matrix(NA,1,finT-debT+1)
# beta_hat_chaine[1:length(chaine[[1]][chaine[[1]][,1]>0, 2])] = chaine[[1]][chaine[[1]][,1]>0,2] 
print("CS")
CS=CS_estimPeriod_Boot(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp,bstrap,biters,STRATE,pond_RD,treated)
beta_hat_CS      = matrix(NA,1,finT-debT+1)
beta_hat_CS[1, 1:length(CS[[1]][CS[[1]][,1]>0, 2])] = CS[[1]][CS[[1]][,1]>0,2]
# print("longDID")
# longDID=longDID_estimPeriod_Boot(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp,bstrap,biters,STRATE,pond_RD,treated)
# beta_hat_longDID = matrix(NA,1,finT-debT+1)
# beta_hat_longDID[1, 1:length(longDID[[1]][longDID[[1]][,1]>0, 2])] = longDID[[1]][longDID[[1]][,1]>0,2]
# print("GMM")
# #specificite de l'estimateur GMM : il faut les outcome de l'estimateur chaine et les poids de l'estimateur CS
# GMM=GMM_estimPeriod_Boot(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp,bstrap,biters,STRATE,pond_RD,treated)
# beta_hat_GMM = matrix(NA,1,finT-debT+1)
# beta_hat_GMM[1, 1:length(GMM[[1]][GMM[[1]][,1]>0, 2])] = GMM[[1]][GMM[[1]][,1]>0,2]
# print("IC")
# # Calcul des ecarts-types avec wild bootstrap
# if (isTRUE(bootStrap)){
# IC<-wild_bootstrap(chaine[[1]],chaine[[2]],chaine[[3]],nom_outcome=out_chaine,biters=1000,nivtest=0.05,seedvec=NULL)
# IC_inf_chaine[1, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
# IC_sup_chaine[1, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
# # wild bootstrap estimateur Cross Section
# IC<-wild_bootstrap(CS[[1]],CS[[2]],CS[[3]],nom_outcome=out_CS,biters=1000,nivtest=0.05,seedvec=NULL)
# IC_inf_CS[1, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
# IC_sup_CS[1, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
# # wild bootstrap estimateur long DID
# IC<-wild_bootstrap(longDID[[1]],longDID[[2]],longDID[[3]],nom_outcome=out_longDID,biters=1000,nivtest=0.05,seedvec=NULL)
# IC_inf_longDID[1, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
# IC_sup_longDID[1, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
# # wild bootstrap estimateur GMM
# IC<-wild_bootstrap(GMM[[1]],GMM[[2]],GMM[[3]],nom_outcome=out_chaine,biters=1000,nivtest=0.05,seedvec=NULL)
# IC_inf_GMM[1, 1:length(IC[[2]][IC[[2]][,1]>0, 2])] = IC[[2]][IC[[2]][,1]>0,2]
# IC_sup_GMM[1, 1:length(IC[[3]][IC[[3]][,1]>0, 2])] = IC[[3]][IC[[3]][,1]>0,2]
#     }
out=list(beta_hat_chaine,beta_hat_CS,beta_hat_longDID,beta_hat_GMM,IC_inf_chaine,IC_sup_chaine,IC_inf_CS,IC_sup_CS,IC_inf_longDID,IC_sup_longDID,IC_inf_GMM,IC_sup_GMM) 
# out=list(beta_hat_CS)
}


results=group_att_estimators(yname=c("Y1_chaine","Y2_chaine"),tname="annee",idname="id",gname="annee_G",xformla=c("X"),data=data,debT=3,finT=8,deb=1,fin=8,select='select',weightsname=c("P_Y1_chaine","P_Y2_chaine"),alp=0.05,bstrap=TRUE,biters=1000,STRATE='strate',pond_RD=NULL,treated='traite_G')


# results=group_att_estimators(yname=c("Y1_CS","Y2_CS"),tname="annee",idname="id",gname="annee_G",xformla=c("X"),data=data,debT=3,finT=8,deb=1,fin=8,select='select',weightsname=c("P_Y1_CS","P_Y2_CS"),alp=0.05,bstrap=FALSE,biters=1000,STRATE='strate',pond_RD=NULL,bootStrap=FALSE,treated='traite_G')
# results=group_att_estimators(yname=c("Y1_longDID","Y2_longDID"),tname="annee",idname="id",gname="annee_G",xformla=c("X"),data=data,debT=3,finT=8,deb=1,fin=8,select='select',weightsname=c("P_Y1_longDID","P_Y2_longDID"),alp=0.05,bstrap=FALSE,biters=1000,STRATE='strate',pond_RD=NULL,bootStrap=FALSE,treated='traite_G')
# results=group_att_estimators(yname=c("Y1_chaine","Y2_chaine"),tname="annee",idname="id",gname="annee_G",xformla=c("X"),data=data,debT=3,finT=8,deb=1,fin=8,select='select',weightsname=c("P_Y1_CS","P_Y2_CS"),alp=0.05,bstrap=FALSE,biters=1000,STRATE='strate',pond_RD=NULL,bootStrap=FALSE,treated='traite_G')

# # outputs : liste de 9 matrices
debT=3
finT=8
if (!file.exists("Simulations")) {dir.create("Simulations")}
repertoire <- file.path(getwd(), "R/chained/Simulations")
save(results,file=paste0(repertoire,"/simu_attri.RData"))

# Mettre en forme les r�sultats pour exportation
# aa_attri<-load(paste0(repertoire,"/simu_attri.RData"))
# result_sim_attri = data.frame()
# result_sim_attri_IC = data.frame()

nb_estimateur=4 #nombre de model; (chaine, CS, longDID, GMM)
result_sim_attri <- matrix(NA, nrow = (finT - debT + 1) * 2, ncol = nb_estimateur)
result_sim_attri_IC <- matrix(NA, nrow = (finT - debT + 1) * 2, ncol = nb_estimateur * 2)

for (j in 1:nb_estimateur){
  for (i in 1:finT-debT+1){ 
    result_sim_attri[(i*2-1),j]=round(mean(results[[j]][,i]),digits=3)
    result_sim_attri[i*2,j]=paste0(paste0("(",round(sd(results[[j]][,i]),digits=3)),")")
  }
}

result_sim_attri<-data.frame(result_sim_attri)
names(result_sim_attri)=c("Chaine","CS","longDID","GMM")

nombre_de_beta <- finT - debT + 1
row_names <- c(sapply(1:nombre_de_beta, function(i) c(paste0("beta_", i), paste0("sd(beta_", i, ")"))))
row.names(result_sim_attri) <- row_names

for (j in 1:nb_estimateur){
  for (i in 1:finT-debT+1){
    result_sim_attri_IC[(i*2-1),(j*2-1)]=round(mean(results[[(j*2-1+nb_estimateur)]][,i]),digits=3)
    result_sim_attri_IC[i*2,(j*2-1)]=paste0(paste0("(",round(sd(results[[(j*2-1+nb_estimateur)]][,i]),digits=3)),")")
    result_sim_attri_IC[(i*2-1),j*2]=round(mean(results[[(j*2+nb_estimateur)]][,i]),digits=3)
    result_sim_attri_IC[i*2,j*2]=paste0(paste0("(",round(sd(results[[(j*2+nb_estimateur)]][,i]),digits=3)),")")
  }
}
result_sim_attri_IC<-data.frame(result_sim_attri_IC)
names(result_sim_attri_IC)=c("Chaine_Inf","Chaine_Sup","CS_Inf","CS_Sup","longDID_Inf","longDID_Sup","GMM_Inf","GMM_Sup")
row_names <- c(sapply(1:nombre_de_beta, function(i) c(paste0("beta_", i), paste0("sd(beta_", i, ")"))))
row.names(result_sim_attri_IC) <- row_names

View(result_sim_attri)




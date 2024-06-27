#' @noRd
library(devtools)
library(openxlsx)
library(data.table) 
library(BMisc)
library(did)
library(jsonlite)
library(devtools)
library(tidyr)
set.seed(123)

# load_all("C:/Users/cuerr/Documents/cdid/R/")
# ls(pos = "package:cdid")
# devtools::document("C:/Users/cuerr/Documents/cdid/R/")

source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition.R")
source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition_nofe.R")
source("C:/Users/cuerr/Documents/cdid/R/fonctions_estimation_Boot.R")
source("C:/Users/cuerr/Documents/cdid/R/mp_spatt_Boot.R")
source("C:/Users/cuerr/Documents/cdid/R/compute_mp_spatt_Boot_alt.R")
source("C:/Users/cuerr/Documents/cdid/R/panelDiffV.R")
source("C:/Users/cuerr/Documents/cdid/R/gg.R")
source("C:/Users/cuerr/Documents/cdid/R/agregat.R")
source("C:/Users/cuerr/Documents/cdid/R/process_attgt.R")
source("C:/Users/cuerr/Documents/cdid/R/pre_process_did.R")
source("C:/Users/cuerr/Documents/cdid/R/DIDparams.R")
source("C:/Users/cuerr/Documents/cdid/R/mboot.R")
source("C:/Users/cuerr/Documents/cdid/R/MP.R")
source("C:/Users/cuerr/Documents/cdid/R/chained.R")
source("C:/Users/cuerr/Documents/cdid/R/compute.aggte.R")
source("C:/Users/cuerr/Documents/cdid/R/aggte.R")
data=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
# data=fonction_simu_attrition_nofe(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)


results=chained(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption="missing_trends",
                bstrap=TRUE, #si je met false je recoit un message erreur pour aggte dynamic. Dans le code did cette erreure est controlée avec un message perso.
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



#prochains travaux futurs: anticipations? voir https://bcallaway11.github.io/did/articles/extensions.html

# Note: MP a été modifié pour rejetter les traités < debT. C'est pourquoi on a 36 combinaisons de g,t et non 42.
# J'uniformise l'output de chained pour que inffunc soit compatibles avec les fonctions de did. 
# Pour se faire, la matrice doit être transformée en dgc matrix. @p correspond au groupe g,t. 


unique_id <- unique(results$DIDparams$data$id)
unique_group <- unique(results$group)
unique_t <- unique(results$t)

inffunc_data <- expand.grid(id = unique_id,
                             group = unique_group,
                             t = unique_t)

inffunc_data$gt <- as.integer(factor(paste0(inffunc_data$group, inffunc_data$t), levels = unique(paste0(inffunc_data$group, inffunc_data$t))))
inffunc_data$x=results$inffunc

sparse_matrix <- sparseMatrix(
  i = inffunc_data$id,
  j = inffunc_data$gt,
  x = inffunc_data$x,
#   dims = c(6597, 36),
  dimnames = list(NULL, NULL)
)

results$inffunc=sparse_matrix



















results


did::ggdid(results)
agg.simple=aggte(MP = results, type = "simple")
summary(agg.simple)

agg.es=aggte(MP = results, type = "dynamic") 
ggdid(agg.es) #voir gt ==0, drole de resultat.

agg.gs <- aggte(MP = results, type = "group")
summary(agg.gs)
# # ggdid(agg.gs)

# # agg.ct <- aggte(MP = results, type = "calendar")
# # summary(agg.ct)
# # ggdid(agg.ct)








# # #  $ V_analytical: num [1:36, 1:36] 130.87 -42.3 4.07 12.73 -3.44 ...
# # # au lieu de 
# # # $ V_analytical:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots  
# # #   .. ..@ i       : int [1:144] 0 1 2 3 4 5 6 7 8 9 ...
# # #   .. ..@ p       : int [1:13] 0 12 24 36 48 60 72 84 96 108 ...

# # #Prochaines étapes
# # # réparer group_ATT_estimators (valider les resulttats)
# # # (voir function gg) et les betas de 0.25,50,.75 qui sont hard coded
# # #améliorer l'efficacité fonction compute
# # #mettre clean
# # #je crois dans les fonctions d'aggregation de did, ils utilisent les poids pg du papier de callaway. voir les fonctions compute.aggte.R (les codes pour faire rouller aggte.)
# #Vérifier les fonctions de compute.aggte.R et les résultats:
#         # dim(data)#lambda probleme de selection
#         #enlever les effets fixes juste garder betas pour generer y att pour valider les betas .
#         # cs pour valider 

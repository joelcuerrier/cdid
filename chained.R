

# install.packages("openxlsx")
# install.packages("data.table")
# install.packages("jsonlite")
# install.packages("devtools")
#' @noRd
library(devtools)
library(openxlsx)
library(data.table) 
library(BMisc)
library(did)
library(jsonlite)
# install.packages("did")

set.seed(123)

library(devtools)
# load_all("C:/Users/cuerr/Documents/cdid/R/")
# ls(pos = "package:cdid")
# devtools::document("C:/Users/cuerr/Documents/cdid/R/")


source("C:/Users/cuerr/Documents/cdid/R/base_de_donnees.R")
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

data=data_sim=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
results=chained(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X,
                data=data,      
                anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption="missing_trends",
                link="logit",
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
                clustervars=NULL)


results

#Prochaines étapes
# réparer group_ATT_estimators (valider les resulttats)
#implémenter les poids différents selon les hypothèses

#retirer les hard code et rendre les fonctions clean pour etre integres dans did (ex les poids de la sim)

# (voir function gg) et les betas de 0.25,50,.75 qui sont hard coded
#retirer le truc des strates, et de pond_RD

#améliorer l'efficacité fonction compute
#mettre clean
#publier
#vérifier les hypothèses
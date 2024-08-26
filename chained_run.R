# devtools::build()
# devtools::install_local("cdid_0.1.0.tar.gz")

# install.packages("openxlsx")
# install.packages("data.table")
# install.packages("jsonlite")
# install.packages("devtools")
# install.packages("Matrix")

#' @noRd
library(devtools)
library(openxlsx)
library(data.table) 
library(BMisc)
library(did)
library(jsonlite)
library(devtools)
library(tidyr)
library(Matrix)

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
source("C:/Users/cuerr/Documents/cdid/R/compute.aggte.R")

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
                # weight_assumption=NULL,
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

agg.es=aggte(MP = results, type = "dynamic") 
agg.es
ggdid(agg.es) 

agg.gs <- aggte(MP = results, type = "group")
summary(agg.gs)
ggdid(agg.gs)

agg.ct <- aggte(MP = results, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)


#did pour comparer les results
library(did)
res=did::att_gt( yname="Y1_chaine",
            tname="annee",
            panel=FALSE,
            idname="id",
            gname="annee_G",
            xformla=~X, 
            data=data,      
            bstrap=FALSE, 
            biters=1000)
res
agg1=aggte(MP = res, type = "dynamic") 
ggdid(agg1)

        


                
                
                
                
                
                
                
                
                
                
                




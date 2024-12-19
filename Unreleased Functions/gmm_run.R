

# data=fonction_simu_attrition_nofe(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
data=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)

source("C:/Users/cuerr/Documents/cdid/R/gmm.R")
results=gmm(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                anticipation=0,      
                weightsname=c("P_Y1_chaine"), 
                # weight_assumption="missing_trends", #Not yet implemented.
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
          

source("C:/Users/cuerr/Documents/cdid/R/compute.aggte.R")
source("C:/Users/cuerr/Documents/cdid/R/aggte.R")

agg.gs <- aggte(MP = results, type = "group")
summary(agg.gs)












agg.ct <- aggte(MP = results, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)



# Group Effects:
#  Group Estimate Std. Error [95% Pointwise  Conf. Band]
#      3   0.4954     0.0533          0.3703      0.6206 *
#      4   0.8004     0.0782          0.6168      0.9841 *
#      5   0.9125     0.1616          0.5329      1.2921 *
#      6   1.1042     0.1026          0.8632      1.3452 *
#      7   1.6720     0.1926          1.2196      2.1244 *
#      8   1.6703     0.1487          1.3208      2.0197 *



          



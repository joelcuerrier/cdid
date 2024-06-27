


#   library(devtools)
#   library(openxlsx)
#   library(data.table) 
#   library(BMisc)
#   library(did)
#   library(jsonlite)
#   library(devtools)
#   library(tidyr)
#   set.seed(123)

#   # load_all("C:/Users/cuerr/Documents/cdid/R/")
#   # ls(pos = "package:cdid")
#   # devtools::document("C:/Users/cuerr/Documents/cdid/R/")

#   source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition.R")
#   source("C:/Users/cuerr/Documents/cdid/R/fonction_simu_attrition_nofe.R")
#   source("C:/Users/cuerr/Documents/cdid/R/fonctions_estimation_Boot.R")
#   source("C:/Users/cuerr/Documents/cdid/R/mp_spatt_Boot.R")
#   source("C:/Users/cuerr/Documents/cdid/R/compute_mp_spatt_Boot_alt.R")
#   source("C:/Users/cuerr/Documents/cdid/R/panelDiffV.R")
#   source("C:/Users/cuerr/Documents/cdid/R/gg.R")
#   source("C:/Users/cuerr/Documents/cdid/R/agregat.R")
#   source("C:/Users/cuerr/Documents/cdid/R/process_attgt.R")
#   source("C:/Users/cuerr/Documents/cdid/R/pre_process_did.R")
#   source("C:/Users/cuerr/Documents/cdid/R/DIDparams.R")
#   source("C:/Users/cuerr/Documents/cdid/R/mboot.R")
#   source("C:/Users/cuerr/Documents/cdid/R/MP.R")
#   source("C:/Users/cuerr/Documents/cdid/R/chained.R")
#   source("C:/Users/cuerr/Documents/cdid/R/compute.aggte.R")
#   source("C:/Users/cuerr/Documents/cdid/R/aggte.R")
#   # data=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
#   # data=fonction_simu_attrition_nofe(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
#   data=fonction_simu_attrition(nbsimu = 1, theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
 
#   gmm <-function(yname,
#                     tname,
#                     idname=NULL,
#                     gname,
#                     xformla=NULL,
#                     propensityformla,
#                     data,
#                     panel=TRUE,  #
#                     allow_unbalanced_panel=TRUE, #
                    
#                     control_group=c("nevertreated","notyettreated"), #
#                     anticipation=0, #
#                     weightsname=NULL,
#                     weight_assumption=NULL,
                    
#                     alp=0.05,
#                     bstrap=TRUE,
#                     cband=TRUE, #
#                     biters=1000,
#                     clustervars=NULL, #
#                     est_method="chained", #
#                     base_period="varying",#
#                     print_details=FALSE, #
#                     pl=FALSE, #
#                     cores=1, #
#                     debT,
#                     finT,
#                     deb,
#                     fin,
#                     select,
#                     treated){

#   dp=pre_process_did(yname=yname,
#                     tname=tname,
#                     idname=idname,
#                     gname=gname,
#                     xformla=xformla,
#                     data=data,
#                     panel=panel,
#                     allow_unbalanced_panel=allow_unbalanced_panel,
#                     control_group=control_group,
#                     anticipation=anticipation,
#                     weightsname=weightsname,
#                     alp=alp,
#                     bstrap=bstrap,
#                     biters=biters,
#                     clustervars=clustervars,
#                     cband=cband,
#                     est_method=est_method,
#                     base_period=base_period,
#                     print_details=print_details,
#                     pl=pl,
#                     cores=cores,
#                     call=match.call()
#                     ,treated=treated
#                     )

    
#     #-----------------------------------------------------------------------------
#     # Compute all ATT(g,t)
#     #-----------------------------------------------------------------------------

#   resultat=GMM_estimPeriod_Boot(yname=yname,
#                             tname=tname,
#                             idname=idname,
#                             gname=gname,
#                             xformla=xformla,
#                             propensityformla=propensityformla,  
#                             data=data,
#                             debT=debT,
#                             finT=finT,
#                             deb=deb,
#                             fin=fin,
#                             select=select,
#                             weightsname=weightsname,
#                             weight_assumption=weight_assumption,
                            
#                             cband=cband,
#                             alp=0.05,
#                             bstrap=bstrap,
#                             biters=biters,
#                             treated=treated)    

# #delete once finished begin{                   
#                     }  
                    
# results=gmm(
#                 yname="Y1_chaine",
#                 tname="annee",
#                 idname="id",
#                 gname="annee_G",
#                 xformla=~X, 
#                 propensityformla=c("X"), 
#                 data=data,      
#                 anticipation=0,      
#                 weightsname=c("P_Y1_chaine"), #St   
#                 # weight_assumption="missing_trends",
#                 weight_assumption=NULL,
#                 bstrap=FALSE, 
#                 biters=1000,
#                 debT=3,
#                 finT=8,
#                 deb=1,
#                 fin=8,
#                 select='select',
#                 treated='traite_G',
#                 pl=FALSE,
#                 cores=1,
#                 cband=TRUE,
#                 clustervars=NULL)
          
#                 #definition hard coded to run function
#                 yname="Y1_chaine"
#                 tname="annee"
#                 idname="id"
#                 gname="annee_G"
#                 xformla=~X
#                 propensityformla=c("X")
#                 alp=0.05
#                 anticipation=0
#                 weightsname=c("P_Y1_chaine")
#                 weight_assumption=NULL
#                 bstrap=FALSE
#                 biters=1000
#                 debT=3
#                 finT=8
#                 deb=1
#                 fin=8
#                 select='select'
#                 treated='traite_G'
#                 pl=FALSE
#                 cores=1
#                 cband=TRUE
#                 clustervars=NULL
          
          
#           dp=pre_process_did(yname=yname,
#                     tname=tname,
#                     idname=idname,
#                     gname=gname,
#                     xformla=xformla,
#                     data=data,
#                     # panel=panel,
#                     allow_unbalanced_panel=TRUE,
#                     # control_group=control_group,
#                     anticipation=anticipation,
#                     weightsname=weightsname,
#                     alp=alp,
#                     bstrap=bstrap,
#                     biters=biters,
#                     clustervars=clustervars,
#                     cband=cband,
#                     # est_method='est_method',
#                     # base_period=base_period,
#                     print_details=FALSE,
#                     pl=pl,
#                     cores=cores,
#                     call=match.call()
#                     ,treated=treated
#                     )
          
          
          

#           #delete end}
# # #Note: 
# # #À discuter avec David.
# # #En print print(results[[1]][1]) on voit que les valeurs de nC et de nG ne sont pas bonnes. Il devrait y avoir une variation dans le temps et en fonction d'alpha. Les valeurs ATT sont bonnes. Il est possible de le valider en utilisant la base de données nofe. Est-ce que ça vaut la peine de passer plus de temps à essayer de trouver la source de l'erreur pour nG et nC si le reste est bon? L'utilisateur ne voit pas print(results[[1]][1]) de toute façon.

# #Les codes suivants sont basés sur les codes de la fonction att_gt.R du package did
# #https://cran.r-project.org/web/packages/did/did.pdf



# resultat=results

# #À discuter avec David:
# #attgt a une dimension de 168,6. On a donc une valeur ATT pour chaque combinaison de g,t et tbis. 
# #6 valeurs pour g, t a des valeurs qui vont de 3 à 8, 4 à 8, etc. et tbis est la différence entre g et t.



# attgt_list=resultat[[1]][1]
# attgt_list=as.data.frame(attgt_list)

# attgt.list <- lapply(1:nrow(attgt_list), function(i) {
# list(
#   att = as.numeric(attgt_list[i, yname]),#42
#   group = as.numeric(attgt_list[i, gname]), #42
#   year = as.numeric(attgt_list[i, tname]),
#   post = ifelse(attgt_list[i, gname] > attgt_list[i, tname], 1, 0))})
# # dim(attgt_list) #168,6
# # class(attgt_list) #data.frame

# #Il y a un ajustement des résultats GMM pour les rendre compatibles avec les résultats de cDiD.
# #L'ajustement consiste à ajouter la dimension tbis. On a donc un ATT associé à chaque g,t et tbis. Sans quoi on a plusieurs ATT pour un seul g,t. 168 au lieu de 42.


# #Pour générer tbis (visuellement plus facile, mais inutile)
# # library(purrr)
# # attgt.list=map_dfr(attgt.list, ~as.data.frame(.x))
# # attgt.list$tbis=attgt.list$group-attgt.list$year

# # inffunc=resultat[[1]][2]
# inffunc=resultat[[2]]
# inffunc=as.data.frame(inffunc)
# inffunc=as.matrix(inffunc)
# # View(inffunc)
# # dim(inffunc) #6597 x 336(?)/ok maintenant 6597x42
# # class(inffunc) #matrix / array



# source("C:/Users/cuerr/Documents/cdid/R/process_attgt_gmm.R")
# # attgt.results=process_attgt_gmm(attgt.list)
# attgt.results=process_attgt(attgt.list)



# group=attgt.results$group
# att=attgt.results$att
# tt=attgt.results$tt

# # analytical standard errors
# # estimate variance
# # this is analogous to cluster robust standard errors that
# # are clustered at the unit level

# # note to self: this def. won't work with unbalanced panel,
# # same with clustered standard errors
# # but it is always ignored b/c bstrap has to be true in that case

# n <- length(unique(data[,idname]))
# V <- Matrix::t(inffunc)%*%inffunc/n
# # View(V)
# se <- sqrt(Matrix::diag(V)/n)  
# se[se <= sqrt(.Machine$double.eps)*10] <- NA
# # View(se)

# # #muted for now because no clustvars argument in the current function  
# # # # if clustering along another dimension...we require using the
# # # # bootstrap (in principle, could come up with analytical standard
# # # # errors here though)
# # if ( (length(clustervars) > 0) & !bstrap) {
# #   warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
# # }

# # Identify entries of main diagonal V that are zero or NA
# zero_na_sd_entry <- unique(which(is.na(se))) 

# # bootstrap variance matrix
# if (bstrap) {
#   bout <- mboot(inffunc, DIDparams=dp, pl=pl, cores=cores)
#   bres <- bout$bres
#   if(length(zero_na_sd_entry)>0) {
#     se[-zero_na_sd_entry] <- bout$se[-zero_na_sd_entry]
#   } else {
#     se <- bout$se
#   }
# }
# # Zero standard error replaced by NA
# se[se <= sqrt(.Machine$double.eps)*10] <- NA

# #-----------------------------------------------------------------------------
# # compute Wald pre-test
# #-----------------------------------------------------------------------------

# # select which periods are pre-treatment
# pre <- which(group > tt)

# #En utilisant process_attgt_gmm il y a 168 ATT. Il y a 42 groupes et 4 périodes. J'ajoute une patch pour retirer les valeurs supérieures à 42.
# pre <- pre[pre<= length(unique(group))*length(unique(tt))]

# # Drop group-periods that have variance equal to zero (singularity problems)
# if(length(zero_na_sd_entry)>0){
#   pre <- pre[!(pre %in% zero_na_sd_entry)]
# }
# # pseudo-atts in pre-treatment periods
# preatt <- as.matrix(att[pre])

# # covariance matrix of pre-treatment atts
# preV <- as.matrix(V[pre,pre])

# # check if there are actually any pre-treatment periods
# if (length(preV) == 0) {
#   message("No pre-treatment periods to test")
#   W  <- NULL
#   Wpval <- NULL
# } else if(sum(is.na(preV))) {
#   warning("Not returning pre-test Wald statistic due to NA pre-treatment values")
#   W <- NULL
#   Wpval <- NULL
# } else if (rcond(preV) <= .Machine$double.eps) {
#   # singluar covariance matrix for pre-treatment periods
#   warning("Not returning pre-test Wald statistic due to singular covariance matrix")
#   W <- NULL
#   Wpval <- NULL
# } else {
  
#   W <- n*t(preatt)%*%solve(preV)%*%preatt
#   q <- length(pre) # number of restrictions
#   Wpval <- round(1-pchisq(W,q),10)
# }


# #-----------------------------------------------------------------------------
# # compute confidence intervals / bands
# #-----------------------------------------------------------------------------

# # critical value from N(0,1), for pointwise
# cval <- qnorm(1-alp/2)

# # in order to get uniform confidence bands
# # HAVE to use the bootstrap
# if (bstrap){
#   if (cband) {
#     # for uniform confidence band
#     # compute new critical value
#     # see paper for details
#     bSigma <- apply(bres, 2,
#                     function(b) (quantile(b, .75, type=1, na.rm = T) -
#                                     quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))

#     bSigma[bSigma <= sqrt(.Machine$double.eps)*10] <- NA

#     # sup-t confidence band
#     bT <- apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = TRUE))
#     cval <- quantile(bT, 1-alp, type=1, na.rm = T)
#     if(cval >= 7){
#       warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
#     }
#   }
# }
# # length(group)
# # length(tt)
# # length(att)
# # length(V)
# # length(se)
# # length(cval)
# # length(inffunc)
# # dim(inffunc)
# # class(inffunc)
# source("C:/Users/cuerr/Documents/cdid/R/MP.R")
# MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp,debT=debT)
# return(MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp,debT=debT))  
# # }

# results=gmm(
#                 yname="Y1_chaine",
#                 tname="annee",
#                 idname="id",
#                 gname="annee_G",
#                 xformla=~X, 
#                 propensityformla=c("X"), 
#                 data=data,      
#                 anticipation=0,      
#                 weightsname=c("P_Y1_chaine"), #St   
#                 # weight_assumption="missing_trends",
#                 weight_assumption=NULL,
#                 bstrap=FALSE, 
#                 biters=1000,
#                 debT=3,
#                 finT=8,
#                 deb=1,
#                 fin=8,
#                 select='select',
#                 treated='traite_G',
#                 pl=FALSE,
#                 cores=1,
#                 cband=TRUE,
#                 clustervars=NULL)

# # str(results)
# View(results[[1]][1])
 
# unique_id <- unique(results$DIDparams$data$id)
# unique_group <- unique(results$group)
# unique_t <- unique(results$t)


# inffunc_data <- expand.grid(id = unique_id,
#                              group = unique_group,
#                              t = unique_t)
# dim(inffunc_data) #à confirmer si ok...
# inffunc_data$gt <- as.integer(factor(paste0(inffunc_data$group, inffunc_data$t), levels = unique(paste0(inffunc_data$group, inffunc_data$t))))
# inffunc_data$x=results$inffunc

# sparse_matrix <- sparseMatrix(
#   i = inffunc_data$id,
#   j = inffunc_data$gt,
#   x = inffunc_data$x,
#   dimnames = list(NULL, NULL)
# )

# results$inffunc=sparse_matrix


# # did::ggdid(results)

# # agg.es=aggte(MP = results, type = "dynamic") 
# # summary(agg.simple)

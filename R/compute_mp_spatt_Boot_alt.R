



##### ---------------------------------------------------------------------------------------------------------------------------
##### CHAINED-----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------

chained.compute.mp.spatt.Boot <- function ( nom_outcome
                                 , nom_traitement
                                 , flen  ## nombre de cohortes dans le traitement 
                                 , tlen  ## nombre d'annee
                                 , flist ## Le vecteur des vrais traitements (sans 0 !!)
                                 , tlist ## vecteur ordonn? des dates
                                 , data  ## Les données
                                 , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donné
                                 , xformla ## le modèle de score ~X
                                 , propensityformla
                                 , tname ## la variable de date
                                 , w ## ??  NULL par défaut
                                 , idname ## identifiant individu
                                 , method ## par défaut "logit"
                                 , seedvec ## par défaut NULL
                                 , se ## TRUE par défaut
                                 , pl ## FALSE par défaut
                                 , cores ## par défaut = 2
                                 , printdetails ## TRUE par defaut
                                 , selection ## indicatrice observations valides
                                 , ponderation ## poids utilise  
                                 , weight_assumption ## Hypothèses sur les poids (missing_trends, missing_outcomes, missing_trends_outcomes)
                                 , debT
                                 , finT
                                 
)
{

  
  

  ################################
  # initialisation des matrices  #
  ################################
  nbligne=flen*(tlen - 1)
  
  rempli<-seq(1:nbligne)
  fatt <- data.frame(compteur=rempli)
  fatt[,first.treat.name]<-rempli
  fatt[,tname]<-rempli
  fatt$nobsG<-rempli
  fatt$nobsC<-rempli 
  for(i in 1:length(nom_outcome)){
    fatt[,nom_outcome[i]]<-rempli}
  
  ### Liste des individus impliqués dans l'estimation pour le wild bootstrap
  indiv<-unique(data[,c(first.treat.name,idname)]) #liste des observations uniques (id) et leur première année de traitement
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
  counter <- 1
  
  #Matrices pour les poids
  ssit <-matrix(0, nrow = nrow(indiv), ncol = tlen-1) 
  st1 <- matrix(0, nrow = nrow(indiv), ncol = tlen-1) 
  st2 <- matrix(0, nrow = nrow(indiv), ncol = tlen-1) 
  cit <- matrix(0, nrow = nrow(indiv), ncol = tlen-1)
  flag=FALSE #In the calculations of reweights, Gg matrix does not have to be calculated every g,t. I use this flag to calculate it only once.
  
##########################################
# Boucles sur les cohortes et les dates  #
##########################################
  for (f in 1:flen) {  
    ## boucle sur la date : t, [1 2 3 4 5 6 7 8]. tlist[t+1] prend les valeurs 2, 3, 4, 5, 6, 7, 8. tlist[t] prend les valeurs 1, 2, 3, 4, 5, 6, 7.
    for (t in 1:(tlen - 1)) {  
        pret <- flist[f]-1 # cohorte G-1, flist est inclut dans [3 4 5 6 7 8], pret  est inclut dans [2,3,4,5,6,7].
        
        
        

        #################################################
        # Mesure des Pscores et Prédictions des Pscores #
        #################################################

        disdat <- data[(data[, tname] == pret), ]
        disdat <- droplevels(disdat) 
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) 
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        if (is.null(xformla)) {
          xformla <- ~1}
        pformla <- xformla
        
        ### On supprime les observations pour lesquelles le modèle de score ne peut être défini
        LesX<-BMisc::rhs.vars(pformla)
        bbb<-length(LesX)
        for (jj in 1:bbb){
          disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,] }
        
        ## finalement définition du modèle de score
        pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla)) #G~X
        base_score<-subset(disdat, C + G == 1) 
        base_score=base_score[base_score$select==1,] 
        pscore.reg <- glm(pformla, family = binomial(link = "logit"),
                          data = base_score)
        thet <- coef(pscore.reg) #juste pour regarder si le modèle est bien estimé
        #On garde les observation avec "annee"==G-1 ou t+1 ou t (pour T et conterfactual)
        
        
        if (any(is.na(thet))) {
          warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
          flist[f], " at time period: ", tlist[t +1]))}
          #disdat contient deja seulement tname ==pret. Donc disdat contient toutes les observations pour tname=pret. on a donc toutes les valeurs de t, pas seulement t et t+1
        disdat <- data[(data[, tname] == pret |data[, tname] == tlist[t + 1] | data[, tname] == tlist[t]), ]#}
        disdat <- makeBalancedPanel(disdat, idname, tname)##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1]). Doc: https://bcallaway11.github.io/BMisc/reference/makeBalancedPanel.html
         ##This function drops observations from data.frame that are not part of balanced panel data set.
        ##disdat à ce point contient les observations aux années flist[f]-1, tlist[t], tlist[t+1]. Si f=2 et t=6, on a les observations années 3,6,7.
        # #dim
        # [1] 19791    28
        # [1] 13194    28
        # ...


        disdat$C <- 1 * (disdat[, first.treat.name] == 0) 
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f]) 
       

        ############################
        # Assumptions on Weights ###
        ############################
          
        if(!is.null(weight_assumption)){
            
            ################################
            # Assumption 7: Missing trends #
            ################################
          
            if (weight_assumption=="missing_trends"){
                
                #######################
                # Calcule de matrices #
                #######################
                
                
                if (t>=debT){ #Pour exclure les valeurs avant le début du traitement. 
                if (flag==FALSE){ # The values are calculated only once.
                    unique_ids <- unique(disdat[[idname]])
                    
                    
                    #Matrice Gg
                    id_gg <- data.frame(id = unique(disdat[[idname]]), annee_G = disdat[[first.treat.name]][match(unique_ids, disdat[[idname]])])
                    ggi <- matrix(0, nrow = nrow(id_gg), ncol = finT)
                    for (i in 1:nrow(id_gg)) {
                      ggi[i, id_gg[i, "annee_G"]] <- 1}
                    ggi <- as.data.frame(ggi)   
                    names(ggi) <- 1:finT
                    ggi <- ggi[, 1:finT-1, drop = FALSE] 
                    ggi = as.matrix(ggi) # dim 6594x8 (id x finT)
                    # ggi <- ggi[, debT:ncol(ggi)]
                    
                    # Matrice cit
                    cit <- disdat[disdat[[tname]]==tlist[debT],]$C
                    cit = as.matrix(cit) # dim 6594x1 (id x 1)

                    # Matrice pour xit et propensityformla
                    lenX=length(propensityformla) #"X"
                    #xit porte a confusion ajouter la 3e dim dans le nom de la var.
                    # xit <- array(0, dim = c(6597, finT, lenX)) #On a donc une dimension pour chacune des éléments de propensityformla
                    # nrow(cit) for the number of rows
                    xit <- array(0, dim = c(nrow(cit), finT, lenX)) #On a donc une dimension pour chacune des éléments de propensityformla

                    for (i in 1:lenX) {
                         value <- pivot_wider(data = data[,c(idname,tname,propensityformla[i])], 
                         
                         id_cols = id, 
                         names_from = annee, 
                         values_from = propensityformla[i])    
                        
                         value <- value[,  -1] #remove id column                
                         value=value[, order(colnames(value))] #NxT  
                         value=as.matrix(value) #NxT
                         xit[, , i] <-value
                         }
                    
                    flag=TRUE 
                    }

                st1[,t] = disdat[disdat[[tname]]==tlist[t],][,c(ponderation)]  # Easier to troubleshoot with [,t]
                st2[,t] = disdat[disdat[[tname]]==tlist[t+1],][,c(ponderation)]
                ssit[,t] <- matrix(unlist(st1[,t]*st2[,t]), ncol = 1) #dim 6594x1
                
                ##############
                # numérateur #
                ##############
                # repeated_ssi <- matrix(rep(ssit[,t], ncol(ggi)), ncol = ncol(ggi), byrow = TRUE) # nom=mean(repeated_ssi * ggi) 
                nom=mean(ssit[,t] * ggi[,t]) # 1x1. mean(list). Je prend [,t], sinon le calcul n'est pas intuitif. Voir ci-dessous.
                #colMeans

                #Une autre option est de conserver le code initial pour le nominateur. (similaire a christophe)
                # gg  = disdat[disdat[[tname]]==tlist[t],][,c(nom_traitement)]
                # nom = mean(gg*st1*st2)

                ################
                # denominateur #
                ################

                # Le calcule des poids ne se fait qu'avec logit pour le moment. 
                # Sinon, il faut ajuster les calculs pour les fonctions d'influence.

                # Espérance de Gg
                mean_gg=mean(ggi[,tlist[t]]) # à g -> 1x1

                #Denominateur G                           
                #desing_matrix ya tout le monde. 
                design_matrix <- cbind(ggi[, debT:ncol(ggi)], xit[,,]) #Je garde seulement t apres debT. Ils servent à rien, et les betas sont 0.
                model <- glm(ssit[,t] ~ ., data = data.frame(design_matrix), family = binomial(link = "logit"))
                betas <- coef(model) #dim Gx1. #pour le moment les betas sur Xit sont les mêmes car X est constant dans le temps pour chq i.
                
                if (any(is.na(betas))) {
                  warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",   flist[f], " at time period: ", tlist[t +1]))}
                pred_ssit_g <- predict(model, newdata = data.frame(design_matrix), type = "response") #Nx1
                denom_g=mean(pred_ssit_g)*mean_gg #Rappel: mean_gg=mean(ggi[,tlist[t]]) # à g -> 1x1. À valider avec David.
                
                # #Denominateur C (Tu peux enlver ca, l'espérance devrait être la meme pour G et C) Utiliser le même denom
                #retirer design matrix (David spécification. pour les Betas cest diff que si tu mets design matrix ici. Tu peux prendre ligne 206, sinon set les Ggi à 0)
                #Si j'ai ggi xi avec des ggi 0 partout, ça veux dire que tu entraine les betas sur des groupe de controle donc ca serait mieux d'utiliser le design matrix d'avant, reset les variables.
                design_matrix <- cbind(cit, xit[,,]) #used in logit  
                model <- glm(ssit[,t] ~ ., data = data.frame(design_matrix), family = binomial(link = "logit"))
                betas <- coef(model) #dim Gx1
                if (any(is.na(betas))) {
                  warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",   flist[f], " at time period: ", tlist[t +1]))}
                pred_ssit_c <- predict(model, newdata = data.frame(design_matrix), type = "response") #Nx1
                denom_c=mean(pred_ssit_c)*mean_gg #pred-> un nombre (car tu fais la moyenne d'une liste) x un chiffre
                
                # # Calcul des poids pour i,g,t.
                git=1-cit
                pp2 = (nom/denom_g) * git 
                # pp1= (nom/denom_g) * cit #unmute ca et valider ca devrait donner la meme affaire pour les c
                pp1= (nom/denom_c) * cit #nom/denom 1x1 et ggi NxT.Cohérant, tu estime StSt+1 pour chq i et tu fais la moyenne.
                pp0=pp2+pp1 # Nx1, le poids est le meme pour chaque i, cohérant puisqu'on fait des moyenne à t.
                
                #normal que les poids soient les memes pour chaque individu (+ ou -). oui St depend pas de X.
                }

            }
          

          else {pp0 <- matrix(1, nrow = unique(disdat[[idname]]), ncol = 1)}
        }
          else {pp0 <- matrix(1, nrow = unique(disdat[[idname]]), ncol = 1)}

        #Delta Y
        
        disdat <- panelDiffV(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],tlist[t+1],selection) #deltaY
        #Masking pour disdat[année==pret,]$delta_y=yt+1-yt. Disdat contient seulement les obs pour l'année pret.
        # View(disdat)

        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
        pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
        ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
        disdat$pscore<-pscore
        disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
        disdat[,c("pscore")]=list(NULL)
        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur

        #####################
        # Calcul des ATTs   #
        #####################

        Dnom_outcome<-paste0("D",nom_outcome) 
        dy <- disdat[,Dnom_outcome] #delta y
        #On refait les mêmes manipulations sur disdat en raison de l'overwriting de disdat
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) 
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        NN<-dim(disdat)[1] ## d�finition de N pour boot
        disdat$pp=disdat[[ponderation[1]]] #hard coded...meme chose pour gmm... Ok, puisque la simulation utilisait 2 colonnes de pondération et elles sont identiques. 
        
        # Calcul de matrices pour calcul strate 
        traite<-cbind(disdat$G,disdat$C)
        traiteS<- traite #G,C
        devant<-(traiteS[,1]-traiteS[,2]) 
        devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2] #1 ou pscore/(1-pscore) selon si G ou C. Partie des equations Wg ou Wc.
        
        traite_boot<-traiteS
        
        dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,ponderation]))
        colnames(dommk)<- ponderation #renommer la colonne, sinon pp_noDenom aura des NaN.
        Denom <- 1/(dommk) 
        colnames(Denom)<- ponderation #renommer la colonne, sinon pp_noDenom aura des NaN.
        
        pp<-devant*(traiteS%*%Denom)*(devant2*disdat[,ponderation]) 
        
        
        
        
        
        #if pp0 exists, we reweight the ATT
        
        if(!is.null(weight_assumption)){
          if (weight_assumption=="missing_trends"){
            
          if (exists("pp0")) {
            pp <- pp * pp0
        } else {
        
        }
  # Handle the case where pp0 doesn't exist, e.g., assign a default value
  # pp <- pp * some_default_value
      
        }}



        # pp=pp*pp0 }
        
        
        
        pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,ponderation])       

        ### ATT ###
        att<- colSums(pp * dy) 

        ########################################################
        ### Bootstrap  : calcul de la fonction d'influence de la brique
        # matrice de covariance pour calculer ATT        
        ########################################################
        esperance_dy<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * dy))
        esperance_dy<- (-1) * esperance_dy * disdat$C + 1 * esperance_dy * disdat$G
        
        psi<-NN*pp *(dy-esperance_dy) #dim 6594x1 (il y a 6594 indiv)
        

        colnames(psi)=nom_outcome
        x <- model.matrix(xformla, data = disdat)
        n <- nrow(disdat)

        for (i in 1:length(nom_outcome)){
          G_i <- disdat$G*disdat[[ponderation[i]]]
          C_i <- disdat$C*disdat[[ponderation[i]]]
          wc1_i <- -1*C_i * pp_noDenom 
          wc_i <- wc1_i/mean(wc1_i) 
          dy_i <- disdat[[Dnom_outcome[i]]]

          M_i <- as.matrix(apply(as.matrix((C_i/(1 - pscore))^2 * gg(x, thet) * (dy_i - mean(wc_i * dy_i)) * x), 2,mean)/mean(wc1_i))
          
          #gg c'est une fonction qui retourne gval <- 1/((1 + exp(x %*% thet))^2). 
          A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
          A1_i <- (t(A1_i * x) %*% x/n)
          

          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x
          A_i <- A2_i %*% MASS::ginv(A1_i)
          
          correction_i<-A_i %*% M_i
          correction_i<-correction_i*disdat$C

          psi[,nom_outcome[i]]<-subset(psi,select=nom_outcome[i])-correction_i
        }






        
        ### mettre les siren dans psi
        psi<-as.data.frame(psi)
        colnames(psi)<-nom_outcome
        # View(psi)
        # psi[,idname]<-disdat[,..idname]
        psi[,idname]<-disdat[,idname] #dim 6594x1
        
        ### apparier psi avec indiv
        psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
        for (i in 1:length(nom_outcome)){
          psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
        }#dim 6594x1

      
        
        
  
        nG<-round(mean(colSums(as.matrix(disdat$G * 1 *(disdat[,ponderation]>0)))))
        nC<-round(mean(colSums(as.matrix(disdat$C * 1 *(disdat[,ponderation]>0)))))

        #fin de la patch
        # else  {
        #   att<- colSums(0 * dy) #filled with 0
        #   # Nombre d'observations utilis�es pour l'estimation de chaque brique
        #   nG<-mean(colSums(0 * dy))
        #   nC<-mean(colSums(0 * dy)) 
        #   psi<-as.data.frame(indiv[,idname])
        #   for (i in 1:length(nom_outcome)){
        #     psi[,nom_outcome[i]]<-0 
        #   }
          
        #   }

        
        fatt[counter,first.treat.name]<-flist[f] 
        fatt[counter,tname]<-tlist[(t + 1)]        
        fatt[counter,nom_outcome]<-att
        fatt[counter,"nobsG"]<-nG
        fatt[counter,"nobsC"]<-nC
        ### Remplissage de la matrice des fonctions d'influence
        
        
        for (i in 1:dim(psi)[2]){ #i=1,2,1,...
          
          mat_influence[,i,counter]=psi[,i]   #counter C'est g,t = 42
        }  #6597x1x42
        #donc on loop pour tous les individus, i et tout counter=42
    
    

    counter <- counter + 1
    
    


    } ### fin boucle sur l'annee



  
  } ## fin boucle sur le traitement f
  
  


  
  # View(mat_influence)
  list(fatt,mat_influence,indiv)
}



























##### ---------------------------------------------------------------------------------------------------------------------------
##### compute.mp.spatt.Boot-----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------

# Fonction compute.mp.spatt.Boot o� on enl�ve les "not yet treated" du groupe de contr�le !
compute.mp.spatt.Boot <- function ( nom_outcome
                                 , nom_traitement
                                 , flen  ## nombre de cohortes dans le traitement 
                                 , tlen  ## nombre d'annee
                                 , flist ## Le vecteur des vrais traitements (sans 0 !!) (traitements G commence 3:8)
                                 , tlist ## vecteur ordonn? des dates
                                 , data ### le panel balanc?
                                 , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donné
                                 , xformla ## le mod?le de score
                                 , tname ## la variable de date
                                 , w ## ??  NULL par d?faut
                                 , idname ## identifiant individu
                                 , method ## par d?faut "logit"
                                 , seedvec ## par d?faut NULL
                                 , se ## TRUE par d?faut
                                 , pl ## FALSE par d?faut
                                 , cores ## par d?faut = 2
                                 , printdetails ## TRUE par defaut
                                 , selection ## indicatrice observations valides
                                 , ponderation ## poids utilise  
                                  #poid du papier david
                                 , strate  ## variable qui donne les classes d'observations, a utiliser
                                 , POND_RD=NULL ## POND_RD : � compl�ter �ventuellement si on utilise une pond�ration suppl�mentaire (enqu�te RD) 
)
{

  


  nbligne=flen*(tlen - 1)
  rempli<-seq(1:nbligne)
  fatt <- data.frame(compteur=rempli)
  fatt[,first.treat.name]<-rempli
  fatt[,tname]<-rempli
  fatt$nobsG<-rempli
  fatt$nobsC<-rempli 

  for(i in 1:length(nom_outcome)){
    fatt[,nom_outcome[i]]<-rempli
  }
  
  ### Liste des individus impliqués dans l'estimation pour le wild bootstrap
  indiv<-unique(data[,c(first.treat.name,idname)]) #liste des observations uniques (id) et leur première année de traitement
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
  counter <- 1

  for (f in 1:flen) {  ### boucle sur les cohortes dans le traitement : f
    for (t in 1:(tlen - 1)) {  ## boucle sur la date : t
        pret <- flist[f]-1
        
        #### Sélection de la base de calcul avec éventuelle prise en compte des pondérations de l'enquete RD 
        if (is.null(POND_RD)){
          disdat <- data[(data[, tname] == pret), ]
        } else {
          nom_pondRD<-paste0(POND_RD,as.character(pret))
          
          data2<-data
          data2[,ponderation]<-data2[,nom_pondRD]*data2[,ponderation]
          disdat <- data2[(data2[, tname] == pret), ] 
        }
        
        disdat <- droplevels(disdat) 
        ### Controle de l'existence de la variable de strate 
        disdat<-disdat[is.na(disdat[,strate])==FALSE,]
        ##Indicatrice de Traitement et de contrefactuel (non traité et non encore traités)
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) 
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        if (is.null(xformla)) {
          xformla <- ~1
        }
        pformla <- xformla
        ### On supprime les observations pour lesquelles le modèle de score ne peut être défini
        LesX<-BMisc::rhs.vars(pformla)
        
        bbb<-length(LesX)
        for (jj in 1:bbb){
          disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,] #keep only if no missing values
        }
        
        ## finlament définitio ndu modèle de score
        pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla)) #G~X

        base_score<-subset(disdat, C + G == 1) 
        base_score=base_score[base_score$select==1,] 

        pscore.reg <- glm(pformla, family = binomial(link = "logit"),
                          data = base_score)
        
        thet <- coef(pscore.reg)
        if (any(is.na(thet))) {
          warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
                         flist[f], " at time period: ", tlist[t +1]))
        }
        if (is.null(POND_RD)){
          disdat <- data[(data[, tname] == pret |data[, tname] == tlist[t + 1] | data[, tname] == tlist[t]), ]
        } else{
          disdat <- data2[(data2[, tname] == pret |data2[, tname] == tlist[t + 1] | data2[, tname] == tlist[t]),]
        }
        
        
        disdat <- makeBalancedPanel(disdat, idname, tname)##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1])
        
        ##This function drops observations from data.frame that are not part of balanced panel data set.
        ##on élimine les observations pour lesquelles on a pas des valeurs pour chaque période de temps.
        disdat <- panelDiffV(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],tlist[t+1],selection)
         
        #muted
        # disdat$P_Y1_chaine <- disdat$P_Y1_chaine[[1]]
        # disdat$P_Y2_chaine <- disdat$P_Y2_chaine[[1]]

        # disdat=data.table(disdat)
        disdat<-disdat[is.na(disdat[,strate])==FALSE,]
        
        ### base de calcul
        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
        pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
        
        ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
        disdat$pscore<-pscore
        disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
        disdat[,c("pscore")]=list(NULL)
        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur

        Dnom_outcome<-paste0("D",nom_outcome) 
        dy <- disdat[,Dnom_outcome] #delta y
        
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t+1])
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        ## d�finition de N pour boot
        NN<-dim(disdat)[1]
    
        ### Calcul par strate 
        ### on ne conserve que les strates qui ont un poid nul et un score non nul
        disdat$pp=disdat[[ponderation[1]]]
        #pondération = Pond_outcome : vecteur des variabels permettnat de pondérer les observations / Pout_chaine=c("P_Y1_chaine","P_Y2_chaine") créé à partir de St_chaine
        for(i in 1:length(ponderation)){
          disdat$pp<-pmin(disdat[,ponderation[i]],disdat[,"pp"])
        }
        
        
        # strate_list=unique(disdat[get(first.treat.name) == flist[f] & pp > 0, .(strate)])
        strate_list<-unique(disdat[disdat[,first.treat.name]==flist[f] & disdat$pp>0,strate])
        #patch pour reproduire le meme résultat...
        strate_list <- data.frame(strate = strate_list)
        strate_list<-sort(strate_list$strate) 
        strate_list<-strate_list[is.na(strate_list)==FALSE]
        strate_nb<-length(strate_list)
        
        ### Patch pour pr�venir les cas o� aucune strate n'est active
        if (strate_nb>0){
        strate_poids<-seq(1:strate_nb)
        for(riri in 1:strate_nb){
          #lu = disdat$G
          lu<-1*(disdat[,strate]==strate_list[riri])*disdat$G #on ne compte que les observations traitées
          strate_poids[riri]<-sum(lu,na.rm=TRUE) #nombre d'obs dans la strate
          
        }
        strate_poids<-strate_poids/sum(strate_poids,na.rm=TRUE) #poids de la strate sur toutes les strates.

        # Calcul de matrices pour calcul strate 
        traite<-cbind(disdat$G,disdat$C)
        for(riri in 1:strate_nb){
          if (riri==1){

            traiteS<- 1*(disdat[,strate]==strate_list[riri])*traite #G,C

            devant<-strate_poids[1]*(traiteS[,1]-traiteS[,2]) 
            
            devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2] #1 ou pscore/(1-pscore) selon si G ou C. Partie des equations Wg ou Wc.
            traite_boot<-traiteS
          
          } else {
            
            ttt<- 1*(disdat[,strate]==strate_list[riri])*traite 
            traiteS<- cbind(traiteS,ttt)  
            ddd<-strate_poids[riri]*(ttt[,1]-ttt[,2])
            
            devant<-devant+ddd
            ddd2<-ttt[,1]+(pscore/(1 - pscore))*ttt[,2]
            devant2<-devant2+ddd2
            traite_boot<-cbind(traite_boot,traiteS) #traite_boot = matrice de traiteS
          }}

         
         #   #garder en tete que les probabilités les formules sont hard codés par la simultation
        #   # 1/exp... = 1/(1+exp(-x))
      
        dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,ponderation]))
        colnames(dommk)<- ponderation #renommer la colonne, sinon pp_noDenom aura des NaN.
        Denom <- 1/(dommk) 
        colnames(Denom)<- ponderation #renommer la colonne, sinon pp_noDenom aura des NaN.
        pp<-devant*(traiteS%*%Denom)*(devant2*disdat[,ponderation]) 
        pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,ponderation]) 
        
        ### ATT
        att<- colSums(pp * dy) 
        ### Bootstrap  : calcul de la fonction d'influence de la brique
        # matrice de covariance pour calculer ATT        
        esperance_dy<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * dy))
        esperance_dy<- (-1) * esperance_dy * disdat$C + 1 * esperance_dy * disdat$G
        
        psi<-NN*pp *(dy-esperance_dy)

        colnames(psi)=nom_outcome
        x <- model.matrix(xformla, data = disdat)
        n <- nrow(disdat)

        for (i in 1:length(nom_outcome)){
          G_i <- disdat$G*disdat[[ponderation[i]]]
          C_i <- disdat$C*disdat[[ponderation[i]]]
          wc1_i <- -1*C_i * pp_noDenom 
          wc_i <- wc1_i/mean(wc1_i) 
          dy_i <- disdat[[Dnom_outcome[i]]]

          M_i <- as.matrix(apply(as.matrix((C_i/(1 - pscore))^2 * gg(x, thet) * (dy_i - mean(wc_i * dy_i)) * x), 2,mean)/mean(wc1_i))
          #gg c'est une fonction qui retourne gval <- 1/((1 + exp(x %*% thet))^2)
          A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
          A1_i <- (t(A1_i * x) %*% x/n)
          

          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x
          A_i <- A2_i %*% MASS::ginv(A1_i)
          


         
          # View(M_i)
          correction_i<-A_i %*% M_i
    
          correction_i<-correction_i*disdat$C

          
          
          psi[,nom_outcome[i]]<-subset(psi,select=nom_outcome[i])-correction_i
        }






        
        ### mettre les siren dans psi
        psi<-as.data.frame(psi)
        colnames(psi)<-nom_outcome

        # psi[,idname]<-disdat[,..idname]
        psi[,idname]<-disdat[,idname]
        ### apparier psi avec indiv
        psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
        for (i in 1:length(nom_outcome)){
          psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
        }

        # Compter le nombre d'observations utilis�es pour l'estimation de chaque brique (moyenne sur les diff�rents outcomes �valu�s)
        # nG<-round(mean(colSums(disdat$G * 1 *(disdat[,..ponderation]>0))))
        # nC<-round(mean(colSums(disdat$C * 1 *(disdat[,..ponderation]>0))))
        
 
        
        
        
        nG<-round(mean(colSums(as.matrix(disdat$G * 1 *(disdat[,ponderation]>0)))))
        nC<-round(mean(colSums(as.matrix(disdat$C * 1 *(disdat[,ponderation]>0)))))
        #fin de la patch
      
        } else  {
          att<- colSums(0 * dy) #filled with 0
          # Nombre d'observations utilis�es pour l'estimation de chaque brique
          nG<-mean(colSums(0 * dy))
          nC<-mean(colSums(0 * dy)) 
          psi<-as.data.frame(indiv[,idname])
          for (i in 1:length(nom_outcome)){
            psi[,nom_outcome[i]]<-0 
          }
          
          }

        
        fatt[counter,first.treat.name]<-flist[f] 
        fatt[counter,tname]<-tlist[(t + 1)]        
        fatt[counter,nom_outcome]<-att
        fatt[counter,"nobsG"]<-nG
        fatt[counter,"nobsC"]<-nC
        ### Remplissage de la matrice des fonctions d'influence
        for (i in 1:dim(psi)[2]){
          mat_influence[,i,counter]=psi[,i]  
        } 
    
    counter <- counter + 1
    
    


    } ### fin boucle sur l'annee



  
  } ## fin boucle sur le traitement f
  
  
  

  list(fatt,mat_influence,indiv)
}











##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR GMM -----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------
compute.mp.spatt.GMM <- function (  nom_outcome
                                    , nom_traitement
                                    , flen  ## nombre de cohortes dans le traitement 
                                    , tlen  ## nombre d'annee
                                    , flist ## Le vecteur des vrais traitements (sans 0 !!)
                                    , tlist ## vecteur ordonn? des dates
                                    , data ### le panel balanc?
                                    , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donn�
                                    , xformla ## le mod?le de score
                                    , propensityformla
                                    , tname ## la variable de date
                                    , w ## ??  NULL par d?faut
                                    , idname ## identifiant individu
                                    , method ## par d?faut "logit"
                                    , seedvec ## par d?faut NULL
                                    , se ## TRUE par d?faut
                                    , pl ## FALSE par d?faut
                                    , cores ## par d?faut = 2
                                    , printdetails ## TRUE par defaut
                                    , selection ## indicatrice observations valides
                                    , ponderation ## poids utilise
                                    , weight_assumption
                                    , debT  ## variable qui donne les classes d'observations, a utiliser 
                                    , finT=NULL ## POND_RD : � compl�ter �ventuellement si on utilise une pond�ration suppl�mentaire (enqu�te RD) 
)


{
  ################################
  # initialisation des matrices  #
  ################################
  # remarque : le nombre possible de Delta ATT pour un g donné est de T(T-1)/2.
  nbligne=flen*(tlen - 1)*tlen/2 #168
  rempli<-seq(1:nbligne)

  

  fatt <- data.frame(compteur=rempli)
  fatt[,first.treat.name]<-rempli
  fatt[,tname]<-rempli
  fatt$nobsG<-rempli
  fatt$nobsC<-rempli 
  for(i in 1:length(nom_outcome)){
    fatt[,nom_outcome[i]]<-rempli}
  
  # Matrix W containing the weights to go from ATT(t) to \Delta(ATT). The nb of possible ATT(g,t) that can be identified is flen*(tlen-1)
  # but we fill in flen*tlen ATT(g,t), given the extra columns are deleted later on
  mat_W <-as.data.frame(matrix(0,nbligne,flen*tlen)) #168x42
  for (f in 1:flen) { 
    for (t in 1:tlen) {
      colnames(mat_W)[(f-1)*tlen+t]<-paste0("ATT(",flist[f],",",tlist[t],")")
    }
  }
  
  ### Liste des individus impliqués dans l'estimation pour le wild bootstrap
  indiv<-unique(data[,c(first.treat.name,idname)])
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
  counter <- 1
  
  
  
  
  
  
  for (f in 1:flen) {  ### boucle sur les cohortes dans le traitement : f
    pret <- flist[f]-1
    for (t in 1:(tlen - 1)) {  ## boucle sur la date : t
      tbis = tlen - t #detlta att_k (les distances possibles)
      
      for (k in 1:tbis) {
        
        #weights for ATT
        ### Fill the matrix W for ATT with positive weight (1), ATT with negative weight (-1), or ATT with zero weight (by default)
        #Possible d'améliorer avec un calcul matricielle. (à venir)
        att_pos_w = (f-1)*tlen+(t+k)
        att_neg_w = (f-1)*tlen+(t)
        mat_W[counter,att_pos_w]<-  (1)
        mat_W[counter,att_neg_w]<- (-1)
        
        disdat <- data[(data[, tname] == pret), ]
        disdat <- droplevels(disdat)  ##suppression d'observations manquantes ??
        
        ##Indicatrice de Traitement et de contrefactuel (non traité et non encore traités)
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t+k])
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        if (is.null(xformla)) {
          xformla <- ~1
        }
        pformla <- xformla
        ### On supprime les observations pour lesquelles le mod�le de score ne peut �tre d�fini
        LesX<-BMisc::rhs.vars(pformla)
        bbb<-length(LesX)
        for (jj in 1:bbb){
          disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,]
        }
        
        ## finlament d�finitio ndu mod�le de score
        pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))
        base_score<-subset(disdat, C + G == 1)
        base_score=base_score[base_score$select==1,]
        
        pscore.reg <- glm(pformla, family = binomial(link = "logit"),
                          data = base_score)
        
        thet <- coef(pscore.reg)
        
        if (any(is.na(thet))) {
          warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
                         flist[f], " at time period: ", tlist[t +k]))
        }
        
        disdat <- data[(data[, tname] == pret |data[, tname] == tlist[t + k] | data[, tname] == tlist[t]), ]
        ##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1])
        disdat <- makeBalancedPanel(disdat, idname, tname)
        disdat <- panelDiffV_GMM(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],tlist[t+k],selection)
        # disdat=as.data.table(disdat)
        
        ### base de calcul
        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
        pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
        
        ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
        disdat$pscore<-pscore
        disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
        disdat[,c("pscore")]=list(NULL)
        pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur
        

        #####################
        # Calcul des ATTs   #
        #####################
        #appliquer le corrolaire 1 à ce contexte.

        Dnom_outcome<-paste0("D",nom_outcome)
        dy <- disdat[,Dnom_outcome]
        
        ##Indicatrice de Traitement et de contrefactuel
        disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t+k])
        disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
        NN<-dim(disdat)[1]
        disdat$pp<-disdat[[ponderation[1]]]

        for(i in 1:length(ponderation)){
          disdat$pp<-pmin(disdat[,ponderation[i]],disdat[,"pp"])}
      
        traite<-cbind(disdat$G,disdat$C)
        traiteS<-traite
        devant<-(traiteS[,1]-traiteS[,2])
        devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2]
        traite_boot<-traiteS

        dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,ponderation]))
        Denom <- 1/(dommk)
        pp<- devant*(traiteS%*%Denom)*(devant2*disdat[,ponderation])
        pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,ponderation])
        att<- colSums(pp * dy)
          
        
        ########################################################
        ### Bootstrap  : calcul de la fonction d'influence de la brique
        # matrice de covariance pour calculer ATT        
        ########################################################
        esperance_dy<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * dy))
        esperance_dy<- (-1) * esperance_dy * disdat$C + 1 * esperance_dy * disdat$G
        psi<-NN*pp *(dy-esperance_dy)
        colnames(psi)=nom_outcome
        x <- model.matrix(xformla, data = disdat)
        n <- nrow(disdat)
        
        for (i in 1:length(nom_outcome)){
          
          # G_i <- disdat$G*(subset(disdat,select=ponderation[i]))
          # C_i <- disdat$C*(subset(disdat,select=ponderation[i]))

          G_i <- disdat$G*disdat[[ponderation[i]]]
          C_i <- disdat$C*disdat[[ponderation[i]]]
        
          # G_i=as.vector(G_i)
          # C_i=as.vector(C_i)
          # G_i=unlist(G_i)
          # C_i=unlist(C_i)

          # wc1_i <- (-1)*C_i * pp_noDenom[[i]]
          wc1_i <- (-1)*C_i * pp_noDenom
          wc_i <- wc1_i/mean(wc1_i)

          # dy_i = subset(disdat,select=Dnom_outcome[i])
          # dy_i=as.vector(dy_i)
          # dy_i=unlist(dy_i)
          dy_i <- disdat[[Dnom_outcome[i]]]

          M_i <- as.matrix(apply(as.matrix(x * (C_i/(1 - pscore))^2 * gg(x, thet) * (dy_i - mean(wc_i * dy_i))), 2,mean)/mean(wc1_i))
          A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
          A1_i <- (t(A1_i * x) %*% x/n)
          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x
          A_i <- A2_i %*% MASS::ginv(A1_i)

          correction_i<-A_i %*% M_i
          correction_i<-correction_i*disdat$C
          
          
          
          psi[,nom_outcome[i]]<-subset(psi,select=nom_outcome[i])-correction_i
        }
          
          ### mettre les siren dans psi
          psi<-as.data.frame(psi)
          colnames(psi)<-nom_outcome
          psi[,idname]<-disdat[,idname]
          
          ### apparier psi avec indiv
          psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
          
          for (i in 1:length(nom_outcome)){
            psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
          }
          
          # Compter le nombre d'observations utilis�es pour l'estimation de chaque brique (moyenne sur les diff�rents outcomes �valu�s)
          # nG<-round(mean(colSums(disdat$G * 1 *(disdat[,ponderation]>0))))
          # nC<-round(mean(colSums(disdat$C * 1 *(disdat[,ponderation]>0))))
        



        #################################
        ###Probleme
        ################################
        #Probleme ici, nG et nC sont constant... pourtant la définition de G et de C dans la simulation depend de t et i.
        #Les ATT sont bons.
        nG<-round(mean(colSums(as.matrix(disdat$G * 1 *(disdat[,ponderation]>0)))))
        nC<-round(mean(colSums(as.matrix(disdat$C * 1 *(disdat[,ponderation]>0)))))
        # View(cbind(disdat$C,disdat[,ponderation]>0,disdat$C*disdat[,ponderation]>0,colSums(as.matrix(disdat$C*disdat[,ponderation]>0))))
        #Est ce que disdat à le bon subset sur t et sur pret? oui...

        ### le resultat de base est de type liste: le coefficient, le numero du traitement, la date, avant ou post
        fatt[counter,first.treat.name]<-flist[f]
        fatt[counter,tname]<-tlist[(t + k)]
        fatt[counter,nom_outcome]<-att
        fatt[counter,"nobsG"]<-nG
        fatt[counter,"nobsC"]<-nC
        
        # View(fatt)
        ### Remplissage de la matrice des fonctions d'influence
        for (i in 1:dim(psi)[2]){
          mat_influence[,i,counter]=psi[,i]  
        }
        counter <- counter + 1
        
      } ## fin boucle sur l'annee bis
    } ### fin boucle sur l'annee
  } ## fin boucle sur le traitement f


  list(fatt,mat_influence,indiv,mat_W,flist,tlist)
  
}






# ##### ---------------------------------------------------------------------------------------------------------------------------
# ##### BACK UP sand modifs POUR L'ESTIMATEUR GMM -----------------------------------------------------------------------------------------------------
# ##### ---------------------------------------------------------------------------------------------------------------------------
# compute.mp.spatt.GMM <- function (  nom_outcome
#                                     , nom_traitement
#                                     , flen  ## nombre de cohortes dans le traitement 
#                                     , tlen  ## nombre d'annee
#                                     , flist ## Le vecteur des vrais traitements (sans 0 !!)
#                                     , tlist ## vecteur ordonn? des dates
#                                     , data ### le panel balanc?
#                                     , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donn�
#                                     , xformla ## le mod?le de score
#                                     , tname ## la variable de date
#                                     , w ## ??  NULL par d?faut
#                                     , idname ## identifiant individu
#                                     , method ## par d?faut "logit"
#                                     , seedvec ## par d?faut NULL
#                                     , se ## TRUE par d?faut
#                                     , pl ## FALSE par d?faut
#                                     , cores ## par d?faut = 2
#                                     , printdetails ## TRUE par defaut
#                                     , selection ## indicatrice observations valides
#                                     , ponderation ## poids utilise
#                                     , strate  ## variable qui donne les classes d'observations, a utiliser 
#                                     , POND_RD=NULL ## POND_RD : � compl�ter �ventuellement si on utilise une pond�ration suppl�mentaire (enqu�te RD) 
# )
# {
#   ################################
#   # initialisation des matrices  #
#   ################################
#   # remarque : le nombre possible de Delta ATT pour un g donné est de T(T-1)/2.
#   nbligne=flen*(tlen - 1)*tlen/2
#   rempli<-seq(1:nbligne)
#   fatt <- data.frame(compteur=rempli)
#   fatt[,first.treat.name]<-rempli
#   fatt[,tname]<-rempli
#   fatt$nobsG<-rempli
#   fatt$nobsC<-rempli 
#   for(i in 1:length(nom_outcome)){
#     fatt[,nom_outcome[i]]<-rempli}
  
#   # Matrix W containing the weights to go from ATT(t) to \Delta(ATT). The nb of possible ATT(g,t) that can be identified is flen*(tlen-1)
#   # but we fill in flen*tlen ATT(g,t), given the extra columns are deleted later on
#   mat_W <-as.data.frame(matrix(0,nbligne,flen*tlen))
#   for (f in 1:flen) { 
#     for (t in 1:tlen) {
      
#       colnames(mat_W)[(f-1)*tlen+t]<-paste0("ATT(",flist[f],",",tlist[t],")")
      
#     }
#   }
  
#   ### Liste des individus impliqu�s dans l'estimation pour le wild bootstrap
#   indiv<-unique(data[,c(first.treat.name,idname)])
  
#   mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
#   # print(dim(mat_influence))
#   counter <- 1
  
#   for (f in 1:flen) {  ### boucle sur les cohortes dans le traitement : f
    
#     pret <- flist[f]-1
    
#     for (t in 1:(tlen - 1)) {  ## boucle sur la date : t
      
#       tbis = tlen - t
      
#       for (k in 1:tbis) {
        
#         #weights for ATT
#         ### Fill the matrix W for ATT with positive weight (1), ATT with negative weight (-1), or ATT with zero weight (by default)
#         att_pos_w = (f-1)*tlen+(t+k)
#         att_neg_w = (f-1)*tlen+(t)
#         mat_W[counter,att_pos_w]<-  (1)
#         mat_W[counter,att_neg_w]<- (-1)
        
#         #### S�lection de la base de calcul avec �ventuelle prise en compte des pond�rations de l'enqu�te RD 
#         if (is.null(POND_RD)){
#           disdat <- data[(data[, tname] == pret), ]
#         } else {
#           nom_pondRD<-paste0(POND_RD,as.character(pret))
#           data2<-data
#           data2[,ponderation]<-data2[,nom_pondRD]*data2[,ponderation]
#           disdat <- data2[(data2[, tname] == pret), ] 
#         }
        
#         disdat <- droplevels(disdat)  ##suppression d'observations manquantes ??
        
#         ### Contr�le de l'existence de la variable de strate 
        
#         disdat<-disdat[is.na(disdat[,strate])==FALSE,]
        
#         ##Indicatrice de Traitement et de contrefactuel (non trait� et non encore trait�s)
#         disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t+k])
#         disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
        
#         if (is.null(xformla)) {
#           xformla <- ~1
#         }
#         pformla <- xformla
#         ### On supprime les observations pour lesquelles le mod�le de score ne peut �tre d�fini
#         LesX<-BMisc::rhs.vars(pformla)
#         # print(LesX)
#         bbb<-length(LesX)
        
#         for (jj in 1:bbb){
#           disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,]
#         }
        
#         # print(dim(disdat))
        
#         ## finlament d�finitio ndu mod�le de score
#         pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))
      
#         base_score<-subset(disdat, C + G == 1)
#         base_score=base_score[base_score$select==1,]
        
#         pscore.reg <- glm(pformla, family = binomial(link = "logit"),
#                           data = base_score)
#         # print(summary(pscore.reg ))
#         thet <- coef(pscore.reg)
#         # print(thet)
#         if (any(is.na(thet))) {
#           warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
#                          flist[f], " at time period: ", tlist[t +k]))
#         }
        
#         if (is.null(POND_RD)){
#           disdat <- data[(data[, tname] == pret |data[, tname] == tlist[t + k] | data[, tname] == tlist[t]), ]
#         } else{
#           disdat <- data2[(data2[, tname] == pret |data2[, tname] == tlist[t + k] | data2[, tname] == tlist[t]),]
#         }
        
        
        
#         #pourquoi on cylindre pour GMM?
#         ##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1])
#         disdat <- makeBalancedPanel(disdat, idname, tname)
#         disdat <- panelDiffV_GMM(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],tlist[t+k],selection)
#         # print(table(disdat[, first.treat.name]))
        
#         ### Contr�le pr�alable de l'existence de la variable de strate 
        
#         # View(disdat[,strate])
#         # print(class(disdat[,strate]))
#         disdat=as.data.table(disdat)

#         disdat<-disdat[is.na(disdat[,strate])==FALSE,]
        
#         ### base de calcul
#         pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
#         pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
        
#         ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
#         disdat$pscore<-pscore
#         disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
#         disdat[,c("pscore")]=list(NULL)
#         pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur
        
#         Dnom_outcome<-paste0("D",nom_outcome)
#         dy <- disdat[,..Dnom_outcome]
        
#         ##Indicatrice de Traitement et de contrefactuel
#         disdat$C <- 1 * (disdat[, ..first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t+k])
#         disdat$G <- 1 * (disdat[, ..first.treat.name] == flist[f])
        
#         ## d�finition de N pour boot
#         NN<-dim(disdat)[1]
        
#         ### Calcul par strate 
#         ### on ne conserve que les strates qui ont un poid nul et un score non nul
#         # disdat$pp<-disdat[,ponderation[1]]
#         disdat$pp<-disdat[[ponderation[1]]]
#         for(i in 1:length(ponderation)){
#           disdat$pp<-pmin(disdat[,ponderation[i]],disdat[,"pp"])
#         }
#         # strate_list<-unique(disdat[disdat[,first.treat.name]==flist[f] & disdat$pp>0,strate])
#         strate_list <- unique(disdat[get(first.treat.name) == flist[f] & pp > 0 & !is.na(strate), .(strate)])
        
#         # strate_list<-sort(strate_list)
#         # strate_list<-strate_list[is.na(strate_list)==FALSE]
#         strate_list <- strate_list[order(strate)]
#         strate_nb<-length(strate_list)
        
        

#         ### Patch pour pr�venir les cas o� aucune strate n'est active
#         if (strate_nb>0){
#           strate_poids<-seq(1:strate_nb)
#           for(riri in 1:strate_nb){
#             # lu<-1*(disdat[,strate]==strate_list[riri])*disdat$G
#             lu <- 1 * as.vector(disdat[[strate]] == strate_list[riri]) * disdat$G
#             strate_poids[riri]<-sum(lu,na.rm=TRUE)
#           }
          
#           strate_poids<-strate_poids/sum(strate_poids,na.rm=TRUE)
          
#           # Calcul de matrices pour calcul strate 
          
#           traite<-cbind(disdat$G,disdat$C)
#           for(riri in 1:strate_nb){
#             if (riri==1){
#               # traiteS<- 1*(disdat[,strate]==strate_list[riri])*traite
#               traiteS<- as.numeric(1*(disdat[,strate])==as.numeric(strate_list[riri]))*traite
#               devant<-strate_poids[1]*(traiteS[,1]-traiteS[,2])
#               devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2]
#               traite_boot<-traiteS
              
#             } else {
#               ttt<- 1*(disdat[,strate]==strate_list[riri])*traite
            
#               traiteS<- cbind(traiteS,ttt)
            
#               ddd<-strate_poids[riri]*(ttt[,1]-ttt[,2])
            
#               devant<-devant+ddd
            
            
#               ddd2<-ttt[,1]+(pscore/(1 - pscore))*ttt[,2]
            
#               devant2<-devant2+ddd2
            
#               traite_boot<-cbind(traite_boot,traiteS)
#             }}
          
#           dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,..ponderation]))
#           Denom <- 1/(dommk)
#           pp<- devant*(traiteS%*%Denom)*(devant2*disdat[,..ponderation])
#           pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,..ponderation])
#           att<- colSums(pp * dy)
          
#           ### Bootstrap  : calcul de la fonction d'influence de la brique
#           esperance_dy<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * dy))
#           esperance_dy<- (-1) * esperance_dy * disdat$C + 1 * esperance_dy * disdat$G
#           psi<-NN*pp *(dy-esperance_dy)
#           colnames(psi)=nom_outcome
#           x <- model.matrix(xformla, data = disdat)
#           n <- nrow(disdat)
#           for (i in 1:length(nom_outcome)){
#             # G_i <- disdat$G*disdat[,ponderation[i]]
#             # C_i <- disdat$C*disdat[,ponderation[i]]
#             G_i <- disdat$G*(subset(disdat,select=ponderation[i]))
#             C_i <- disdat$C*(subset(disdat,select=ponderation[i]))
          
#             G_i=as.vector(G_i)
#             C_i=as.vector(C_i)
#             G_i=unlist(G_i)
#             C_i=unlist(C_i)

#             # wc1_i <- (-1)*C_i * pp_noDenom[,i]
#             wc1_i <- (-1)*C_i * pp_noDenom[[i]]
#             wc_i <- wc1_i/mean(wc1_i)
#             # dy_i <- disdat[,Dnom_outcome[i]]
#             dy_i = subset(disdat,select=Dnom_outcome[i])
#             dy_i=as.vector(dy_i)
#             dy_i=unlist(dy_i)

#             M_i <- as.matrix(apply(as.matrix(x * (C_i/(1 - pscore))^2 * gg(x, thet) * (dy_i - mean(wc_i * dy_i))), 2,mean)/mean(wc1_i))
#             A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
#             A1_i <- (t(A1_i * x) %*% x/n)
#             A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x
#             A_i <- A2_i %*% MASS::ginv(A1_i)

#             correction_i<-A_i %*% M_i
#             correction_i<-correction_i*disdat$C
#             # psi[,nom_outcome[i]]<-psi[,nom_outcome[i]]-correction_i #est-ce que c'est la correction que David m'a présenté, pour le la sous-population soit repr.sentative de la population (Em)
#             psi[,nom_outcome[i]]<-subset(psi,select=nom_outcome[i])-correction_i
#           }
          
#           ### mettre les siren dans psi
#           psi<-as.data.frame(psi)
#           colnames(psi)<-nom_outcome
#           psi[,idname]<-disdat[,..idname]
          
#           ### apparier psi avec indiv
#           psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
          
#           for (i in 1:length(nom_outcome)){
#             psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
#           }
          
#           # Compter le nombre d'observations utilis�es pour l'estimation de chaque brique (moyenne sur les diff�rents outcomes �valu�s)
#           nG<-round(mean(colSums(disdat$G * 1 *(disdat[,..ponderation]>0))))
#           nC<-round(mean(colSums(disdat$C * 1 *(disdat[,..ponderation]>0))))
          
#         } else  {
          
#           att<- colSums(0 * dy)
#           # Nombre d'observations utilis�es pour l'estimation de chaque brique
#           nG<-mean(colSums(0 * dy))
#           nC<-mean(colSums(0 * dy))
      
#           psi<-as.data.frame(indiv[,idname])
#           for (i in 1:length(nom_outcome)){
#             psi[,nom_outcome[i]]<-0
#           }
          
#         }  
        
#         ### le resultat de base est de type liste: le coefficient, le numero du traitement, la date, avant ou post
#         fatt[counter,first.treat.name]<-flist[f]
#         fatt[counter,tname]<-tlist[(t + k)]
#         fatt[counter,nom_outcome]<-att
#         fatt[counter,"nobsG"]<-nG
#         fatt[counter,"nobsC"]<-nC
#         ### Remplissage de la matrice des fonctions d'influence
#         for (i in 1:dim(psi)[2]){
#           mat_influence[,i,counter]=psi[,i]  
#         }
#         counter <- counter + 1
        
#       } ## fin boucle sur l'annee bis
#     } ### fin boucle sur l'annee
#   } ## fin boucle sur le traitement f
  
#   list(fatt,mat_influence,indiv,mat_W,flist,tlist)
  
# }














































##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR EN CROSS SECTION ----------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------
compute.mp.spatt.CS_Boot <- function (  nom_outcome
                                      , nom_traitement
                                      , flen  ## nombre de cohortes dans le traitement 
                                      , tlen  ## nombre d'annee
                                      , flist ## Le vecteur des vrais traitements (sans 0 !!)
                                      , tlist ## vecteur ordonn? des dates
                                      , data ### le panel balanc?
                                      , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donn�
                                      , xformla ## le mod?le de score
                                      , tname ## la variable de date
                                      , w ## ??  NULL par d?faut
                                      , idname ## identifiant individu
                                      , method ## par d?faut "logit"
                                      , seedvec ## par d?faut NULL
                                      , se ## TRUE par d?faut
                                      , pl ## FALSE par d?faut
                                      , cores ## par d?faut = 2
                                      , printdetails ## TRUE par defaut
                                      , selection ## indicatrice observations valides
                                      , ponderation ## poids utilise
                                      , strate  ## variable qui donne les classes d'observations, a utiliser 
                                      , POND_RD=NULL ## POND_RD : � compl�ter �ventuellement si on utilise une pond�ration suppl�mentaire (enqu�te RD) 
)
{
  ### remplissage de l'objet en sortie : toujours utile ??
  nbligne=flen*tlen
  rempli<-seq(1:nbligne)
  fatt <- data.frame(compteur=rempli)
  fatt[,first.treat.name]<-rempli
  fatt[,tname]<-rempli
  fatt$nobsG<-rempli
  fatt$nobsC<-rempli 
  for(i in 1:length(nom_outcome)){
    fatt[,nom_outcome[i]]<-rempli
  }
  
  ### Liste des individus impliqu�s dans l'estimation pour le wild bootstrap + initialisation matrice d'influence
  indiv<-unique(data[,c(first.treat.name,idname)])
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
  
  counter <- 1
  
  for (f in 1:flen) {  ### boucle sur les cohortes dans le traitement : f
    
    for (t in 1:(tlen)) {  ## boucle sur la date : t
      pret <- flist[f]-1
      
      ## la base de calcul
      if (is.null(POND_RD)){
        disdat <- data[(data[, tname] == pret), ]
      } else {
        nom_pondRD<-paste0(POND_RD,as.character(pret))
        data2<-data
        data2[,ponderation]<-data2[,nom_pondRD]*data2[,ponderation]
        disdat <- data2[(data2[, tname] == pret), ] 
      }
      
      disdat <- droplevels(disdat)  ##suppression d'observations manquantes ??
      
      ### Contr�le de l'existence de la variable de strate 
      disdat<-disdat[is.na(disdat[,strate])==FALSE,]

      # on calcule le score sur toutes les non encore traitees une fois pour toute
      disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+  1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t])
      disdat$G <- 1 * (disdat[, first.treat.name] == flist[f]) ## les traitees f
      
      if (is.null(xformla)) {
        xformla <- ~1
      }
      pformla <- xformla
      ### On supprime les observations pour lesquelles le mod�le de score ne peut �tre d�fini
      LesX<-BMisc::rhs.vars(pformla)
      bbb<-length(LesX)
      for (jj in 1:bbb){
        disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,]
      }
      
      ## finlament d�finitio ndu mod�le de score
      pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))
      base_score<-subset(disdat, C + G == 1)
      base_score=base_score[base_score$select==1,]
      pscore.reg <- glm(pformla, family = binomial(link = "logit"),
                        data = base_score)
      thet <- coef(pscore.reg)
      if (any(is.na(thet))) {
        warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
                       flist[f], " at time period: ", tlist[t +1]))
      }
      
      
      if (is.null(POND_RD)){
        disdat <- data[(data[, tname] == pret | data[, tname] == tlist[t]), ]
      } else{
        disdat <- data2[(data2[, tname] == pret | data2[, tname] == tlist[t]),]
      }
 
      # ou l'on prend l'annee de pretraitement  et l'annee observee
      disdat <- makeBalancedPanel(disdat, idname, tname)##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1])
      disdat <- panelDiffV_CS(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],selection)
      # write.xlsx(disdat,"disdat.xlsx")
      ### Contr�le pr�alable de l'existence de la variable de strate 
  
      #J'ai ajouté cette ligne car il y avait un problème avec la variable Dnom_outcome[i] qui n'était pas reconnue
      disdat=as.data.table(disdat)
      




      
      disdat<-disdat[is.na(disdat[,strate])==FALSE,]

      
      
      ### base de calcul
      pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
      pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
      
      ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
      disdat$pscore<-pscore
      disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
      disdat[,c("pscore")]=list(NULL)
      pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur
      
      
      Lnom_outcome<-paste0("L",nom_outcome)
      # ly <- disdat[,Lnom_outcome]
      ly <- disdat[,..Lnom_outcome]
      ##Indicatrice de Traitement et de contrefactuel
      # disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t])
      disdat$C <- 1 * (disdat[, ..first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t])
      # disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
      disdat$G <- 1 * (disdat[, ..first.treat.name] == flist[f])
      
      ## d�finition de N pour boot
      NN<-dim(disdat)[1]
      




      
      ### Calcul par strate 
      ### on ne conserve que les strates qui ont un poid non nul et un score non nul
      # disdat$pp<-disdat[,ponderation[1]]
      disdat$pp=disdat[[ponderation[1]]]
      for(i in 1:length(ponderation)){
        disdat$pp<-pmin(disdat[,ponderation[i]],disdat[,"pp"])
      }
        
      
      # strate_list<-unique(disdat[disdat[,first.treat.name]==flist[f] & disdat$pp>0,strate])
      # strate_list=unique(disdat[get(first.treat.name) == flist[f] & pp > 0, .(strate)])
      strate_list <- unique(disdat[get(first.treat.name) == flist[f] & pp > 0 & !is.na(strate), .(strate)])
      

      # strate_list<-strate_list[is.na(strate_list)==FALSE]
      strate_nb<-length(strate_list)
      ### Patch pour pr�venir les cas o� aucune strate n'est active
      if (strate_nb>0){

        strate_poids<-seq(1:strate_nb)
        for(riri in 1:strate_nb){


          #Probleme dans la facon de storer les données, cette ligne le règle.
          strate_list=strate_list$strate
        #             strate
        #     <num>
        # 1:      1
        # [1] 0 0 0 0 0 0
        #    strate
        #     <num>
        # 1:      0
        # [1] 1 1 1 1 1 1
     

          # lu<-1*(disdat[,strate]==strate_list[riri])*disdat$G
          lu<-1*(disdat[[strate]]==strate_list[riri])*disdat$G

            strate_poids[riri]<-sum(lu,na.rm=TRUE)
        }







        strate_poids<-strate_poids/sum(strate_poids,na.rm=TRUE)
        # Calcul de matrices pour calcul strate 
        traite<-cbind(disdat$G,disdat$C)
        for(riri in 1:strate_nb){
          if (riri==1){
            traiteS<- 1*(disdat[,strate]==strate_list[riri])*traite
            devant<-strate_poids[1]*(traiteS[,1]-traiteS[,2])
            devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2]
            traite_boot<-traiteS
            
          } else {
            ttt<- 1*(disdat[,strate]==strate_list[riri])*traite
            traiteS<- cbind(traiteS,ttt)
            ddd<-strate_poids[riri]*(ttt[,1]-ttt[,2])
            devant<-devant+ddd
            ddd2<-ttt[,1]+(pscore/(1 - pscore))*ttt[,2]
            devant2<-devant2+ddd2
            traite_boot<-cbind(traite_boot,traiteS)
            
          }}
        
        
        
        
        

        #Reformat devant2...
        devant2_df <- data.frame(values = as.numeric(devant2))
        row.names(devant2_df) <- names(devant2)
        devant2=devant2_df

        devant2=as.matrix(devant2)
        # dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,ponderation]))
        dommk<-t(as.matrix(traiteS))%*%as.matrix((as.matrix(devant2)*subset(disdat,select=ponderation)))
        
        Denom <- 1/(dommk)
        # pp<- devant*(traiteS%*%Denom)*(devant2*disdat[,ponderation])
        pp<- devant*(traiteS%*%Denom)*(devant2*subset(disdat,select=ponderation))
        # pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,ponderation])
        pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*subset(disdat,select=ponderation))


        
        att<- colSums(pp * ly)
        
        ### Bootstrap  : calcul de la fonction d'influence de la brique
        esperance_ly<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * ly))












        esperance_ly<- (-1) * esperance_ly * disdat$C + 1 * esperance_ly * disdat$G
        psi<-NN*pp *(ly-esperance_ly)
        colnames(psi)=nom_outcome
        x <- model.matrix(xformla, data = disdat)
        n <- nrow(disdat)
        for (i in 1:length(nom_outcome)){
          # G_i <- disdat$G*disdat[,ponderation[i]]
          G_i <- disdat$G*disdat[[ponderation[i]]]
          # C_i <- disdat$C*disdat[,ponderation[i]]
          C_i <- disdat$C*disdat[[ponderation[i]]]


          
          # wc1_i <- -1*C_i * pp_noDenom[,i]
          wc1_i <- -1*C_i * pp_noDenom[[i]]
          wc_i <- wc1_i/mean(wc1_i)
          # ly_i <- disdat[,Lnom_outcome[i]]
          ly_i <- disdat[[Lnom_outcome[i]]]
          
          
          
          
          
          
          
        
          
          M_i <- as.matrix(apply(as.matrix((C_i/(1 - pscore))^2 * gg(x, thet) * (ly_i - mean(wc_i * ly_i)) * x), 2,mean)/mean(wc1_i))
          A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
          A1_i <- (t(A1_i * x) %*% x/n)
          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x
          A_i <- A2_i %*% MASS::ginv(A1_i)
          correction_i<-A_i %*% M_i
          correction_i<-correction_i*disdat$C
          # psi[,nom_outcome[i]]<-psi[,nom_outcome[i]]-correction_i
          psi[,nom_outcome[i]]<-psi[[nom_outcome[i]]]-correction_i
        }
        
        ### mettre les siren dans psi
        psi<-as.data.frame(psi)
        colnames(psi)<-nom_outcome
        # psi[,idname]<-disdat[,idname]
        psi[,idname]<-disdat[[idname]]
        
        ### apparier psi avec indiv
        psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
        for (i in 1:length(nom_outcome)){
          psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
        }
        
        # Compter le nombre d'observations utilis�es pour l'estimation de chaque brique (moyenne sur les diff�rents outcomes �valu�s)
        # nG<-round(mean(colSums(disdat$G * 1 *(disdat[,ponderation]>0))))
        # nC<-round(mean(colSums(disdat$C * 1 *(disdat[,ponderation]>0))))
        nG<-round(mean(colSums(disdat$G * 1 *(disdat[,..ponderation]>0))))
        nC<-round(mean(colSums(disdat$C * 1 *(disdat[,..ponderation]>0))))

      } else  {
        
        att<- colSums(0 * ly)
        # Nombre d'observations utilis�es pour l'estimation de chaque brique
        nG<-mean(colSums(0 * ly))
        nC<-mean(colSums(0 * ly))
        psi<-as.data.frame(indiv[,idname])
        for (i in 1:length(nom_outcome)){
          psi[,nom_outcome[i]]<-0
        }
        
      }  
      ### le resultat de base est de type liste: le coefficient, le numero du traitement, la date, avant ou post
      fatt[counter,first.treat.name]<-flist[f]
      fatt[counter,tname]<-tlist[t]
      fatt[counter,nom_outcome]<-att
      fatt[counter,"nobsG"]<-nG
      fatt[counter,"nobsC"]<-nC
      ### Remplissage de la matrice des fonctions d'influence
      for (i in 1:dim(psi)[2]){
        mat_influence[,i,counter]=psi[,i]  
      }
    counter <- counter + 1
    
    } ### fin boucle sur l'annee
  
  } ## fin boucle sur le traitement f
  
  list(fatt,mat_influence,indiv)
}




##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR EN LONG DID ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------

compute.mp.spatt.Boot_longDID <- function (  nom_outcome
                                           , nom_traitement
                                           , flen  ## nombre de cohortes dans le traitement 
                                           , tlen  ## nombre d'annee
                                           , flist ## Le vecteur des vrais traitements (sans 0 !!)
                                           , tlist ## vecteur ordonn? des dates
                                           , data ### le panel balanc?
                                           , first.treat.name ## la variable indiquant l'annee de traitement par un dispositif donn�
                                           , xformla ## le mod?le de score
                                           , tname ## la variable de date
                                           , w ## ??  NULL par d?faut
                                           , idname ## identifiant individu
                                           , method ## par d?faut "logit"
                                           , seedvec ## par d?faut NULL
                                           , se ## TRUE par d?faut
                                           , pl ## FALSE par d?faut
                                           , cores ## par d?faut = 2
                                           , printdetails ## TRUE par defaut
                                           , selection ## indicatrice observations valides
                                           , ponderation ## poids utilise
                                           , strate  ## variable qui donne les classes d'observations, a utiliser 
)
{
  ### remplissage de l'objet en sortie : toujours utile ??
  nbligne=flen*tlen
  rempli<-seq(1:nbligne)
  fatt <- data.frame(compteur=rempli)
  fatt[,first.treat.name]<-rempli
  fatt[,tname]<-rempli
  fatt$nobsG<-rempli
  fatt$nobsC<-rempli 
  for(i in 1:length(nom_outcome)){
    fatt[,nom_outcome[i]]<-rempli
  }
  ### Liste des individus impliqu�s dans l'estimation pour le wild bootstrap
  indiv<-unique(data[,c(first.treat.name,idname)])
  ### Matrice d'influence initialis�e � 0
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(nom_outcome)),nbligne))
  
  counter <- 1
  
  for (f in 1:flen) {  ### boucle sur les cohortes dans le traitement : f
    
    for (t in 1:(tlen)) {  ## boucle sur la date : t
      ## calcul de la formule du score une fois pour toute
      pret <- flist[f]-1
      disdat <- data[(data[, tname] == pret), ]  ## la base de calcul
      
      disdat <- droplevels(disdat)  ##suppression d'observations manquantes ??
      
      ### Contr�le de l'existence de la variable de strate 
      disdat<-disdat[is.na(disdat[,strate])==FALSE,]
      
      # on calcule le score sur toutes les non encore traitees
      
      #disdat<-disdat[ (disdat[,first.treat.name]==0) | (disdat[,first.treat.name]==flist[f ]) | (disdat[,first.treat.name]>tlist[t]) ,]
      disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] >flist[f] & disdat[, first.treat.name] > tlist[t]) ### les contrefactuels
      disdat$G <- 1 * (disdat[, first.treat.name] == flist[f]) ## les traitees f
      if (is.null(xformla)) {
        xformla <- ~1
      }
      pformla <- xformla
      ### On supprime les observations pour lesquelles le mod�le de score ne peut �tre d�fini
      LesX<-BMisc::rhs.vars(pformla)
      
      bbb<-length(LesX)
      for (jj in 1:bbb){
        disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,]
      }
      
      ## finlament d�finitio ndu mod�le de score
      pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla))
      base_score<-subset(disdat, C + G == 1)
      base_score=base_score[base_score$select==1,]
      pscore.reg <- glm(pformla, family = binomial(link = "logit"),
                        data = base_score)
      
      thet <- coef(pscore.reg)
      
      if (any(is.na(thet))) {
        warning(paste0("Problems estimating propensity score...likely perfectly predicting treatment for group: ",
                       flist[f], " at time period: ", tlist[t +1]))
      }
      
      disdat <- data[(data[, tname] == pret | data[, tname] == tlist[t]), ]  ## la base de calcul
      # ou l'on prend l'annee de pretraitement  et l'annee observee
      disdat <- makeBalancedPanel(disdat, idname, tname)##on cylindre sur les trois dates utiles au calcul (pret,tlist[t],tlist[t+1])
      disdat <- panelDiffV_longDID(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],selection)
      
      #disdat<-disdat[ (disdat[,first.treat.name]==0) | (disdat[,first.treat.name]==flist[f ]) | (disdat[,first.treat.name]>tlist[t]) ,]
      
      ### Contr�le pr�alable de l'existence de la variable de strate 
      
      
      
      
      
      
      
      #J'ai ajouté cette ligne car il y avait un problème avec la variable Dnom_outcome[i] qui n'était pas reconnue
      disdat=as.data.table(disdat)
      disdat<-disdat[is.na(disdat[,strate])==FALSE,]
      
      ### base de calcul
      pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on calcule un score pour tout disdat
      pscore[is.na(pscore)]<-0   ## on met un score nul pour les obseravation pour lesquelles on ne peut pas calculer de score
      
      ### on supprime les scores = 0 | = 1 pour le bootstrap (et m�me pour l'estimateur)
      disdat$pscore<-pscore
      disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,]
      disdat[,c("pscore")]=list(NULL)
      pscore <- predict(pscore.reg, newdata = disdat,type = "response")  ### on recalcule un score pour la nouvelle base disdat pour que pscore ait la bonne longueur
      
      Dnom_outcome<-paste0("D",nom_outcome)
      # dy <- disdat[,Dnom_outcome]
      dy <- disdat[,..Dnom_outcome]
      ##Indicatrice de Traitement et de contrefactuel
      
      # disdat$C <- 1 * (disdat[, first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t])
      disdat$C <- 1 * (disdat[, ..first.treat.name] == 0) #+ 1 * (disdat[, first.treat.name] > flist[f] & disdat[, first.treat.name] > tlist[t])
      # disdat$G <- 1 * (disdat[, first.treat.name] == flist[f])
      disdat$G <- 1 * (disdat[, ..first.treat.name] == flist[f])
      
      ## d�finition de N pour boot
      NN<-dim(disdat)[1]
      
      ### Calcul par strate 
      ### on ne conserve que les strates qui ont un poid nul et un score non nul
      
      # disdat$pp<-disdat[,ponderation[1]]
      disdat$pp=disdat[[ponderation[1]]]
      for(i in 1:length(ponderation)){
        disdat$pp<-pmin(disdat[,ponderation[i]],disdat[,"pp"])
      }
      
      # strate_list<-unique(disdat[disdat[,first.treat.name]==flist[f] & disdat$pp>0,strate])
      # strate_list<-strate_list[is.na(strate_list)==FALSE]
      strate_list <- unique(disdat[get(first.treat.name) == flist[f] & pp > 0 & !is.na(strate), .(strate)])
      
      strate_nb<-length(strate_list)

      

      ### Patch pour pr�venir les cas o� aucune strate n'est active
      if (strate_nb>0){
        strate_poids<-seq(1:strate_nb)
        for(riri in 1:strate_nb){
          # lu<-1*(disdat[,strate]==strate_list[riri])*disdat$G
          # lu<-1*(disdat[[strate]]==strate_list[riri])*disdat$G
          lu <- 1 * as.vector(disdat[[strate]] == strate_list[riri]) * disdat$G

          

          
          
          
          strate_poids[riri]<-sum(lu,na.rm=TRUE)
        }
        
        strate_poids<-strate_poids/sum(strate_poids,na.rm=TRUE)
        # Calcul de matrices pour calcul strate 
        traite<-cbind(disdat$G,disdat$C)
        for(riri in 1:strate_nb){
          if (riri==1){
            
            
            
          
            traiteS<- as.numeric(1*(disdat[,strate])==as.numeric(strate_list[riri]))*traite
            
            
            







            
            devant<-strate_poids[1]*(traiteS[,1]-traiteS[,2])
            devant2<-traiteS[,1]+(pscore/(1 - pscore))*traiteS[,2]
            traite_boot<-traiteS
            
          } else {
            
            ttt<- 1*(disdat[,strate]==strate_list[riri])*traite
            
            traiteS<- cbind(traiteS,ttt)
            ddd<-strate_poids[riri]*(ttt[,1]-ttt[,2])
            devant<-devant+ddd
            ddd2<-ttt[,1]+(pscore/(1 - pscore))*ttt[,2]
            devant2<-devant2+ddd2
            traite_boot<-cbind(traite_boot,traiteS)
            
          }}
        
        # dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,ponderation]))
        dommk<-t(as.matrix(traiteS))%*%as.matrix((devant2*disdat[,..ponderation]))
        Denom <- 1/(dommk)
        pp<- devant*(traiteS%*%Denom)*(devant2*disdat[,..ponderation])
        pp_noDenom<- devant*(traiteS%*%(dommk*Denom))*(devant2*disdat[,..ponderation])
        att<- colSums(pp * dy)
        
        ### Bootstrap  : calcul de la fonction d'influence de la brique
        esperance_dy<- traite_boot%*% t(traite_boot) %*%as.matrix((pp * dy))
        esperance_dy<- (-1) * esperance_dy * disdat$C + 1 * esperance_dy * disdat$G
        
        psi<-NN*pp *(dy-esperance_dy)
        colnames(psi)=nom_outcome
        x <- model.matrix(xformla, data = disdat)
        n <- nrow(disdat)
      
      
        for (i in 1:length(nom_outcome)){
          
          G_i <- disdat$G*(subset(disdat,select=ponderation[i]))
          C_i <- disdat$C*(subset(disdat,select=ponderation[i]))
          
          G_i=as.vector(G_i)
          C_i=as.vector(C_i)
          G_i=unlist(G_i)
          C_i=unlist(C_i)

          wc1_i <- -1*C_i * pp_noDenom[[i]]
          wc_i <- wc1_i/mean(wc1_i)
          # wc_i <- as.numeric(wc_i)
          dy_i = subset(disdat,select=Dnom_outcome[i])
          dy_i=as.vector(dy_i)
          dy_i=unlist(dy_i)



    
    

          M_i <- as.matrix(apply(as.matrix((C_i/(1 - pscore))^2 * gg(x, thet) * (dy_i - mean(wc_i * dy_i)) * x), 2,mean)/mean(wc1_i))
          
          A1_i <- (G_i + C_i) * gg(x, thet)^2/(pscore * (1 -   pscore))
          A1_i <- (t(A1_i * x) %*% x/n)
          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg(x, thet)/(pscore *  (1 - pscore))) * x       
          A_i <- A2_i %*% MASS::ginv((A1_i))
          

        

  




          
          correction_i<-A_i %*% M_i
          correction_i<-correction_i*disdat$C


          # psi[,nom_outcome[i]]<-psi[,nom_outcome[i]]-correction_i
          psi[,nom_outcome[i]]<-subset(psi,select=nom_outcome[i])-correction_i
        }
        
        
        ### mettre les siren dans psi
        psi<-as.data.frame(psi)
        colnames(psi)<-nom_outcome
        
        # psi[,idname]<-disdat[,idname]
        psi[,idname]<-disdat[,..idname]
        
        ### apparier psi avec indiv
        psi<-merge(indiv[idname], psi, by.x=idname, by.y=idname, all.x=TRUE)
        for (i in 1:length(nom_outcome)){
          psi[is.na(psi[,nom_outcome[i]]),nom_outcome[i]]<-0
        }
        

        # Compter le nombre d'observations utilis�es pour l'estimation de chaque brique (moyenne sur les diff�rents outcomes �valu�s)
        # nG<-round(mean(colSums(disdat$G * 1 *(disdat[,ponderation]>0))))
        # nC<-round(mean(colSums(disdat$C * 1 *(disdat[,ponderation]>0))))
        nG<-round(mean(colSums(disdat$G * 1 *(disdat[,..ponderation]>0))))
        nC<-round(mean(colSums(disdat$C * 1 *(disdat[,..ponderation]>0))))
        
      } else  {
        
        att<- colSums(0 * dy)
        # Nombre d'observations utilis�es pour l'estimation de chaque brique
        nG<-mean(colSums(0 * dy))
        nC<-mean(colSums(0 * dy))
        psi<-as.data.frame(indiv[,idname])
        for (i in 1:length(nom_outcome)){
        psi[,nom_outcome[i]]<-0
        }  
      }  

      ### le resultat de base est de type liste: le coefficient, le numero du traitement, la date, avant ou post
      fatt[counter,first.treat.name]<-flist[f]
      fatt[counter,tname]<-tlist[t]
      #??? a verifier= ou bien :   fatt[counter,tname]<-tlist[(t + 1)] ???
      fatt[counter,nom_outcome]<-att
      fatt[counter,"nobsG"]<-nG
      fatt[counter,"nobsC"]<-nC
      ### Remplissage de la matrice des fonctions d'influence
      for (i in 1:dim(psi)[2]){
        mat_influence[,i,counter]=psi[,i]  
      }
      counter <- counter + 1
      
    } ### fin boucle sur l'annee
  } ## fin boucle sur le traitement f

  list(fatt,mat_influence,indiv)
}


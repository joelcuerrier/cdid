#Fonctions d'estimation des effets de traitement

# ###### ----------------------------------------------------------------------------------------------------------------------------------------
# ###### Basic EstimPeriod_Boot ------------------------------------------------------------------------------------------------------------------
# ###### ----------------------------------------------------------------------------------------------------------------------------------------
# estimPeriod_Boot<-function(yname,
#                            tname,
#                            idname,
#                            gname,
#                            xformla,
#                            data,
#                            debT,
#                            finT,
#                            deb,
#                            fin,
#                            select,
#                            weightsname,
#                            alp=0.05,
#                            bstrap,
#                            biters=1000,
#                            STRATE,
#                            pond_RD=NULL,
#                            treated){

#   ### Créer un siren numérique (id)
#   list_id <-as.data.frame(unique(data[,c(idname)]))
#   list_id$iden_num<-1:dim(list_id)[1] 
#   colnames(list_id)[1]=idname
#   data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)
  
#   ### selection des données 
#   bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
#   bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]
  
#   ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
#   ## puis on bascule entièrement dans le  contrefactuel les traités après finT
#   for(i in 1:length(weightsname)){
#   bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
#   bebe[bebe[,gname]>finT,gname]<-0
  
#   ### definition du modèele utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
    
#   # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
#   resultat<- mp.spatt.Boot(nom_outcome=yname,nom_traitement=treated,xformla=xxF,data=bebe,
#                          first.treat.name = gname,
#                          idname = idname, tname=tname,
#                          bstrap = FALSE,se=TRUE,cband =FALSE
#                          ,selection=select,ponderation=weightsname,debT2=debT,finT2=finT,strate=STRATE,POND_RD=pond_RD)
  







  
#   # Poids pour aggréger les effets des différentes cohortes de traitement 
#   list_id<-merge(list_id,resultat[[3]])
#   list_id_poids=resultat[[3]]
#   list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
#   list_id_poids<-unique(list_id_poids)
#   dum<-as.data.frame(table(list_id_poids[,gname]))
  
#   result<-agregatChris(tab=resultat[[1]],nom_outcome=yname,tname=tname,first.treat.name=gname,dum)
#   influ<-agregat_influence(tab=resultat[[1]],array_inf=resultat[[2]],listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   list(result,influ,list_id)
# }


###### ----------------------------------------------------------------------------------------------------------------------------------------
###### ESTIMATEUR EN CHAINES ------------------------------------------------------------------------------------------------------------------
###### ----------------------------------------------------------------------------------------------------------------------------------------
chained_estimPeriod_Boot<-function(yname,
                                   tname,
                                   idname,
                                   gname,
                                   xformla,
                                   propensityformla,
                                   data,
                                   debT,
                                   finT,
                                   deb,
                                   fin,
                                   select,
                                   weightsname, #St
                                   weight_assumption,
                                   
                                   cband=cband,
                                   alp=0.05,
                                   bstrap,
                                   biters=1000,
                                   treated){
  set.seed(123)
  
  ### Créer un siren numérique (id)
  list_id <-as.data.frame(unique(data[,c(idname)]))
  list_id$iden_num<-1:dim(list_id)[1] 
  colnames(list_id)[1]=idname
  data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)
  
  ### selection des données 
  bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
  bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]
  
  ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
  ## puis on bascule entièrement dans le  contrefactuel les traités après finT
  for(i in 1:length(weightsname)){
  bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
  bebe[bebe[,gname]>finT,gname]<-0
  
  ### definition du modèele utilisant le score 
  xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
  # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
  resultat<- chained.mp.spatt.Boot(nom_outcome=yname,nom_traitement=treated,xformla=xxF,propensityformla=propensityformla,data=bebe,
                         first.treat.name = gname,
                         idname = idname, tname=tname,
                         bstrap = FALSE,se=TRUE,cband =FALSE
                         ,selection=select,ponderation=weightsname,weight_assumption=weight_assumption,debT=debT,finT=finT)
  View(resultat[[1]])
  # # Ajouté le 4 décembre pour reproduire l'aggrégation inititiale réalisée par Christophe.
  # # Poids pour aggr�ger les effets des diff�rentes cohortes de traitement 
  list_id_poids=resultat[[3]]
  list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
  list_id_poids<-unique(list_id_poids)
  dum<-as.data.frame(table(list_id_poids[,gname]))

  # #intialement
  # result<-agregatChris(tab=resultat[[1]],nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum)
  result<-agregatChris(tab=resultat[[1]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
# influ<-agregat_influence(tab=resultat[[1]],array_inf=resultat[[2]],listG=resultat[[3]],nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,poids=dum)
  influ<-agregat_influence(tab=resultat[[1]],array_inf=resultat[[2]],listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
  list(result,influ,list_id)
}


wild_bootstrap<-function(result,influ,list_id,nom_outcome,biters,nivtest,seedvec=NULL){
  NN=dim(list_id)[1]
  borne_inf<-result
  borne_sup<-result
  mat_alea = matrix(sample(c(-1, 1),NN*biters , replace = T),ncol=biters)
  mat_boot = array(0,dim=c(biters,length(nom_outcome),dim(result)[1]))
  mat_boot_diff = mat_boot
  sdboot = result
  quant_005 = result
  quant_010 = result
  quant_025 = result
  quant_050 = result
  quant_100 = result
  quant_900 = result
  quant_950 = result
  quant_975 = result
  quant_990 = result
  quant_995 = result
  waldtest = result
  sdalt = result
  sdalt2 = result
  SigmaBoot = result
  ccc = result
    
  for (i in 1:length(nom_outcome)){
    bout <- lapply(1:biters, FUN = function(b) {
      Ub <- mat_alea[,b]
      Rb <- sqrt(NN) * (apply(Ub * (influ[,(i+1),]), 2, mean))
      Rb
    })
    for (b in 1:biters){
      Ub <- mat_alea[,b]
      Rb <- apply(Ub * (influ[,(i+1),]), 2, mean)
      mat_boot[b,i,]<-as.vector(Rb) + as.vector(result[,i+1])
      mat_boot_diff[b,i,]<-as.vector(abs( mat_boot[b,i,])) - as.vector(abs(result[,i+1]))
    }  
    bres <- t(simplify2array(bout))
    V <- cov(bres)
    cval <- qnorm(1 - nivtest/2)
    bSigma <- apply(bres, 2, function(b) 
               (quantile(b, 0.75, type = 1) - quantile(b, 0.25, type = 1))/(qnorm(0.75) - qnorm(0.25)))
    
    # Erreur possible ???
    #bSigma_moins_un_demi<-bSigma^(-0.5)
    bSigma_moins_un_demi<-bSigma^(-1)
    bSigma_moins_un_demi[is.infinite(bSigma_moins_un_demi)]<-0
    bT <- apply(bres, 1, function(b) max(abs(b) * bSigma_moins_un_demi))
    cval <- quantile(bT, 1 - nivtest, type = 1)
    # V <- diag(bSigma)
    #borne_inf[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #borne_sup[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    borne_inf[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    borne_sup[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    sdboot[,nom_outcome[i]]<-apply(mat_boot[,i,],2,sd)
    CoVar_boot=cov(mat_boot[,i,])
    #Autre méthode pour calculter la matrice de covariance des estimateurs (pas robuste à l'hétéroscédasticité)
    #CoVar_boot=cov(influ[,i+1,])/dim(influ)[1]
    # quant_005[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.005),na.rm=TRUE)})
    # quant_010[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.010),na.rm=TRUE)})
    # quant_025[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.025),na.rm=TRUE)})
    # quant_050[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.050),na.rm=TRUE)})
    # quant_100[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.100),na.rm=TRUE)})
    # quant_900[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.900),na.rm=TRUE)})
    # quant_950[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.950),na.rm=TRUE)})
    # quant_975[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.975),na.rm=TRUE)})
    # quant_990[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.990),na.rm=TRUE)})
    # quant_995[,nom_outcome[i]]<-apply(mat_boot[,i,],2, FUN = function(vec) {quantile(vec,probs=c(0.995),na.rm=TRUE)})
    ## Correction de la méthode qui est apparement fausse ;: adaptation de la méthode de l'intervalle de confiance (VDB - 19 aout)
    cval <- quantile(bT,0.99, type = 1)
    #quant_005[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #quant_995[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    quant_005[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    quant_995[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    cval <- quantile(bT,0.98, type = 1)
    #quant_010[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #quant_990[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    quant_010[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    quant_990[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    cval <- quantile(bT,0.95, type = 1)
    #quant_025[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #quant_975[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    quant_025[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    quant_975[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    
    ## Standard Deviation : calcul alternatif
    sdalt[,nom_outcome[i]]<- (quant_975[,nom_outcome[i]] - quant_025[,nom_outcome[i]]) / 3.92 
    sdalt2[,nom_outcome[i]]<- (quant_975[,nom_outcome[i]] - quant_025[,nom_outcome[i]]) / (2*cval) 
    SigmaBoot[,nom_outcome[i]] <- bSigma
    ccc[,nom_outcome[i]] <- cval / sqrt(NN) 
    
    cval <- quantile(bT,0.9, type = 1)
    #quant_050[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #quant_950[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    quant_050[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    quant_950[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    cval <- quantile(bT,0.8, type = 1)
    #quant_100[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma^(0.5)
    #quant_900[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma^(0.5)
    quant_100[,nom_outcome[i]]<-result[,nom_outcome[i]]-(1/sqrt(NN))*cval*bSigma
    quant_900[,nom_outcome[i]]<-result[,nom_outcome[i]]+(1/sqrt(NN))*cval*bSigma
    ## Test de Wald pour la pré-trend
    V<-CoVar_boot
    V=cbind(result[,1],V)
    V=rbind(matrix(-1,1,dim(V)[2]),V)
    V[1,(2:dim(V)[2])]<-result[,1]
    result_befor<-result[result[,1]<=(-1),]
    result_after<-result[result[,1]>0,]
    V_befor<-V[V[,1]<=(-1),V[1,]<=(-1)]
    V_after<-V[V[,1]>0,V[1,]>0]
    waldtest_befor<-waldtest[waldtest[,1]<=(-1),]
    waldtest_after<-waldtest[waldtest[,1]>0,]
    for (j in 1:(max(result[,1]))){
      theta=result_after[(1:j),(i+1)]
      R=diag(j)
      r=matrix(0,j,1)
      CoVar=V_after[(1:j),(1:j)]
      stat= t(R %*% theta - r) %*% solve(R %*% (CoVar/length(r)) %*% t(R)) %*% (R %*% theta - r)
      pval=1-pchisq(stat,length(r))
      waldtest_after[j,(i+1)]=pval
    }
    for (j in 1:(abs(min(result[,1])))){
      k= abs(min(result[,1])) - (j-1)
      theta=result_befor[(k:(abs(min(result[,1])))),(i+1)]
      R=diag(j)
      r=matrix(0,j,1)
      CoVar=V_befor[(k:(abs(min(result[,1])))),(k:(abs(min(result[,1]))))]
      stat= t(R %*% theta - r) %*% solve(R %*% (CoVar/length(r)) %*% t(R)) %*% (R %*% theta - r)
      pval=1-pchisq(stat,length(r))
      waldtest_befor[k,(i+1)]=pval
    }
    waldtest_befor[abs(min(result[,1]))+1,]<-0
    waldtest=rbind(waldtest_befor,waldtest_after)
    waldtest[,1]=result[,1]
  }
  
  list(result,borne_inf,borne_sup,sdboot,quant_005,quant_010,quant_025,quant_050,quant_100,quant_900,quant_950,quant_975,quant_990,quant_995,waldtest,sdalt,sdalt2,SigmaBoot,ccc)
}

# ## Fonctions pour calcul bootstrap 
# precizionS<-function(NOM,CHEMI,baz,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,Les_outcome,
#                      Pond_outcome,NN,strate){
  
#   sisiz<-unique(baz[,c("siren",anneeT)])
#   Martu<-function(samp){
#     sampi<-subset(samp,select=siren)
#     nn<-dim(sampi)[1]
#     sampi$identif<-c(1:nn)
#     PAPA<-merge(baz,sampi,by="siren")
#     lele<-estimPeriod6(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident="identif",Les_outcome,Pond_outcome,STRATE=strate)
#     Les_outcome=c(Les_outcome,"nobsG","nobsC")
#     resu<-lele[,Les_outcome[1]]
#     nnn<-length(Les_outcome)
#     for (i in 2:nnn){
#       resu<-c(resu,lele[,Les_outcome[i]])
#     }
#     resu
#   }
#   estu<-function(tab,ind){
#     print(compteur) ; compteur<-compteur+1
#     Martu(tab[ind,])
#   }
#   compteur<-1
#   vv<-boot(sisiz,estu,NN,sim="ordinary",stype="i",sisiz[,anneeT])
#   nom_tt<-paste0(NOM,"_tt")
#   sauve(vv$t,nom_tt,CHEMI)
#   nom_t0<-paste0(NOM,"_t0")
#   sauve(vv$t0,nom_t0,CHEMI)
#   vv}


# #### Pour les donnnées pondérées de l'enquête R&D (un argument supplémentaire)
# estimPeriodRD<-function(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident,Les_outcome
#                         ,Pond_outcome,STRATE,popo){
#   ### selection des données 
#   bebe<-PAPA[PAPA$annee>=deb & PAPA$annee <= fin ,]
#   bebe<-bebe[bebe[,anneeT]>=debT | bebe[,anneeT]==0,]### cette condition reste n�cessaire
#   ## Pour les entreprises traitées aprés finT on met une pond�ration nulle pour les observation apr�s finT
#   for(i in 1:length(Pond_outcome)){
#     bebe[(bebe[,anneeT]>finT)&(bebe[,"annee"]>=finT),Pond_outcome[i]]<-0
#   }
#   ## puis on bascule enti�rement dans le  contrefactuel les trait�s apr�s finT
#   bebe[bebe[,anneeT]>finT,anneeT]<-0
#   #bebe<-bebe[bebe[,anneeT]<=finT | bebe[,anneeT]==0,]### cette conditions supprime du contrefactuel les traitements suivants !! 
#   ### mais indispensable pour limiter la liste des traitements pris en compte : du coup on arr�te des faire de estimations par sous-p�riode
#   ### definition du mod�le utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(var_score, collapse=" + "),sep=""))
#   baba<-bebe
#   resultat<- mp.spatt.RD(nom_outcome=Les_outcome,nom_taitement=TT,xformla=xxF,data=baba,
#                          first.treat.name = anneeT,
#                          idname = nom_ident, tname="annee",
#                          bstrap = FALSE,se=TRUE,cband =FALSE
#                          ,selection=nom_select,ponderation=Pond_outcome,
#                          debT2=debT,finT2=finT,strate=STRATE,PondRD=popo)
#   # resultat
#   dum<-as.data.frame(table(bebe[bebe[,anneeT]>=debT & bebe[,anneeT]<=finT,anneeT]))
#   # agregat(tab=resultat,nom_outcome=outC,tname="annee",first.treat.name=anneeT) 
#   agregatChris(tab=resultat,nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum) 
# }

# ## Fonctions pour calcul bootstrap 
# precizionRD<-function(NOM,CHEMI,baz,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,
#                       Les_outcome,Pond_outcome,NN,strate,PP){
#   sisiz<-unique(baz[,c("siren",anneeT)])
#   Martu<-function(samp){
#     sampi<-subset(samp,select=siren)
#     nn<-dim(sampi)[1]
#     sampi$identif<-c(1:nn)
#     PAPA<-merge(baz,sampi,by="siren")
#     lele<-estimPeriodRD(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident="identif",Les_outcome,Pond_outcome
#                         ,STRATE=strate,popo=PP)
#     Les_outcome=c(Les_outcome,"nobsG","nobsC")
#     resu<-lele[,Les_outcome[1]]
#     nnn<-length(Les_outcome)
#     for (i in 2:nnn){
#       resu<-c(resu,lele[,Les_outcome[i]])
#     }
#     resu
#   }
#   estu<-function(tab,ind){
#     print(compteur) ; compteur<-compteur+1
#     Martu(tab[ind,])
#   }
#   compteur<-1
#   vv<-boot(sisiz,estu,NN,sim="ordinary",stype="i",sisiz[,anneeT])
#   nom_tt<-paste0(NOM,"_tt")
#   sauve(vv$t,nom_tt,CHEMI)
#   nom_t0<-paste0(NOM,"_t0")
#   sauve(vv$t0,nom_t0,CHEMI)
#   vv}


###### ----------------------------------------------------------------------------------------------------------------------------------------
###### ESTIMATEUR GMM -------------------------------------------------------------------------------------------------------------------------
###### Version du 5 décembre 2024
###### ----------------------------------------------------------------------------------------------------------------------------------------

GMM_estimPeriod_Boot<-function(yname,
                                 tname,
                                 idname,
                                 gname,
                                 xformla,
                                 propensityformla,
                                 data,
                                 debT,
                                 finT,
                                 deb,
                                 fin,
                                 select,
                                 weightsname, #ST
                                 weight_assumption,
                                 cband=cband,
                                 alp=0.05,
                                 bstrap=FALSE,
                                 biters=1000,
                                 
                                 
                                 treated){
 
  set.seed(123)
  
  ### Créer un siren numérique (id)
  list_id <-as.data.frame(unique(data[,c(idname)]))
  list_id$iden_num<-1:dim(list_id)[1]
  colnames(list_id)[1]=idname
  data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)

  ### selection des données 
  bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
  bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]

  ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
  ## puis on bascule entièrement dans le  contrefactuel les traités après finT
  for(i in 1:length(weightsname)){
  bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
  bebe[bebe[,gname]>finT,gname]<-0
  
  ### definition du modèle utilisant le score 
  xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
  # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
  resultat<- mp.spatt.GMM(nom_outcome=yname,nom_traitement=treated,xformla=xxF,propensityformla=propensityformla,data=bebe,
                           first.treat.name = gname,
                           idname = idname, tname=tname,
                           bstrap = FALSE,se=TRUE,cband =FALSE
                           ,selection=select,ponderation=weightsname,weight_assumption=weight_assumption,debT=debT,finT=finT)
  
  
  View(resultat[[1]])
  
  # Poids pour aggréger les effets des différentes cohortes de traitement 
  list_id<-merge(list_id,resultat[[3]]) #indiv
  list_id_poids=resultat[[3]] #id, gname
  list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
  list_id_poids<-unique(list_id_poids)
  dum<-as.data.frame(table(list_id_poids[,gname]))

  flist<-resultat[[5]]
  tlist<-resultat[[6]]
  flen <- length(flist)
  tlen <- length(tlist)

  # Covariance matrix of delta ATT estimator 
  #VERIFIER QUE J'AI LE BON NOMBRE D'OBSERVATIONS POUR NN ??
  omega_deltaATT <- array(0,dim=c(length(yname),dim(resultat[[1]])[1],dim(resultat[[1]])[1]))
  for(i in 1:length(yname)){
    omega_deltaATT[i,,] <- (1/dim(resultat[[2]])[1]) * t(resultat[[2]][,i+1,]) %*% resultat[[2]][,i+1,]   
  }
  # Remove columns from mat_W 
  mat_W <- resultat[[4]]
  remov_col=c()
  for (f in 1:flen) { 
    for (t in 1:tlen) {
      if (flist[f]-1 == tlist[t]){remov_col<-c(remov_col,(f-1)*tlen+t)}
    }
  }
  mat_W <- subset( mat_W, select = -c(remov_col ) )
  mat_W <- as.matrix(mat_W)
  # Covariance Matrix of ATTgt
  Sigma_ATTgt <- array(0,dim=c(length(yname),dim(mat_W)[2],dim(mat_W)[2]))
  for(i in 1:length(yname)){
    Sigma_ATTgt[i,,] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)
  }
  # Optimal estimator of ATTgt
  delta_ATT <- as.matrix(resultat[[1]][,yname])
  ATTgt <- matrix(NA,ncol(mat_W),length(yname))
  colnames(ATTgt) <- yname
  rownames(ATTgt) <- colnames(mat_W)
  for(i in 1:length(yname)){
    #donc laTT est défini ici et non dans compute. Donc att dépend de l'inffluence.
    ATTgt[,i] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%delta_ATT[,i]
  }
  #inf pour inference
  # Mise en forme de la matrice ATTgt pour pouvoir l'agreger a des effets dynamiques
  ATTgt<-as.data.frame(ATTgt)
  ncol_ATTgt<-ncol(ATTgt)
  for (f in 1:flen) { 
    for (t in 1:(tlen-1)) {
      
      ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+1]<-flist[f]
      if(flist[f]-1>tlist[t]){ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]}else{ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]+1}
    }
  }
  colnames(ATTgt)[ncol_ATTgt+1]<-gname
  colnames(ATTgt)[ncol_ATTgt+2]<-tname
  # Influence function of ATTgt : PHI
  PHI_influence_ATTgt=array(0,dim=c(dim(resultat[[2]])[1],length(yname),dim(ATTgt)[1])) #6000,1,attgt
  for(i in 1:length(yname)){
    #Equation 50 fonction influence ATTgt, pas delta attgt.
    PHI_influence_ATTgt[,i,] <- t(MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%t(resultat[[2]][,i+1,]))
  }
  
  # return(list(ATTgt,PHI_influence_ATTgt,resultat,list_id,dum))
  
  # Agregation des effets dynamiques à partir de ATTgt
  result<-agregatChris_GMM(tab=ATTgt,nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
  # Calculer la fonction d'influence des effets agreges dynamiques
  influ<-agregat_influence_GMM(tab=ATTgt,array_influ=resultat[[2]],agreg_influence=PHI_influence_ATTgt,listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
  #influ calcule les effets dynamiques (?) #sum pondéré donc moins de dim
  # Liste des sorties utiles
  # list(result,influ,list_id)
  # View(as.data.frame(influ))
  # print(rownames(influ))
  


  # list(resultat,PHI_influence_ATTgt)
  list(result,influ,list_id)
}
###### ----------------------------------------------------------------------------------------------------------------------------------------
###### ESTIMATEUR GMM -------------------------------------------------------------------------------------------------------------------------
###### Version avant le 5 décembre 2024
###### ----------------------------------------------------------------------------------------------------------------------------------------

#   GMM_estimPeriod_Boot<-function(yname,
#                                  tname,
#                                  idname,
#                                  gname,
#                                  xformla,
#                                  propensityformla,
#                                  data,
#                                  debT,
#                                  finT,
#                                  deb,
#                                  fin,
#                                  select,
#                                  weightsname, #ST
#                                  weight_assumption,
#                                  cband=cband,
#                                  alp=0.05,
#                                  bstrap=FALSE,
#                                  biters=1000,
                                 
                                 
#                                  treated){
 
#   ### Créer un siren numérique (id)
#   list_id <-as.data.frame(unique(data[,c(idname)]))
#   list_id$iden_num<-1:dim(list_id)[1]
#   colnames(list_id)[1]=idname
#   data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)

#   ### selection des données 
#   bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
#   bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]

#   ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
#   ## puis on bascule entièrement dans le  contrefactuel les traités après finT
#   for(i in 1:length(weightsname)){
#   bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
#   bebe[bebe[,gname]>finT,gname]<-0
  
#   ### definition du modèle utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
#   # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
#   resultat<- mp.spatt.GMM(nom_outcome=yname,nom_traitement=treated,xformla=xxF,propensityformla=propensityformla,data=bebe,
#                            first.treat.name = gname,
#                            idname = idname, tname=tname,
#                            bstrap = FALSE,se=TRUE,cband =FALSE
#                            ,selection=select,ponderation=weightsname,weight_assumption=weight_assumption,debT=debT,finT=finT)
  
  
#   # str(resultat)
#   # Poids pour aggréger les effets des différentes cohortes de traitement 
#   list_id<-merge(list_id,resultat[[3]]) #indiv
#   list_id_poids=resultat[[3]] #id, gname
#   list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
#   list_id_poids<-unique(list_id_poids)
#   dum<-as.data.frame(table(list_id_poids[,gname]))

#   flist<-resultat[[5]]
#   tlist<-resultat[[6]]
#   flen <- length(flist)
#   tlen <- length(tlist)

#   # Covariance matrix of delta ATT estimator 
#   #VERIFIER QUE J'AI LE BON NOMBRE D'OBSERVATIONS POUR NN ??
#   omega_deltaATT <- array(0,dim=c(length(yname),dim(resultat[[1]])[1],dim(resultat[[1]])[1]))
#   for(i in 1:length(yname)){
#     omega_deltaATT[i,,] <- (1/dim(resultat[[2]])[1]) * t(resultat[[2]][,i+1,]) %*% resultat[[2]][,i+1,]   
#   }
#   # Remove columns from mat_W 
#   mat_W <- resultat[[4]]
#   remov_col=c()
#   for (f in 1:flen) { 
#     for (t in 1:tlen) {
#       if (flist[f]-1 == tlist[t]){remov_col<-c(remov_col,(f-1)*tlen+t)}
#     }
#   }
#   mat_W <- subset( mat_W, select = -c(remov_col ) )
#   mat_W <- as.matrix(mat_W)
#   # Covariance Matrix of ATTgt
#   Sigma_ATTgt <- array(0,dim=c(length(yname),dim(mat_W)[2],dim(mat_W)[2]))
#   for(i in 1:length(yname)){
#     Sigma_ATTgt[i,,] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)
#   }
#   # Optimal estimator of ATTgt
#   delta_ATT <- as.matrix(resultat[[1]][,yname])
#   ATTgt <- matrix(NA,ncol(mat_W),length(yname))
#   colnames(ATTgt) <- yname
#   rownames(ATTgt) <- colnames(mat_W)
#   for(i in 1:length(yname)){
#     #donc laTT est défini ici et non dans compute. Donc att dépend de l'inffluence.
#     ATTgt[,i] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%delta_ATT[,i]
#   }
#   #inf pour inference
#   # Mise en forme de la matrice ATTgt pour pouvoir l'agreger a des effets dynamiques
#   ATTgt<-as.data.frame(ATTgt)
#   ncol_ATTgt<-ncol(ATTgt)
#   for (f in 1:flen) { 
#     for (t in 1:(tlen-1)) {
      
#       ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+1]<-flist[f]
#       if(flist[f]-1>tlist[t]){ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]}else{ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]+1}
#     }
#   }
#   colnames(ATTgt)[ncol_ATTgt+1]<-gname
#   colnames(ATTgt)[ncol_ATTgt+2]<-tname
#   # Influence function of ATTgt : PHI
#   PHI_influence_ATTgt=array(0,dim=c(dim(resultat[[2]])[1],length(yname),dim(ATTgt)[1])) #6000,1,attgt
#   for(i in 1:length(yname)){
#     #Equation 50 fonction influence ATTgt, pas delta attgt.
#     PHI_influence_ATTgt[,i,] <- t(MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%t(resultat[[2]][,i+1,]))
#   }
    
  
#   # Agregation des effets dynamiques à partir de ATTgt
#   # result<-agregatChris_GMM(tab=ATTgt,nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   # Calculer la fonction d'influence des effets agreges dynamiques

#   # influ<-agregat_influence_GMM(tab=ATTgt,array_influ=resultat[[2]],agreg_influence=PHI_influence_ATTgt,listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   #influ calcule les effets dynamiques (?) #sum pondéré donc moins de dim
#   # Liste des sorties utiles
#   # list(result,influ,list_id)
#   # View(as.data.frame(influ))
#   # print(rownames(influ))
  


#   list(resultat,PHI_influence_ATTgt)
# }

# ###### ----------------------------------------------------------------------------------------------------------------------------------------
# ###### back up (sans modifications) ESTIMATEUR GMM -------------------------------------------------------------------------------------------------------------------------
# ###### ----------------------------------------------------------------------------------------------------------------------------------------

#   GMM_estimPeriod_Boot<-function(yname,
#                                  tname,
#                                  idname,
#                                  gname,
#                                  xformla,
#                                  propensityformla
#                                  data,
#                                  debT,
#                                  finT,
#                                  deb,
#                                  fin,
#                                  select,
#                                  weightsname, #ST
#                                  weight_assumption,
#                                  cband=cband,
#                                  alp=0.05,
#                                  bstrap=FALSE,
#                                  biters=1000,
                                 
                                 
#                                  treated){
 
#   ### Créer un siren numérique (id)
#   list_id <-as.data.frame(unique(data[,c(idname)]))
#   list_id$iden_num<-1:dim(list_id)[1]
#   colnames(list_id)[1]=idname
#   data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)

#   ### selection des données 
#   bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
#   bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]

#   ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
#   ## puis on bascule entièrement dans le  contrefactuel les traités après finT
#   for(i in 1:length(weightsname)){
#   bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
#   bebe[bebe[,gname]>finT,gname]<-0
  
#   ### definition du modèle utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
#   # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
#   resultat<- mp.spatt.GMM(nom_outcome=yname,nom_traitement=treated,xformla=xxF,data=bebe,
#                            first.treat.name = gname,
#                            idname = idname, tname=tname,
#                            bstrap = FALSE,se=TRUE,cband =FALSE
#                            ,selection=select,ponderation=weightsname,debT2=debT,finT2=finT,strate=STRATE,POND_RD=pond_RD)
  
  
#   # Poids pour aggréger les effets des différentes cohortes de traitement 
#   list_id<-merge(list_id,resultat[[3]])
#   list_id_poids=resultat[[3]]
#   list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
#   list_id_poids<-unique(list_id_poids)
#   dum<-as.data.frame(table(list_id_poids[,gname]))

#   flist<-resultat[[5]]
#   tlist<-resultat[[6]]
#   flen <- length(flist)
#   tlen <- length(tlist)

#   # Covariance matrix of delta ATT estimator 
#   #VERIFIER QUE J'AI LE BON NOMBRE D'OBSERVATIONS POUR NN ??
#   omega_deltaATT <- array(0,dim=c(length(yname),dim(resultat[[1]])[1],dim(resultat[[1]])[1]))
#   for(i in 1:length(yname)){
#     omega_deltaATT[i,,] <- (1/dim(resultat[[2]])[1]) * t(resultat[[2]][,i+1,]) %*% resultat[[2]][,i+1,]   
#   }
#   # Remove columns from mat_W 
#   mat_W <- resultat[[4]]
#   remov_col=c()
#   for (f in 1:flen) { 
#     for (t in 1:tlen) {
#       if (flist[f]-1 == tlist[t]){remov_col<-c(remov_col,(f-1)*tlen+t)}
#     }
#   }
#   mat_W <- subset( mat_W, select = -c(remov_col ) )
#   mat_W <- as.matrix(mat_W)
#   # Covariance Matrix of ATTgt
#   Sigma_ATTgt <- array(0,dim=c(length(yname),dim(mat_W)[2],dim(mat_W)[2]))
#   for(i in 1:length(yname)){
#     Sigma_ATTgt[i,,] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)
#   }
#   # Optimal estimator of ATTgt
#   delta_ATT <- as.matrix(resultat[[1]][,yname])
#   ATTgt <- matrix(NA,ncol(mat_W),length(yname))
#   colnames(ATTgt) <- yname
#   rownames(ATTgt) <- colnames(mat_W)
#   for(i in 1:length(yname)){
#     ATTgt[,i] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%delta_ATT[,i]
#   }
#   # Mise en forme de la matrice ATTgt pour pouvoir l'agreger a des effets dynamiques
#   ATTgt<-as.data.frame(ATTgt)
#   ncol_ATTgt<-ncol(ATTgt)
#   for (f in 1:flen) { 
#     for (t in 1:(tlen-1)) {
      
#       ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+1]<-flist[f]
#       if(flist[f]-1>tlist[t]){ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]}else{ATTgt[(f-1)*(tlen-1)+t,ncol_ATTgt+2]<-tlist[t]+1}
#     }
#   }
#   colnames(ATTgt)[ncol_ATTgt+1]<-gname
#   colnames(ATTgt)[ncol_ATTgt+2]<-tname
#   # Influence function of ATTgt : PHI
#   PHI_influence_ATTgt=array(0,dim=c(dim(resultat[[2]])[1],length(yname),dim(ATTgt)[1]))
#   for(i in 1:length(yname)){
#     PHI_influence_ATTgt[,i,] <- t(MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[i,,])%*%t(resultat[[2]][,i+1,]))
#   }
    
  
#   # Agregation des effets dynamiques à partir de ATTgt
#   result<-agregatChris_GMM(tab=ATTgt,nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   # Calculer la fonction d'influence des effets agreges dynamiques
#   influ<-agregat_influence_GMM(tab=ATTgt,array_influ=resultat[[2]],agreg_influence=PHI_influence_ATTgt,listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   # Liste des sorties utiles
#   list(result,influ,list_id)
# }



###### ----------------------------------------------------------------------------------------------------------------------------------------
###### ESTIMATEUR EN CROSS-SECTION ------------------------------------------------------------------------------------------------------------
###### ----------------------------------------------------------------------------------------------------------------------------------------

CS_estimPeriod_Boot<-function(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp=0.05,bstrap=FALSE,biters=1000,STRATE,pond_RD=NULL,treated){

  ### Créer un siren numérique (id)
  list_id <-as.data.frame(unique(data[,c(idname)]))
  list_id$iden_num<-1:dim(list_id)[1]
  colnames(list_id)[1]=idname
  data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)

  ### selection des données 
  bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
  bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]
  
  ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
  ## puis on bascule entièrement dans le  contrefactuel les traités après finT
  for(i in 1:length(weightsname)){
  bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
  bebe[bebe[,gname]>finT,gname]<-0
  
  ### definition du mod'ele utilisant le score 
  xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
  # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
  resultat<- mp.spatt.CS_Boot(nom_outcome=yname,nom_traitement=treated,xformla=xxF,data=bebe,
                                   first.treat.name = gname,idname = idname, tname=tname,
                                   bstrap = FALSE,se=TRUE,cband =FALSE,selection=select,
                                   ponderation=weightsname,debT2=debT,finT2=finT,strate=STRATE,POND_RD=pond_RD)
  
  # Poids pour aggréger les effets des différentes cohortes de traitement 
  list_id<-merge(list_id,resultat[[3]])
  list_id_poids=resultat[[3]]
  list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
  list_id_poids<-unique(list_id_poids)
  dum<-as.data.frame(table(list_id_poids[,gname]))

  result<-agregatChris_CS(tab=resultat[[1]],nom_outcome=yname,tname=tname,first.treat.name=gname,dum)
  
  influ<-agregat_influence_CS(tab=resultat[[1]],array_inf=resultat[[2]],listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
  list(result,influ,list_id)
}



# ### SUR DONNEES DE L'ENQUETE RD
# estimPeriodRD_CS<-function(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident,Les_outcome
#                             ,Pond_outcome,STRATE,popo){
#   ### selection des donn�es 
#   bebe<-PAPA[PAPA$annee>=deb & PAPA$annee <= fin ,]
#   bebe<-bebe[bebe[,anneeT]>=debT | bebe[,anneeT]==0,]### cette condition reste n�cessaire
#   ## Pour les entreprises trait�es apr�s finT on met une pond�ration nulle pour les observation apr�s finT
#   for(i in 1:length(Pond_outcome)){
#     bebe[(bebe[,anneeT]>finT)&(bebe[,"annee"]>=finT),Pond_outcome[i]]<-0
#   }
#   ## puis on bascule enti�rement dans le  contrefactuel les trait�s apr�s finT
#   bebe[bebe[,anneeT]>finT,anneeT]<-0
#   #bebe<-bebe[bebe[,anneeT]<=finT | bebe[,anneeT]==0,]### cette conditions supprime du contrefactuel les traitements suivants !! 
#   ### mais indispensable pour limiter la liste des traitements pris en compte : du coup on arr�te des faire de estimations par sous-p�riode
#   ### definition du mod�le utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(var_score, collapse=" + "),sep=""))
#   baba<-bebe
#   resultat<- mp.spatt.RD_CS(nom_outcome=Les_outcome,nom_taitement=TT,xformla=xxF,data=baba,
#                          first.treat.name = anneeT,
#                          idname = nom_ident, tname="annee",
#                          bstrap = FALSE,se=TRUE,cband =FALSE
#                          ,selection=nom_select,ponderation=Pond_outcome,
#                          debT2=debT,finT2=finT,strate=STRATE,PondRD=popo)
#   # resultat
#   dum<-as.data.frame(table(bebe[bebe[,anneeT]>=debT & bebe[,anneeT]<=finT,anneeT]))
#   # agregat(tab=resultat,nom_outcome=outC,tname="annee",first.treat.name=anneeT) 
#   agregatChris_CS(tab=resultat,nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum) 
# }

# ## Fonctions pour calcul bootstrap 
# precizionRD_CS<-function(NOM,CHEMI,baz,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,
#                       Les_outcome,Pond_outcome,NN,strate,PP){
#   sisiz<-unique(baz[,c("siren",anneeT)])
#   Martu<-function(samp){
#     sampi<-subset(samp,select=siren)
#     nn<-dim(sampi)[1]
#     sampi$identif<-c(1:nn)
#     PAPA<-merge(baz,sampi,by="siren")
#     lele<-estimPeriodRD_CS(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident="identif",Les_outcome,Pond_outcome
#                         ,STRATE=strate,popo=PP)
#     Les_outcome=c(Les_outcome,"nobsG","nobsC")
#     resu<-lele[,Les_outcome[1]]
#     nnn<-length(Les_outcome)
#     for (i in 2:nnn){
#       resu<-c(resu,lele[,Les_outcome[i]])
#     }
#     resu
#   }
#   estu<-function(tab,ind){
#     print(compteur) ; compteur<-compteur+1
#     Martu(tab[ind,])
#   }
#   compteur<-1
#   vv<-boot(sisiz,estu,NN,sim="ordinary",stype="i",sisiz[,anneeT])
#   nom_tt<-paste0(NOM,"_tt")
#   sauve(vv$t,nom_tt,CHEMI)
#   nom_t0<-paste0(NOM,"_t0")
#   sauve(vv$t0,nom_t0,CHEMI)
#   vv}



# ### SUR DONNEES EXHAUSTIVES
# estimPeriod6_CS<-function(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident,
#                           Les_outcome,Pond_outcome,STRATE){
#   ## PAPA : la base de calcul
#   # EXCLU : la variable donnant la condtiion d'exclusion 
#   # TT : 
#   # var_score : variable(s) permettant de d�finir le score
#   # anneeT : la varaible donnant la g�n�ratio nde traitement 
#   # debT : premi�re g�n�ration de traitement consid�r�e
#   # finT : derni�re g�n�ratio nde traitement consid�r�e
#   # deb : premi�re ann�e de la base de calcul
#   # fin : erni�re ann�e de la base d calcul 
#   # nom_select : nom de la variable de s�lection 
#   # nom_ident : nom de la variable permettant d'identifier les individus
#   # Les_outcome : vecteur des variables outcome
#   # Pond_outcome : vecteur des variabels permettnat de pond�rer les observations
#   ### selection des donn�es 
#   bebe<-PAPA[PAPA$annee>=deb & PAPA$annee <= fin ,]
#   bebe<-bebe[bebe[,anneeT]>=debT | bebe[,anneeT]==0,]### cette condition reste n�cessaire
#   ## Pour les entreprises trait�es apr�s finT on met une pond�ration nulle pour les observation apr�s finT
#   for(i in 1:length(Pond_outcome)){
#     bebe[(bebe[,anneeT]>finT)&(bebe[,"annee"]>=finT),Pond_outcome[i]]<-0
#   }
#   ## puis on bascule enti�rement dans le  contrefactuel les trait�s apr�s finT
#   bebe[bebe[,anneeT]>finT,anneeT]<-0
  
#   #bebe<-bebe[bebe[,anneeT]<=finT | bebe[,anneeT]==0,]### cette conditions supprime du contrefactuel les traitements suivants !! 
#   ### mais indispensable pour limiter la liste des traitements pris en compte : du coup on arr�te des faire de estimations par sous-p�riode
#   ### definition du mod�le utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(var_score, collapse=" + "),sep=""))
#   baba<-bebe
#   resultat<- mp.spatt.V6_CS(nom_outcome=Les_outcome,nom_taitement=TT,xformla=xxF,data=baba,
#                          first.treat.name = anneeT,
#                          idname = nom_ident, tname="annee",
#                          bstrap = FALSE,se=TRUE,cband =FALSE
#                          ,selection=nom_select,ponderation=Pond_outcome,debT2=debT,finT2=finT,strate=STRATE)
#   # resultat
#   #dum<-table(bebe[bebe[,anneeT]>=debT & bebe[,anneeT]<=finT,anneeT])
#   dum<-as.data.frame(table(bebe[bebe[,anneeT]>=debT & bebe[,anneeT]<=finT,anneeT]))
  
#   # agregat(tab=resultat,nom_outcome=outC,tname="annee",first.treat.name=anneeT) 
#   #agregatP(tab=resultat,nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum)
#   agregatChris_CS(tab=resultat,nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum)
# }

# ## Fonctions pour calcul bootstrap 
# precizionS_CS<-function(NOM,CHEMI,baz,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,Les_outcome,
#                      Pond_outcome,NN,strate){
#   sisiz<-unique(baz[,c("siren",anneeT)])
#   Martu<-function(samp){
#     sampi<-subset(samp,select=siren)
#     nn<-dim(sampi)[1]
#     sampi$identif<-c(1:nn)
#     PAPA<-merge(baz,sampi,by="siren")
#     lele<-estimPeriod6_CS(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident="identif",Les_outcome,Pond_outcome,STRATE=strate)
#     Les_outcome=c(Les_outcome,"nobsG","nobsC")
#     resu<-lele[,Les_outcome[1]]
#     nnn<-length(Les_outcome)
#     for (i in 2:nnn){
#       resu<-c(resu,lele[,Les_outcome[i]])
#     }
#     resu
#   }
#   estu<-function(tab,ind){
#     print(compteur) ; compteur<-compteur+1
#     Martu(tab[ind,])
#   }
#   compteur<-1
#   vv<-boot(sisiz,estu,NN,sim="ordinary",stype="i",sisiz[,anneeT])
#   nom_tt<-paste0(NOM,"_tt")
#   sauve(vv$t,nom_tt,CHEMI)
#   nom_t0<-paste0(NOM,"_t0")
#   sauve(vv$t0,nom_t0,CHEMI)
#   vv}

# ###### ----------------------------------------------------------------------------------------------------------------------------------------
# ###### ESTIMATEUR EN LONG DID  ----------------------------------------------------------------------------------------------------------------
# ###### ----------------------------------------------------------------------------------------------------------------------------------------
# longDID_estimPeriod_Boot<-function(yname,tname,idname,gname,xformla,data,debT,finT,deb,fin,select,weightsname,alp=0.05,bstrap,biters=1000,STRATE,pond_RD=NULL,treated){
  
  
#   ### Créer un siren numérique (id)
#   list_id <-as.data.frame(unique(data[,c(idname)]))
#   list_id$iden_num<-1:dim(list_id)[1]
#   colnames(list_id)[1]=idname
#   data=merge(data,list_id,by.x=idname,by.y=idname,all.x=TRUE)
  
#   ### selection des données 
#   bebe<-data[data[[tname]]>=deb & data[[tname]] <= fin ,]
#   bebe<-bebe[bebe[,gname]>=debT | bebe[,gname]==0,]
  
#   ## Pour les observations traitées après finT on met une pondération nulle pour les observation après finT
#   ## puis on bascule entièrement dans le  contrefactuel les traités après finT
#   for(i in 1:length(weightsname)){
#   bebe[(bebe[,gname]>finT)&(bebe[[tname]]>=finT),weightsname[i]]<-0}
#   bebe[bebe[,gname]>finT,gname]<-0
  
#   ### definition du modèle utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(xformla, collapse=" + "),sep=""))
  
#   # resultat: objet liste avec trois éléments dans la liste: (1) att, (2) mat_influence et (3) indiv la liste des individus 
#   resultat<- mp.spatt.Boot_longDID(nom_outcome=yname,nom_traitement=treated,xformla=xxF,data=bebe,
#                             first.treat.name = gname,
#                             idname = idname, tname=tname,
#                             bstrap = FALSE,se=TRUE,cband =FALSE
#                             ,selection=select,ponderation=weightsname,debT2=debT,finT2=finT,strate=STRATE)
  
#   # Poids pour aggréger les effets des différentes cohortes de traitement 
#   list_id<-merge(list_id,resultat[[3]])
#   list_id_poids=resultat[[3]]
#   list_id_poids=list_id_poids[list_id_poids[,gname]>0,]
#   list_id_poids<-unique(list_id_poids)
#   dum<-as.data.frame(table(list_id_poids[,gname]))
  
#   result<-agregatChris_longDID(tab=resultat[[1]],nom_outcome=yname,tname=tname,first.treat.name=gname,dum)
#   influ<-agregat_influence_longDID(tab=resultat[[1]],array_inf=resultat[[2]],listG=resultat[[3]],nom_outcome=yname,tname=tname,first.treat.name=gname,poids=dum)
#   list(result,influ,list_id)
# }



# ### SUR DONNEES EXHAUSTIVES (FORCEMENT, PAS POSSIBLE SUR DONNEES ENQUETE RD)
# estimPeriod6_longDID<-function(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,idname,
#                                    Les_outcome,Pond_outcome,STRATE){
#   ## PAPA : la base de calcul
#   # EXCLU : la variable donnant la condtiion d'exclusion 
#   # TT : 
#   # var_score : variable(s) permettant de d�finir le score
#   # anneeT : la varaible donnant la g�n�ratio nde traitement 
#   # debT : premi�re g�n�ration de traitement consid�r�e
#   # finT : derni�re g�n�ratio nde traitement consid�r�e
#   # deb : premi�re ann�e de la base de calcul
#   # fin : erni�re ann�e de la base d calcul 
#   # nom_select : nom de la variable de s�lection 
#   # nom_ident : nom de la variable permettant d'identifier les individus
#   # Les_outcome : vecteur des variables outcome
#   # Pond_outcome : vecteur des variabels permettnat de pond�rer les observations
#   ### selection des donn�es 
#   bebe<-PAPA[PAPA$annee>=deb & PAPA$annee <= fin ,]
#   bebe<-bebe[bebe[,anneeT]>=debT | bebe[,anneeT]==0,]### cette condition reste n�cessaire
#   ## Pour les entreprises trait�es apr�s finT on met une pond�ration nulle pour les observation apr�s finT
#   for(i in 1:length(Pond_outcome)){
#     bebe[(bebe[,anneeT]>finT)&(bebe[,"annee"]>=finT),Pond_outcome[i]]<-0
#   }
#   ## puis on bascule enti�rement dans le  contrefactuel les trait�s apr�s finT
#   bebe[bebe[,anneeT]>finT,anneeT]<-0
  
#   #bebe<-bebe[bebe[,anneeT]<=finT | bebe[,anneeT]==0,]### cette conditions supprime du contrefactuel les traitements suivants !! 
#   ### mais indispensable pour limiter la liste des traitements pris en compte : du coup on arr�te des faire de estimations par sous-p�riode
#   ### definition du mod�le utilisant le score 
#   xxF<-as.formula(paste(" ~ ",paste(var_score, collapse=" + "),sep=""))
#   baba<-bebe
#   resultat<- mp.spatt.V6_longDID(nom_outcome=Les_outcome,nom_taitement=TT,xformla=xxF,data=baba,
#                                    first.treat.name = anneeT,
#                                    idname = nom_ident, tname="annee",
#                                    bstrap = FALSE,se=TRUE,cband =FALSE
#                                    ,selection=nom_select,ponderation=Pond_outcome,debT2=debT,finT2=finT,strate=STRATE)
#   dum<-as.data.frame(table(bebe[bebe[,anneeT]>=debT & bebe[,anneeT]<=finT,anneeT]))
#   # agregat(tab=resultat,nom_outcome=outC,tname="annee",first.treat.name=anneeT) 
#   agregatChris_longDID(tab=resultat,nom_outcome=Les_outcome,tname="annee",first.treat.name=anneeT,dum)
# }



# ## Fonctions pour calcul bootstrap 
# precizionS_longDID<-function(NOM,CHEMI,baz,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,Les_outcome,
#                         Pond_outcome,NN,strate){
#   sisiz<-unique(baz[,c("siren",anneeT)])
#   Martu<-function(samp){
#     sampi<-subset(samp,select=siren)
#     nn<-dim(sampi)[1]
#     sampi$identif<-c(1:nn)
#     PAPA<-merge(baz,sampi,by="siren")
#     lele<-estimPeriod6_longDID(PAPA,TT,var_score,anneeT,debT,finT,deb,fin,nom_select,nom_ident="identif",Les_outcome,Pond_outcome,STRATE=strate)
#     Les_outcome=c(Les_outcome,"nobsG","nobsC")
#     resu<-lele[,Les_outcome[1]]
#     nnn<-length(Les_outcome)
#     for (i in 2:nnn){
#       resu<-c(resu,lele[,Les_outcome[i]])
#     }
#     resu
#   }
#   estu<-function(tab,ind){
#     print(compteur) ; compteur<-compteur+1
#     Martu(tab[ind,])
#   }
#   compteur<-1
#   vv<-boot(sisiz,estu,NN,sim="ordinary",stype="i",sisiz[,anneeT])
#   nom_tt<-paste0(NOM,"_tt")
#   sauve(vv$t,nom_tt,CHEMI)
#   nom_t0<-paste0(NOM,"_t0")
#   sauve(vv$t0,nom_t0,CHEMI)
#   vv}

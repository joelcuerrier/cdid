panelDiffV<-function (data, timevars,pondera,idname, tname,t0,t1,t2,selecta)
{
####### panelDiffV(disdat,nom_outcome,ponderation,idname, tname,pret,tlist[t],tlist[t+1],selection)

  

  ## cette fonction part d'un panele cylindr� (ordonn� sur les dates !!!) sur trois dates t0,t1,t2 et fabrique une base avec une seule
  # ligne par individu contenant des informations sur t0 (en particuler le score), les �volutions entre t1 et t2 des outcome et la pond�ration de 
  # cette �volution
  # t0 : la date de r�f�rence (typiquement celle pr�c�dant l'entree dans le traitement)
  # t1,t2 : on calcule des �volution entre les deux dates t1 et t2 (typiquement t et t+1)

  #poids
  ppp = data[data[[tname]] == t1,]
  ppp=subset(ppp, select = c(idname, pondera))
  ppp<-ppp[order(ppp[, idname]),]
  
  retdat <- data[data[[tname]] == t0, ]
  retdat<-retdat[order(retdat[, idname]),]
  retdat[,pondera] <- retdat[,selecta]*ppp[,pondera]
    
  dt1=data[data[[tname]]==t1,]
  dt1=subset(dt1, select = c(timevars, idname)) #timevars: c'est la variable Y

  dt2=data[data[[tname]]==t2,]
  dt2=subset(dt2, select = c(timevars, idname))

  dt1<-dt1[order(dt1[, idname]),] ; dt2<-dt2[order(dt2[, idname]),]

  Dtimevars<-paste0("D",timevars) #"DY1_chaine" "DY2_chaine"
  retdat[,Dtimevars] <- dt2[,timevars] - dt1[,timevars] #Donc tu fais DY= Yt+1-Yt, ça devrait pas être Yt-Yt-1?
    
  return(retdat)
}


## Pour l'estimateur en cross section
panelDiffV_CS<-function (data, timevars,pondera,idname, tname,t0,t1,selecta)
{
  ## cette fonction part d'un panele cylindr� (ordonn� sur les dates !!!) sur deux dates t0,t1 et fabrique une base avec une seule
  # ligne par individu contenant des informations sur t0 (en particuler le score)
  # t0 : la date de r�f�rence (typiquement celle pr�c�dant l'entree dans le traitement)
  # t1 : typiquement t

  # ppp<-data[data[, tname] == t1,c(idname,pondera)]
  ppp = data[data[[tname]] == t1,]
  ppp=subset(ppp, select = c(idname, pondera))
  ppp<-ppp[order(ppp[, idname]),]
  # ppp<-ppp[order(ppp[, ..idname]),]

  # retdat <- subset(data, data[, tname] == t0)
  retdat <- data[data[[tname]] == t0, ]
  retdat<-retdat[order(retdat[, idname]),]
  # retdat<-retdat[order(retdat[, ..idname]),]
  retdat[,pondera] <- retdat[,selecta]*ppp[,pondera]

#   for (i in seq_along(pondera)) {
#     retdat[[pondera[i]]] <- retdat[,..selecta]*ppp[[pondera[i]]]
# }


  # dt1<-data[data[, tname] == t1,c(timevars, idname)]
  # dt1<-dt1[order(dt1[, idname]),] 
  
  dt1=data[data[[tname]]==t1,]
  dt1=subset(dt1, select = c(timevars, idname))
  # dt1<-dt1[order(dt1[, ..idname]),]
  dt1<-dt1[order(dt1[, idname]),]
  
  Ltimevars<-paste0("L",timevars)
  # retdat[,Ltimevars] <- dt1[,timevars]
  retdat[,Ltimevars] <- subset(dt1,select=timevars)
  return(retdat)
}


## Pour l'estimateur en long Diff-in-Diff
panelDiffV_longDID<-function (data, timevars,pondera,idname, tname,t0,t1,selecta)
{
  ## cette fonction part d'un panele cylindr� (ordonn� sur les dates !!!) sur trois dates t0,t1 et fabrique une base avec une seule
  # ligne par individu contenant des informations sur t0 (en particuler le score)
  # t0 : la date de r�f�rence g-1 (typiquement celle pr�c�dant l'entree dans le traitement)
  # t1 : typiquement t
  # ppp<-data[data[, tname] == t1,c(idname,pondera)]
  
  
  ppp = data[data[[tname]] == t1,]
  ppp=subset(ppp, select = c(idname, pondera))  
  ppp<-ppp[order(ppp[, idname]),]
  # ppp<-ppp[order(ppp[, ..idname]),]

  # retdat <- subset(data, data[, tname] == t0)
  retdat <- data[data[[tname]] == t0, ]
  # retdat<-retdat[order(retdat[, ..idname]),]
  retdat<-retdat[order(retdat[, idname]),]
  # Pond�ration pour �tre sur que l'entreprise est bien pr�sente en g-1 et en t
  retdat[,pondera] <- retdat[,selecta]*retdat[,pondera]*ppp[,pondera]
#   for (i in seq_along(pondera)) {
#     retdat[[pondera[i]]] <- retdat[,..selecta]*ppp[[pondera[i]]] *retdat[[pondera[[i]]]]
#                             # retdat[,selecta]  *ppp[,pondera]     *retdat[,pondera]
# }

  # dt1<-data[data[, tname] == t1,c(timevars, idname)]
  dt1=data[data[[tname]]==t1,]
  dt1=subset(dt1, select = c(timevars, idname))
  
  # dt0<-data[data[, tname] == t0,c(timevars, idname)]
  dt0=data[data[[tname]]==t0,]
  dt0=subset(dt0, select = c(timevars, idname))
  
  # dt1<-dt1[order(dt1[, ..idname]),] 
  # dt0<-dt0[order(dt0[, ..idname]),]
    dt1<-dt1[order(dt1[, idname]),] 
  dt0<-dt0[order(dt0[, idname]),]
  
  Dtimevars<-paste0("D",timevars)
  
  # retdat[,Dtimevars] <- dt1[,timevars]-dt0[,timevars]
  # retdat[,Dtimevars] <- dt1[,..timevars]-dt0[,..timevars]
  retdat[,Dtimevars] <- dt1[,timevars]-dt0[,timevars]
  return(retdat)
}


## Pour l'estimateur en long Diff-in-Diff : pour simulation, pour travailler sur un sous-�chantillon
panelDiffV_longDID_alt<-function (data, timevars,pondera,idname, tname,t0,t1,selecta)
{
  ## cette fonction part d'un panele cylindr� (ordonn� sur les dates !!!) sur trois dates t0,t1 et fabrique une base avec une seule
  # ligne par individu contenant des informations sur t0 (en particuler le score)
  # t0 : la date de r�f�rence g-1 (typiquement celle pr�c�dant l'entree dans le traitement)
  # t1 : typiquement t
  ppp<-data[data[, tname] == t1,c(idname,pondera)]
  ppp<-ppp[order(ppp[, idname]),]
  retdat <- subset(data, data[, tname] == t0)
  retdat<-retdat[order(retdat[, idname]),]
  # Pond�ration pour �tre sur que l'entreprise est bien pr�sente en g-1 et en t
  retdat[,pondera] <- retdat[,selecta]*retdat[,pondera]
  #commande � ex�cuter normalement
  #retdat[,pondera] <- retdat[,selecta]*retdat[,pondera]*ppp[,pondera]
  dt1<-data[data[, tname] == t1,c(timevars, idname)]
  dt0<-data[data[, tname] == t0,c(timevars, idname)]
  dt1<-dt1[order(dt1[, idname]),] 
  dt0<-dt0[order(dt0[, idname]),]
  Dtimevars<-paste0("D",timevars)
  retdat[,Dtimevars] <- dt1[,timevars]-dt0[,timevars]
  return(retdat)
}

# Pour l'estimateur GMM avec les delta k ATT
panelDiffV_GMM<-function (data, timevars,pondera,idname, tname,t0,t1,t2,selecta)
{
  ## cette fonction part d'un panele cylindr� (ordonn� sur les dates !!!) sur trois dates t0,t1,t2 et fabrique une base avec une seule
  # ligne par individu contenant des informations sur t0 (en particuler le score), les �volutions entre t1 et t2 des outcome et la pond�ration de 
  # cette �volution
  # t0 : la date de r�f�rence (typiquement celle pr�c�dant l'entree dans le traitement)
  # t1,t2 : on calcule des �volution entre les deux dates t1 et t2 (typiquement t et t+k)
  # on applique la ponderation St*St+k (i.e., le fait d'etre echantillonne en t et en t+k)
  
  # ppp1<-data[data[, tname] == t1,c(idname,pondera)]
  # ppp1<-ppp1[order(ppp1[, idname]),]
  ppp1 = data[data[[tname]] == t1,]
  ppp1=subset(ppp1, select = c(idname, pondera))
  # ppp1<-ppp1[order(ppp1[, ..idname]),]
  ppp1<-ppp1[order(ppp1[, idname]),]

  # ppp2<-data[data[, tname] == t2,c(idname,pondera)]
  # ppp2<-ppp2[order(ppp2[, idname]),]
  ppp2 = data[data[[tname]] == t2,]
  ppp2=subset(ppp2, select = c(idname, pondera))
  # ppp2<-ppp2[order(ppp2[, ..idname]),]
  ppp2<-ppp2[order(ppp2[, idname]),]



  # retdat <- subset(data, data[, tname] == t0)
  # retdat<-retdat[order(retdat[, idname]),]
  retdat <- data[data[[tname]] == t0, ]
  # retdat<-retdat[order(retdat[, ..idname]),]
  retdat<-retdat[order(retdat[, idname]),]

  

  retdat[,pondera] <- retdat[,selecta]*ppp1[,pondera]*ppp2[,pondera]
    # for (i in seq_along(pondera)) {
    # retdat[[pondera[i]]] <- retdat[,..selecta]*ppp1[[pondera[i]]]*ppp2[[pondera[i]]]
  # }



  # dt1<-data[data[, tname] == t1,c(timevars, idname)]
  dt1=data[data[[tname]]==t1,]
  dt1=subset(dt1, select = c(timevars, idname))
  # dt1<-dt1[order(dt1[, ..idname]),]
  dt1<-dt1[order(dt1[, idname]),]

  # dt2<-data[data[, tname] == t2,c(timevars, idname)]
  dt2<-data[data[[tname]] == t2,]
  dt2=subset(dt2,select=c(timevars, idname))
  # dt2<-dt2[order(dt2[, ..idname]),]
  dt2<-dt2[order(dt2[, idname]),]

  Dtimevars<-paste0("D",timevars)
  # retdat[,Dtimevars] <- dt2[,..timevars] - dt1[,..timevars]
  retdat[,Dtimevars] <- dt2[,timevars] - dt1[,timevars]
  return(retdat)
}



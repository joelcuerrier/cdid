## fonction d'agr�gation 
library(Matrix)
#Voir l'équation de la page 17. Tu vas catch aussi pouruqoi il fait cumsum.
agregatChris<-function(tab,nom_outcome,tname,first.treat.name,poids){
  
  
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  ATT=data.frame()
  
  for (i in 1:NN){
    tab_cohort=tab[(nn*(i-1)+1):(nn*i),]
    after=tab_cohort[tab_cohort[,tname]>=an_cohort[i],]
    befor=tab_cohort[tab_cohort[,tname]<an_cohort[i],]
    befor=befor[order(-befor[,tname]),]
    ATT_after=cumsum(after[,nom_outcome])
    ATT_befor=cumsum(-befor[,nom_outcome])
    ATT_after=cbind(after[,2:5],ATT_after)
    ATT_befor=cbind(befor[,2:5],ATT_befor)
    ATT_after$time_to_treatment=ATT_after[,tname]-ATT_after[,first.treat.name]+1
    ATT_befor$time_to_treatment=ATT_befor[,tname]-ATT_befor[,first.treat.name]
    ATT_temp=rbind(ATT_after,ATT_befor)
    ATT_temp=merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT=rbind(ATT,ATT_temp)
  }
  num=ATT[,nom_outcome]*ATT$Freq
  num$time_to_treatment=ATT$time_to_treatment
  num$nobsG=ATT$nobsG
  num$nobsC=ATT$nobsC
  num=aggregate(num[,c(nom_outcome,"nobsG","nobsC")],by=list(num$time_to_treatment),FUN=sum) 
  denom=aggregate(ATT$Freq,by=list(ATT$time_to_treatment),FUN=sum)
  denom=1/denom$x
  result=num 
  
  result[,nom_outcome]=result[,nom_outcome]*denom
  
  result=rbind(result,c(0,rep(0,length(nom_outcome)+2)))
  
  
  result[result$Group.1==-1,"nobsG"]=result[result$Group.1==-2,"nobsG"]
  result[result$Group.1==-1,"nobsC"]=result[result$Group.1==-2,"nobsC"]
  
  result=result[order(result[,"Group.1"]),]
  result
}


agregatChris_GMM<-function(tab,nom_outcome,tname,first.treat.name,poids){

  # View(tab)
  tab$time_to_treatment=tab[,tname]-(tab[,first.treat.name]-1)
  tab=merge(tab,poids,by.x=first.treat.name, by.y="Var1",all.x)
  num=tab[,nom_outcome]*tab$Freq
  num$time_to_treatment=tab$time_to_treatment
  num=aggregate(num[,c(nom_outcome)],by=list(num$time_to_treatment),FUN=sum) 
  denom=aggregate(tab$Freq,by=list(tab$time_to_treatment),FUN=sum)
  denom=1/denom$x
  result=num 
  result[,nom_outcome]=result[,nom_outcome]*denom
  
  
  
  # print(c(0,rep(0,length(nom_outcome)+2))) #Voir la matrice results, elle est differente avec le 3.6.1. Dans celle ci, il y am oins de colonne donc j'enleve le +2 pour ne pas avoir le message d'avertissement...
  result=rbind(result,c(0,rep(0,length(nom_outcome))))
  
  

  result=result[order(result[,"Group.1"]),]
  result
}

agregatChris_CS<-function(tab,nom_outcome,tname,first.treat.name,poids){
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  ATT=data.frame()
  for (i in 1:NN){
    tab_cohort=tab[(nn*(i-1)+1):(nn*i),]
    ref=matrix(rep(as.numeric(tab_cohort[tab_cohort[,tname]==(an_cohort[i]-1),nom_outcome]),nrow(tab_cohort)),nrow=nrow(tab_cohort),byrow=TRUE)
    tab_cohort[,nom_outcome]=as.matrix(tab_cohort[,nom_outcome])-ref
    ATT_temp=tab_cohort[,2:ncol(tab_cohort)]
    ATT_temp$time_to_treatment=tab_cohort[,tname]-(tab_cohort[,first.treat.name]-1)
    ATT_temp=merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT=rbind(ATT,ATT_temp)
  }
  num=ATT[,nom_outcome]*ATT$Freq
  num$time_to_treatment=ATT$time_to_treatment
  num$nobsG=ATT$nobsG
  num$nobsC=ATT$nobsC
  num=aggregate(num[,c(nom_outcome,"nobsG","nobsC")],by=list(num$time_to_treatment),FUN=sum) 
  denom=aggregate(ATT$Freq,by=list(ATT$time_to_treatment),FUN=sum)
  denom=1/denom$x
  result=num 
  result[,nom_outcome]=result[,nom_outcome]*denom
  result
}


agregatChris_longDID<-function(tab,nom_outcome,tname,first.treat.name,poids){
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  ATT=data.frame()
  for (i in 1:NN){
    tab_cohort=tab[(nn*(i-1)+1):(nn*i),]
    ATT_temp=tab_cohort[,2:ncol(tab_cohort)]
    ATT_temp$time_to_treatment=tab_cohort[,tname]-(tab_cohort[,first.treat.name]-1)
    ATT_temp=merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT=rbind(ATT,ATT_temp)
  }
  num=ATT[,nom_outcome]*ATT$Freq
  num$time_to_treatment=ATT$time_to_treatment
  num$nobsG=ATT$nobsG
  num$nobsC=ATT$nobsC
  num=aggregate(num[,c(nom_outcome,"nobsG","nobsC")],by=list(num$time_to_treatment),FUN=sum) 
  denom=aggregate(ATT$Freq,by=list(ATT$time_to_treatment),FUN=sum)
  denom=1/denom$x
  result=num 
  result[,nom_outcome]=result[,nom_outcome]*denom
  result
}


agregat_influence<-function(tab,array_inf,listG,nom_outcome,tname,first.treat.name,poids){
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  agreg_influence = array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),nn*NN))
  cumsum_after<-function(array_f){
    for (k in 2:dim(array_f)[3]){
      array_f[,,k]=array_f[,,k]+array_f[,,k-1]
    }
    array_f
  }
  cumsum_befor<-function(array_f){
    len = dim(array_f)[3]
    for (k in 1:(len-1)){
      j=len-k
      array_f[,,j]=array_f[,,j]+array_f[,,j+1]
    }
    array_f[,,1:len]=(-1)*array_f[,,1:len]
  }
  sum_array<-function(array_f){
    len = dim(array_f)[3]
    mat_result=array_f[,,1]
    if (len>1){
    for (k in 2:len){
      mat_result=mat_result+array_f[,,k]
      }
    }
    mat_result
  }
  ATT=data.frame()
  for (i in 1:NN){
    tab_cohort=tab[(nn*(i-1)+1):(nn*i),]
    after=tab_cohort[tab_cohort[,tname]>=an_cohort[i],]
    befor=tab_cohort[tab_cohort[,tname]<an_cohort[i],]
    after$time_to_treatment=after[,tname]-after[,first.treat.name]
    befor$time_to_treatment=befor[,tname]-befor[,first.treat.name]-1
    after$compteur
    befor$compteur
    influence_cohort_after=array_inf[,,after$compteur]
    influence_cohort_befor=array_inf[,,befor$compteur]
    if (is.na(dim(influence_cohort_after)[3])==TRUE){
      influence_cohort_after=array_inf[,,after$compteur]
    } else if (dim(influence_cohort_after)[3]==0) {
      influence_cohort_after=array_inf[,,after$compteur]
    } else {
      influence_cohort_after=cumsum_after(influence_cohort_after)
      # Modif
      influence_cohort_after[,1,]=array_inf[,1,after$compteur]
    }
    if (is.na(dim(influence_cohort_befor)[3])==TRUE){
      influence_cohort_befor=array_inf[,,befor$compteur]
    } else if (dim(influence_cohort_befor)[3]==0) {
      influence_cohort_befor=array_inf[,,befor$compteur]
    } else {
      influence_cohort_befor=cumsum_befor(influence_cohort_befor)
      # Modif
      influence_cohort_befor[,1,]=array_inf[,1,befor$compteur]
    }
    agreg_influence[,,after$compteur]<-influence_cohort_after
    agreg_influence[,,befor$compteur]<-influence_cohort_befor
    ATT_temp=rbind(befor,after)
    ATT_temp=merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT=rbind(ATT,ATT_temp)
  }
  ttt=unique(sort(ATT$time_to_treatment))
  agreg_influence_final = array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),length(ttt)+1))
  for (i in 1:length(ttt)){
    ATTi=ATT[ATT$time_to_treatment==ttt[i],]
    ATTi$Freq=ATTi$Freq/sum(ATTi$Freq)
    agreg_inf_i=agreg_influence[,,ATTi$compteur]
    if (is.na(dim(agreg_inf_i)[3])==TRUE){
      agreg_inf_i=array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),1))
      agreg_inf_i[,,1]=agreg_influence[,,ATTi$compteur]
    } else {
        agreg_inf_i=agreg_influence[,,ATTi$compteur]
        for (j in 1:dim(ATTi)[1]){
          agreg_inf_i[,,j]=ATTi$Freq[j]*agreg_inf_i[,,j]
        }
        # Modif
        agreg_inf_i[,1,]=agreg_influence[,1,ATTi$compteur]
    }
    # D�but modif
    # Comme les poids sont des fr�quences empiriques, ils sont al�atoires, 
    # il faut ajouter un terme suppl�mentaire dans la fonction d'influence pour ne pas sous-�valuer la variance
    for (j in 1:dim(ATTi)[1]){
      traitement_j <- listG
      traitement_j[,1][traitement_j[,1]!= ATTi[j,1]] <- 0
      traitement_j[,1][traitement_j[,1]== ATTi[j,1]] <- 1
      traitement_j[,1] <- (traitement_j[,1]-mean(traitement_j[,1])) * 1/(dim(ATTi)[1])
      # Doute sur le fait que listG & agreg_inf_i ait le m�me nb d'observations
      # mat_gamma <- as.data.frame(agreg_inf_i[,1,j])
      # colnames(mat_gamma) <- c("var1")
      # mat_gamma$gamma <- merge(mat_gamma,traitement_j,by.x="var1", by.y="iden_num",all.x)
      # agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ mat_gamma$gamma%*%ATTi[j,nom_outcome]
      agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ t(t(traitement_j[,1]))%*%as.matrix(ATTi[j,nom_outcome])
    }
    # End modif
    if (ttt[i]< -1 ) {
      agreg_influence_final[,,i]=sum_array(agreg_inf_i)
    } else {
      agreg_influence_final[,,i+1]=sum_array(agreg_inf_i)
    }
  }
  # Remettre les identifiants des individus
  for (i in 1:(length(ttt)+1)){
    agreg_influence_final[,1,i]=array_inf[,1,1]
  }
  agreg_influence_final
}  



agregat_influence_CS<-function(tab,array_inf,listG,nom_outcome,tname,first.treat.name,poids){
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  sum_array<-function(array_f){
    len = dim(array_f)[3]
    mat_result=array_f[,,1]
    if (len>1){
      for (k in 2:len){
        mat_result=mat_result+array_f[,,k]
      }
    }
    mat_result
  }
  ATT=data.frame()
  array_inf_did = array(0,dim=c(dim(array_inf)[1],dim(array_inf)[2],dim(array_inf)[3]))
  for (i in 1:NN){
    # Diff de diff pour les ATT
    tab_cohort = tab[(nn*(i-1)+1):(nn*i),]
    ref = matrix(rep(as.numeric(tab_cohort[tab_cohort[,tname]==(an_cohort[i]-1),nom_outcome]),nrow(tab_cohort)),nrow=nrow(tab_cohort),byrow=TRUE)
    tab_cohort[,nom_outcome] = as.matrix(tab_cohort[,nom_outcome])-ref
    ATT_temp = tab_cohort[,1:ncol(tab_cohort)]
    ATT_temp$time_to_treatment = tab_cohort[,tname]-tab_cohort[,first.treat.name]
    ATT_temp = merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT = rbind(ATT,ATT_temp)
    # Diff de diff pour les fonction d'influence
    num_compteur = tab_cohort[tab_cohort[,tname]==(an_cohort[i]-1),1]
    inf_cohort = array_inf[,,(nn*(i-1)+1):(nn*i)]
    inf_ref = as.matrix(array_inf[,2:dim(array_inf)[2],num_compteur])
    inf_ref_array = array(inf_ref,c(dim(inf_ref)[1],dim(inf_ref)[2],dim(inf_cohort)[3]))
    inf_cohort[,2:dim(array_inf)[2],] = inf_cohort[,2:dim(array_inf)[2],] - inf_ref_array
    array_inf_did[,,(nn*(i-1)+1):(nn*i)] = inf_cohort 
  }
  ttt=unique(sort(ATT$time_to_treatment))
  agreg_influence_final = array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),length(ttt)))
  # Agr�gation de la fonction d'influence au niveau "time to treatment"
  for (i in 1:length(ttt)){
    ATTi=ATT[ATT$time_to_treatment==ttt[i],]
    ATTi$Freq=ATTi$Freq/sum(ATTi$Freq)
    agreg_inf_i=array_inf_did[,,ATTi$compteur]
    if (is.na(dim(agreg_inf_i)[3])==TRUE){
      agreg_inf_i=array(0,dim=c(dim(array_inf_did)[1],(1+length(nom_outcome)),1))
      agreg_inf_i[,,1]=array_inf_did[,,ATTi$compteur]
    } else {
      agreg_inf_i=array_inf_did[,,ATTi$compteur]
      for (j in 1:dim(ATTi)[1]){
        agreg_inf_i[,,j]=ATTi$Freq[j]*agreg_inf_i[,,j]
      }
      # Modif
      agreg_inf_i[,1,]=array_inf_did[,1,ATTi$compteur]
    }
    # D�but modif
    # Comme les poids sont des fr�quences empiriques, ils sont al�atoires, 
    # il faut ajouter un terme suppl�mentaire dans la fonction d'influence pour ne pas sous-�valuer la variance
    for (j in 1:dim(ATTi)[1]){
      traitement_j <- listG
      traitement_j[,1][traitement_j[,1]!= ATTi[j,1]] <- 0
      traitement_j[,1][traitement_j[,1]== ATTi[j,1]] <- 1
      traitement_j[,1] <- (traitement_j[,1]-mean(traitement_j[,1])) * 1/(dim(ATTi)[1])
      # Doute sur le fait que listG & agreg_inf_i ait le m�me nb d'observations
      # mat_gamma <- as.data.frame(agreg_inf_i[,1,j])
      # colnames(mat_gamma) <- c("var1")
      # mat_gamma$gamma <- merge(mat_gamma,traitement_j,by.x="var1", by.y="iden_num",all.x)
      # agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ mat_gamma$gamma%*%ATTi[j,nom_outcome]
      agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ t(t(traitement_j[,1]))%*%as.matrix(ATTi[j,nom_outcome])
    }
    # End modif
    agreg_influence_final[,,i]=sum_array(agreg_inf_i)
  }
  # Remettre les identifiants des individus
  for (i in 1:(length(ttt))){
    agreg_influence_final[,1,i]=array_inf_did[,1,1]
  }
  agreg_influence_final
}  



agregat_influence_longDID<-function(tab,array_inf,listG,nom_outcome,tname,first.treat.name,poids){
  nb_an<-unique(tab[,tname])
  nn<-length(nb_an)
  nb_an_long<-c(min(nb_an)-1,nb_an)
  an_cohort<-unique(tab[,first.treat.name])
  NN<-length(an_cohort)
  agreg_influence = array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),nn*NN))
  ATT=data.frame()
  sum_array<-function(array_f){
    len = dim(array_f)[3]
    mat_result=array_f[,,1]
    if (len>1){
      for (k in 2:len){
        mat_result=mat_result+array_f[,,k]
      }
    }
    mat_result
  }
  for (i in 1:NN){
    tab_cohort=tab[(nn*(i-1)+1):(nn*i),]
    ATT_temp=tab_cohort[,1:ncol(tab_cohort)]
    ATT_temp$time_to_treatment=tab_cohort[,tname]-tab_cohort[,first.treat.name]
    ATT_temp=merge(ATT_temp,poids,by.x=first.treat.name, by.y="Var1",all.x)
    ATT=rbind(ATT,ATT_temp)
  }
  ttt=unique(sort(ATT$time_to_treatment))
  agreg_influence_final = array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),length(ttt)))
  for (i in 1:length(ttt)){
    ATTi=ATT[ATT$time_to_treatment==ttt[i],]
    ATTi$Freq=ATTi$Freq/sum(ATTi$Freq)
    agreg_inf_i=array_inf[,,ATTi$compteur]
    if (is.na(dim(agreg_inf_i)[3])==TRUE){
      agreg_inf_i=array(0,dim=c(dim(array_inf)[1],(1+length(nom_outcome)),1))
      agreg_inf_i[,,1]=array_inf[,,ATTi$compteur]
    } else {
      agreg_inf_i=array_inf[,,ATTi$compteur]
      for (j in 1:dim(ATTi)[1]){
        agreg_inf_i[,,j]=ATTi$Freq[j]*agreg_inf_i[,,j]
      }
      # Modif
      agreg_inf_i[,1,]=array_inf[,1,ATTi$compteur]
    }
    # D�but modif
    # Comme les poids sont des fr�quences empiriques, ils sont al�atoires, 
    # il faut ajouter un terme suppl�mentaire dans la fonction d'influence pour ne pas sous-�valuer la variance
    for (j in 1:dim(ATTi)[1]){
      traitement_j <- listG
      traitement_j[,1][traitement_j[,1]!= ATTi[j,1]] <- 0
      traitement_j[,1][traitement_j[,1]== ATTi[j,1]] <- 1
      traitement_j[,1] <- (traitement_j[,1]-mean(traitement_j[,1])) * 1/(dim(ATTi)[1])
      # Doute sur le fait que listG & agreg_inf_i ait le m�me nb d'observations
      # mat_gamma <- as.data.frame(agreg_inf_i[,1,j])
      # colnames(mat_gamma) <- c("var1")
      # mat_gamma$gamma <- merge(mat_gamma,traitement_j,by.x="var1", by.y="iden_num",all.x)
      # agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ mat_gamma$gamma%*%ATTi[j,nom_outcome]
      agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(2:dim(agreg_inf_i)[2]),j]+ t(t(traitement_j[,1]))%*%as.matrix(ATTi[j,nom_outcome])
    }
    # End modif
    agreg_influence_final[,,i]=sum_array(agreg_inf_i)
  }
  # Remettre les identifiants des individus
  for (i in 1:(length(ttt))){
    agreg_influence_final[,1,i]=array_inf[,1,1]
  }
  agreg_influence_final
}  


agregat_influence_GMM<-function(tab,array_influ,agreg_influence,listG,nom_outcome,tname,first.treat.name,poids){
  # fonction pour agreger les fonctions d'influence PHI_influence_ATTgt
  sum_array<-function(array_f){
    len = dim(array_f)[3]
    mat_result=array_f[,,1]
    if (len>1){
      for (k in 2:len){
        mat_result=mat_result+array_f[,,k]
      }
    }
    mat_result
  }
  tab$time_to_treatment<-tab[,tname]-(tab[,first.treat.name]-1)
  tab$compteur<-1:nrow(tab)
  tab=merge(tab,poids,by.x=first.treat.name, by.y="Var1",all.x)
  ttt=unique(sort(tab$time_to_treatment))
  agreg_influence_final <- array(0,dim=c(dim(agreg_influence)[1],length(nom_outcome)+1,length(ttt)+1))
  for (i in 1:length(ttt)){
    ATTgti=tab[tab$time_to_treatment==ttt[i],]
    ATTgti$Freq=ATTgti$Freq/sum(ATTgti$Freq)
    agreg_inf_i=agreg_influence[,,ATTgti$compteur]
    if (is.na(dim(agreg_inf_i)[3])==TRUE){
      agreg_inf_i=array(0,dim=c(dim(agreg_influence)[1],length(nom_outcome),1))
      agreg_inf_i[,,1]=agreg_influence[,,ATTgti$compteur]
    } else {
      agreg_inf_i=agreg_influence[,,ATTgti$compteur]
      for (j in 1:dim(ATTgti)[1]){
        agreg_inf_i[,,j]=ATTgti$Freq[j]*agreg_inf_i[,,j]
      }
    }
    # D�but modif
    # Comme les poids sont des fr�quences empiriques, ils sont al�atoires, 
    # il faut ajouter un terme suppl�mentaire dans la fonction d'influence pour ne pas sous-�valuer la variance
    for (j in 1:dim(ATTgti)[1]){
      traitement_j <- listG
      traitement_j[,1][traitement_j[,1]!= ATTgti[j,1]] <- 0
      traitement_j[,1][traitement_j[,1]== ATTgti[j,1]] <- 1
      traitement_j[,1] <- (traitement_j[,1]-mean(traitement_j[,1])) * 1/(dim(ATTgti)[1])
      agreg_inf_i[,(1:dim(agreg_inf_i)[2]),j] <- agreg_inf_i[,(1:dim(agreg_inf_i)[2]),j]+ t(t(traitement_j[,1]))%*%as.matrix(ATTgti[j,nom_outcome])
    }
    # End modif
    if (ttt[i]< 0 ) {
      agreg_influence_final[,2:dim(agreg_influence_final)[2],i]=sum_array(agreg_inf_i)
    } else {
      agreg_influence_final[,2:dim(agreg_influence_final)[2],i+1]=sum_array(agreg_inf_i)
    }
  }
  # Remettre les identifiants des individus
  for (i in 1:(length(ttt)+1)){
    agreg_influence_final[,1,i]=array_influ[,1,1]
  }
  agreg_influence_final
}



# Fonction qui calcule les pvalues � partir des bootstraps
pvalue<-function(tab,agreg_inf,nom_outcome){
  for (i in 1:length(nom_outcome)){
    for (j in 1:dim(tab)[1])
    aa=tab[j,(i+1)]
    # View(influ[,i,j])
  }
}


agg1<-function(vv,typ){aggregate(vv,by=list(vv[,typ]),FUN=mean)}

agregat<-function(tab,nom_outcome,tname,first.treat.name){
  loo<-unique(tab[,tname])
  nn<-length(loo)
  loo_long<-c(min(loo)-1,loo)
  # print(loo_long)
  laa<-unique(tab[,first.treat.name])
  NN<-length(laa)
  # print(loo);print(laa)
  list1<-list() ; list2<-list() ; an_traitV<-c()
  for(i in 1:NN){
    vv<-1*(loo<laa[i])
    # print(vv) ; print(nn)
    list1[[i]]<-matrix(vv,nrow=nn+1,ncol=nn,byrow=TRUE) 
    jj<-matrix(0,nn+1,nn)
    jj[lower.tri(jj)]<-1
    list2[[i]]<-jj
    vvv<- loo_long-laa[i]
    an_traitV<-c(an_traitV,vvv)
    }
mat1<-bdiag(list1)
mat2<-bdiag(list2)
mama<-as.matrix(tab[,nom_outcome])
mum<-mat2%*%mama-mat1%*%mama
mo<-data.frame(an_trait=an_traitV)
mo[,nom_outcome]<-mum
#mo
agg1(mo,"an_trait")
} 

agg2<-function(vv,typ){aggregate(vv,by=list(vv[,typ]),FUN=sum)}            
agregatP<-function(tab,nom_outcome,tname,first.treat.name,poids){
  loo<-unique(tab[,tname])
  nn<-length(loo)
  loo_long<-c(min(loo)-1,loo)
  laa<-unique(tab[,first.treat.name])
  NN<-length(laa)
  # print(loo);print(laa)
  list1<-list() ; list2<-list() ; an_traitV<-c(); pp<-c()
  for(i in 1:NN){
    vv<-1*(loo<laa[i])
    hh<-rep(poids[i],nn+1)
    # print(vv) ; print(nn)
    list1[[i]]<-matrix(vv,nrow=nn+1,ncol=nn,byrow=TRUE) 
    jj<-matrix(0,nn+1,nn)
    jj[lower.tri(jj)]<-1
    list2[[i]]<-jj
    vvv<- loo_long-laa[i]
    an_traitV<-c(an_traitV,vvv)
    pp<-c(pp,hh)
  }
  mat1<-bdiag(list1)
  mat2<-bdiag(list2)
  mama<-as.matrix(tab[,nom_outcome])
  mum<-mat2%*%mama-mat1%*%mama
  mo<-data.frame(an_trait=an_traitV,poids=pp)
  mo[,nom_outcome]<-mum*pp
  mo<-agg2(mo,"an_trait")
  ma<-mo
  ma[,nom_outcome]<-mo[,nom_outcome]/mo[,"poids"]
  ma
} 


  
# tete<-essai1[[1]]
# agregat(tab=tete,nom_outcome=outcomeF,tname="annee",first.treat.name="annee_FUI") 
# agregatP(tab=tete,nom_outcome=outcomeF,tname="annee",first.treat.name="annee_FUI",c(100,100,200)) 
# 
# tete<-essai1[[2]]
# agregat(tab=tete,nom_outcome=outcomeD,tname="annee",first.treat.name="annee_FUI") 
# agregatP(tab=tete,nom_outcome=outcomeD,tname="annee",first.treat.name="annee_FUI",c(100,100,200)) 
# 
# 
# tete<-essai1[[3]]
# agregat(tab=tete,nom_outcome=outcomeDS,tname="annee",first.treat.name="annee_FUI")   
# agregatP(tab=tete,nom_outcome=outcomeDS,tname="annee",first.treat.name="annee_FUI",c(100,100,200)) 
# 
# tete<-essai2[[2]]
# agregat(tab=tete,nom_outcome=outcomeD,tname="annee",first.treat.name="annee_FUI")
# 
# ujj<-c(4634,6986,6608,6020,5292,4284,4438,4102,1890,3318,2688)
# length(ujj)
# agregatP(tab=tete,nom_outcome=outcomeD,tname="annee",first.treat.name="annee_FUI",poids=ujj)

# 
# 
# 
# lk<-agregat(tab=tete,nom_outcome=outcomeF,tname="annee",first.treat.name="annee_FUI")   
# 
# 
# 
# mean2<-function(hh){
#   weighted.mean(hh, w)
# }
# 
# 
# tete<-essai1[[2]]
# agregat(tab=tete,nom_outcome=outcomeD,tname="annee",first.treat.name="annee_FUI")   
# tete<-essai1[[3]]
# agregat(tab=tete,nom_outcome=outcomeDS,tname="annee",first.treat.name="annee_FUI")   
# 
# 
# #   
# # mlist<-list(matrix(1:8,2,2),matrix(1:20,4,5),matrix(1:5,1,5))
# # lo<-bdiag(mlist)
# # mm<-matrix(1,1,7)
# # mm%*%lo
# # matrix(c(1,1,0,0,0),nrow=10,ncol=5,byrow=TRUE)  
# # jj<-matrix(0,length(looo)+1,length(looo))
# # jj[lower.tri(jj)]<-1
# # jj
# # tete<-essai1[[1]]
# # looo<-unique(tete$annee)
# # laaa<-unique(tete$annee_FUI)
# # looo
# # laaa
# # ## annee-annee_trai : donne la bonnne position
# # ## des matrcies carr�e de longuer le nombre total d'annn�e, c'est � dire annee +1



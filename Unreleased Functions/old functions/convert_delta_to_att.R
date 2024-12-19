convert_delta_to_att <- function(dp){

# We extract the necessary variables from the DIDparams object
data <- dp$data
yname = dp$yname
tname = dp$tname
gname = dp$gname
idname = dp$idname
deb = min(dp$tlist)
fin = max(dp$tlist)
debT = min(dp$glist)
finT = max(dp$glist)
tlist <- dp$tlist
glist <- dp$glist
weightsname <- dp$weightsname
xformla <- dp$xformla

result <- dp$delta.att.influ #delta att and influence matrix from compute_delta_att.

tlen <- length(tlist) #how many dates
glen <- length(glist) #how many treatment cohorts


# 3. Post-estimation aggregation step (to be cleaned)

# #useless?
# # Poids pour aggr�ger les effets des diff�rentes cohortes de traitement 
list_id_w=result[[3]]
list_id_w=list_id_w[list_id_w[,gname]>0,]
list_id_w<-unique(list_id_w)
weightmat<-as.data.frame(table(list_id_w[,gname]))

# 3.1. Aggregate Delta ATT into ATT
# resultat<-agregatChris(tab=result[[1]],yname =yname,tname=tname,gname=gname,poids=weightmat)

tab <- result[[1]]
tlist_long<-c(min(tlist)-1,tlist)
ATT=data.frame()

# Loop over each cohort in glist
for (i in glist) {
  # Filter rows for the current cohort
  tab_cohort <- tab[tab[[gname]] == i, ]
  
  # Separate rows into "after" and "before" treatment groups
  after <- tab_cohort[tab_cohort[[tname]] >= i, ]
  befor <- tab_cohort[tab_cohort[[tname]] < i, ]
  befor <- befor[order(-befor[[tname]]), ]  # Sort 'befor' by descending tname
  
  # Calculate cumulative sums
  after[[yname]] <- cumsum(after[[yname]])
  befor[[yname]] <- cumsum(-befor[[yname]])
  
  # Add time_to_treatment
  after$time_to_treatment <- after[[tname]] - after[[gname]] + 1
  befor$time_to_treatment <- befor[[tname]] - befor[[gname]]
  
  # Combine and rename columns in one step
  ATT_temp <- rbind(
    befor[, c(names(befor)[2:5], yname, "time_to_treatment")],
    after[, c(names(after)[2:5], yname, "time_to_treatment")]
  )
  
  # Append to ATT
  ATT <- rbind(ATT, ATT_temp)
}

# 3.2. Influence function for ATT's
tab <- result[[1]]
array_inf <- result[[2]]
listG <- result[[3]]
tlist_long<-c(min(tlist)-1,tlist)
agreg_influence = array(0,dim=c(dim(array_inf)[1],(1+length(yname)),(tlen-1)*glen))

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

for (i in glist){
  # Filter rows for the current cohort
  tab_cohort <- tab[tab[[gname]] == i, ]
  
  after <- tab_cohort[tab_cohort[[tname]] >= i, ]
  befor <- tab_cohort[tab_cohort[[tname]] < i, ]
  
  after$time_to_treatment=after[,tname]-after[,gname] + 1
  befor$time_to_treatment=befor[,tname]-befor[,gname]
  
  influence_cohort_after=array_inf[,,after$attgt_id]
  influence_cohort_befor=array_inf[,,befor$attgt_id]
  
  if (is.na(dim(influence_cohort_after)[3])==TRUE){
    influence_cohort_after=array_inf[,,after$attgt_id]
  } else if (dim(influence_cohort_after)[3]==0) {
    influence_cohort_after=array_inf[,,after$attgt_id]
  } else {
    influence_cohort_after=cumsum_after(influence_cohort_after)
    influence_cohort_after[,1,]=array_inf[,1,after$attgt_id] #replace for first dim which is not phi, but used to keep track, see mat_influence
  }
  if (is.na(dim(influence_cohort_befor)[3])==TRUE){
    influence_cohort_befor=array_inf[,,befor$attgt_id]
  } else if (dim(influence_cohort_befor)[3]==0) {
    influence_cohort_befor=array_inf[,,befor$attgt_id]
  } else {
    influence_cohort_befor=cumsum_befor(influence_cohort_befor)
    influence_cohort_befor[,1,]=array_inf[,1,befor$attgt_id] #replace for first dim which is not phi, but used to keep track, see mat_influence
  }
  agreg_influence[,,after$attgt_id]<-influence_cohort_after
  agreg_influence[,,befor$attgt_id]<-influence_cohort_befor
}

# ATT: contains all ATT(g,t)'s, sorted as desired
# agreg_influence: influence matrix is dim (nbindiv x 2 x nb of ATT(g,t)), influence values are in agreg_influence[,2,]
dp$att.influ <-list(ATT,agreg_influence)
dp
}
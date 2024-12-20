#' @title convert_delta_to_att
#'
#' @description Function to process arguments passed to the main methods in the
#'  `cdid` package to compute ATT.
#' 
#' @return a [`DIDparams`] object
#'
#' @export

gmm_convert_delta_to_att <- function(dp){



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
chained <- dp$chained #TRUE for chained, FALSE for GMM.


result <- dp$delta.att.influ #delta att and influence matrix from compute_delta_att.
tlen <- length(tlist) #how many dates
glen <- length(glist) #how many treatment cohorts

  # 3. Post-estimation aggregation step
  # 3.1. Aggregate Delta ATT into ATT
  result[[1]][] <- lapply(result[[1]], function(x) as.numeric(as.character(x)))
  result[[1]][is.na(result[[1]])] <- 0
  
  # Covariance matrix of delta ATT estimator 
  omega_deltaATT <- array(0,dim=c(length(yname),dim(result[[1]])[1],dim(result[[1]])[1]))
  omega_deltaATT[1,,] <- (1/dim(result[[2]])[1]) * t(result[[2]][,2,]) %*% result[[2]][,2,]   
  
  # Remove columns from mat_W 
  mat_W <- result[[4]]
  remov_col=c()
  for (f in 1:glen) { 
    for (t in 1:tlen) {
      if (glist[f]-1 == tlist[t]){remov_col<-c(remov_col,(f-1)*tlen+t)}
    }
  }
  mat_W <- subset( mat_W, select = -c(remov_col ) )
  mat_W <- as.matrix(mat_W)
  
  # Covariance Matrix of ATTgt
  Sigma_ATTgt <- array(0,dim=c(length(yname),dim(mat_W)[2],dim(mat_W)[2]))
  Sigma_ATTgt[1,,] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[1,,])%*%mat_W)
  
  # Optimal estimator of ATTgt
  delta_ATT <- as.matrix(result[[1]][,yname])
  ATT <- matrix(NA,ncol(mat_W),length(yname))
  colnames(ATT) <- yname
  rownames(ATT) <- colnames(mat_W)
  

  if (chained==FALSE){
  ATT[,1] <- MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[1,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[1,,])%*%delta_ATT[,1]
  }else{
  ATT[,1] <- MASS::ginv(t(mat_W)%*%diag(dim(omega_deltaATT[1,,])[1])%*%mat_W)%*%t(mat_W)%*%diag(dim(omega_deltaATT[1,,])[1])%*%delta_ATT[,1] #chained
  }


  # Improve formatting
  ATT<-as.data.frame(ATT)
  ncol_ATT<-ncol(ATT)
  for (f in 1:glen) { 
    for (t in 1:(tlen-1)) {
      
      ATT[(f-1)*(tlen-1)+t,ncol_ATT+1]<-glist[f]
      if(glist[f]-1>tlist[t]){ATT[(f-1)*(tlen-1)+t,ncol_ATT+2]<-tlist[t]}else{ATT[(f-1)*(tlen-1)+t,ncol_ATT+2]<-tlist[t]+1}
    }
  }
  colnames(ATT)[ncol_ATT+1]<-gname
  colnames(ATT)[ncol_ATT+2]<-tname
  
  # Influence function of ATT : PHI
  agreg_influence=array(0,dim=c(dim(result[[2]])[1],length(yname),dim(ATT)[1])) 
  
  if (chained==FALSE){
  agreg_influence[,1,] <- t(MASS::ginv(t(mat_W)%*%MASS::ginv(omega_deltaATT[1,,])%*%mat_W)%*%t(mat_W)%*%MASS::ginv(omega_deltaATT[1,,])%*%t(result[[2]][,2,]))
  }else{
  agreg_influence[,1,] <- t(MASS::ginv(t(mat_W)%*%diag(dim(omega_deltaATT[1,,])[1])%*%mat_W)%*%t(mat_W)%*%diag(dim(omega_deltaATT[1,,])[1])%*%t(result[[2]][,2,])) #chained
  }

  # We add the result to dp.
  # ATT: contains all ATT(g,t)'s, sorted as desired
  # agreg_influence: influence matrix is dim (nbindiv x 2 x nb of ATT(g,t)), influence values are in agreg_influence[,2,]
  dp$att.influ <-list(ATT,agreg_influence)
  

  dp
}

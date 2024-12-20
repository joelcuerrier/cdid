#' @title GMM_compute_delta_att
#'
#' @description Function to process arguments passed to the main methods in the
#'  `cdid` package to compute Delta ATT as well as conducting some tests to make sure
#'  data is in proper format. 
#' 
#' @param call Function call to gg
#'
#' @return a [`DIDparams`] object
#'
#' @export
#' 
#' 
gmm_compute_delta_att <-function(dp) {

  # all parameters are taken from pre-process dp object
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
  control_group <- dp$control_group
  chained <- dp$chained 


  #gmm.R calls GMM_estimPeriod_Boot dans fonctions_estimation_Boot.R
  ########################################
  ## Part 1. Prelim checks (previously done with mp.spatt.GMM)
  #Make sure we have a dataframe
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  
  #enforce start/end date #useless if data well constructed
  data<-data[(data[[tname]]>=deb & data[[tname]] <= fin) & (data[,gname]>=debT | data[,gname]==0),] 
  #treated after enddate are considered as controls, probably useless in most cases
  data[data[,gname]>finT,gname]<-0
  
  #Condition on sample size
  gsize <- aggregate(data[, gname], by = list(data[,gname]), function(x) length(x)/length(tlist))
  gsize <- subset(gsize, x < (length(rhs.vars(xformla)) + 5)) #check groups not satisfying condition
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("Some groups are probably too small: \n", gpaste))
  }
  
  tlen <- length(tlist) #how many dates
  glen <- length(glist) #how many treatment cohorts

    
  ########################################
  ## Part 2. Compute Delta ATT (previously done with compute.mp.spatt.GMM)
  
  # 2.1 Prelim
  # Delta_k ATT(g,t) matrix init  
  nbobjects=glen*(tlen - 1)*tlen/2 
  fatt <- data.frame(
    attgt_id = seq_len(nbobjects),   # Generate sequence for attgt_id to track object
    nobsG = seq_len(nbobjects),      # Same sequence for nobsG
    nobsC = seq_len(nbobjects)      # Same sequence for nobsC
  )
  
  # Initialize matrix W containing the weights to go from ATT(g,t) to \Delta_k ATT(g,t). The nb of possible ATT(g,t) that can be identified is glen*(tlen-1)
  # but we fill in glen*tlen ATT(g,t), the extra columns are deleted later on
  mat_W <-as.data.frame(matrix(0,nbobjects,glen*tlen)) #168x42
  for (f in 1:glen) { 
    for (t in 1:tlen) {
      colnames(mat_W)[(f-1)*tlen+t]<-paste0("ATT(",glist[f],",",tlist[t],")")
    }
  }
  
  # Influence matrix
  indiv<-unique(data[,c(gname,idname)]) #unique list of indiv with their cohort
  mat_influence = array(0,dim=c(nrow(indiv),(1+length(yname)),nbobjects)) #influence matrix init
  counter <- 1
  
  # 2.2 Loop over cohorts and dates
  for (f in 1:glen) { #cohort loop (outer loop)
    for (t in 1:(tlen - 1)) { #date loop (inner loop)
      tbis = tlen - t #possible distances 
      
      for (k in 1:tbis) {
        
        # 2.2.1 Propensity scores
        # Here we estimate propensity scores using data on all indivs in current treatment cohort and control cohort.
        # This can vary if we use the not-yet-treated, so it's more convenient to keep it in the inner loop
        if (control_group == "nevertreated" & k==1) {
          disdat <- data[!duplicated(data[["id"]]), ]
          disdat <- data[data[[gname]] == 0 | data[[gname]] == glist[f], ]
          disdat <- droplevels(disdat) 
          
          if (control_group == "nevertreated") {
            disdat$C <- 1 * (disdat[, gname] == 0) 
          } else { # + notyettreated
            disdat$C <- 1 * (disdat[, gname] == 0) + 1 * (disdat[, gname] > glist[f] & disdat[, gname] > tlist[t+k])
          }
          
          disdat$G <- 1 * (disdat[, gname] == glist[f])
          
          if (is.null(xformla)) {
            xformla <- ~1}
          pformla <- xformla
          
          # Discard observations for which scores cannot be computed
          LesX<-BMisc::rhs.vars(pformla)
          bbb<-length(LesX)
          for (jj in 1:bbb){
            disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,] }
          
          # Propensity score estimation
          pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla)) #G~X
          base_score<-subset(disdat, C + G == 1) #keep only control + current cohort G
          pscore.reg <- glm(pformla, family = binomial(link = "logit"), data = base_score)
          thet <- coef(pscore.reg) #can be used to check coeffs
          
          if (any(is.na(thet))) {
            warning(paste0("Problems estimating propensity score... likely perfectly predicting treatment for group: ",
                           glist[f], " at time period: ", tlist[t +k]))}
        }
        if (control_group == "notyettreated") {
          disdat <- data[!duplicated(data[["id"]]), ]
          disdat <- data[data[[gname]] == 0 | data[[gname]] == glist[f], ]
          disdat <- droplevels(disdat) 
          
          if (control_group == "nevertreated") {
            disdat$C <- 1 * (disdat[, gname] == 0) 
          } else { # + notyettreated
            disdat$C <- 1 * (disdat[, gname] == 0) + 1 * (disdat[, gname] > glist[f] & disdat[, gname] > tlist[t+k])
          }
          
          disdat$G <- 1 * (disdat[, gname] == glist[f])
          
          if (is.null(xformla)) {
            xformla <- ~1}
          pformla <- xformla
          
          # Discard observations for which scores cannot be computed
          LesX<-BMisc::rhs.vars(pformla)
          bbb<-length(LesX)
          for (jj in 1:bbb){
            disdat<-disdat[is.na(disdat[,LesX[jj]])==FALSE,] }
          
          # Propensity score estimation
          pformla <- BMisc::toformula("G", BMisc::rhs.vars(pformla)) #G~X
          base_score<-subset(disdat, C + G == 1) #keep only control + current cohort G
          pscore.reg <- glm(pformla, family = binomial(link = "logit"), data = base_score)
          thet <- coef(pscore.reg) #can be used to check coeffs
          
          if (any(is.na(thet))) {
            warning(paste0("Problems estimating propensity score... likely perfectly predicting treatment for group: ",
                           glist[f], " at time period: ", tlist[t +k]))}
        }
        ### Fill W
        # Calculate indices for positive and negative weights
        att_pos_index = (f - 1) * tlen + (t + k)
        att_neg_index = (f - 1) * tlen + t
        
        # Update the matrix with positive and negative weights
        mat_W[counter, att_pos_index] = 1  # Assign positive weight
        mat_W[counter, att_neg_index] = -1  # Assign negative weight
        
        
        # 2.2.2 Prep data for estimation !! Dec 14: it comes from this part.
        # Make balanced panel for (tlist[t],tlist[t+k]). https://bcallaway11.github.io/BMisc/reference/makeBalancedPanel.html
        disdat <- data[(data[, tname] == tlist[t + k] | data[, tname] == tlist[t]), ]
        disdat <- makeBalancedPanel(disdat, idname, tname) #This function drops observations from data.frame that are not part of balanced panel data set.
        
        if (dim(disdat)[1]>=1) { #if at least one indiv observed twice
          if (control_group == "nevertreated") {
            disdat$C <- 1 * (disdat[, gname] == 0) 
          } else { # + notyettreated
            disdat$C <- 1 * (disdat[, gname] == 0) + 1 * (disdat[, gname] > glist[f] & disdat[, gname] > tlist[t+k])
          }
          disdat$G <- 1 * (disdat[, gname] == glist[f]) 
          
          # Calculate Delta Y: here we only compute one-period differences: disdat$delta_y=yt+k-yt
          # (Formerly done with panelDiffV.R)
          t1 <- tlist[t]
          t2 <- tlist[t+k]
          
          # Select and order weights for year t1
          ppp1 <- disdat[disdat[[tname]] == t1, c(idname, weightsname)]
          ppp1 <- ppp1[order(ppp1[[idname]]), ]
          ppp2 <- disdat[disdat[[tname]] == t2, c(idname, weightsname)]
          ppp2 <- ppp2[order(ppp2[[idname]]), ]
          
          # Filter data for year t1 and assign weights
          retdat <- disdat[disdat[[tname]] == t1, ]
          retdat <- retdat[order(retdat[[idname]]), ]
          retdat[[weightsname]] <- ppp1[[weightsname]]*ppp2[[weightsname]]
          
          # Select and compute differences in y for t1 and t2
          dt1 <- disdat[disdat[[tname]] == t1, c(idname, yname)]
          dt2 <- disdat[disdat[[tname]] == t2, c(idname, yname)]
          dt1 <- dt1[order(dt1[[idname]]), ]
          dt2 <- dt2[order(dt2[[idname]]), ]
          
          retdat[,paste0("D",yname)] <- dt2[,yname] - dt1[,yname] 
          disdat <- retdat
          
          # 2.2.3 Predict propensity scores
          pscore <- predict(pscore.reg, newdata = disdat,type = "response")  
          pscore[is.na(pscore)]<-0  
          disdat$pscore<-pscore
          disdat<-disdat[disdat$pscore<1 & disdat$pscore>0,] ### discard obs with scores = 1 or 0 to prevent NA/Inf
          disdat[,c("pscore")]=list(NULL)
          pscore <- predict(pscore.reg, newdata = disdat,type = "response")  # recompute scores after discarding so it has the correct length
          
          # 2.2.4 Compute DeltaATT
          dy <- disdat[,paste0("D",yname)] #delta y
          
          #re-do because overwriting
          if (control_group == "nevertreated") {
            disdat$C <- 1 * (disdat[, gname] == 0) 
          } else { # + notyettreated
            disdat$C <- 1 * (disdat[, gname] == 0) + 1 * (disdat[, gname] > glist[f] & disdat[, gname] > tlist[t+k])
          }
          disdat$G <- 1 * (disdat[, gname] == glist[f])
          
          # Calculate Wg and Wc
          GC<- cbind(disdat$G,disdat$C)
          devant<-(GC[,1]-GC[,2]) #= (G-C)
          devant2<-GC[,1]+(pscore/(1 - pscore))*GC[,2] #= G + C*pscore/(1-pscore). 
          
          # dommk has 2 elements: [sum G*S, sum C*S*(pscore/(1 - pscore))]
          dommk<-t(as.matrix(GC))%*%as.matrix((devant2*disdat[,weightsname]))
          
          # Denom are the denominators of Wg and Wc
          Denom <- 1/(dommk) 
          
          # Wg-Wc
          pp<-devant*(GC%*%Denom)*(devant2*disdat[,weightsname]) 
          
          # Delta ATT
          att<- t(as.numeric(pp))%*%as.numeric(dy)
          
          # A useful step for the influence matrix
          pp_noDenom<- devant*(GC%*%matrix(c(1, 1), ncol = 1))*(devant2*disdat[,weightsname])   
          
          
          # 2.2.5 Influence matrix 
          n <- length(unique(data$id))# nrow(disdat) 
          ## nb units (in total in the sample, not just the observed cohort!)
          #this was an error in our previous version, where we also counted unobserved indiv in the simulations
          
          # Calculate E(dy)
          esperance_dy <- GC %*% t(GC) %*% as.matrix(pp * dy)
          esperance_dy <- esperance_dy * (disdat$G - disdat$C)
          
          # Calculate psi
          psi <- n* pp *(dy-esperance_dy) #dim nbunitsx1
          colnames(psi) <- yname
          
          # Generate model matrix
          x <- model.matrix(xformla, data = disdat)
          
          # Precompute weights
          G_i <- disdat$G * disdat[[weightsname]]
          C_i <- disdat$C * disdat[[weightsname]]
          wc1_i <- -C_i * pp_noDenom
          mean_wc1_i <- mean(wc1_i)
          wc_i <- wc1_i / mean_wc1_i
          
          # Function gg(x, thet) returns 1 / ((1 + exp(x %*% thet))^2)
          gg_val <- gg(x, thet)
          
          # Compute M_i
          adjusted_dy_i <- dy - mean(wc_i * dy)
          M_i <- colMeans(((C_i / (1 - pscore))^2 * gg_val * adjusted_dy_i) * x) / mean_wc1_i
          #M_i <- as.matrix(apply(as.matrix((C_i/(1 - pscore))^2 * gg(x, thet) * (dy - mean(wc_i * dy)) * x), 2,mean)/mean(wc1_i))
          
          # Compute A1_i
          A1_i <- t((G_i + C_i) * gg_val^2 / (pscore * (1 - pscore)) * x) %*% x / n
          
          # Compute A2_i
          A2_i <- ((G_i + C_i) * (G_i - pscore) * gg_val / (pscore * (1 - pscore))) * x
          
          # Calculate A_i and correction
          A_i <- A2_i %*% MASS::ginv(A1_i)
          correction_i <- A_i %*% M_i * disdat$C
          
          # Update psi with correction
          psi[, yname] <- psi[, yname] - correction_i
          
          # Convert psi to a data frame and set column names
          psi <- as.data.frame(psi)
          colnames(psi) <- yname
          
          # Add unit id to psi
          psi[[idname]] <- disdat[[idname]]
          
          # Merge with unit ids, filling missing values with 0
          psi <- merge(indiv[idname], psi, by = idname, all.x = TRUE)
          
          # Replace NA values in all yname columns with 0
          psi[yname] <- lapply(psi[yname], function(col) ifelse(is.na(col), 0, col))
          
          
          # 2.2.6 Update `fatt` values
          nG<-round(mean(colSums(as.matrix(disdat$G  *(disdat[,weightsname]>0)))))
          nC<-round(mean(colSums(as.matrix(disdat$C  *(disdat[,weightsname]>0)))))
          
          fatt[counter, c(gname, tname, yname, "nobsG", "nobsC")] <- list(
            glist[f],
            tlist[t + k],
            att,
            nG,
            nC
          )
          
          for (i in 1:dim(psi)[2]){
            mat_influence[,i,counter]=psi[,i]   #counter keeps track of (g,t) pair
          }
        } else { #if sample size for this delta_k att(g,t)
          mat_W[counter, att_pos_index] = 0  # Assign 0 weight
          mat_W[counter, att_neg_index] = 0  # Assign 0 weight
        }
        
        counter <- counter + 1
        
      }## distance loop end 
    }## date loop end
  } ## cohort loop end
  
  # Results:
  result <- list(fatt,mat_influence,indiv,mat_W)  
  # We add the result to dp.
dp$delta.att.influ <- result
dp$chained <-chained
dp
}
#   colnames(mat_W[1,colSums(abs(mat_W),na.rm=TRUE)>0])

  
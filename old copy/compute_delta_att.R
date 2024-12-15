
# compute_delta_att <-function(
#                     data,
#                     yname,
#                     tname,
#                     gname,
#                     idname,
#                     deb,
#                     fin,
#                     debT,
#                     finT,
#                     tlist,
#                     glist,
#                     weightsname,
#                     xformla) {


compute_delta_att <-function(dp) {


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

########################################
## Part 1. Prelim checks (previously done with chained.mp.spatt.Boot)
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
## Part 2. Compute Delta ATT (previously done with chained.compute.mp.spatt.Boot)


# 2.1 Prelim
# ATT(g,t) matrix init  
nbobjects=glen*(tlen - 1)
fatt <- data.frame(
  attgt_id = seq_len(nbobjects),   # Generate sequence for attgt_id
  nobsG = seq_len(nbobjects),      # Same sequence for nobsG
  nobsC = seq_len(nbobjects)      # Same sequence for nobsC
)

# Influence matrix
indiv<-unique(data[,c(gname,idname)]) #unique list of indiv with their cohort
mat_influence = array(0,dim=c(nrow(indiv),(1+length(yname)),nbobjects)) #influence matrix init
counter <- 1

# 2.2 Loop over cohorts and dates
for (f in 1:glen) { #cohort loop (outer loop)
  for (t in 1:(tlen - 1)) { #date loop (inner loop)
    # 2.2.1 Propensity scores
    # Here we estimate propensity scores using data on all indivs in current treatment cohort and control cohort.
    # This can vary if we use the not-yet-treated, so it's more convenient to keep it in the inner loop (not yet implemented though)
    disdat <- data[!duplicated(data[["id"]]), ]
    disdat <- data[data[[gname]] == 0 | data[[gname]] == glist[f], ]
    disdat <- droplevels(disdat) 
    disdat$C <- 1 * (disdat[, gname] == 0) 
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
                     glist[f], " at time period: ", tlist[t +1]))}
    
    # 2.2.2 Prep data for estimation
    # Make balanced panel for (tlist[t],tlist[t+1]). https://bcallaway11.github.io/BMisc/reference/makeBalancedPanel.html
    disdat <- data[(data[, tname] == tlist[t + 1] | data[, tname] == tlist[t]), ]
    disdat <- makeBalancedPanel(disdat, idname, tname) #This function drops observations from data.frame that are not part of balanced panel data set.
    
    disdat$C <- 1 * (disdat[, gname] == 0) 
    disdat$G <- 1 * (disdat[, gname] == glist[f]) 
    
    # Calculate Delta Y: here we only compute one-period differences: disdat$delta_y=yt+1-yt
    # (Formerly done with panelDiffV.R)
    t1 <- tlist[t]
    t2 <- tlist[t+1]
    
    # Select and order weights for year t1
    ppp <- disdat[disdat[[tname]] == t1, c(idname, weightsname)]
    ppp <- ppp[order(ppp[[idname]]), ]
    
    # Filter data for year t1 and assign weights
    retdat <- disdat[disdat[[tname]] == t1, ]
    retdat <- retdat[order(retdat[[idname]]), ]
    retdat[[weightsname]] <- ppp[[weightsname]]
    
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
    disdat$C <- 1 * (disdat[, gname] == 0) 
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
    n <- nrow(disdat) ## nb units
    
    # Calculate E(dy)
    esperance_dy <- GC %*% t(GC) %*% as.matrix(disdat[,weightsname] * dy)
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
      tlist[t + 1],
      att,
      nG,
      nC
    )

    for (i in 1:dim(psi)[2]){
      mat_influence[,i,counter]=psi[,i]   #counter keeps track of (g,t) pair
    } 
    
    counter <- counter + 1

  } ### date loop end

} ## cohort loop end

# Results:
result <- list(fatt,mat_influence,indiv)

# We add the result to dp.
dp$delta.att.influ <- result

dp
}
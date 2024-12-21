#' @title convert_result
#'
#' @description Result must be converted to be used in the aggte function.  . 
#'
#' @return a [`DIDparams`] object
#'
#' @export
    
    gmm_convert_result <- function(dp,type){

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
      cband <- dp$cband
      
      clustervars <- dp$clustervars
      bstrap <- dp$bstrap
      alp <- dp$alp
      pl <- dp$pl
      cores <- dp$cores
      
      ##
      results <- dp$att.influ
     
      attgt_list=results[[1]]
      attgt_list=as.data.frame(attgt_list)
      
      attgt.list <- lapply(1:nrow(attgt_list), function(i) {
        list(
          att = as.numeric(attgt_list[i, type]), #2-step==1
          group = as.numeric(attgt_list[i, gname]),
          year = as.numeric(attgt_list[i, tname]),
          post = ifelse(attgt_list[i, gname] >= attgt_list[i, tname], 1, 0))})
      
      inffunc=results[[2]][,type,] #2-step==1
      inffunc=as.data.frame(inffunc)
      inffunc=as.matrix(inffunc)
      
      attgt.results=process_attgt_gmm(attgt.list)
      #print(attgt.results)
      group=attgt.results$group
      att=attgt.results$att
      tt=attgt.results$tt
      
      # analytical standard errors
      # estimate variance
      # this is analogous to cluster robust standard errors that
      # are clustered at the unit level
      
      # note to self: this def. won't work with unbalanced panel,
      # same with clustered standard errors
      # but it is always ignored b/c bstrap has to be true in that case
      
      n <- length(unique(data[,idname]))
      V <- Matrix::t(inffunc)%*%inffunc/n
      se <- sqrt(Matrix::diag(V)/n)  
      se[se <= sqrt(.Machine$double.eps)*10] <- NA
      
      # #muted for now because no clustvars argument in the current function  
      # # # if clustering along another dimension...we require using the
      # # # bootstrap (in principle, could come up with analytical standard
      # # # errors here though)
      # if ( (length(clustervars) > 0) & !bstrap) {
      #   warning("clustering the standard errors requires using the bootstrap, resulting standard errors are NOT accounting for clustering")
      # }
      
      # Identify entries of main diagonal V that are zero or NA
      zero_na_sd_entry <- unique(which(is.na(se))) 
      
      # bootstrap variance matrix
      if (bstrap) {
        bout <- mboot(inffunc, DIDparams=dp, pl=pl, cores=cores)
        bres <- bout$bres
        if(length(zero_na_sd_entry)>0) {
          se[-zero_na_sd_entry] <- bout$se[-zero_na_sd_entry]
        } else {
          se <- bout$se
        }
      }
      # Zero standard error replaced by NA
      se[se <= sqrt(.Machine$double.eps)*10] <- NA
      
      #-----------------------------------------------------------------------------
      # compute Wald pre-test
      #-----------------------------------------------------------------------------
      
      # select which periods are pre-treatment
      
      pre <- which(group > tt)
      
      # Drop group-periods that have variance equal to zero (singularity problems)
      if(length(zero_na_sd_entry)>0){
        pre <- pre[!(pre %in% zero_na_sd_entry)]
      }
      
      
      # pseudo-atts in pre-treatment periods
      preatt <- as.matrix(att[pre])
      
      # covariance matrix of pre-treatment atts
      preV <- as.matrix(V[pre,pre])
      
      # check if there are actually any pre-treatment periods
      if (length(preV) == 0) {
        message("No pre-treatment periods to test")
        W  <- NULL
        Wpval <- NULL
      } else if(sum(is.na(preV))) {
        warning("Not returning pre-test Wald statistic due to NA pre-treatment values")
        W <- NULL
        Wpval <- NULL
      } else if (rcond(preV) <= .Machine$double.eps) {
        # singluar covariance matrix for pre-treatment periods
        warning("Not returning pre-test Wald statistic due to singular covariance matrix")
        W <- NULL
        Wpval <- NULL
      } else {
        
        W <- n*t(preatt)%*%solve(preV)%*%preatt
        q <- length(pre) # number of restrictions
        Wpval <- round(1-pchisq(W,q),10)
      }
      
      
      #-----------------------------------------------------------------------------
      # compute confidence intervals / bands
      #-----------------------------------------------------------------------------
      
      # critical value from N(0,1), for pointwise
      cval <- qnorm(1-alp/2)
      
      # in order to get uniform confidence bands
      # HAVE to use the bootstrap
      if (bstrap){
        if (cband) {
          # for uniform confidence band
          # compute new critical value
          # see paper for details
          bSigma <- apply(bres, 2,
                          function(b) (quantile(b, .75, type=1, na.rm = T) -
                                         quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
          
          bSigma[bSigma <= sqrt(.Machine$double.eps)*10] <- NA
          
          # sup-t confidence band
          bT <- apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = TRUE))
          cval <- quantile(bT, 1-alp, type=1, na.rm = T)
          if(cval >= 7){
            warning("Simultaneous critical value is arguably `too large' to be realible. This usually happens when number of observations per group is small and/or there is no much variation in outcomes.")
          }
        }
      }
      
      
      #THIS IS GIVING MY ATT(g,t) with CI.
      results=MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp,debT=debT)
      
      #
      
      unique_id <- (1:length(unique(data[,idname]))) #unique(results$DIDparams$data$id)
      unique_gt <- paste0("(",results$group,",",results$t,")")
      
      inffunc_data <- expand.grid(id = unique_id,
                                  gt = unique_gt)
      
      inffunc_data$x <- as.vector(results$inffunc)
      
      #results$inffunc[inffunc_data[20891,"id"]==unique_id, inffunc_data[20891,"gt"]==paste0("(",results$group,",",results$t,")")]
      #inffunc_data[20891,]
      
      sparse_matrix <- Matrix::sparseMatrix(
        i = inffunc_data$id,
        j = inffunc_data$gt,
        x = inffunc_data$x,
        dimnames = list(NULL, NULL)
      )
      
      results$inffunc=sparse_matrix
      # results$inffunc=inffunc
      results$t[results$group>results$t] <- results$t[results$group>results$t]+1
      print(results)


    return(results)
    }  
            
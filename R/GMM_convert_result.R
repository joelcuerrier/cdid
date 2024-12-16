
    
    gmm_convert_result <- function(dp){


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

    resultat <- dp$att.influ
    
    attgt_list=resultat[[1]]
    attgt_list=as.data.frame(attgt_list)
    

    attgt.list <- lapply(1:nrow(attgt_list), function(i) { #168,6
    list(
    att = as.numeric(attgt_list[i, yname]),#42
    group = as.numeric(attgt_list[i, gname]), #42
    
    year = as.numeric(attgt_list[i, tname]),
    post = ifelse(attgt_list[i, gname] > attgt_list[i, tname], 1, 0))})
    
    #Il y a un ajustement des résultats GMM pour les rendre compatibles avec les résultats de cDiD.
    #L'ajustement consiste à ajouter la dimension tbis. On a donc un ATT associé à chaque g,t et tbis. Sans quoi on a plusieurs ATT pour un seul g,t. 168 au lieu de 42.

    inffunc=resultat[[2]] #6597x42
    inffunc=as.data.frame(inffunc)
    inffunc=as.matrix(inffunc)

    # View(attgt.list) #ok

    attgt.results=process_attgt_gmm(attgt.list)
    # View(attgt.results) #ok
    
    #This patch is used to aggregate the results by g,t. It allows to have 8 groups instead of 2 JC 30-06-2024
    # attgt.results<- aggregate(att ~ group + tt, data = as.data.frame(attgt.results), FUN = mean)
      
    
    group=attgt.results$group
    att=attgt.results$att
    tt=attgt.results$tt

    
    #Cette partie du code provient de la librarie cdid. L'output de cdid est adapté. 
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

    #En utilisant process_attgt_gmm il y a 168 ATT. Il y a 42 groupes et 4 périodes. J'ajoute une patch pour retirer les valeurs supérieures à 42.
    pre <- pre[pre<= length(unique(group))*length(unique(tt))]

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

  
    

    results=MP(group=group, t=tt, att=att, V_analytical=V, se=se, c=cval, inffunc=inffunc, n=n, W=W, Wpval=Wpval, alp = alp, DIDparams=dp,debT=debT)
    
    
    unique_id <- unique(results$DIDparams$data$id)
    unique_group <- unique(results$group)
    unique_t <- unique(results$t)

    inffunc_data <- expand.grid(id = unique_id,
                                group = unique_group,
                                t = unique_t)

    inffunc_data$gt <- as.integer(factor(paste0(inffunc_data$group, inffunc_data$t), levels = unique(paste0(inffunc_data$group, inffunc_data$t))))
    
    # inffunc_data$x=results$inffunc
    inffunc_data$x <- if (length(results$inffunc) >= nrow(inffunc_data)) results$inffunc[1:nrow(inffunc_data)] else c(results$inffunc, rep(NA, nrow(inffunc_data) - length(results$inffunc)))
    
    sparse_matrix <- sparseMatrix(
    i = inffunc_data$id,
    j = inffunc_data$gt,
    x = inffunc_data$x,
    dimnames = list(NULL, NULL)
    )
    print(str(results))
    results$inffunc=sparse_matrix
    # results$inffunc=inffunc
    return(results)
    }  
            
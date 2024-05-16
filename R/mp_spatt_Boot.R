

##### ---------------------------------------------------------------------------------------------------------------------------
##### CHAINED-----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------

chained.mp.spatt.Boot <- function (nom_outcome,nom_traitement, xformla = NULL,propensityformla,data, tname, aggte = TRUE,
                         w = NULL, idname = NULL, first.treat.name,
                         alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                         biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                         seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                         ,selection , ponderation,weight_assumption,debT,finT)  ## de nombreux arguments ne sont plus utilis�s


{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  
  data2<-data[,c(first.treat.name,tname)]


  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))] #liste des traitements
  

#finT2 c'est  dernière génératio nde traitement considérée finT
  data2<-data[data[,first.treat.name]<=finT,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  
  
  flist <- flist[flist > 0]
  
  
  if (is.null(xformla)) {
    xformla <- ~1
  }

  ### tests sur le nombre de traitees par generation
  #group size (proportion % par rapport au nombre de périodes t)
  #on regarde si les groupes sont assez grands 
  
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  

  reqsize <- length(rhs.vars(xformla)) + 5 #right hand side variables (value=6). Le +5 est arbitraire?
  gsize <- subset(gsize, x < reqsize) #x is the column x in the gsize. we keep only the ones respecting the condition.
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }

tlen <- length(tlist) #tlist c'est la liste des traitements. Formé par tlist (unique traitements) ordonné
flen <- length(flist) # flist c'est la liste des cohortes. Formé par flist (first.treatment.name unique traitements) ordonné

results <- chained.compute.mp.spatt.Boot(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, propensityformla,tname, w,
                                 idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,weight_assumption,debT, finT)



return(results)
}






##### ---------------------------------------------------------------------------------------------------------------------------
##### mp.spatt.Boot-----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------



mp.spatt.Boot <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                         w = NULL, idname = NULL, first.treat.name,
                         alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                         biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                         seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                         ,selection ,debT2,finT2,strate,POND_RD=NULL)  ## de nombreux arguments ne sont plus utilis�s


{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  
  data2<-data[,c(first.treat.name,tname)]


  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  

#finT2 c'est  dernière génératio nde traitement considérée finT
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  
  
  flist <- flist[flist > 0]
  
  
  if (is.null(xformla)) {
    xformla <- ~1
  }

  ### tests sur le nombre de traitees par generation
  #group size (proportion % par rapport au nombre de périodes t)
  #on regarde si les groupes sont assez grands 
  
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  

  reqsize <- length(rhs.vars(xformla)) + 5 #right hand side variables (value=6). Le +5 est arbitraire?
  gsize <- subset(gsize, x < reqsize) #x is the column x in the gsize. we keep only the ones respecting the condition.
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }

tlen <- length(tlist)
flen <- length(flist)

results <- compute.mp.spatt.Boot(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                 idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate,POND_RD)



return(results)
}



##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR GMM -----------------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------

mp.spatt.GMM <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                           w = NULL, idname = NULL, first.treat.name,
                           alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                           biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                           seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                           ,selection , ponderation,debT2,finT2,strate,POND_RD=NULL)  ## de nombreux arguments ne sont plus utilis�s
{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  
  data2<-data[,c(first.treat.name,tname)]
  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  flist <- flist[flist > 0]
  
  if (is.null(xformla)) {
    xformla <- ~1
  }
  ### tests sur le nombre de traitees par generation
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  reqsize <- length(rhs.vars(xformla)) + 5
  gsize <- subset(gsize, x < reqsize)
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }
  
  tlen <- length(tlist)
  flen <- length(flist)
  
  results <- compute.mp.spatt.GMM(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                   idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate,POND_RD)
  
  
  return(results)
}



##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR EN CROSS SECTION ----------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------
mp.spatt.CS_Boot <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                            w = NULL, idname = NULL, first.treat.name,
                            alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                            biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                            seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                            ,selection , ponderation,debT2,finT2,strate,POND_RD=NULL)  ## de nombreux arguments ne sont plus utilis�s
{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  data2<-data[,c(first.treat.name,tname)]
  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  flist <- flist[flist > 0]
  if (is.null(xformla)) {
    xformla <- ~1
  }
  ### tests sur le nombre de traitees par generation
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  reqsize <- length(rhs.vars(xformla)) + 5
  gsize <- subset(gsize, x < reqsize)
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }
  tlen <- length(tlist)
  flen <- length(flist)
  results <- compute.mp.spatt.CS_Boot(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                    idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate,POND_RD)
  
  return(results)
}



mp.spatt.V6_CS <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                         w = NULL, idname = NULL, first.treat.name,
                         alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                         biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                         seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                         ,selection , ponderation,debT2,finT2,strate)  ## de nombreux arguments ne sont plus utilis�s
{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  data2<-data[,c(first.treat.name,tname)]
  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  flist <- flist[flist > 0]
  if (is.null(xformla)) {
    xformla <- ~1
  }
  ### tests sur le nombre de traitees par generation
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  reqsize <- length(rhs.vars(xformla)) + 5
  gsize <- subset(gsize, x < reqsize)
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }
  tlen <- length(tlist)
  flen <- length(flist)
  results <- compute.mp.spatt.V6_CS(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                 idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate)
  return(results)
}


##### ---------------------------------------------------------------------------------------------------------------------------
##### POUR L'ESTIMATEUR EN LONG DID ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------------------------------------
mp.spatt.Boot_longDID <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                                 w = NULL, idname = NULL, first.treat.name,
                                 alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                                 biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                                 seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                                 ,selection , ponderation,debT2,finT2,strate)  ## de nombreux arguments ne sont plus utilis�s
{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  data2<-data[,c(first.treat.name,tname)]
  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  ## liste des traitementsjusqu'� finT
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  flist <- flist[flist > 0]
  if (is.null(xformla)) {
    xformla <- ~1
  }
  ### tests sur le nombre de traitees par generation
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  reqsize <- length(rhs.vars(xformla)) + 5
  gsize <- subset(gsize, x < reqsize)
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }
  tlen <- length(tlist)
  flen <- length(flist)
  results <- compute.mp.spatt.Boot_longDID(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                         idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate)
  return(results)
}


mp.spatt.V6_longDID <- function (nom_outcome,nom_traitement, xformla = NULL, data, tname, aggte = TRUE,
                            w = NULL, idname = NULL, first.treat.name,
                            alp = 0.05, method = "logit", se = TRUE, bstrap = FALSE,
                            biters = 100, clustervars = NULL, cband = FALSE, citers = 100,
                            seedvec = NULL, pl = FALSE, cores = 2, printdetails = TRUE
                            ,selection , ponderation,debT2,finT2,strate)  ## de nombreux arguments ne sont plus utilis�s
{
  if (!all(class(data) == "data.frame")) {
    warning("class of data object was not data.frame; converting...")
    data <- as.data.frame(data)
  }
  data2<-data[,c(first.treat.name,tname)]
  tlist <- unique(data2[, tname])[order(unique(data2[, tname]))]
  data2<-data[data[,first.treat.name]<=finT2,c(first.treat.name,tname)]
  flist <- unique(data2[, first.treat.name])[order(unique(data2[,first.treat.name]))]
  flist <- flist[flist > 0]
  if (is.null(xformla)) {
    xformla <- ~1
  }
  ### tests sur le nombre de traitees par generation
  gsize <- aggregate(data[, first.treat.name], by = list(data[,first.treat.name]), function(x) length(x)/length(tlist))
  reqsize <- length(rhs.vars(xformla)) + 5
  gsize <- subset(gsize, x < reqsize)
  if (nrow(gsize) > 0) {
    gpaste <- paste(gsize[, 1], collapse = ",")
    warning(paste0("There are some very small groups in your dataset...\n  This is a very common source of bugs...\n  Check groups: ",
                   gpaste, "\n  and consider dropping these..."))
  }
  tlen <- length(tlist)
  flen <- length(flist)
  results <- compute.mp.spatt.V6_longDID(nom_outcome,nom_traitement,flen, tlen, flist, tlist, data, first.treat.name, xformla, tname, w,
                                    idname, method, seedvec, se, pl, cores, printdetails,selection,ponderation,strate)
  return(results)
}
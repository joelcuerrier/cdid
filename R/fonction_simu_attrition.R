#' Cross-sectional data simulation
#' 
#' @docType data
#' @usage data=data_sim=fonction_simu_attrition(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5)
#' @format data: a data frame
#' @description This function generates a simulated dataset with attrition. The dataset is used to estimate the treatment effect using the chained method. See the simulation in "The Chained Difference-in-Differences" paper.
#' @keywords datasets
#' @examples 
#' data=data_sim=fonction_simu_attrition(theta2_alpha_Gg=0.2, lambda1_alpha_St=0.2, sigma_alpha=2, sigma_epsilon=0.5,alpha_percentile=0.75)
#' @export

fonction_simu_attrition <- function(theta2_alpha_Gg, lambda1_alpha_St, sigma_alpha, sigma_epsilon, alpha_percentile){
# Remarque : pour l'estimateur long DID, on l'estime sur un panel cylindr� qui drop les observations manquantes
# c-à-d, on utilise pour cet estimateur la pond�ration utilis�e pour l'estimateur en Cross Section !!!
set.seed(123)
# Settings
T = 9           # periods 
N0 = 150        # individuals (per two periods!), so N0*2 per periods 
nsims = 1  # simulations hard coded to 1
  
### Uncertainty parameters
mu_a = 1  # Moyenne de alpha_i
sig_a = sigma_alpha # Ecart-type de alpha_i
mu_d = 1  # Moyenne de delta_t
sig_d = 1 # Ecart-type de delta_t
sig_e = sigma_epsilon # Ecart-type de epsilon_it (le terme d'erreur)
mu_x = 1  # Moyenne de X_ig-1
sig_x = 1 # Ecart-type de X_ig-1
# Paramètres qui génèrent le score, i.e., la proba d'être traité, et donc, G
theta0 = -1  # Constante
theta1 = 0.4 # Param�tre X_i
theta2 = theta2_alpha_Gg # 0.0 ou 0.01 # Param�tre alpha_i*t
# Param�tres qui g�n�rent kt (le nombre d'observations) et S_{t,t+1}, le fait d'�tre observ� 2 p�riodes de suite
lambda0 = -1   # Constante
lambda1 = lambda1_alpha_St  # 0.0 ou 0.1 # Param�tre alpha_i*t 

for (simu_i in 1:nsims){
  ALPHA = mu_a+sig_a*matrix(rnorm(1e5), ncol=1) 
  prob = 1/(1+exp(lambda0 + lambda1*ALPHA%*%c(0:(T+1))))
  kt = 1/colMeans(prob)
  rm(ALPHA, prob) 
  
  # For each period: 0 to T, draw individuals conditionally on being treated
  # at t+1 or being in the control group
  N = round(2*max(N0*kt)) # population among which we will draw N0/2 individuals per period
  
  # Unobservable heterogeneity : N = 410 * 8 = 3280 individus
  ALPHA = mu_a+sig_a*matrix(rnorm(N*T),ncol=T) # individual heterogeneity: N indivs * T+2 periods
  DELTA = mu_d+sig_d*matrix(rnorm(1*T),ncol=T) # period-specific heretogeneity
  
  # Treatment selection
  X = mu_x + sig_x*matrix(rnorm(N*T),ncol=T) # observables for propensity scores
  mat_t = matrix(1,nrow=N,ncol=T)
  for (t in 1:T){
    mat_t[,t]=mat_t[,t]*t # C'est juste une matrice pour le temps. Voir la definition du propensity score. mat_t correspond à g. C'est pour que la prob d'être traité dépende de t.
  }

  Score = 1/(1+exp((theta0 + theta1*X + theta2*ALPHA*mat_t))) # propensity score qui d�pend de alpha * t
  #Score = 1/(1+exp((theta0 + theta1*X + theta2*ALPHA))) # propensity score : qui ne d�pend que de alpha
  G = 1*(Score>matrix(runif(N*T),ncol=T))
  ##runif(N * T): This generates a matrix of random numbers drawn from a uniform distribution between 0 and 1. The matrix has dimensions N rows and T columns.
  ## Score > matrix(runif(N * T), ncol = T): This creates a logical matrix by comparing each element of the Score variable with the corresponding element in the matrix of random numbers. Each element of the resulting matrix is TRUE if the corresponding Score is greater than the random number, and FALSE otherwise.
  G[,(1:3)] = 0 # earliest treatment date is t=2 (so t=0,1 no effect)
  colMeans(G)
  
  D = data.frame() # D correspond à si la personne est traitée ou pas dans le passé.
  for (d in 0:(T-1)){
    D1 = cbind( matrix(0,nrow=N,ncol=d) , matrix(G[,(1+d)],nrow=N,ncol=(T-d)) )
    D = rbind(D,D1)
  }
  D = data.matrix(D)
  
  # Terme d'erreur
  EPSI = sig_e*matrix(rnorm(N*T*T), ncol=T)
  # Treatment effect : l'effet d�pend du temps relatif � la date de traitement !
  alpha = matrix(as.vector(ALPHA[,1:T]), nrow=length(as.vector(ALPHA[,1:T])), ncol=T) # stack up ALPHA's

  delta = t(matrix(DELTA, nrow= T, ncol=N*T )) # stack up DELTA's
  beta = c(seq(T-2,1,-1),0,0)/4


  mat_beta = data.frame()
  for (d in 1:T){
    beta_ok = beta[1:(length(beta)-(d-1))]  #Retrait de la derniere valeur de beta
    beta_ok=c(rep(0,(d-1)),beta_ok) 
    
    beta_temp = t(matrix(beta_ok,nrow=T,ncol=N))
    mat_beta = rbind(mat_beta,beta_temp)
  }
  
  effet = t(t(D * mat_beta))
  Y = effet + alpha + delta + EPSI
  
  # Sampling process:
  # Once we have drawn the full population, sample individuals based on the
  # following rule. Here sampling is correlated with time and
  # individual-specific unobservable heterogeneity \alpha_i
  # However, any given individual can only be observed for two consecutive
  # periods. Take the first N0 individuals in each period.
  prob = matrix(0,dim(alpha)[1],dim(alpha)[2]) #fill with zeros
  for (t in 1:dim(prob)[2]){
    prob[,t] = 1/(1+exp((lambda0 + lambda1*alpha[,t]*t)))
  } #replace with probas
  St = 1*(prob>matrix(runif(N*T*T),ncol=T)) # observed at t and t+1

  idx0 = t(1:dim(St)[1]) # full population (15,000 obs)
  idx1 = matrix(0,N0,T) # sampled individuals (observations index) (150 obs)
  idx_drop = matrix(, nrow =0, ncol = 0) #(13,500 obs)
  for (t in 1:T){
    idx = which(St[,t] %in% c(1)) # individuals observed at t and t+1
    idx = setdiff(idx,idx_drop) # tirage sans remise
    perm1 = sample(length(idx))
    idx1[,t] = idx[perm1[1:N0]]
    idx_drop = c(idx_drop,idx1[,t])
  }
  length(unique(as.vector(idx1[,1:T])))==length(idx1[,1:T]) #check : le check 
  
  # matrice de pond�ration li�e � la probabilit� d'�tre observ� deux ann�es de suite
  # pour l'estimateur en chaine
  St_chaine = matrix(0,dim(St)[1],dim(St)[2]) #Matrice vide
  for (t in 1:T){
    St_chaine[idx1[,t],t]=1 #Remplissage de la matrice avec les individus observés
  } #sample des individus observés
  # on g�re les cas o� on s'observe uniquement une entreprise pendant la derni�re p�riode et pas pendant l'avant derni�re p�riode
  # i.e., on met le poids � 0 dans ce cas
  St_chaine[,T]=0
  # Mixture sur la probabilit� d'�tre observ� deux ann�es de suite:
  # Les individus qui ont un alpha_i sup�rieur au quantile � 75% sont tout le temps observ�
  for (t in 1:T){
    St_chaine[alpha[,t] > quantile(alpha[,1], prob=c(alpha_percentile)), t] = 1
  }
  
  # matrice de pondération liée à la probabilité d'�tre observ� une ann�e donn�e t (i.e., P(S_t,t+1)=1 et P(S_t-1,t)=1 )
  # pour l'estimateur cross section
  
  St_CS = matrix(0,dim(St)[1],dim(St)[2])
  for (t in 1:T){
    St_CS[idx1[,t],t]=1
  }
  St_bis = cbind(rep(0,N*T),St_CS[,1:(T-1)])
  St_CS = St_CS + St_bis
  # Mixture sur la probabilit� d'�tre observ� deux ann�es de suite:
  # Les individus qui ont un alpha_i sup�rieur au quantile � 75% sont tout le temps observ�
  for (t in 1:T){
    St_CS[alpha[,t] > quantile(alpha[,1], prob=c(alpha_percentile)), t] = 1
  }
  
  ID = matrix(seq.int(nrow(Y)),nrow=(N*T), ncol=T)
  time_var = as.vector(0:8)
  time_var = t(matrix(time_var,nrow=T,ncol=(N*T)))
  
  XX = as.vector(X)
  XX = matrix(XX,nrow=(N*T),ncol=T) 
  
  GG = as.vector(G*kronecker(t(c(0:(T-1))),rep(1,dim(G)[1])))
  GG = matrix(GG,nrow=(N*T),ncol=T) 
  
  # Cr�ation des variables dans les donn�es finales : il faut cr�er au moins deux Y car le code ne fonctionne qu'avec au moins deux variables d�pendantes
  data_sim = data.frame(id = as.vector(ID) )
  data_sim$annee = as.vector(time_var)
  data_sim$Y1_chaine = as.vector(Y)
  data_sim$Y2_chaine = as.vector(Y)
  data_sim$Y1_CS = as.vector(Y)
  data_sim$Y2_CS = as.vector(Y)
  data_sim$Y1_longDID = as.vector(Y)
  data_sim$Y2_longDID = as.vector(Y)
  data_sim$X = as.vector(XX)
  data_sim$annee_G = as.vector(GG)
  data_sim$P_Y1_chaine = as.vector(St_chaine)
  data_sim$P_Y2_chaine = as.vector(St_chaine)
  data_sim$P_Y1_CS = as.vector(St_CS)
  data_sim$P_Y2_CS = as.vector(St_CS)
  data_sim$P_Y1_longDID = as.vector(St_CS)
  data_sim$P_Y2_longDID = as.vector(St_CS)
  data_sim$traite_G = 1*(data_sim$annee_G>0)
  data_sim$select = 1
  for (i in 1:8){
    pondRDi<-paste0("pondRD",as.character(i))
    data_sim[,pondRDi] = 1
  }
  data_sim$strate = 1
  # Supprimer la premi�re observation t=0 (car on a cr�� cette date uniquement pour initialiser P(S_t,t+1))
  data_sim = subset(data_sim, annee!= 0)
  for (i in 1:8){
  pondRDi<-paste0("pondRD",as.character(i))
  data_sim[,pondRDi] = 1
  }

  # Added 5th od December
  # We are fixing the missing observations.
  # We remove Y2.
  
  data_sim$Y1_longDID[data_sim$P_Y1_longDID!=1] <- NA
  data_sim$Y1_CS[data_sim$P_Y1_longDID!=1] <- NA
  
  data_sim <- subset(data_sim, select = -Y2_chaine)
  data_sim <- subset(data_sim, select = -Y2_CS)
  data_sim <- subset(data_sim, select = -Y2_longDID)
  
  data_sim <- subset(data_sim, select = -P_Y2_chaine)
  data_sim <- subset(data_sim, select = -P_Y2_CS)
  data_sim <- subset(data_sim, select = -P_Y2_longDID)
  
  




  
  
  return(data_sim)
    } 
  }
  
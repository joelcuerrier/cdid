#' Simulate Unbalanced Panel Data
#'
#' @description This function generates a simulated dataset with treatment assignment,
#' individual-level heterogeneity, and time-varying effects. It incorporates attrition
#' based on individual characteristics and time periods. For more details on the methodology, see:
#' Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences",
#' Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2023.11.002.
#'
#' @importFrom stats rnorm runif
#'
#' @param N Number of units
#' @param TT Number of periods
#' @param theta2_alpha_Gg Coefficient for interaction between individual heterogeneity and time in the propensity score.
#' @param lambda1_alpha_St Coefficient for individual heterogeneity in the propensity score.
#' @param sigma_alpha Standard deviation of individual heterogeneity (alpha).
#' @param sigma_epsilon Standard deviation of the error term (epsilon).
#' @param tprob Probability target to get approximately N*TT*tprob observations
#' @return A data frame containing simulated data.
#' @examples
#' data_sim <- fonction_simu_attrition(N=150,TT=9,theta2_alpha_Gg = 0.01,
#' lambda1_alpha_St = 0.5, sigma_alpha = 2, sigma_epsilon = 0.5, tprob=0.5)
#' @export

fonction_simu_attrition <- function(N,TT,theta2_alpha_Gg, lambda1_alpha_St, sigma_alpha, sigma_epsilon, tprob){

# Settings
TT = TT+1           # periods
N = floor(N/TT)        # unique individuals
k = 1
tprob = tprob
# Parameters
mu_a <- 1        # Mean of individual heterogeneity (alpha_i)
sig_a <- sigma_alpha  # Standard deviation of alpha_i
mu_d <- 1        # Mean of time-specific heterogeneity (delta_t)
sig_d <- 1       # Standard deviation of delta_t
sig_e <- sigma_epsilon # Standard deviation of the error term (epsilon_it)
mu_x <- 1        # Mean of covariate X_i
sig_x <- 1       # Standard deviation of X_i

# Parameters for propensity score (probability of treatment assignment)
theta0 <- -1     # Constant
theta1 <- 0.4    # Coefficient for covariate X_i
theta2 <- 0.2 #theta2_alpha_Gg # Interaction coefficient for alpha_i * t

# Sampling parameters
lambda0 = -1   # Constante
lambda1 = 0.1 #lambda1_alpha_St # Parameter assigned to alpha_i*t

# Generate individual heterogeneity (alpha)
ALPHA <- mu_a + sig_a * matrix(rnorm(N * TT), ncol = TT)  # Individual heterogeneity for N individuals over TT periods
DELTA <- mu_d + sig_d * matrix(rnorm(TT), ncol = TT)      # Time-specific heterogeneity for TT periods

# Generate covariates
X <- mu_x + sig_x * matrix(rnorm(N * TT), ncol = TT)      # Covariates for propensity scores

# Create a time matrix for propensity score calculation
  mat_t = matrix(1,nrow=N,ncol=TT)
  for (t in 1:TT){
    mat_t[,t]=mat_t[,t]*t
  }

  Score = 1/(1+exp((theta0 + theta1*X + theta2*ALPHA*mat_t))) # propensity score qui d�pend de alpha * t
  G = 1*(Score>matrix(runif(N*TT),ncol=TT))

  # Enforce no treatment before period 2 (t = 2)
  G[,(1:3)] = 0
  colMeans(G)

  # Create a matrix indicating past treatment (D)
  D = data.frame()
  for (d in 0:(TT-1)){
    D1 = cbind( matrix(0,nrow=N,ncol=d) , matrix(G[,(1+d)],nrow=N,ncol=(TT-d)) )
    D = rbind(D,D1)
  }
  D = data.matrix(D)

  # Generate errors
  EPSI = sig_e*matrix(rnorm(N*TT*TT), ncol=TT)

  # Treatment effects
  alpha = matrix(as.vector(ALPHA[,1:TT]), nrow=length(as.vector(ALPHA[,1:TT])), ncol=TT) # stack up ALPHA's
  delta = t(matrix(DELTA, nrow= TT, ncol=N*TT )) # stack up DELTA's
  beta = c(seq(TT-2,1,-1),0,0)/4

  mat_beta = data.frame()
  for (d in 1:TT){
    beta_ok = beta[1:(length(beta)-(d-1))]  #Retrait de la derniere valeur de beta
    beta_ok=c(rep(0,(d-1)),beta_ok)

    beta_temp = t(matrix(beta_ok,nrow=TT,ncol=N))
    mat_beta = rbind(mat_beta,beta_temp)
  }

  effet = t(t(D * mat_beta))
  Y = effet + alpha + delta + EPSI

  # Create time and ID variables
  ID = matrix(seq.int(nrow(Y)),nrow=(N*TT), ncol=TT)
  time_var = as.vector(0:8)
  time_var = t(matrix(time_var,nrow=TT,ncol=(N*TT)))

  XX = as.vector(X)
  XX = matrix(XX,nrow=(N*TT),ncol=TT)

  GG = as.vector(G*kronecker(t(c(0:(TT-1))),rep(1,dim(G)[1])))
  GG = matrix(GG,nrow=(N*TT),ncol=TT)

  # Assemble the balanced panel dataset
  data_sim = data.frame(id = as.vector(ID) )
  data_sim$date = as.vector(time_var)
  data_sim$Y = as.vector(Y)
  data_sim$X = as.vector(XX)
  data_sim$date_G = as.vector(GG)
  data_sim$S = 1

  # Sampling process
  prob = matrix(0,dim(alpha)[1],dim(alpha)[2]) #fill with zeros
  for (t in 1:dim(prob)[2]){
    prob[,t] = 1/(1+exp((lambda0 + lambda1*alpha[,t]*t - 0.0*X)))
  } #replace with probas
  prob <- tprob*prob/mean(colMeans(prob))
  St = 1*(prob>matrix(runif(N*TT*TT),ncol=TT)) # observed at t and t+1
  data_sim$S <- as.vector(St)

  data_sim <- data_sim[data_sim$date!=0,]

  return(data_sim)

  }

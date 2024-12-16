library(MASS)
library(simsurv)
library(mvtnorm)
library(cubature)
library(pracma)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(survival)
library(timereg)
library(mice)
library(lattice)
library(hexbin)
library(tidyverse)
library(ggridges)
library(patchwork)
library(viridis)
library(hrbrthemes)
library(gapminder)
library(kableExtra)
library(patchwork)
library(tayloRswift)


#### SIMULATE DATA ####

sim_data_cox <- function(n, mu, true){
  
  X1 <- rbinom(n, size = 1, prob = 0.3)
  X2 <- rbinom(n, size = 1, prob = 0.7)
  X3 <- runif(n, min = 0, max = 3)
  
  X <- cbind(X1, X2, X3)
  A <- rbinom(n, size = 1, prob = mu)
  
  covariates_coefficients <- c(log(0.7), log(1.3), log(0.8))
  linear_combination <- X %*% covariates_coefficients
  
  U <- runif(n, min = 0, max = 1)
  
  Tstar <- -log(U)/(0.1*exp(true*A + linear_combination))
  
  C <- runif(n, min = 0, max = 50)
  C = apply(cbind(10,C),1,min)
  TT <- pmin(Tstar,C)
  status <- as.integer(Tstar <= C)
  Cstatus <- as.integer(C <= Tstar)
  
  dat <- data.frame(eventtime = TT, status = status, Cstatus = Cstatus, trt = A, covariates = X)
  
  return(dat)
}


sim_data_additive <- function(n, mu, true){
  X1 <- rbinom(n, size = 1, prob = 0.3)
  X2 <- rbinom(n, size = 1, prob = 0.7)
  X3 <- runif(n, min = 0, max = 3)
  
  X <- cbind(X1, X2, X3)
  A <- rbinom(n, size = 1, prob = mu)
  
  covariates_coefficients <- c(0.6, 0.5, 0.3)
  linear_combination <- X %*% covariates_coefficients
  
  U <- runif(n, min = 0, max = 1)
  
  Tstar <- -log(U)/(0.1 + true*A + linear_combination)
  
  C <- runif(n, min = 0, max = 5)
  C <- apply(cbind(0.7, C),1,min)
  TT <- pmin(Tstar,C)
  status <- as.integer(Tstar <= C)
  Cstatus <- as.integer(C <= Tstar)
  
  dat <- data.frame(eventtime = TT, status = status, Cstatus = Cstatus, trt = A, covariates = X)
  
  return(dat)
}


#### Plot to show distribution of variables #####


X1 <- rbinom(500, size = 1, prob = 0.3)
X2 <- rbinom(500, size = 1, prob = 0.7)
X3 <- runif(500, min = 0, max = 3)
X <- cbind(X1, X2, X3)
A <- rbinom(500, size = 1, prob = 0.5)
U <- runif(500, min = 0, max = 1)

covariates_coefficients <- c(log(0.7), log(1.3), log(0.8))
linear_combination <- X %*% covariates_coefficients

Tstar <- -log(U)/(0.1*exp(log(2)*A + linear_combination))

C <- runif(500, min = 0, max = 50)
C = apply(cbind(10,C),1,min)
TT <- pmin(Tstar,C)


time_data <- data.frame(
  Tstar = Tstar,
  C = C,
  Tobs = TT
)
sim_long <- pivot_longer(time_data, cols = everything(), names_to = "variable", values_to = "value")

plot_cox <- ggplot(sim_long, aes(x = value, fill = factor(variable, levels = c("Tstar", "C", "Tobs")))) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  geom_density(alpha = 0.7) +
  labs(
    x = "Time",
    y = "Density",
    fill = "Variable"
  ) +
  xlim(0, 32) +
  scale_fill_manual(
    values = c("Tstar" = "#F77F00", "C" = "#00B4D8", "Tobs" = "#7209B7"),
    labels = c(expression(T^"*"), "C", "T")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.75, 0.75),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )


X1 <- rbinom(500, size = 1, prob = 0.3)
X2 <- rbinom(500, size = 1, prob = 0.7)
X3 <- runif(500, min = 0, max = 3)
X <- cbind(X1, X2, X3)
A <- rbinom(500, size = 1, prob = 0.5)
U <- runif(500, min = 0, max = 1)

covariates_coefficients <- c(0.6, 0.5, 0.3)
linear_combination <- X %*% covariates_coefficients

U <- runif(500, min = 0, max = 1)

Tstar <- -log(U)/(0.1 + 0*A + linear_combination)

quantile(Tstar, 0.6)
C <- runif(500, min = 0, max = 5)
C <- apply(cbind(0.7, C),1,min)
TT <- pmin(Tstar,C)

time_data <- data.frame(
  Tstar = Tstar,
  C = C,
  Tobs = TT
)
sim_long <- pivot_longer(time_data, cols = everything(), names_to = "variable", values_to = "value")


plot_additive <- ggplot(sim_long, aes(x = value, fill = factor(variable, levels = c("Tstar", "C", "Tobs")))) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  geom_density(alpha = 0.7) +
  labs(
    x = "Time",
    y = "Density",
    fill = "Variable"
  ) +
  xlim(0, 3) +
  scale_fill_manual(
    values = c("Tstar" = "#F77F00", "C" = "#00B4D8", "Tobs" = "#7209B7"),
    labels = c(expression(T^"*"), "C", "T") 
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )

combined_plot <- plot_additive + plot_cox
  plot_layout(ncol = 2)


  

#### OBTAIN ESTIMATES ####

One_Step_estimate_cox <- function(DATA, tau1, tau2, mu){
  
  g <- function(x) {log(-log(x))}
  gdiff <- function(x) {-1 / (x * log(x))}
  
  covariates <- DATA[, !colnames(DATA) %in% c("eventtime", "status", "Cstatus", "trt")]
  covariate_terms <- paste0("prop(", colnames(covariates), ")", collapse = " + ")
  formula_str_temp <- paste("Surv(eventtime, status) ~ prop(trt) + ", covariate_terms) 
  formula_str <- as.formula(formula_str_temp)
  
  fit <- cox.aalen(formula_str, data = DATA) 
  
  M = dim(fit$cum)[1] - 1 # Number of events for T*
  jump_times = fit$cum[,1][-1] # Jumping times
  Cum_base_haz = fit$cum[,2] # Cumulative baseline hazard in jumping times
  betahat = fit$gamma
  Ahat_cond = (exp(as.matrix(cbind(DATA$trt,covariates))%*%c(betahat))%*%c(Cum_base_haz))[,2:(M + 1)]
  Shat_cond = round(exp(-Ahat_cond), 6)
  Shat_cond0 = exp(-exp(as.matrix(cbind(0, covariates))%*%c(betahat))%*%c(Cum_base_haz))[,2:(M + 1)] 
  Shat_cond1 = exp(-exp(as.matrix(cbind(1, covariates))%*%c(betahat))%*%c(Cum_base_haz))[,2:(M + 1)]
  dAhat_cond = cbind(Ahat_cond[,1],Ahat_cond[,2:M]-Ahat_cond[,-M])
  

  Ahat <- numeric(M)  # Population-level cumulative hazard
  for (j in 1:M) {
    t <- jump_times[j]
    risk_set <- which(DATA$eventtime >= t) # Finding risk set: individuals i with observed event time >= t
    sum_A_cond_AL <- sum(Ahat_cond[risk_set, j]) # Compute numerator
    if (length(risk_set) > 0) {
      Ahat[j] <- sum_A_cond_AL / length(risk_set)
    } else {
      Ahat[j] <- ifelse(j > 1, Ahat[j - 1], 0)  # Use the previous value or 0 if j == 1
    }
  }
  dAhat = c(Ahat[1], diff(Ahat))
  
  f_Tdt <- exp(-Ahat)*dAhat
  
  ## FOR CENSORING ##
  formula_str_C <- paste("Surv(eventtime, Cstatus) ~  prop(trt) + ", covariate_terms)
  formula_obj_C <- as.formula(formula_str_C)
  fitC = cox.aalen(formula_obj_C, data = DATA)
  M_C = dim(fitC$cum)[1] - 1 # Number of events for C
  t_C = fitC$cum[,1][-1] # Jumping times for C
  Cum_base_hazC = fitC$cum[,2] # Cumulative baseline hazard in jumping times
  betahatC = fitC$gamma
  AhatC_cond = (exp(as.matrix(cbind(DATA$trt,covariates))%*%c(betahatC))%*%c(Cum_base_hazC))[,2:(M_C + 1)]
  Khat_cond = exp(-AhatC_cond)
  dAhatC_cond = cbind(AhatC_cond[,1],AhatC_cond[,2:M_C]-AhatC_cond[,-M_C])
  ##finding K in jump_times. Taking the column from Khat corresponding to the largest jumping time for C before the i'th jumping time for T^*
  ##since Khat is constant between the jumping times this is the same as evaluating Khat in the jumping times for S
  K_in_jump_times = matrix(NA, nrow = nrow(DATA), ncol = M)
  for(j in 1:M) {
    valid_indices <- which(t_C <= jump_times[j])
    if (length(valid_indices) > 0) {
      K_in_jump_times[,j] <- Khat_cond[,max(valid_indices)]
    } else {
      K_in_jump_times[,j] <- 1  # Default to 1
    }
  }
  
  
  S_in_t_C <- matrix(NA, nrow = nrow(DATA), ncol = M_C)
  for(i in 1:M_C) {
    valid_indices <- which(jump_times <= t_C[i])
    if (length(valid_indices) > 0) {
      S_in_t_C[,i] <- Shat_cond[,max(valid_indices)]
    } else {
      S_in_t_C[,i] <- 1  # Default to 1
    }
  }

  
  
  ## OBTAIN S(t|A) (as two vectors)
  Shat_condA0 <- numeric(M)
  Shat_condA1 <- numeric(M)
  for (j in 1:M){
    Shat_condA0[j] <- mean(Shat_cond[DATA$trt == 0, j])
    Shat_condA1[j] <- mean(Shat_cond[DATA$trt == 1, j])
  }
  
  g_Shat_condA0 <- numeric(M)
  g_Shat_condA1 <- numeric(M)
  for (j in 1:M){
    g_Shat_condA0[j] <- g(Shat_condA0[j])
    g_Shat_condA1[j] <- g(Shat_condA1[j])
  }
  
  g_difference <- g_Shat_condA1 - g_Shat_condA0  
  
  
  ######## STOCHASTIC INTEGRALS ######
  at.risk = matrix(NA, nrow = nrow(DATA), ncol = M)
  dN = matrix(NA, nrow = nrow(DATA), ncol = M)
  dM = matrix(NA, nrow = nrow(DATA), ncol = M)
  for (j in 1:M){
    for(i in 1:nrow(DATA)){
      at.risk[i,j] = as.numeric(DATA$eventtime[i]>=jump_times[j])
      dN[i,j] = as.numeric(DATA$eventtime[i]==jump_times[j])*as.numeric(DATA$status[i]==1)
      dM[i,j] = dN[i,j]-at.risk[i,j]*dAhat_cond[i,j]
    }
  }
  
  at.risk_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  dN_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  dM_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  for (j in 1:M_C){
    for(i in 1:nrow(DATA)){
      at.risk_C[i,j] = as.numeric(DATA$eventtime[i]>=t_C[j])
      dN_C[i,j] = as.numeric(DATA$eventtime[i]==t_C[j])*as.numeric(DATA$status[i]==0)
      dM_C[i,j] = dN_C[i,j]-at.risk_C[i,j]*dAhatC_cond[i,j]
    }
  }
  
  ## integral with M ##
  
  first_integral <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    first_integral[i] <- 0
    
    for (j in 1:M) {
      t_j <- jump_times[j]
      
      # Check if the current event time is within the desired range
      if (t_j >= tau1 && t_j <= tau2) {
        inner_sum <- 0
        
        for (m in 1:M) {
          t_m <- jump_times[m]
          if (t_m <= t_j) {
            if (!is.na(Shat_cond[i, m]) && Shat_cond[i, m] > 1e-10 &&
                !is.na(K_in_jump_times[i, m]) && K_in_jump_times[i, m] > 1e-10) {
              inner_sum <- inner_sum + (1 / (Shat_cond[i, m] * K_in_jump_times[i, m])) * dM[i, m]
            }
          }
        }
        # Outer sum: add the contribution from t_j
        first_integral[i] <- first_integral[i] + 
          gdiff(if(DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] }) *(Shat_cond[i, j] * (1 - inner_sum) - (if(DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] })) * f_Tdt[j]
      }
    }
    
  }
  
  
  ## integral with M_C ##
  f_T_cond_dt <- matrix(NA, nrow = nrow(DATA), ncol = M)
  for (i in 1:nrow(DATA)){
    for (j in 1:M){
      f_T_cond_dt[i,j] <- exp(-Ahat_cond[i,j])*dAhat_cond[i,j]
    }
  }
  
  
  integral_with_M_C_result <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    integral_with_M_C_result[i] <- 0
    
    for (j in 1:M) {
      t_j <- jump_times[j]
      
      # Check if the current event time is within [tau1, tau2]
      if (t_j >= tau1 && t_j <= tau2) {
        inner_sum <- 0 # Inner sum over censoring times
        
        for (m in 1:M_C) {
          u_m <- t_C[m]
          if (u_m <= t_j &&
              !is.na(S_in_t_C[i, m]) && S_in_t_C[i, m] > 1e-10 &&
              !is.na(Khat_cond[i, m]) && Khat_cond[i, m] > 1e-10) {
            inner_sum <- inner_sum + (1 / (S_in_t_C[i, m] * Khat_cond[i, m])) * dM_C[i, m]
          }
        }
        
        # Outer sum: add the contribution from t_j
        integral_with_M_C_result[i] <- integral_with_M_C_result[i] + g_difference[j] * inner_sum * f_T_cond_dt[i, j]
      }
    }

    }
    
  
  #### TERMS ####
  
  first_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    first_term_for_i[i] <- (DATA$trt[i] - mu)*first_integral[i]
  }
  
  #### ####
  
  second_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)) {
    terms_in_second_integral <- numeric(M)  # Initialize for each individual
    for (j in 1:M) {
      t_j <- jump_times[j]
      # Only include terms where t_j is within [tau1, tau2]
      if (tau1 <= t_j && t_j <= tau2) {
        S_value <- if (DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] }
        terms_in_second_integral[j] <- (g(S_value) - (mu * g_Shat_condA1[j] + (1 - mu) * g_Shat_condA0[j])) * f_Tdt[j]
      } else {
        terms_in_second_integral[j] <- 0  # Set to 0 for out-of-range t_j
      }
    }
    second_term_for_i[i] <- (DATA$trt[i] - mu) * sum(terms_in_second_integral)
  }
  
  #### ####
  
  terms <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    indicator <- ifelse(tau1 <= DATA$eventtime[i] & DATA$eventtime[i] <= tau2, 1, 0)
    if (DATA$status[i] == 1){
      closest_index <- max(which(jump_times <= DATA$eventtime[i]))
      terms[i] <- (g_difference[closest_index]*indicator)/(K_in_jump_times[i, closest_index])
    } else {terms[i] <- 0}
  }
  third_term_for_i <- mu*(1-mu)*terms
  
  
  #### ####
  
  fourth_term_for_i <- mu*(1-mu)*integral_with_M_C_result
  
  #### ####
  
  
  #### Additional terms after projecting onto tangent space ####
  
fifth_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)) {
    terms_in_int <- numeric(M)
    for (j in 1:M) {
      t_j <- jump_times[j]  # Get the time corresponding to index j
      if (tau1 <= t_j && t_j <= tau2) {  # Include only if within the range
        terms_in_int[j] <- (mu * ((1 - mu)^2 * gdiff(Shat_condA1[j]) * (Shat_cond1[i, j] - Shat_condA1[j])) + 
                              (1 - mu) * (mu^2 * gdiff(Shat_condA0[j]) * (Shat_cond0[i, j] - Shat_condA0[j]))) * f_Tdt[j]
      } else {
        terms_in_int[j] <- 0  # Set to 0 if outside the range
      }
    }
    fifth_term_for_i[i] <- ((DATA$trt[i] - mu) / (mu * (1 - mu))) * sum(terms_in_int)
  }
  
  
  #### ####
  
  sixth_term_for_i <- numeric(nrow(DATA))
  
  
  for (i in 1:nrow(DATA)){
    terms_in_int <- numeric(M)
    for (j in 1:M){
      t_j <- jump_times[j]
      if (tau1 <= t_j && t_j <= tau2){
        terms_in_int[j] <- (DATA$trt[i]-mu)*(1-2*mu)*(g_difference[j])*f_Tdt[j]
      } else { terms_in_int[j] <- 0 }
    }
    sixth_term_for_i[i] <- sum(terms_in_int)
  }
  
  #### ####
  
  new_variable <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    indicator <- ifelse(tau1 <= DATA$eventtime[i] & DATA$eventtime[i] <= tau2, 1, 0)
    if (DATA$status[i] == 1){
      closest_index <- max(which(jump_times <= DATA$eventtime[i]))
      new_variable[i] <- (DATA$trt[i] - mu)*g_difference[closest_index]*indicator/K_in_jump_times[i, closest_index]
    }
    else {new_variable[i] <- 0}
  }
  
  data_for_fit <- covariates
  data_for_fit$Y <- new_variable
  covariate_terms <- paste0(colnames(covariates), collapse = " + ")
  formula_fit <- as.formula(paste("Y ~", covariate_terms))
  fit_regression <- lm(formula_fit, data_for_fit)
  
  expectation_for_i <- predict(fit_regression, newdata = covariates)
  
  
  seventh_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    seventh_term_for_i[i] <- (DATA$trt[i]-mu)/(mu*(1-mu))*expectation_for_i[i]
  }
  
  
  #### TERMS FOR BETA2 ####
  
  closest_index_tau1 <- ifelse(any(jump_times <= tau1), max(which(jump_times <= tau1)), NA)
  closest_index_tau2 <- max(which(jump_times <= tau2))
  
  if (is.na(closest_index_tau1)) {S_tau1 <- 1} else {S_tau1 <- exp(-Ahat[closest_index_tau1])}
  S_tau2 <- exp(-Ahat[closest_index_tau2])
  
  first_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    first_term_for_i_2[i] <- (DATA$trt[i] - mu)^2*(S_tau1 - S_tau2)
  }
  
  #### ####
  
  
  integral_tau1 <- numeric(nrow(DATA))
  integral_tau2 <- numeric(nrow(DATA))
  
  M_tau1 <- length(which(jump_times <= tau1))
  M_tau2 <- length(which(jump_times <= tau2))
  
  for (i in 1:nrow(DATA)){
    if(M_tau1 == 0){
      integral_tau1[i] <- 0
    } else {
      terms_in_int1 <- numeric(M_tau1)
      for (j in 1:M_tau1) {
        terms_in_int1[j] <- (1 / (Shat_cond[i, j] * K_in_jump_times[i,j])) * dM[i, j]
      }
      integral_tau1[i] <- sum(terms_in_int1)
    }
    
    
    if(M_tau2 == 0){
      integral_tau2[i] <- 0
    } else {
      terms_in_int2 <- numeric(M_tau2)
      for (j in 1:M_tau2) {
        terms_in_int2[j] <- (1 / (Shat_cond[i, j] * K_in_jump_times[i,j])) * dM[i, j]
      }
      integral_tau2[i] <- sum(terms_in_int2)
    }
  }
  
  S_tau1_cond <- numeric(nrow(DATA))
  S_tau2_cond <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    if (is.na(closest_index_tau1)) {S_tau1_cond[i] <- 1} else {S_tau1_cond[i] <- Shat_cond[i, closest_index_tau1]}
    S_tau2_cond[i] <- Shat_cond[i, closest_index_tau2]
  }
  
  
  
  
  second_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    second_term_for_i_2[i] <- mu*(1-mu)*(S_tau1_cond[i]*(1 - integral_tau1[i]) - S_tau2_cond[i]*(1 - integral_tau2[i]))
  }
  
  #### ####
  
  third_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    third_term_for_i_2[i] <- (DATA$trt[i]-mu)*(1-2*mu)*(S_tau1 - S_tau2)
  }
  
  #### ####
  
  S_tau1_cond0 <- numeric(nrow(DATA))
  S_tau2_cond0 <- numeric(nrow(DATA))
  S_tau1_cond1 <- numeric(nrow(DATA))
  S_tau2_cond1 <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    if (is.na(closest_index_tau1)) {S_tau1_cond0[i] <- 1} else {S_tau1_cond0[i] <- Shat_cond0[i, closest_index_tau1]}
    S_tau2_cond0[i] <- Shat_cond0[i, closest_index_tau2]
    
    if (is.na(closest_index_tau1)) {S_tau1_cond1[i] <- 1} else {S_tau1_cond1[i] <- Shat_cond1[i, closest_index_tau1]}
    S_tau2_cond1[i] <- Shat_cond1[i, closest_index_tau2]
  }
  
  fourth_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    fourth_term_for_i_2[i] <- (DATA$trt[i]-mu)*(mu*(1-mu)*(S_tau1_cond1[i] - S_tau2_cond1[i]) + (1-mu)*(-mu)*(S_tau1_cond0[i] - S_tau2_cond0[i]))
  }
  
  #### ####
  
  
  #### Obtain final estimates for beta1 and beta2

  
  beta1_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    beta1_for_i[i] <- first_term_for_i[i] + second_term_for_i[i] + third_term_for_i[i] + fourth_term_for_i[i] - fifth_term_for_i[i] - sixth_term_for_i[i] - seventh_term_for_i[i]
  }
  
  
  beta2_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    beta2_for_i[i] <- first_term_for_i_2[i] + second_term_for_i_2[i] - third_term_for_i_2[i] - fourth_term_for_i_2[i]
  }
  
  beta_1_est <- 1/2*mean(beta1_for_i, na.rm = TRUE)
  beta_2_est <- 1/2*mean(beta2_for_i, na.rm = TRUE)
  
  OS_est <- beta_1_est/beta_2_est
  
  ## Variance estimation ##
  
  varphi_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    varphi_for_i[i] <- (((beta1_for_i[i] - 2*beta_1_est) - OS_est*(beta2_for_i[i] - 2*beta_2_est))/(beta_2_est))^2
  }
  
  variance_est <- (1/nrow(DATA))*mean(varphi_for_i)
  
  CI_OS <- c(OS_est - 1.96*sqrt(variance_est), OS_est + 1.96*sqrt(variance_est))
  
  ## Plug-in estimator ##
  
  PI_terms_in_int <- numeric(M)
  for (j in 1:M){
    t_j <- jump_times[j]
    if (tau1 <= t_j && t_j <= tau2){
      PI_terms_in_int[j] <- PI_terms_in_int[j] <- g_difference[j]*f_Tdt[j]/(S_tau1 - S_tau2)
    } else { terms_in_int[j] <- 0 }
  }
  PI_estimate <- sum(PI_terms_in_int)
  
  out <- list(OS_est = OS_est, OS_var = variance_est, CI_OS = CI_OS, PI_estimate = PI_estimate)
  
  return(out)
}


One_Step_estimate_additive <- function(DATA, tau1, tau2, mu){
  
  g <- function(x) {log(-log(x))} 
  gdiff <- function(x) {-1 / (x * log(x))}
  
  covariates <- DATA[, !colnames(DATA) %in% c("eventtime", "status", "Cstatus", "trt")]
  covariate_terms <- paste0("const(", colnames(covariates), ")", collapse = " + ")
  formula_str_temp <- paste("Surv(eventtime, status) ~ const(trt) + ", covariate_terms)
  formula_str <- as.formula(formula_str_temp)
  
  fit <- aalen(formula_str, data = DATA) 
  
  M = dim(fit$cum)[1] - 1 # Number of events for T*
  jump_times = fit$cum[,1][-1] # Jumping times
  Cum_base_haz = fit$cum[,2][-1] # Cumulative baseline hazard in jumping times
  betahat = fit$gamma
  baseline_matrix = matrix(Cum_base_haz, nrow = nrow(DATA), ncol = M, byrow = TRUE)
  time_matrix = matrix(jump_times, nrow = nrow(DATA), ncol = M, byrow = TRUE)
  cov_and_trt = DATA[, !colnames(DATA) %in% c("eventtime", "status", "Cstatus")]
  covariate_contributions = as.matrix(cov_and_trt) %*% betahat
  time_scaled_contributions = sapply(jump_times, function(t) covariate_contributions*t)
  Ahat_cond = baseline_matrix + time_scaled_contributions
  Shat_cond = round(exp(-Ahat_cond), 6)
  dAhat_cond = cbind(Ahat_cond[,1],Ahat_cond[,2:M]-Ahat_cond[,-M])
  
  cov_and_trt$trt <- 1
  covariate_contributions1 = as.matrix(cov_and_trt) %*% betahat
  time_scaled_contributions1 = sapply(jump_times, function(t) covariate_contributions1*t)
  Ahat_cond1 = baseline_matrix + time_scaled_contributions1
  Shat_cond1 = round(exp(-Ahat_cond1), 6)
  
  cov_and_trt$trt <- 0
  covariate_contributions0 = as.matrix(cov_and_trt) %*% betahat
  time_scaled_contributions0 = sapply(jump_times, function(t) covariate_contributions0*t)
  Ahat_cond0 = baseline_matrix + time_scaled_contributions0
  Shat_cond0 = round(exp(-Ahat_cond0), 6)
  
  Ahat <- numeric(M)  # Initialize population-level cumulative hazard
  for (j in 1:M) {
    t <- jump_times[j]
    risk_set <- which(DATA$eventtime >= t) # Finding risk set: individuals i with observed eventtime >= t
    sum_A_cond_AL <- sum(Ahat_cond[risk_set, j]) # Compute numerator
    if (length(risk_set) > 0) {
      Ahat[j] <- sum_A_cond_AL / length(risk_set)
    } else {
      Ahat[j] <- ifelse(j > 1, Ahat[j - 1], 0)  # Use the previous value or 0 if j == 1
    }
  }
  dAhat = c(Ahat[1], diff(Ahat))
  
  f_Tdt <- exp(-Ahat)*dAhat
  
  ## FOR CENSORING ##
  formula_str_C <- paste("Surv(eventtime, Cstatus) ~  const(trt) + ", covariate_terms)
  formula_obj_C <- as.formula(formula_str_C)
  fitC = aalen(formula_obj_C, data = DATA)
  M_C = dim(fitC$cum)[1] - 1 # Number of events for C
  t_C = fitC$cum[,1][-1] # Jumping times for C
  Cum_base_hazC = fitC$cum[,2][-1] # Cumulative baseline hazard in jumping times
  betahatC = fitC$gamma
  baseline_matrixC = matrix(Cum_base_hazC, nrow = nrow(DATA), ncol = M_C, byrow = TRUE)
  time_matrixC = matrix(t_C, nrow = nrow(DATA), ncol = M_C, byrow = TRUE)
  cov_and_trt = DATA[, !colnames(DATA) %in% c("eventtime", "status", "Cstatus")]
  covariate_contributionsC = as.matrix(cov_and_trt) %*% betahatC
  time_scaled_contributionsC = sapply(t_C, function(t) covariate_contributionsC*t)
  AhatC_cond = baseline_matrixC + time_scaled_contributionsC
  Khat_cond = round(exp(-AhatC_cond), 6)
  dAhatC_cond = cbind(AhatC_cond[,1],AhatC_cond[,2:M_C]-AhatC_cond[,-M_C])

  K_in_jump_times = matrix(NA, nrow = nrow(DATA), ncol = M)
  for(j in 1:M) {
    valid_indices <- which(t_C <= jump_times[j])
    if (length(valid_indices) > 0) {
      K_in_jump_times[,j] <- Khat_cond[,max(valid_indices)]
    } else {
      K_in_jump_times[,j] <- 1  # Default to 1
    }
  }
  
  
  S_in_t_C <- matrix(NA, nrow = nrow(DATA), ncol = M_C)
  for(i in 1:M_C) {
    valid_indices <- which(jump_times <= t_C[i])
    if (length(valid_indices) > 0) {
      S_in_t_C[,i] <- Shat_cond[,max(valid_indices)]
    } else {
      S_in_t_C[,i] <- 1  # Default to 1
    }
  }
  
  
  
  ## OBTAIN S(t|A) (as two vectors)
  Shat_condA0 <- numeric(M)
  Shat_condA1 <- numeric(M)
  for (j in 1:M){
    Shat_condA0[j] <- mean(Shat_cond[DATA$trt == 0, j])
    Shat_condA1[j] <- mean(Shat_cond[DATA$trt == 1, j])
  }
  
  g_Shat_condA0 <- numeric(M)
  g_Shat_condA1 <- numeric(M)
  for (j in 1:M){
    g_Shat_condA0[j] <- g(Shat_condA0[j])
    g_Shat_condA1[j] <- g(Shat_condA1[j])
  }
  
  g_difference <- g_Shat_condA1 - g_Shat_condA0  
  
  
  ######## STOCHASTIC INTEGRALS ######
  at.risk = matrix(NA, nrow = nrow(DATA), ncol = M)
  dN = matrix(NA, nrow = nrow(DATA), ncol = M)
  dM = matrix(NA, nrow = nrow(DATA), ncol = M)
  for (j in 1:M){
    for(i in 1:nrow(DATA)){
      at.risk[i,j] = as.numeric(DATA$eventtime[i]>=jump_times[j])
      dN[i,j] = as.numeric(DATA$eventtime[i]==jump_times[j])*as.numeric(DATA$status[i]==1)
      dM[i,j] = dN[i,j]-at.risk[i,j]*dAhat_cond[i,j]
    }
  }
  
  at.risk_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  dN_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  dM_C = matrix(NA, nrow = nrow(DATA), ncol = M_C)
  for (j in 1:M_C){
    for(i in 1:nrow(DATA)){
      at.risk_C[i,j] = as.numeric(DATA$eventtime[i]>=t_C[j])
      dN_C[i,j] = as.numeric(DATA$eventtime[i]==t_C[j])*as.numeric(DATA$status[i]==0)
      dM_C[i,j] = dN_C[i,j]-at.risk_C[i,j]*dAhatC_cond[i,j]
    }
  }
  
  ## integral with M ##
  
  first_integral <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    first_integral[i] <- 0
    
    for (j in 1:M) {
      t_j <- jump_times[j]
      
      # Check if the current event time is within the desired range
      if (t_j >= tau1 && t_j <= tau2) {
        inner_sum <- 0
        
        for (m in 1:M) {
          t_m <- jump_times[m]
          if (t_m <= t_j) {
            if (!is.na(Shat_cond[i, m]) && Shat_cond[i, m] > 1e-10 &&
                !is.na(K_in_jump_times[i, m]) && K_in_jump_times[i, m] > 1e-10) {
              inner_sum <- inner_sum + (1 / (Shat_cond[i, m] * K_in_jump_times[i, m])) * dM[i, m]
            }
          }
        }
        # Outer sum: add the contribution from t_j
        first_integral[i] <- first_integral[i] + 
          gdiff(if(DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] }) *(Shat_cond[i, j] * (1 - inner_sum) - (if(DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] })) * f_Tdt[j]
      }
    }
    
  }
  
  
  
  ## integral with M_C ##
  f_T_cond_dt <- matrix(NA, nrow = nrow(DATA), ncol = M)
  for (i in 1:nrow(DATA)){
    for (j in 1:M){
      f_T_cond_dt[i,j] <- exp(-Ahat_cond[i,j])*dAhat_cond[i,j]
    }
  }
  
  
  integral_with_M_C_result <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    integral_with_M_C_result[i] <- 0
    
    for (j in 1:M) {
      t_j <- jump_times[j]
      
      # Check if the current event time is within (tau1, tau2)
      if (t_j >= tau1 && t_j <= tau2) {
        inner_sum <- 0 # Inner sum over censoring times
        
        for (m in 1:M_C) {
          u_m <- t_C[m]
          if (u_m <= t_j &&
              !is.na(S_in_t_C[i, m]) && S_in_t_C[i, m] > 1e-10 &&
              !is.na(Khat_cond[i, m]) && Khat_cond[i, m] > 1e-10) {
            inner_sum <- inner_sum + (1 / (S_in_t_C[i, m] * Khat_cond[i, m])) * dM_C[i, m]
          }
        }
        
        # Outer sum: add the contribution from t_j
        integral_with_M_C_result[i] <- integral_with_M_C_result[i] + g_difference[j] * inner_sum * f_T_cond_dt[i, j]
      }
    }
    
  }
  
  
  #### TERMS ####
  
  first_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    first_term_for_i[i] <- (DATA$trt[i] - mu)*first_integral[i]
  }
  
  #### ####
  
  second_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)) {
    terms_in_second_integral <- numeric(M)  # Initialize for each individual
    for (j in 1:M) {
      t_j <- jump_times[j]
      # Only include terms where t_j is within [tau1, tau2]
      if (tau1 <= t_j && t_j <= tau2) {
        S_value <- if (DATA$trt[i] == 1) { Shat_condA1[j] } else { Shat_condA0[j] }
        terms_in_second_integral[j] <- (g(S_value) - (mu * g_Shat_condA1[j] + (1 - mu) * g_Shat_condA0[j])) * f_Tdt[j]
      } else {
        terms_in_second_integral[j] <- 0  # Set to 0 for out-of-range t_j
      }
    }
    second_term_for_i[i] <- (DATA$trt[i] - mu) * sum(terms_in_second_integral)
  }
  
  #### ####
  
  terms <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    indicator <- ifelse(tau1 <= DATA$eventtime[i] & DATA$eventtime[i] <= tau2, 1, 0)
    if (DATA$status[i] == 1){
      closest_index <- max(which(jump_times <= DATA$eventtime[i]))
      terms[i] <- (g_difference[closest_index]*indicator)/(K_in_jump_times[i, closest_index])
    } else {terms[i] <- 0}
  }
  third_term_for_i <- mu*(1-mu)*terms
  
  
  #### ####
  
  fourth_term_for_i <- mu*(1-mu)*integral_with_M_C_result
  
  #### ####
  
  
  #### Terms after projecting onto the tangent space ####
  
  fifth_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)) {
    terms_in_int <- numeric(M)
    for (j in 1:M) {
      t_j <- jump_times[j]  # Get the time corresponding to index j
      if (tau1 <= t_j && t_j <= tau2) {  # Include only if within the range
        terms_in_int[j] <- (mu * ((1 - mu)^2 * gdiff(Shat_condA1[j]) * (Shat_cond1[i, j] - Shat_condA1[j])) + 
                              (1 - mu) * (mu^2 * gdiff(Shat_condA0[j]) * (Shat_cond0[i, j] - Shat_condA0[j]))) * f_Tdt[j]
      } else {
        terms_in_int[j] <- 0  # Set to 0 if outside the range
      }
    }
    fifth_term_for_i[i] <- ((DATA$trt[i] - mu) / (mu * (1 - mu))) * sum(terms_in_int)
  }
  
  
  #### ####
  
  sixth_term_for_i <- numeric(nrow(DATA))
  
  
  for (i in 1:nrow(DATA)){
    terms_in_int <- numeric(M)
    for (j in 1:M){
      t_j <- jump_times[j]
      if (tau1 <= t_j && t_j <= tau2){
        terms_in_int[j] <- (DATA$trt[i]-mu)*(1-2*mu)*(g_difference[j])*f_Tdt[j]
      } else { terms_in_int[j] <- 0 }
    }
    sixth_term_for_i[i] <- sum(terms_in_int)
  }
  
  #### ####
  
  new_variable <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    indicator <- ifelse(tau1 <= DATA$eventtime[i] & DATA$eventtime[i] <= tau2, 1, 0)
    if (DATA$status[i] == 1){
      closest_index <- max(which(jump_times <= DATA$eventtime[i]))
      new_variable[i] <- (DATA$trt[i] - mu)*g_difference[closest_index]*indicator/K_in_jump_times[i, closest_index]
    }
    else {new_variable[i] <- 0}
  }
  
  data_for_fit <- covariates
  data_for_fit$Y <- new_variable
  covariate_terms <- paste0(colnames(covariates), collapse = " + ")
  formula_fit <- as.formula(paste("Y ~", covariate_terms))
  fit_regression <- lm(formula_fit, data_for_fit)
  
  expectation_for_i <- predict(fit_regression, newdata = covariates)
  
  
  seventh_term_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    seventh_term_for_i[i] <- (DATA$trt[i]-mu)/(mu*(1-mu))*expectation_for_i[i]
  }
  
  
  #### TERMS FOR BETA2 ####
  
  closest_index_tau1 <- ifelse(any(jump_times <= tau1), max(which(jump_times <= tau1)), NA)
  closest_index_tau2 <- max(which(jump_times <= tau2))
  
  if (is.na(closest_index_tau1)) {S_tau1 <- 1} else {S_tau1 <- exp(-Ahat[closest_index_tau1])}
  S_tau2 <- exp(-Ahat[closest_index_tau2])
  
  first_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    first_term_for_i_2[i] <- (DATA$trt[i] - mu)^2*(S_tau1 - S_tau2)
  }
  
  #### ####
  
  
  integral_tau1 <- numeric(nrow(DATA))
  integral_tau2 <- numeric(nrow(DATA))
  
  M_tau1 <- length(which(jump_times <= tau1))
  M_tau2 <- length(which(jump_times <= tau2))
  
  for (i in 1:nrow(DATA)){
    if(M_tau1 == 0){
      integral_tau1[i] <- 0
    } else {
      terms_in_int1 <- numeric(M_tau1)
      for (j in 1:M_tau1) {
        terms_in_int1[j] <- (1 / (Shat_cond[i, j] * K_in_jump_times[i,j])) * dM[i, j]
      }
      integral_tau1[i] <- sum(terms_in_int1)
    }
    
    
    if(M_tau2 == 0){
      integral_tau2[i] <- 0
    } else {
      terms_in_int2 <- numeric(M_tau2)
      for (j in 1:M_tau2) {
        terms_in_int2[j] <- (1 / (Shat_cond[i, j] * K_in_jump_times[i,j])) * dM[i, j]
      }
      integral_tau2[i] <- sum(terms_in_int2)
    }
  }
  
  S_tau1_cond <- numeric(nrow(DATA))
  S_tau2_cond <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    if (is.na(closest_index_tau1)) {S_tau1_cond[i] <- 1} else {S_tau1_cond[i] <- Shat_cond[i, closest_index_tau1]}
    S_tau2_cond[i] <- Shat_cond[i, closest_index_tau2]
  }
  
  
  
  
  second_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    second_term_for_i_2[i] <- mu*(1-mu)*(S_tau1_cond[i]*(1 - integral_tau1[i]) - S_tau2_cond[i]*(1 - integral_tau2[i]))
  }
  
  #### ####
  
  third_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    third_term_for_i_2[i] <- (DATA$trt[i]-mu)*(1-2*mu)*(S_tau1 - S_tau2)
  }
  
  #### ####
  
  S_tau1_cond0 <- numeric(nrow(DATA))
  S_tau2_cond0 <- numeric(nrow(DATA))
  S_tau1_cond1 <- numeric(nrow(DATA))
  S_tau2_cond1 <- numeric(nrow(DATA))
  
  for (i in 1:nrow(DATA)){
    if (is.na(closest_index_tau1)) {S_tau1_cond0[i] <- 1} else {S_tau1_cond0[i] <- Shat_cond0[i, closest_index_tau1]}
    S_tau2_cond0[i] <- Shat_cond0[i, closest_index_tau2]
    
    if (is.na(closest_index_tau1)) {S_tau1_cond1[i] <- 1} else {S_tau1_cond1[i] <- Shat_cond1[i, closest_index_tau1]}
    S_tau2_cond1[i] <- Shat_cond1[i, closest_index_tau2]
  }
  
  fourth_term_for_i_2 <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    fourth_term_for_i_2[i] <- (DATA$trt[i]-mu)*(mu*(1-mu)*(S_tau1_cond1[i] - S_tau2_cond1[i]) + (1-mu)*(-mu)*(S_tau1_cond0[i] - S_tau2_cond0[i]))
  }
  
  #### ####
  
  
  
  beta1_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    beta1_for_i[i] <- first_term_for_i[i] + second_term_for_i[i] + third_term_for_i[i] + fourth_term_for_i[i] - fifth_term_for_i[i] - sixth_term_for_i[i] - seventh_term_for_i[i]
  }
  
  
  beta2_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    beta2_for_i[i] <- first_term_for_i_2[i] + second_term_for_i_2[i] - third_term_for_i_2[i] - fourth_term_for_i_2[i]
  }
  
  beta_1_est <- 1/2*mean(beta1_for_i, na.rm = TRUE)
  beta_2_est <- 1/2*mean(beta2_for_i, na.rm = TRUE)
  
  OS_est <- beta_1_est/beta_2_est
  
  ## Variance estimation ##
  
  varphi_for_i <- numeric(nrow(DATA))
  for (i in 1:nrow(DATA)){
    varphi_for_i[i] <- (((beta1_for_i[i] - 2*beta_1_est) - OS_est*(beta2_for_i[i] - 2*beta_2_est))/(beta_2_est))^2
  }
  
  variance_est <- (1/nrow(DATA))*mean(varphi_for_i)
  
  CI_OS <- c(OS_est - 1.96*sqrt(variance_est), OS_est + 1.96*sqrt(variance_est))
  
  ## Plug-in estimator ##
  
  PI_terms_in_int <- numeric(M)
  for (j in 1:M){
    t_j <- jump_times[j]
    if (tau1 <= t_j && t_j <= tau2){
      PI_terms_in_int[j] <- PI_terms_in_int[j] <- g_difference[j]*f_Tdt[j]/(S_tau1 - S_tau2)
    } else { terms_in_int[j] <- 0 }
  }
  PI_estimate <- sum(PI_terms_in_int)
  
  
  out <- list(OS_est = OS_est, OS_var = variance_est, CI_OS = CI_OS, PI_estimate = PI_estimate)
  
  return(out)
}



cox_estimate <- function(DATA){
  covariates <- DATA[, !colnames(DATA) %in% c("eventtime", "status", "Cstatus", "trt")]
  covariate_terms <- paste0(colnames(covariates), collapse = " + ")
  formula_str <- paste("Surv(eventtime, status) ~ trt + ", covariate_terms)
  formula_str <- as.formula(formula_str)
  
  
  cox_fit <- coxph(formula_str, data = DATA)
  cox_est <- coef(summary(cox_fit))["trt", "coef"]
  cox_se_est <- coef(summary(cox_fit))["trt", "se(coef)"]
  
  CI <- c(cox_est - 1.96*cox_se_est, cox_est + 1.96*cox_se_est)
  
  out <- list(cox_est = cox_est, cox_var = cox_se_est^2, CI_cox = CI)
  return(out)
}




############## RUN SIM #############


run_sim <- function(true_data, N, n, mu, tau1, tau2, true){
  
  cox_estimates <- numeric(N)
  variance_estimates_cox <- numeric(N)
  cox_CI <- matrix(NA, nrow = N, ncol = 2)
  
  OS_estimates_cox <- numeric(N)
  variance_estimates_OS_cox <- numeric(N)
  OS_CI_cox <- matrix(NA, nrow = N, ncol = 2)
  
  OS_estimates_additive <- numeric(N)
  variance_estimates_OS_additive <- numeric(N)
  OS_CI_additive <- matrix(NA, nrow = N, ncol = 2)
  
  PI_estimates_cox <- numeric(N)
  PI_estimates_additive <- numeric(N)
  
  for (m in 1:N){
    if (true_data == "cox") {data <- sim_data_cox(n, mu, true)} 
    else if (true_data == "additive") {data <- sim_data_additive(n, mu, true)} 
    else {stop("Invalid true_data type. Choose 'cox' or 'additive'.")}
    
    cox_estimate_return <- cox_estimate(data)
    OS_estimate_return_cox <- One_Step_estimate_cox(data, tau1, tau2, mu)
    OS_estimate_return_additive <- One_Step_estimate_additive(data, tau1, tau2, mu)
    
    cox_estimates[m] <- cox_estimate_return$cox_est
    variance_estimates_cox[m] <- cox_estimate_return$cox_var
    cox_CI[m, ] <- cox_estimate_return$CI_cox
    
    OS_estimates_cox[m] <- OS_estimate_return_cox$OS_est
    variance_estimates_OS_cox[m] <- OS_estimate_return_cox$OS_var
    OS_CI_cox[m, ] <- OS_estimate_return_cox$CI_OS
    
    OS_estimates_additive[m] <- OS_estimate_return_additive$OS_est
    variance_estimates_OS_additive[m] <- OS_estimate_return_additive$OS_var
    OS_CI_additive[m, ] <- OS_estimate_return_additive$CI_OS
    
    PI_estimates_cox[m] <- OS_estimate_return_cox$PI_estimate
    PI_estimates_additive[m] <- OS_estimate_return_additive$PI_estimate
    
  }
  
  out <- list(cox_estimates = cox_estimates, OS_estimates_cox = OS_estimates_cox, OS_estimates_additive = OS_estimates_additive, PI_estimates_cox = PI_estimates_cox, PI_estimates_additive = PI_estimates_additive, 
              variance_estimates_cox = variance_estimates_cox, variance_estimates_OS_cox = variance_estimates_OS_cox, variance_estimates_OS_additive = variance_estimates_OS_additive,
              cox_CI = cox_CI, OS_CI_cox = OS_CI_cox, OS_CI_additive = OS_CI_additive)
  return(out)
}


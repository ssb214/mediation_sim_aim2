#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Install packages 

userLib <-  "/projects/emerge/SSHOFFMAN/R"
.libPaths(userLib)

packages <- c("here",
              "Matrix",
              "lme4",
              "tidyverse",
              "quantreg",
              "MatrixModels",
              "doParallel",
              "clusterGeneration",
              "MASS",
              "hdmed",
              "HIMA")




# for (package in packages) {
#   if (!require(package, character.only=T, quietly=T)) {
#     install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
#   }
# }

for (package in packages) {
  library(package, character.only=T)
}

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("qvalue")

library(qvalue)

RNGkind("L'Ecuyer-CMRG")

num_cores <- as.numeric(args[3])

#### Functions #### 
# sim_id = 1
# sample_size = 100
# cov_n = 2 # number of covariates
# p = 10 # number of mediators 
# p_sig = 2 # the number of p features that are true mediators
# med_beta = 0.1 # the beta value of the metabolites that are mediators for both exposure and outcome 
# independent = "no" # Are the mediators independent or correlated 

data_generation <- function(sim_id, 
                            sample_size, 
                            cov_n,
                            p,
                            p_sig,
                            med_beta,
                            independent = "no"){
  set.seed(sim_id)
  
  library(clusterGeneration)
  library(MASS)
  
  # Sample Size
  n = sample_size
  
  # Covariates
  sigma <- diag(1, nrow = cov_n, ncol = cov_n) # std error 
  
  c <- mvtnorm::rmvnorm(n, mean = rep(0, cov_n), sigma = sigma)
  
  # Exposure ~ N(0,1) 
  beta <- rep(0.2, cov_n)
  e <- rnorm(n, mean = 0 + c%*%beta, sd = 1)
  
  # Mediators
  
  # Mediator means 
  med_means <- rep(0, p)
  
  ## Correlated mediators generation 
  if(independent == "no"){
    # Creating variance covariance matrix 
    R <- genPositiveDefMat(dim = p,
                           covMethod = "unifcorrmat", 
                           rangeVar = rep(1,p))$Sigma
  }
  
  ## independent mediators generation 
  if(independent == "yes"){
    R <- diag(1, nrow = p, ncol = p)
  }
  
  # This creates a matrix of variables using a multivariate normal dist. 
  mediators_corr <- MASS::mvrnorm(n, 
                                  mu = med_means, 
                                  Sigma = R)
  
  
  # Check the correlation matrix and compare this to mediators_ind
  # cor <- cor(mediators_corr)
  
  design_matrix <- matrix(1, nrow = n, ncol = (p-p_sig))
  
  if(p_sig > 0){
    add_cols <- matrix(rep(e, p_sig), ncol = p_sig) # inducing dependence on exposure
    design_matrix <- cbind(add_cols, design_matrix)
  }
  
  ex_med_means <- rep(0, p)
  if(p_sig > 0){
    ex_med_means[1:p_sig] <- med_beta
  }
  
  mediators <- mediators_corr + design_matrix%*%diag(ex_med_means)
  
  # Testing to see if the mediators are actually associated with e as expected 
  # test <- lm(mediators[,3] ~ e)
  
  # Outcome 
  
  # error term 
  e_y <- matrix(rnorm(n, mean = 0, sd = 1), n, 1)
  
  beta_m <- matrix(rep(0, p))
  
  if (p_sig > 0) {
    beta_m[1:p_sig, ] <- med_beta
  }
  
  beta_c <- rep(0.01, cov_n)
  
  y <- 0.8*e + mediators%*%(beta_m) + c%*%(beta_c) + e_y
  
  data <- data.frame(Exposure = e,
                     Outcome = y, 
                     C1 = c[,1],
                     C2 = c[,2],
                     M = mediators)
  return(data)
  
}

# sim_id = 1
# cov_n = 2
# p = 20
# p_sig = 2
# med_beta = 0.1
# independent = "no"
# exposure = 1
# exposure_m = 1

true_value_function <- function(sim_id, 
                                cov_n,
                                p,
                                p_sig,
                                med_beta,
                                independent = "no",
                                exposure, # Force exposure 
                                exposure_m # Force mediator under varying exposure  
){
  set.seed(sim_id)
  
  library(clusterGeneration)
  library(MASS)
  
  # Sample Size
  n = 10000
  
  # Covariates
  sigma <- diag(1, nrow = cov_n, ncol = cov_n) # std error 
  
  c <- mvtnorm::rmvnorm(n, mean = rep(0, cov_n), sigma = sigma)
  
  # Exposure ~ N(0,1) 
  beta <- rep(0.2, cov_n)
  e <- rep(exposure, times = n)
  
  # Mediators
  
  # Mediator means 
  med_means <- rep(0, p) 
  
  
  ## Correlated mediators generation 
  if(independent == "no"){
    # Creating variance covariance matrix 
    R <- genPositiveDefMat(dim = p,
                           covMethod = "unifcorrmat", 
                           rangeVar = rep(1,p))$Sigma
  }
  
  ## independent mediators generation 
  if(independent == "yes"){
    R <- diag(1, nrow = p, ncol = p)
  }
  
  # This creates a matrix of variables using a multivariate normal dist. 
  mediators_corr <- MASS::mvrnorm(n, 
                                  mu = med_means, 
                                  Sigma = R)
  
  
  # Check the correlation matrix and compare this to mediators_ind
  cor <- cor(mediators_corr)
  
  design_matrix <- matrix(1, nrow = n, ncol = (p-p_sig))
  
  if(p_sig > 0){
    e_m <- rep(exposure_m, times = n)
    add_cols <- matrix(rep(e_m, p_sig), ncol = p_sig) # inducing dependence on exposure
    design_matrix <- cbind(add_cols, design_matrix)
  }
  
  ex_med_means <- rep(0, p)
  if(p_sig > 0){
    ex_med_means[1:p_sig] <- med_beta
  }
  
  mediators <- mediators_corr + design_matrix%*%diag(ex_med_means) 
  
  # Outcome 
  
  # error term 
  e_y <- matrix(rnorm(n, mean = 0, sd = 1), n, 1)
  
  beta_m <- matrix(rep(0, p))
  
  if (p_sig > 0) {
    beta_m[1:p_sig, ] <- med_beta
  }
  
  beta_c <- rep(0.01, cov_n)
  
  y <- 0.8*e + mediators%*%(beta_m) + c%*%(beta_c) + e_y
  
  data <- data.frame(Exposure = e,
                     Outcome = y, 
                     C1 = c[,1],
                     C2 = c[,2],
                     M = mediators)
  return(data)
  
}




#### Generating Truths #### 

##### Setting parameters #####

p = as.numeric(args[1]) #### NEEDS TO MATCH DATA GENERATION 

# THESE NEED TO MATCH THE SIMULATION 

parm_data_truth <- expand.grid(
  p_sig = ceiling(c(0.02*p, 0.05*p, 0.1*p)),
  med_beta = c(0.1, 0.3),
  independent = c("yes", "no"),
  cov_n = 2
)



##### Empty dataframe #####

results_truth <- data.frame(matrix(ncol = 7))
colnames(results_truth) <- c("p", "truth","independent", "cov_n", "med_beta", 
                             "TIE", "CIE")

##### Truth Generation #####
truths <- mclapply(1:nrow(parm_data_truth), function(x) {
  # Get parameters 
  cov_n <- parm_data_truth[x, "cov_n"]
  p_sig <- parm_data_truth[x, "p_sig"]
  independent <- parm_data_truth[x, "independent"]
  med_beta <- parm_data_truth[x, "med_beta"]
  p <- p  
  
  ##### Generating the truth #####
  d1 <- true_value_function(sim_id = 1,
                            cov_n = cov_n,
                            p = p,
                            p_sig = p_sig,
                            med_beta = med_beta,
                            exposure = 1,
                            exposure_m = 1,
                            independent = independent)
  d0 <- true_value_function(sim_id = 1,
                            cov_n = cov_n,
                            p = p,
                            p_sig = p_sig,
                            med_beta = med_beta,
                            exposure = 1,
                            exposure_m = 0,
                            independent = independent)
  
  u_e1m1 <- sum(d1$Outcome)/nrow(d1)
  u_e1m0 <- sum(d0$Outcome)/nrow(d0)
  
  tie <- u_e1m1 - u_e1m0
  
  if(independent == "yes") {
    cie <- med_beta*med_beta
  } else {
    cie <- NA
  }
  
  out <- data.frame(
    p =  p,
    truth = p_sig,
    independent = independent, 
    cov_n = cov_n, 
    med_beta = med_beta,
    TIE = tie, 
    CIE = cie
  )
  
  results_truth <- rbind(results_truth, out)
  
  # Return the results
  return(list(results_truth = results_truth))
  
},
mc.cores = num_cores
)

for (res in truths) {
  results_truth <- rbind(results_truth, res$results_truth)
}
results_truth <- results_truth[complete.cases(results_truth$TIE), ]

# save results in csv

output_dir <- "/projects/emerge/SSHOFFMAN/Aim2/output2/"

write.csv(results_truth, file = paste0(output_dir, "truth_", p, ".csv"))


#### Simulation ####

# Need to create one giant function using lapply for all three simulations 

# Need to make sure to capture the simulation parameters in the output 

# parameter specification can occur outside the lapply function 

##### Libraries #####
library(HIMA)
library(qvalue)
library(hdmed)

##### Set seed #####

seed_data <- 1:1000
array_num <- as.numeric(args[2])

# this code works for 5 arrays - the denominator must match the number of arrays you select
start_ <- array_num*(length(seed_data)/5) - ((length(seed_data)/5) - 1)
end_ <- array_num*(length(seed_data))/5

start_
end_

##### Setting parameters #####

p = as.numeric(args[1]) #### NEEDS TO MATCH DATA GENERATION 

## p_sig is 2%, 5%, and 10% of p

parm_data <- expand.grid(
  num_sims = start_:end_,
  n = c(200, 500, 1000),
  p_sig = ceiling(c(0.02*p, 0.05*p, 0.1*p)),
  med_beta = c(0.1, 0.3),
  independent = c("yes", "no"),
  cov_n = 2
)


##### Empty dataframes for results #####

# dataframes for results 

## HIMA
results_HIMA <- data.frame(matrix(ncol = 16))
colnames(results_HIMA) <- c("alpha", "beta", "gamma", "alpha*beta", "Bonferroni.p", 
                            "BH.FDR", "% total effect", "sim_id", "captured.mets",
                            "indirect", "n", "p", "truth", "independent", 
                            "cov_n", "med_beta")

sesp_HIMA <- data.frame(matrix(ncol = 7+p))
colnames(sesp_HIMA) <- c("sim_id", "n", "p", "truth", "independent", 
                         "cov_n", "med_beta",
                         paste0("M.", as.character(1:p)))

## HDMA
results_HDMA <- data.frame(matrix(ncol = 17))
colnames(results_HDMA) <- c("mediator", "alpha", "alpha_pv", "beta", 
                            "beta_pv", "alpha_beta", "ab_pv", "sim_id",
                            "indirect", "direct", "total", 
                            "n", "p", "truth", "independent", "cov_n", "med_beta")

sesp_HDMA <- data.frame(matrix(ncol = 7+p))
colnames(sesp_HDMA) <- c("sim_id", "n", "p", "truth", "independent", 
                         "cov_n", "med_beta",
                         paste0("M.", as.character(1:p)))

## MITM

results_MITM <- data.frame(matrix(ncol = 16))
colnames(results_MITM) <- c("metabolite", "alpha", "alpha.p", "beta", 
                            "beta.p", "gamma", "alphabeta", "per.totaleffect",
                            "tie", "sim_id", "n", "p", "truth", 
                            "independent", "cov_n", "med_beta")

sesp_MITM <- data.frame(matrix(ncol = 7+p))
colnames(sesp_MITM) <- c("sim_id", "n", "p", "truth", "independent", 
                         "cov_n", "med_beta",
                         paste0("M.", as.character(1:p)))

counter <- 0


##### Simulation function #####

simulation_results <- mclapply(1:nrow(parm_data), function(x) {
  
  counter <<- counter + 1 
  
  # Get parameters 
  num_sims <- parm_data[x, "num_sims"]
  n <- parm_data[x, "n"]
  cov_n <- parm_data[x, "cov_n"]
  p_sig <- parm_data[x, "p_sig"]
  independent <- parm_data[x, "independent"]
  med_beta <- parm_data[x, "med_beta"]
  p <- p  # Change to match data generation 
  
  ###### Data generation ######
  my_data <- data_generation(
    sim_id = num_sims,
    sample_size = n,
    cov_n = cov_n,
    p = p,
    p_sig = p_sig,
    med_beta = med_beta,
    independent = independent
  )
  
  range <- 3+cov_n
  # Pulling out mediating set 
  mediators <- my_data[range:ncol(my_data)]
  
  # For HDMA
  covs <- as.matrix(my_data[3:4]) ### CHANGE THIS TO MATCH cov_n ###
  
  # for MITM
  nmets = ncol(mediators)
  
  sim_id = parm_data[x,]$num_sims
  sample_size = parm_data[x,]$n
  
  
  ###### HIMA ######
  out <- hima(
    X = my_data$Exposure,
    Y = my_data$Outcome,
    M = mediators,
    COV.XM = my_data[, c("C1", "C2")],
    COV.MY = my_data[, c("C1", "C2")],
    Y.family = 'gaussian', 
    M.family = 'gaussian',
    penalty = 'MCP',
    parallel = F,
    topN = NULL, 
    scale = T,
    verbose = F
  )
  
  # Need to create a null row of results to r bind to the list and capture the null findings in the output 
  # Or pull the captured metabolites and the TIE
  if(is.null(out)){
    out <- data.frame(matrix(NA, nrow = 1, ncol = 16))
    colnames(out) <- colnames(results_HIMA)
    out$sim_id <- rep(sim_id, nrow(out))
    out$n <- rep(sample_size, nrow(out))
    out$truth <- rep(p_sig, nrow(out))
    out$independent <- rep(independent, nrow(out))
    out$cov_n <- rep(cov_n, nrow(out))
    out$p <- rep(p, nrow(out))
    out$med_beta <- rep(med_beta, nrow(out))
  } else {
    out$sim_id <- rep(sim_id, nrow(out))
    out$captured.mets <- row.names(out)
    out$indirect <- rep(sum(out$`alpha*beta`, na.rm = T), nrow(out))
    out$n <- rep(sample_size, nrow(out))
    out$truth <- rep(p_sig, nrow(out))
    out$independent <- rep(independent, nrow(out))
    out$cov_n <- rep(cov_n, nrow(out))
    out$p <- rep(p, nrow(out))
    out$med_beta <- rep(med_beta, nrow(out))
  }
  
  
  # Storing results 
  results_HIMA <- rbind(results_HIMA, out)
  
  
  
  # Se/Sp results 
  mets <- row.names(out)
  true_mets <- paste0("M.", 1:p_sig)
  
  # Create a matrix to store the values
  out_matrix <- matrix(0, nrow = 1, ncol = 7 + p)
  colnames(out_matrix) <- c("sim_id", "n", "p", "truth", "independent", 
                            "cov_n", "med_beta",
                            paste0("M.", as.character(1:p)))
  
  # Set sim_id
  out_matrix[1, "sim_id"] <- sim_id
  # Set n 
  out_matrix[1, "n"] <- sample_size
  # Set p
  out_matrix[1, "p"] <- p
  # Set truth column
  out_matrix[1, "truth"] <- p_sig
  # Set independence 
  out_matrix[1, "independent"] <- independent
  # Set covs
  out_matrix[1, "cov_n"] <- cov_n
  # set mediator betas 
  out_matrix[1, "med_beta"] <- med_beta
  
  # Set columns for each M. value
  for (i in 1:p) {
    out_matrix[1, paste0("M.", i)] <- as.integer(paste0("M.", i) %in% mets)
  }
  
  # Convert to a data frame
  out_sesp <- as.data.frame(out_matrix, stringsAsFactors = FALSE)
  
  # Storing results 
  sesp_HIMA <- rbind(sesp_HIMA, out_sesp)
  
  ###### HDMA ######
  
  # HDMA 
  out <- mediate_hdma(
    A = my_data$Exposure,
    Y = my_data$Outcome,
    M = mediators,
    C1 = covs, # covariates for the outcome model 
    C2 = covs, # covariates for the mediator model 
    binary_y = F,
    n_include = NULL # ceiling(n/log(n))
  )
  
  out_con <- out$contributions
  out_ef <- out$effects
  
  # Need to create a null row of results to r bind to the list and capture the null findings in the output 
  if(is.null(out_con)){
    out_con <- data.frame(matrix(NA, nrow = 1, ncol = 17))
    colnames(out_con) <- colnames(results_HDMA)
    out_con$sim_id <- rep(sim_id, nrow(out_con))
    out_con$n <- rep(sample_size, nrow(out_con))
    out_con$p <- rep(p, nrow(out_con))
    out_con$truth <- rep(p_sig, nrow(out_con))
    out_con$independent <- rep(independent, nrow(out_con))
    out_con$cov_n <- rep(cov_n, nrow(out_con))
    out_con$med_beta <- rep(med_beta, nrow(out_con))
  } else {
    out_con$sim_id <- rep(sim_id, nrow(out_con))
    out_con$indirect <- rep(out_ef[1,2], nrow(out_con))
    out_con$direct <- rep(out_ef[2,2], nrow(out_con))
    out_con$total <- rep(out_ef[3,2], nrow(out_con))
    out_con$n <- rep(sample_size, nrow(out_con))
    out_con$p <- rep(p, nrow(out_con))
    out_con$truth <- rep(p_sig, nrow(out_con))
    out_con$independent <- rep(independent, nrow(out_con))
    out_con$cov_n <- rep(cov_n, nrow(out_con))
    out_con$med_beta <- rep(med_beta, nrow(out_con))
  }
  
  # Storing results 
  results_HDMA <- rbind(results_HDMA, out_con)
  
  
  # Se/Sp results 
  mets <- out_con$mediator
  true_mets <- paste0("M.", 1:p_sig)
  
  # Create a matrix to store the values
  out_matrix <- matrix(0, nrow = 1, ncol = 7 + p)
  colnames(out_matrix) <- c("sim_id", "n", "p", "truth", "independent", 
                            "cov_n", "med_beta",
                            paste0("M.", as.character(1:p)))
  
  # Set sim_id
  out_matrix[1, "sim_id"] <- sim_id
  # Set n 
  out_matrix[1, "n"] <- sample_size
  # Set p
  out_matrix[1, "p"] <- p
  # Set truth column
  out_matrix[1, "truth"] <- p_sig
  # Set independence 
  out_matrix[1, "independent"] <- independent
  # Set covs
  out_matrix[1, "cov_n"] <- cov_n
  # set mediator betas 
  out_matrix[1, "med_beta"] <- med_beta
  
  # Set columns for each M. value
  for (i in 1:p) {
    out_matrix[1, paste0("M.", i)] <- as.integer(paste0("M.", i) %in% mets)
  }
  
  # Convert to a data frame
  out_sesp <- as.data.frame(out_matrix, stringsAsFactors = FALSE)
  
  # Storing results 
  sesp_HDMA <- rbind(sesp_HDMA, out_sesp)
  
  ###### MITM ######
  
  ### Exposure - metabolites 
  
  dataout <- data.frame(matrix(nrow = nmets, ncol = 5))
  colnames(dataout) <- c("e_met", "e_estimate", "e_stderror", "e_t", "e_p")
  
  for (i in 1:nmets) {
    tryCatch(
      {
        
        lmfit <- glm(my_data[, i+4] ~ my_data$Exposure +
                       my_data$C1 +
                       my_data$C2,
                     family = gaussian())
        
        dataout[i,2]<-summary(lmfit)$coefficients[2,1]
        dataout[i,3]<-summary(lmfit)$coefficients[2,2]
        dataout[i,4]<-summary(lmfit)$coefficients[2,3]
        dataout[i,5]<-summary(lmfit)$coefficients[2,4]
        dataout[i,1] <- names(my_data)[i + 4]
        
        if(i %% 100 == 0){print(i)}
      }, error=function(e){})
  }
  
  met_ex <- dataout
  met_ex$e_fdr <- p.adjust(met_ex$e_p, method = "BH")
  
  ### Outcome - metabolites 
  
  dataout <- data.frame(matrix(nrow = nmets, ncol = 5))
  colnames(dataout) <- c("o_met", "o_estimate", "o_stderror", "o_t", "o_p")
  
  for (i in 1:nmets) {
    tryCatch(
      {
        
        lmfit <- glm(my_data[, i+4] ~ my_data$Outcome +
                       my_data$Exposure +
                       my_data$C1 +
                       my_data$C2,
                     family = gaussian())
        
        dataout[i,2]<-summary(lmfit)$coefficients[2,1]
        dataout[i,3]<-summary(lmfit)$coefficients[2,2]
        dataout[i,4]<-summary(lmfit)$coefficients[2,3]
        dataout[i,5]<-summary(lmfit)$coefficients[2,4]
        dataout[i,1] <- names(my_data)[i + 4]
        
        if(i %% 100 == 0){print(i)}
      }, error=function(e){})
  }
  
  met_out <- dataout
  
  met_out$o_fdr <- p.adjust(met_out$o_p, method = "BH")
  
  ### Exposure - outcome 
  
  ex_out <- glm(my_data$Outcome ~ my_data$Exposure +
                  my_data$C1 +
                  my_data$C2,
                family = gaussian())
  gamma_est <- summary(ex_out)$coefficients[2,1]
  
  ### Finding overlap 
  
  overlapfeats <- cbind(met_ex, met_out)
  overlapfeats <- subset(overlapfeats, 
                         (e_fdr<0.2 & o_fdr<0.2))
  
  ab_est <- overlapfeats$e_estimate*overlapfeats$o_estimate
  te_est <- ab_est/gamma_est*100
  indirect <- sum(ab_est)
  
  if(nrow(overlapfeats) == 0){
    out_mitm = NULL
  } else {
    out_mitm <- data.frame(
      metabolite = overlapfeats$e_met,
      alpha = overlapfeats$e_estimate,
      alpha.p = overlapfeats$e_fdr,
      beta = overlapfeats$o_estimate,
      beta.p = overlapfeats$o_fdr,
      gamma = gamma_est,
      alphabeta = ab_est,
      per.totaleffect = te_est,
      tie = indirect
    )
  }
  
  
  # Need to create a null row of results to r bind to the list and capture the null findings in the output 
  if(is.null(out_mitm)){
    out_mitm <- data.frame(matrix(NA, nrow = 1, ncol = 16))
    colnames(out_mitm) <- colnames(results_MITM)
    out_mitm$sim_id <- rep(sim_id, nrow(out_mitm))
    out_mitm$n <- rep(n, nrow(out_mitm))
    out_mitm$p <- rep(p, nrow(out_mitm))
    out_mitm$truth <- rep(p_sig, nrow(out_mitm))
    out_mitm$independent <- rep(independent, nrow(out_mitm))
    out_mitm$cov_n <- rep(cov_n, nrow(out_mitm))
    out_mitm$med_beta <- rep(med_beta, nrow(out_mitm))
  } else {
    out_mitm$sim_id <- rep(sim_id, nrow(out_mitm))
    out_mitm$n <- rep(n, nrow(out_mitm))
    out_mitm$p <- rep(p, nrow(out_mitm))
    out_mitm$truth <- rep(p_sig, nrow(out_mitm))
    out_mitm$independent <- rep(independent, nrow(out_mitm))
    out_mitm$cov_n <- rep(cov_n, nrow(out_mitm))
    out_mitm$med_beta <- rep(med_beta, nrow(out_mitm))
  }
  
  # Storing results 
  results_MITM <- rbind(results_MITM, out_mitm)
  
  
  
  # Se/Sp results 
  mets <- out_mitm$metabolite
  true_mets <- paste0("M.", 1:p_sig)
  
  # Create a matrix to store the values
  out_matrix <- matrix(0, nrow = 1, ncol = 7 + p)
  colnames(out_matrix) <- c("sim_id", "n", "p", "truth", "independent", 
                            "cov_n", "med_beta",
                            paste0("M.", as.character(1:p)))
  
  # Set sim_id
  out_matrix[1, "sim_id"] <- sim_id
  # Set n 
  out_matrix[1, "n"] <- sample_size
  # Set p
  out_matrix[1, "p"] <- p
  # Set truth column
  out_matrix[1, "truth"] <- p_sig
  # Set independence 
  out_matrix[1, "independent"] <- independent
  # Set covs
  out_matrix[1, "cov_n"] <- cov_n
  # set mediator betas 
  out_matrix[1, "med_beta"] <- med_beta
  
  # Set columns for each M. value
  for (i in 1:p) {
    out_matrix[1, paste0("M.", i)] <- as.integer(paste0("M.", i) %in% mets)
  }
  
  # Convert to a data frame
  out_sesp <- as.data.frame(out_matrix, stringsAsFactors = FALSE)
  
  # Storing results 
  sesp_MITM <- rbind(sesp_MITM, out_sesp)
  
  
  # Return the results
  return(list(results_HDMA = results_HDMA, 
              sesp_HDMA = sesp_HDMA,
              results_HIMA = results_HIMA,
              sesp_HIMA = sesp_HIMA,
              results_MITM = results_MITM,
              sesp_MITM = sesp_MITM
  ))
},
mc.cores = num_cores)

# Pulling together results 

for (res in simulation_results) {
  results_HIMA <- rbind(results_HIMA, res$results_HIMA)
  sesp_HIMA <- rbind(sesp_HIMA, res$sesp_HIMA)
  results_HDMA <- rbind(results_HDMA, res$results_HDMA)
  sesp_HDMA <- rbind(sesp_HDMA, res$sesp_HDMA)
  results_MITM <- rbind(results_MITM, res$results_MITM)
  sesp_MITM <- rbind(sesp_MITM, res$sesp_MITM)
}
results_HIMA <- results_HIMA[complete.cases(results_HIMA$sim_id), ]
sesp_HIMA <- sesp_HIMA[complete.cases(sesp_HIMA$sim_id), ]
results_HDMA <- results_HDMA[complete.cases(results_HDMA$sim_id), ]
sesp_HDMA <- sesp_HDMA[complete.cases(sesp_HDMA$sim_id), ]
results_MITM <- results_MITM[complete.cases(results_MITM$sim_id), ]
sesp_MITM <- sesp_MITM[complete.cases(sesp_MITM$sim_id), ]

# save results in csv

output_dir <- "/projects/emerge/SSHOFFMAN/Aim2/output2/"

## HIMA 

if(file.exists(paste0(output_dir, "HIMA", p, "_res_", array_num, ".csv"))){
  write_csv(results_HIMA, 
            file = paste0(output_dir, "HIMA", p, "_res_", array_num, ".csv"),
            append = T)
} else {
  write_csv(results_HIMA, 
            file = paste0(output_dir, "HIMA", p, "_res_", array_num, ".csv"),
            append = F)
}

if(file.exists(paste0(output_dir, "HIMA", p, "_sesp_", array_num, ".csv"))){
  write_csv(sesp_HIMA, 
            file = paste0(output_dir, "HIMA", p, "_sesp_", array_num, ".csv"),
            append = T)
} else {
  write_csv(sesp_HIMA, 
            file = paste0(output_dir, "HIMA", p, "_sesp_", array_num, ".csv"),
            append = F)
}

## HDMA

if(file.exists(paste0(output_dir, "HDMA", p, "_res_", array_num, ".csv"))){
  write_csv(results_HDMA, 
            file = paste0(output_dir, "HDMA", p, "_res_", array_num, ".csv"),
            append = T)
} else {
  write_csv(results_HDMA, 
            file = paste0(output_dir, "HDMA", p, "_res_", array_num, ".csv"),
            append = F)
}

if(file.exists(paste0(output_dir, "HDMA", p, "_sesp_", array_num, ".csv"))){
  write_csv(sesp_HDMA, 
            file = paste0(output_dir, "HDMA", p, "_sesp_", array_num, ".csv"),
            append = T)
} else {
  write_csv(sesp_HDMA, 
            file = paste0(output_dir, "HDMA", p, "_sesp_", array_num, ".csv"),
            append = F)
}

## MITM

if(file.exists(paste0(output_dir, "MITM", p, "_res_", array_num, ".csv"))){
  write_csv(results_MITM, 
            file = paste0(output_dir, "MITM", p, "_res_", array_num, ".csv"),
            append = T)
} else {
  write_csv(results_MITM, 
            file = paste0(output_dir, "MITM", p, "_res_", array_num, ".csv"),
            append = F)
}

if(file.exists(paste0(output_dir, "MITM", p, "_sesp_", array_num, ".csv"))){
  write_csv(sesp_MITM, 
            file = paste0(output_dir, "MITM", p, "_sesp_", array_num, ".csv"),
            append = T)
} else {
  write_csv(sesp_MITM, 
            file = paste0(output_dir, "MITM", p, "_sesp_", array_num, ".csv"),
            append = F)
}


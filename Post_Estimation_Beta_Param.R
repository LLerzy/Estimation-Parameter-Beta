###########################################################################
#                                                                         #
# Description:                                                            #
# This document contains a set of routines aimed at estimating the shape  #
# parameters of the Beta distribution for the variable X from a Bayesian  #
# perspective. Although a new bivariate prior distribution and its        #
# respective method for simulating random samples are utilized, the main  #
# focus of the code is to develop a simulation study to monitor the       #
# behavior of estimates obtained for these parameters. The routines       #
# include empirical Bayesian and subjective approaches for hyperparameter #
# estimation, using both bootstrap intervals and expert intervals         #
# established with different biases and semi-amplitudes.                  #
#                                                                         #
# The primary procedures in this document include:                        #
# 1. Configuration of empirical and subjective Bayesian approaches for    #
#    hyperparameter estimation.                                           #
# 2. Calculation of prior moments based on the given hyperparameters.     #
# 3. Posterior estimation using importance sampling for the shape         #
#    parameters of the beta distribution.                                 #
# 4. Conducting simulation studies to compare and monitor the performance #
#    of posterior estimates under different scenarios, sample sizes, and  #
#    sets of hyperparameters.                                             #
# 5. Visualization of results and comparison of hyperparameters.          #
#                                                                         #
# Author: Llerzy Torres Ome                                               #
# Creation Date: August 06, 2024                                            #
#                                                                         #
# Notes:                                                                  #
# 1. It is recommended to review and adapt each section according to the  #
#    specific needs of each analysis.                                     #
# 2. Ensure you understand each step before executing it to guarantee     #
#    accurate results and avoid potential errors.                         #
# 3. For any questions or suggestions, contact                            #
#    llerzy.torres@correounivalle.edu.co                                  #
#                                                                         #
###########################################################################

source(file = "requiredfunctions.R")

# Define the three scenarios for the parameters alpha and beta of the beta distribution of the variable X.
Vect_Alpha = c(0.5, 16, 3)
Vect_Beta = c(0.5, 4, 12)
Vect_mu = Vect_Alpha / (Vect_Alpha + Vect_Beta)  # Mean of the beta distribution of X.
Vect_var = Vect_Alpha * Vect_Beta / ((Vect_Alpha + Vect_Beta)^2 * (Vect_Alpha + Vect_Beta + 1))  # Variance of the beta distribution of X.
Vect_cv = sqrt(Vect_var) / Vect_mu  # Coefficient of variation of the beta distribution of X.

# For the subjective approach, the quantile intervals for the expert are specified,
# considering two bias values for the mean and the coefficient of variation, and one value for the semi-amplitude 
# of the mean and variance.
Vect_epsion_mu = c(0.01, 0.05)  # Bias of the midpoint of the expert interval for mu.
Vect_delta_mu = c(0.07)  # Semi-amplitude of the expert interval for mu.
Vect_epsion_cv = c(0.02, 0.06)  # Bias of the midpoint of the expert interval for cv.
Vect_delta_cv = c(0.09)  # Semi-amplitude of the expert interval for cv.

#######
## Settings for the subjective approach.
#######
# Matrix of all combinations between the scenarios of the alpha and beta parameters, biases and amplitudes of the mean and the coefficient of variation.
Param = matrix(data = NA, nrow = 10*3, ncol = 16, dimnames = list(NULL, c("Scenario", "Alpha", "Beta", "EpsionMu", "DeltaMu", "EpsionCv", "DeltaCv", "Va", "Vb", "Vc", "Vd", "Method", "LMu", "UMu", "LCV", "UCV")))
Param = as.data.frame(Param)

# Assigning scenarios and parameter values
Param$Scenario = c(rep(1, 10), rep(2, 10), rep(3, 10))
Param$Alpha = c(rep(Vect_Alpha[1], 10), rep(Vect_Alpha[2], 10), rep(Vect_Alpha[3], 10))
Param$Beta = c(rep(Vect_Beta[1], 10), rep(Vect_Beta[2], 10), rep(Vect_Beta[3], 10))
Param$EpsionMu = c(0, 0, rep(Vect_epsion_mu[1], 4), rep(Vect_epsion_mu[2], 4))
Param$DeltaMu = Vect_delta_mu
Param$EpsionCv = c(0, 0, rep(Vect_epsion_cv[1], 2), rep(Vect_epsion_cv[2], 2), rep(Vect_epsion_cv[1], 2), rep(Vect_epsion_cv[2], 2))
Param$DeltaCv = Vect_delta_cv
Param$Method = c("BM", "BT", "EM1", "ET1", "EM2", "ET2", "EM3", "ET3", "EM4", "ET4")

# Limits of the expert quantile intervals for the mean of X, considering the biases and semi-amplitude.
Param$LMu = Param$Alpha / (Param$Alpha + Param$Beta) + Param$EpsionMu - Param$DeltaMu
Param$UMu = Param$Alpha / (Param$Alpha + Param$Beta) + Param$EpsionMu + Param$DeltaMu

# Limits of the expert quantile intervals for the CV of X, considering the biases and semi-amplitude.
Param$LCV = sqrt(Param$Alpha * Param$Beta / ((Param$Alpha + Param$Beta)^2 * (Param$Alpha + Param$Beta + 1))) / (Param$Alpha / (Param$Alpha + Param$Beta)) + Param$EpsionCv - Param$DeltaCv
Param$UCV = sqrt(Param$Alpha * Param$Beta / ((Param$Alpha + Param$Beta)^2 * (Param$Alpha + Param$Beta + 1))) / (Param$Alpha / (Param$Alpha + Param$Beta)) + Param$EpsionCv + Param$DeltaCv


#######
## Configuration for the Empirical Bayes Approach
#######
# Suppose a sample of size 10 is available, from which the quantile intervals for mu and cv
# are obtained using the subjective approach.

# Initialize a matrix to store beta samples for different scenarios
sim_beta <- matrix(nrow = 10, ncol = length(Vect_Alpha), dimnames = list(NULL, c("Scenario1", "Scenario2", "Scenario3")))
set.seed(123)  # For reproducibility
for (i in 1:length(Vect_Alpha)) {
  sim_beta[, i] <- rbeta(10, shape1 = Vect_Alpha[i], shape2 = Vect_Beta[i])
}


sim_beta = as.data.frame(sim_beta)

#####
## Obtaining hyperparameters for the Empirical Bayes approach and the subjective approach
## considering expert intervals (with epsilon errors and amplitude 2delta).
#####
# The parameter n corresponds to the scenarios for the shape parameters of the beta distribution.
# When n=0, it represents scenario 1. For n=1, it represents scenario 2, and for n=2, it represents scenario 3.
# n_combi: takes values from 0 to 3 and represents the combinations between the epsilon errors of the mean and the cv.
# bound_var: parameter that indicates the way in which the upper limit for the conditional variance of the mean is established.
bound_var = "min"
BoostrapIntervalsMu = NULL
BoostrapIntervalsCv = NULL
BoostrapIntervalsVar = NULL

for (n in 0:2) {
  hip_BM = Hyperparameters(sim_beta[, n + 1], r_boostrap = 100, q_boostrap = c(0.025, 0.975), option_mu = "moments",
                           sig_mu = 0, bound_var = bound_var,
                           sig_var = 0.1, digits = 4, graphs_boot = F)  # Bootstrap interval and Moments method.
  
  hip_BT = Hyperparameters(sim_beta[, n + 1], r_boostrap = 100, q_boostrap = c(0.025, 0.975), option_mu = "tovar", 
                           sig_mu = 0.1, bound_var = bound_var,
                           sig_var = 0.1, digits = 4, graphs_boot = F)  # Bootstrap interval and Tovar method.
  
  BoostrapIntervalsMu = rbind(BoostrapIntervalsMu, hip_BM$Q_mean, hip_BT$Q_mean)
  BoostrapIntervalsCv = rbind(BoostrapIntervalsCv, hip_BM$Q_cv, hip_BT$Q_cv)
  BoostrapIntervalsVar = rbind(BoostrapIntervalsVar, hip_BM$Q_var, hip_BT$Q_var)
  
  Param[10 * n + 1, c(8, 9)] = hip_BM$hiper_mean
  Param[10 * n + 1, c(10, 11)] = c(hip_BM$hiper_var$a, hip_BM$hiper_var$b)  # Store in Param.
  
  Param[10 * n + 2, c(8, 9)] = c(hip_BT$hiper_mean$a, hip_BT$hiper_mean$b)
  Param[10 * n + 2, c(10, 11)] = c(hip_BT$hiper_var$a, hip_BT$hiper_var$b)  # Store in Param.
  
  rm(hip_BM, hip_BT)
  
  for (n_combi in 0:3) {
    hip_EM = Hyperparameters(sim_beta[, (n + 1)], r_boostrap = 0, q_boostrap = c(0.025, 0.975), option_mu = "moments",
                             sig_mu = 0, bound_var = bound_var,
                             sig_var = 0.1, digits = 4, graphs_boot = F, 
                             Q_E_mu = c(Param$LMu[10 * n + 3 + 2 * n_combi], Param$UMu[10 * n + 3 + 2 * n_combi]), 
                             Q_E_cv = c(Param$LCV[10 * n + 3 + 2 * n_combi], Param$UCV[10 * n + 3 + 2 * n_combi]))  # Expert interval and Moments method.
    
    hip_ET = Hyperparameters(sim_beta[, (n + 1)], r_boostrap = 0, q_boostrap = c(0.025, 0.975), option_mu = "tovar",
                             sig_mu = 0.1, bound_var = bound_var,
                             sig_var = 0.1, digits = 4, graphs_boot = F,
                             Q_E_mu = c(Param$LMu[10 * n + 4 + 2 * n_combi], Param$UMu[10 * n + 4 + 2 * n_combi]), 
                             Q_E_cv = c(Param$LCV[10 * n + 4 + 2 * n_combi], Param$UCV[10 * n + 4 + 2 * n_combi]))  # Expert interval and Tovar method.
    
    Param[10 * n + 3 + 2 * n_combi, c(8, 9)] = c(hip_EM$hiper_mean$a, hip_EM$hiper_mean$b)
    Param[10 * n + 3 + 2 * n_combi, c(10, 11)] = c(hip_EM$hiper_var$a, hip_EM$hiper_var$b)  # Store in Param.
    
    Param[10 * n + 4 + 2 * n_combi, c(8, 9)] = c(hip_ET$hiper_mean$a, hip_ET$hiper_mean$b)
    Param[10 * n + 4 + 2 * n_combi, c(10, 11)] = c(hip_ET$hiper_var$a, hip_ET$hiper_var$b)  # Store in Param.
    
    rm(hip_EM, hip_ET)
  }
}

# Verify that in scenario 1 the value obtained for parameter c is greater than 2 to ensure the mean and variance exist.
n = 0
Param[(10 * n + 1):(10 * n + 1 + 9), ]
# Param

# Convert the matrices containing the hyperparameter values obtained from the Empirical Bayes approach to data frames.
BoostrapIntervalsMu = as.data.frame(BoostrapIntervalsMu)
BoostrapIntervalsCv = as.data.frame(BoostrapIntervalsCv)
BoostrapIntervalsVar = as.data.frame(BoostrapIntervalsVar)
BoostrapIntervalsMu = cbind(BoostrapIntervalsMu, BoostrapIntervalsVar, BoostrapIntervalsCv)
names(BoostrapIntervalsMu) = c("LMu", "UMu", "LV", "UV", "LCv", "UCv")
BoostrapIntervalsMu$Method = c("BM", "BT")

###########################
# Obtaining the moments of the prior distributions defined by each set of hyperparameters
###########################
rm(i, Param1)
Measure_Prior_Distribution = NULL

for (i in 1:30) {
  Measure_Prior_Distribution = rbind(Measure_Prior_Distribution, Measure_Analy(Param$Va[i], Param$Vb[i], Param$Vc[i], Param$Vd[i], digits = 9))
}
Param1 = cbind(Param, Measure_Prior_Distribution)

Parameters_Min = Param1
Parameters_Min_boot = BoostrapIntervalsMu

## Repeat the procedure from line 80, modifying bound_var = "mean" and save the results as:
# Parameters_Mean = Param1
# Parameters_Mean_boot = BoostrapIntervalsMu

## Repeat the procedure from line 80, modifying bound_var = "max" and save the results as:
# Parameters_Max = Param1
# Parameters_Max_boot = BoostrapIntervalsMu

## Analyzing the moments of the prior distributions, select those that best represent the simulated scenario.
## Construct the matrix that contains the hyperparameters obtained with bound_var = "min" in scenario 2, hyperparameters generated by bound_var = "mean" in scenario 1, and
## hyperparameters found by bound_var = "max" in scenario 3.
Param = rbind(Parameters_Mean[which(Parameters_Mean$Scenario == 1), ],
              Parameters_Min[which(Parameters_Min$Scenario == 2), ],
              Parameters_Max[which(Parameters_Max$Scenario == 3), ])

BoostrapIntervals = rbind(Parameters_Mean_boot[c(1, 2), ],
                          Parameters_Min_boot[c(3, 4), ],
                          Parameters_Max_boot[c(5, 6), ])

#####
# Create a workbook object to save hyperparameters and moments generated for each scenario and bound_var specification.
##
titulo_xlsx = "ParametersSigMu01SigV01.xlsx"
wb <- createWorkbook()
# Add sheets and write data for each set of parameters and bootstrap intervals
addWorksheet(wb, "ParametersMin")
writeData(wb, sheet = "ParametersMin", x = Parameters_Min, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ParametersMinBoot")
writeData(wb, sheet = "ParametersMinBoot", x = Parameters_Min_boot, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ParametersMean")
writeData(wb, sheet = "ParametersMean", x = Parameters_Mean, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ParametersMeanBoot")
writeData(wb, sheet = "ParametersMeanBoot", x = Parameters_Mean_boot, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ParametersMax")
writeData(wb, sheet = "ParametersMax", x = Parameters_Max, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ParametersMaxBoot")
writeData(wb, sheet = "ParametersMaxBoot", x = Parameters_Max_boot, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ChoosenParameters")
writeData(wb, sheet = "ChoosenParameters", x = Param, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "ChoosenParametersBoot")
writeData(wb, sheet = "ChoosenParametersBoot", x = BoostrapIntervals, colNames = TRUE, rowNames = TRUE)
# Save the workbook with the specified title
saveWorkbook(wb, file = titulo_xlsx, overwrite = TRUE)
#####



###########################
# Using Importance Sampling to obtain
# posterior estimates for the shape parameters of the beta distribution
# using the new bivariate prior distribution, the sets of hyperparameters obtained and stored
# in Parameters_Min, Parameters_Mean, Parameters_Max, Param.
# For each sample size defined in n_sample, the Sim_study function generates a prior distribution sample of size N1,
# generates 1000 samples of the established sample size by n_sample, and with each sample, the posterior estimate
# is obtained from the importance sampling.
###########################

Parall_sim_study=function(n,hyperparameters_values,N_Iter_Sim){
# Number of cores to use
numCores <- detectCores() - 1  # Use all but one to avoid saturating the machine
cl <- makeCluster(numCores)
registerDoParallel(cl, c("Sim_study"))

# Define the variables Method_name, Param, Vect_Alpha, Vect_Beta
# Here you should define the variables used in the code
#Param = Parameters_Min  # Set of hyperparameters established for the prior distribution and obtained by choosing the bound_var parameter.
Param = hyperparameters_values
#n = 0  # Represents the scenario in question; n = 0 corresponds to scenario 1, n = 1 to scenario 2, and n = 2 to scenario 3.

# Initialize the data frames
Results_Scen_Alpha <- data.frame()
Results_Scen_Beta <- data.frame()
Results_Scen_Diff <- data.frame()
Method_name = c("BM", "BT", "EM1", "ET1", "EM2", "ET2", "EM3", "ET3", "EM4", "ET4")

# Parallelize the for loop
Results_list <- foreach(n_hiper = 1:length(Method_name), .combine = 'c',.export = c("Sim_study","Gen_Joint_Dist","Gen_FC_X1_X2","Measure_Diagnostic","Mom_Prior_Dist",
                                                                                    "Vect_Alpha","Vect_Beta"),
                        .packages = c("betafunctions", "coda")) %dopar% {
  Results_Scen <- Sim_study(N = 10^4, N_FC = 2, prop_prec = 3, a = Param$Va[n_hiper + 10 * n], b = Param$Vb[n_hiper + 10 * n], c = Param$Vc[n_hiper + 10 * n], d = Param$Vd[n_hiper + 10 * n],
                                  thin1 = 1, X10_given = "random", dig_tol = 15, thin2 = 1, burnin = 2000, n_sample = c(seq(10, 50, 5), seq(100, 200, 50)),
                                  alpha_real = Vect_Alpha[n + 1], beta_real = Vect_Beta[n + 1], N_Iter_Sim = N_Iter_Sim, sample_size_IS = (10^4-2000)/1)
  
  Results_Scen$Result_Alpha$Method <- Method_name[n_hiper]
  Results_Scen$Result_Beta$Method <- Method_name[n_hiper]
  Results_Scen$Descriptive_Sample_Prior$Differences$Method <- Method_name[n_hiper]
  
  list(
    Alpha = Results_Scen$Result_Alpha,
    Beta = Results_Scen$Result_Beta,
    Diff = Results_Scen$Descriptive_Sample_Prior$Differences
  )
}
# Stop the cluster
stopCluster(cl)
# Combine the results after the foreach loop
for (result in seq_along(Results_list)) {
  Results_Scen_Alpha <- rbind(Results_Scen_Alpha, Results_list[result]$Alpha)
  Results_Scen_Beta <- rbind(Results_Scen_Beta, Results_list[result]$Beta)
  Results_Scen_Diff <- rbind(Results_Scen_Diff, Results_list[result]$Diff)
}


return(list(Results_Scen_Alpha=Results_Scen_Alpha, Results_Scen_Beta=Results_Scen_Beta,Results_Scen_Diff=Results_Scen_Diff))
}

# Obtaining posterior estimates for each scenario, set of hyperparameters, 12 sample sizes and 1000 repetitions.

Parameters_Scenario1=Parall_sim_study(n=0,hyperparameters_values=Param,N_Iter_Sim=1000)
Parameters_Scenario2=Parall_sim_study(n=1,hyperparameters_values=Param,N_Iter_Sim=1000)
Parameters_Scenario3=Parall_sim_study(n=2,hyperparameters_values=Param,N_Iter_Sim=1000)

Min_Scenario1=Parall_sim_study(n=0,hyperparameters_values=Parameters_Min,N_Iter_Sim=1000)
Min_Scenario2=Parall_sim_study(n=1,hyperparameters_values=Parameters_Min,N_Iter_Sim=1000)
Min_Scenario3=Parall_sim_study(n=2,hyperparameters_values=Parameters_Min,N_Iter_Sim=1000)

Mean_Scenario1=Parall_sim_study(n=0,hyperparameters_values=Parameters_Mean,N_Iter_Sim=1000)
Mean_Scenario2=Parall_sim_study(n=1,hyperparameters_values=Parameters_Mean,N_Iter_Sim=1000)
Mean_Scenario3=Parall_sim_study(n=2,hyperparameters_values=Parameters_Mean,N_Iter_Sim=1000)

Max_Scenario1=Parall_sim_study(n=0,hyperparameters_values=Parameters_Max,N_Iter_Sim=1000)
Max_Scenario2=Parall_sim_study(n=1,hyperparameters_values=Parameters_Max,N_Iter_Sim=1000)
Max_Scenario3=Parall_sim_study(n=2,hyperparameters_values=Parameters_Max,N_Iter_Sim=1000)


# The following graphs represent the scenario defined by the parameter n and the Hyperparameters obtained by the choice of bound_var.
Graphs_Comparison_Hyper=function(n,results_sim_study){
# Posterior estimates obtained for the Alpha parameter of the Beta distribution of the variable X.
Comparison_Hyper(results_sim_study$Results_Scen_Alpha, value_real = Vect_Alpha[n + 1], lim_x = c(10, 25, 50, 100, 150, 200),
                title_text = " ", y_Text = "Alpha")
# MinAlphaSigMu01SigV01Scen-3

# Posterior estimates obtained for the Beta parameter of the Beta distribution of the variable X.
Comparison_Hyper(results_sim_study$Results_Scen_Beta, value_real = Vect_Beta[n + 1], lim_x = c(10, 25, 50, 100, 150, 200),
                title_text = " ", y_Text = "Beta")
# MinBetaSigMu01SigV01Scen-3
}

# Graph of the results obtained in the simulation study for each scenario, set of hyperparameters, sample sizes and 1000 repetitions.
Graphs_Comparison_Hyper(n=0,results_sim_study = Parameters_Scenario1)
Graphs_Comparison_Hyper(n=1,results_sim_study = Parameters_Scenario2)
Graphs_Comparison_Hyper(n=2,results_sim_study = Parameters_Scenario3)

Graphs_Comparison_Hyper(n=0,results_sim_study = Min_Scenario1)
Graphs_Comparison_Hyper(n=1,results_sim_study = Min_Scenario2)
Graphs_Comparison_Hyper(n=2,results_sim_study = Min_Scenario3)

Graphs_Comparison_Hyper(n=0,results_sim_study = Mean_Scenario1)
Graphs_Comparison_Hyper(n=1,results_sim_study = Mean_Scenario2)
Graphs_Comparison_Hyper(n=2,results_sim_study = Mean_Scenario3)

Graphs_Comparison_Hyper(n=0,results_sim_study = Max_Scenario1)
Graphs_Comparison_Hyper(n=1,results_sim_study = Max_Scenario2)
Graphs_Comparison_Hyper(n=2,results_sim_study = Max_Scenario3)


# Storage of results in spreadsheets.
Parameters = Max_Scenario3
# Create a workbook object to save Scen results
titulo_xlsx = "MaxSigMu01SigV01Scen-3.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Parameters")
writeData(wb, sheet = "Parameters", x = Param, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "Alpha")
writeData(wb, sheet = "Alpha", x = Parameters$Results_Scen_Alpha, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "Beta")
writeData(wb, sheet = "Beta", x = Parameters$Results_Scen_Beta, colNames = TRUE, rowNames = TRUE)
addWorksheet(wb, "Difference")
writeData(wb, sheet = "Difference", x =Parameters$Results_Scen_Diff, colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, file = titulo_xlsx, overwrite = TRUE)

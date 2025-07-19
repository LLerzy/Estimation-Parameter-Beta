source(file = "requiredfunctions.R")
###########################################################
# Dataset: Quality Control
###########################################################
QC_Data=read.xlsx("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx")
var(QC_Data$Defect)

length(QC_Data$Defect)
# Skewness vs. kurtosis plot to assess whether the MCHC distribution follows a standard Beta distribution
library(fitdistrplus)
descdist(QC_Data$Defect, discrete = FALSE)

# Estimation of hyperparameter values using bootstrap quantile intervals for the mean and variance
hyper_Quality_Control = Hyperparameters(
  QC_Data$Defect,
  r_boostrap = 200, q_boostrap = c(0.025, 0.975), 
  option_mu = "moments", 
  sig_mu = 0.05, 
  bound_var = "min", 
  sig_var = 0.2, 
  digits = 4,
  graphs_boot = T, 
  Q_E_mu = 0, Q_E_cv = 0,language = "Spanish"
  )
hyper_Quality_Control$Q_var
# Extracted hyperparameters: a and b for the mean; c and d for the variance
c(hyper_Quality_Control$hiper_mean[[1]],hyper_Quality_Control$hiper_mean[[2]],hyper_Quality_Control$hiper_var$a[[1]],hyper_Quality_Control$hiper_var$b[[1]])

a2 = hyper_Quality_Control$hiper_mean[[1]]
b2 = hyper_Quality_Control$hiper_mean[[2]]
c2 = hyper_Quality_Control$hiper_var$a[[1]]
d2 = hyper_Quality_Control$hiper_var$b[[1]]

# Importance sampling:
# Posterior numerical estimates for the moments of the prior distribution using the above hyperparameter values
results_Quality_Control = Est_Post(
  ssample = round(QC_Data$Defect,5), 
  N=10^5, 
  N_FC=2, 
  Precision=3, 
  a=a2, b=b2, c=c2, d=d2, 
  thin1=2, thin2=10, burnin=5000, dig_tol=5
  )

round(t(rbind(results_Quality_Control$Descriptivo$Numerical,
              results_Quality_Control$Descriptivo$Analytical,
              results_Quality_Control$Descriptivo$Differences)),3)

# Sampling chains used in the previous importance sampling
results_Quality_Control$SA
results_Quality_Control$SB

# Importance resampling using the previous sampling chains
selected_elements_Quality_Control = sample(
  1:length(results_Quality_Control$SA),
  size = 2000,
  replace = T,
  prob = results_Quality_Control$Zeta/results_Quality_Control$EstPost$Prior_Lik
  )
selected_alpha_Quality_Control = results_Quality_Control$SA[selected_elements_Quality_Control]
selected_beta_Quality_Control = results_Quality_Control$SB[selected_elements_Quality_Control]

# Credible intervals for the mean of alpha and beta
RC_A_CC = quantile(selected_alpha_Quality_Control, probs = c(0.025, 0.975), na.rm = T)
RC_B_CC = quantile(selected_beta_Quality_Control, probs = c(0.025, 0.975), na.rm = T)

# Descriptive measures for the posterior sample obtained via importance resampling
Post_results_QC = Measure_Diagnostic(
  data1 = selected_alpha_Quality_Control, 
  data2 = selected_beta_Quality_Control, 
  var = "original", 
  burnin=0, 
  thin=1, digits = 5, 
  a=a2, b=b2, c=c2, d=d2, 
  cred_level = 0.95, batch_size = 100
  )

round(t(Post_results_QC$Numerical),3)

#########################################################################
# Dengue dataset: Mean Corpuscular Hemoglobin Concentration (MCHC)
#########################################################################
Dengue_Data = read.xlsx("Dengue2.xlsx", sheet = "Ratio",colNames = T)

# round(Dengue_Data$CMHC,4)
# var(Dengue_Data$CMHC)
# prod(1-Dengue_Data$CMHC)

# Skewness vs. kurtosis plot to assess whether the MCHC distribution follows a standard Beta distribution
library(fitdistrplus)
descdist(Dengue_Data$CMHC, discrete = FALSE)

# Estimation of hyperparameter values using bootstrap quantile intervals for the mean and variance
hyper_dengue = Hyperparameters(
  round(Dengue_Data$CMHC,5), 
  r_boostrap = 300, 
  q_boostrap = c(0.025, 0.975), 
  option_mu = "moments",
  sig_mu = 0.05, 
  bound_var = "max", 
  sig_var = 0.05, 
  digits = 4,
  graphs_boot = T, 
  Q_E_mu = 0, Q_E_c = 0,
  language = "Spanish"
  )

# Extracted hyperparameters: a and b for the mean; c and d for the variance
round(c(hyper_dengue$hiper_mean[[1]],hyper_dengue$hiper_mean[[2]],hyper_dengue$hiper_var$a[[1]],hyper_dengue$hiper_var$b[[1]]),3)

a3=hyper_dengue$hiper_mean[[1]]
b3=hyper_dengue$hiper_mean[[2]]
c3=hyper_dengue$hiper_var$a[[1]]
d3=hyper_dengue$hiper_var$b[[1]]

# Importance sampling:
# Posterior numerical estimates for the moments of the prior distribution using the above hyperparameter values
results_dengue_post = Est_Post(
  ssample = round(Dengue_Data$CMHC,5), 
  N=10^5, 
  N_FC=2,
  Precision=3, 
  a=a3, b=b3, c=c3, d=d3, 
  thin1=2, thin2=10, burnin=5000, dig_tol=5
  )

(100000-5000)/10
results_dengue_post$EstPost
round(t(rbind(results_dengue_post$Descriptivo$Numerical,
              results_dengue_post$Descriptivo$Analytical,
              results_dengue_post$Descriptivo$Differences)),3)

# Sampling chains used in the previous importance sampling
results_dengue_post$SA
results_dengue_post$SB

# Importance resampling using the previous sampling chains
selected_elements = sample(
  1:length(results_dengue_post$SA),
  size = 2000,
  replace = T,
  prob = results_dengue_post$Zeta/results_dengue_post$EstPost$Prior_Lik
  )
selected_alpha = results_dengue_post$SA[selected_elements]
selected_beta = results_dengue_post$SB[selected_elements]

# Credible intervals for the mean of alpha and beta
RC_A = quantile(selected_alpha, probs = c(0.025, 0.975), na.rm = T)
RC_B = quantile(selected_beta, probs = c(0.025, 0.975), na.rm = T)

# Descriptive measures for the posterior sample obtained via importance resampling
Post_results_Dengue=Measure_Diagnostic(
  data1 = selected_alpha, data2 = selected_beta, 
  var = "original", 
  burnin=0, thin=1, digits = 5, 
  a=a3, b=b3, c=c3, d=d3, 
  cred_level = 0.95, batch_size = 100
  )
round(t(Post_results_Dengue$Numerical),3)

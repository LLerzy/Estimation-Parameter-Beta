


####
## Momentos de la distribución conjunta para la media y la varianza de la variable aleatoria X.
####

# Momentos marginales de mu.
Mom_Mu=function(k,a,b,c,d){
  beta(a+k,b)/ beta(a,b)
  # gamma(a+k) * gamma(a+b) / (gamma(a) * gamma(a+b+k))
}

# Momentos marginales de sigma.
Mom_Sigma=function(k,a,b,c,d){
  beta(c+k,d) * beta(a+k,b+k) / (beta(a,b)*beta(c,d))
}

# Momentos conjuntos de orden k=k1+k2 para mu y sigma.
Mom_Mu_Sigma=function(k1,k2,a,b,c,d){
  beta(c+k2,d) * beta(a+k1+k2,b+k2) / (beta(a,b)*beta(c,d))
}

# Momentos para los valores de hiperparámetros obtenidos en los escenarios propuestos,
# utilizando información subjetiva (obtenida de una persona considerada experta) e infomración boostrap.
Param_Prior_2=Param
Param_Prior_2$Mean.1 = Mom_Mu_Sigma(1,0,Param$Va,Param$Vb,Param$Vc,Param$Vd)
Param_Prior_2$Var.1 = Mom_Mu_Sigma(2,0,Param$Va,Param$Vb,Param$Vc,Param$Vd)-Mom_Mu_Sigma(1,0,Param$Va,Param$Vb,Param$Vc,Param$Vd)^2
Param_Prior_2$Mean.2 = Mom_Mu_Sigma(0,1,Param$Va,Param$Vb,Param$Vc,Param$Vd)
Param_Prior_2$Var.2 = Mom_Mu_Sigma(0,2,Param$Va,Param$Vb,Param$Vc,Param$Vd) - Mom_Mu_Sigma(0,1,Param$Va,Param$Vb,Param$Vc,Param$Vd)^2
Param_Prior_2$Cov = Mom_Mu_Sigma(1,1,Param$Va,Param$Vb,Param$Vc,Param$Vd)-Mom_Mu_Sigma(0,1,Param$Va,Param$Vb,Param$Vc,Param$Vd)*Mom_Mu_Sigma(1,0,Param$Va,Param$Vb,Param$Vc,Param$Vd)

# Tablas para Latex
library(xtable)
xtable(cbind(Param_Prior_2[,12],Param_Prior_2[,17],format(Param_Prior_2[,18],scientific = TRUE,digits=3),
             Param_Prior_2[,19],format(Param_Prior_2[,c(20,21)],scientific = TRUE,digits=3)),digits = 2)

xtable(cbind(Param_Prior_2[,12],Param_Prior_2[,17],format(Param_Prior_2[,18:21],scientific = TRUE,digits=3)),digits = 2)


##################################################################################################
##################################################################################################


##############
# Estudio de simulación
#############

Measure_Diagnostic = function(data1, data2, var = "original", burnin, thin, digits = 5, a, b, c, d) {
  N = length(data1)
  data1 <- data1[seq((burnin + 1), N, by = thin)]
  data2 <- data2[seq((burnin + 1), N, by = thin)]
  
  if (var == "original") {
    new_names = c("Mean_X1", "Var_X1", "ESS_X1", "Mean_X2", "Var_X2", "ESS_X2", "Cov", "Length")
    
    # Numerical results
    Numerical_results = round(data.frame(
      "mean1" = mean(data1), "var1" = var(data1), "ESS1" = effectiveSize(mcmc(data1))[[1]],
      "mean2" = mean(data2), "var2" = var(data2), "ESS2" = effectiveSize(mcmc(data2))[[1]],
      "cov12" = cov(data1, data2), "length" = length(data1)
    ), digits)
    names(Numerical_results) = new_names
    
    # Analytical results
    # K = Mom_Mu_Sigma(0, 0, a, b, c, d)
    Analytic_results = round(data.frame(
      "Mean.1" = Mom_Mu_Sigma(1, 0, a, b, c, d),
      "Var.1" = Mom_Mu_Sigma(2, 0, a, b, c, d) - (Mom_Mu_Sigma(1, 0, a, b, c, d))^2,
      "ESS.1" = length(data1),
      "Mean.2" = Mom_Mu_Sigma(0, 1, a, b, c, d),
      "Var.2" = Mom_Mu_Sigma(0, 2, a, b, c, d) - (Mom_Mu_Sigma(0, 1, a, b, c, d))^2,
      "ESS.2" = length(data1),
      "Cov" = Mom_Mu_Sigma(1, 1, a, b, c, d) - (Mom_Mu_Sigma(1, 0, a, b, c, d)) * Mom_Mu_Sigma(0, 1, a, b, c, d),
      "length" = length(data1)
    ), digits)
    names(Analytic_results) = new_names
    
    # Differences between analytical and numerical results.
    Differences = round(Analytic_results - Numerical_results, digits)
    return(list(Numerical = Numerical_results, Analytical = Analytic_results, Differences = Differences))
    
  } else if (var == "transform") {
    piece = (data1 * (1 - data1) / data2 - 1)
    new_data1 = data1 * piece
    new_data2 = (1 - data1) * piece
    new_names = c("Mean_Y1", "Var_Y1", "ESS_Y1", "Mean_Y2", "Var_Y2", "ESS_Y2", "Cov", "Length")
    
    # Numerical results
    Numerical_results = round(data.frame(
      "mean1" = mean(new_data1), "var1" = var(new_data1), "ESS1" = effectiveSize(mcmc(new_data1))[[1]],
      "mean2" = mean(new_data2), "var2" = var(new_data2), "ESS2" = effectiveSize(mcmc(new_data2))[[1]],
      "cov12" = cov(new_data1, new_data2), "length" = length(new_data1)
    ), digits)
    names(Numerical_results) = new_names
    
    # Analytical results
    K = Mom_Prior_Dist(0, 0, a, b, c, d)
    Analytic_results = round(data.frame(
      "Mean.1" = Mom_Prior_Dist(1, 0, a, b, c, d) / K,
      "Var.1" = Mom_Prior_Dist(2, 0, a, b, c, d) / K - (Mom_Prior_Dist(1, 0, a, b, c, d) / K)^2,
      "ESS.1" = length(new_data1),
      "Mean.2" = Mom_Prior_Dist(0, 1, a, b, c, d) / K,
      "Var.2" = Mom_Prior_Dist(0, 2, a, b, c, d) / K - (Mom_Prior_Dist(0, 1, a, b, c, d) / K)^2,
      "ESS.2" = length(new_data1),
      "Cov" = Mom_Prior_Dist(1, 1, a, b, c, d) / K - (Mom_Prior_Dist(1, 0, a, b, c, d) / K) * Mom_Prior_Dist(0, 1, a, b, c, d) / K,
      "length" = length(new_data1)
    ), digits)
    names(Analytic_results) = new_names
    
    # Differences between analytical and numerical results.
    Differences = round(Analytic_results - Numerical_results, digits)
    return(list(Numerical = Numerical_results, Analytical = Analytic_results, Differences = Differences))
  }
}




##################
## Función para realiar la simulación: obtener estimaciones posteriores y criterios de convergencia de las estimaciones.
##################
Sim_study_MS = function(N, N_FC, prop_prec, a, b, c, d, thin1, X10_given = "random", dig_tol = 15, thin2, burnin, n_sample, alpha_real, beta_real, N_Iter_Sim,
                     sample_size_IS) {
  # Average length 
  meanlenght=function(Interval){
    lenint=sum(Interval[,2]-Interval[,1])/nrow(Interval)
    return(lenint=lenint)
  }
  
  # Coverage probability
  covprob=function(Interval,param_value){
    indicator=function(j){
      if(Interval[j,1]<=param_value & Interval[j,2]>=param_value){
        return(1)
      }else{return(0)}
    }
    resind=NA
    for(l in 1:nrow(Interval)){
      resind[l]=indicator(l)
    }
    covint=sum(resind)/nrow(Interval)
    
    return(covint=covint)
  }
  
  mu_real=alpha_real/ (alpha_real + beta_real)
  sigma_real=alpha_real * beta_real / ((alpha_real + beta_real)^2 * (alpha_real + beta_real + 1))
  
  # Generate prior sample
  Sample_Prior_Hip = Gen_Joint_Dist(N1 = N, N2 = N_FC, prop_prec = 3, a, b, c, d, thin = thin1, X10_given = "random", dig_tol = dig_tol)
  
  sample_mu = Sample_Prior_Hip$X1[seq((burnin + 1), N, by = thin2)]
  sample_sigma = Sample_Prior_Hip$X2[seq((burnin + 1), N, by = thin2)]
  
  # Descriptive statistics for the prior sample
  Descriptive_Sample_Prior = Measure_Diagnostic(data1 = Sample_Prior_Hip$X1, data2 = Sample_Prior_Hip$X2, var = "original",
                                                digits = 4, burnin, thin = thin2, a, b, c, d)
  
  Iter_Mu = NA
  Iter_Sigma = NA
  VIter_Mu = matrix(data = NA, nrow = 0, ncol = 12, dimnames = list(NULL, c("Min", "Q2.5", "Q50", "Q97.5", "Max", "Mean", "Var", "Bias", "MSE", "SampleSize","Coverage","Length")))
  VIter_Sigma = VIter_Mu
  RC_Mu=matrix(data = NA, nrow = 0, ncol = 2, dimnames = list(NULL, c("Q2.5", "Q97.5")))
  RC_Sig=RC_Mu
  for (j in 1:length(n_sample)) {
    for (i in 1:N_Iter_Sim) {
      sampe_prueba = rbeta(n_sample[j], alpha_real, beta_real)
      x0 = prod(sampe_prueba)
      y0 = prod(1 - sampe_prueba)
      
      sample_alpha= sample_mu * (sample_mu * (1 - sample_mu) / sample_sigma - 1)
      sample_beta = (1 - sample_mu) * (sample_mu * (1 - sample_mu) / sample_sigma - 1)
        
      zeta=exp((sample_alpha-1)*log(x0) + (sample_beta-1)*log(y0)-(n_sample[j])*log(beta(sample_alpha,sample_beta)))
      marginalike=sum(zeta)
      Iter_Mu[i] = sum(zeta*exp(log(sample_mu)-log(marginalike)))
      Iter_Sigma[i] = sum(zeta*exp(log(sample_sigma)-log(marginalike))) 
      
      # Sampling importance resampling
      selected_elements=sample(1:length(sample_mu),size = sample_size_IS,replace = T,prob = zeta/marginalike)
      selected_mu=sample_mu[selected_elements]
      selected_sigma=sample_sigma[selected_elements]
      
      # Regions of credibility
      RC_Mu = rbind(RC_Mu, quantile(selected_mu, probs = c(0.025, 0.975), na.rm = T))
      RC_Sig = rbind(RC_Sig, quantile(selected_sigma, probs = c(0.025, 0.975), na.rm = T))
      
    }
    VIter_Mu = rbind(VIter_Mu, c(quantile(Iter_Mu, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = T), mean(Iter_Mu), var(Iter_Mu), (mu_real - mean(Iter_Mu)),
                                       (var(Iter_Mu) + (mu_real - mean(Iter_Mu))^2), n_sample[j],covprob(RC_Mu,mu_real),meanlenght(RC_Mu)))
    VIter_Sigma = rbind(VIter_Sigma, c(quantile(Iter_Sigma, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = T), mean(Iter_Sigma), var(Iter_Sigma), (sigma_real - mean(Iter_Sigma)),
                                     (var(Iter_Sigma) + (sigma_real - mean(Iter_Sigma))^2), n_sample[j], covprob(RC_Sig,sigma_real), meanlenght(RC_Sig)))
  }
  
  return(list(Result_Mu = as.data.frame(VIter_Mu), Result_Sigma = as.data.frame(VIter_Sigma), Descriptive_Sample_Prior = Descriptive_Sample_Prior))
}


Results_Scen <- Sim_study_MS(N = 10^2, N_FC = 2, prop_prec = 3, a = Param$Va[1], b = Param$Vb[1], c = Param$Vc[1], d = Param$Vd[1],
                          thin1 = 1, X10_given = "random", dig_tol = 15, thin2 = 1, burnin = 20, n_sample = 10,
                          alpha_real = Vect_Alpha[1], beta_real = Vect_Beta[1], N_Iter_Sim = 10, sample_size_IS = (10^2-20)/1)


Results_Scen$Result_Mu




Parall_sim_study_MS=function(n,hyperparameters_values,N_Iter_Sim){
  # Number of cores to use
  numCores <- detectCores() - 1  # Use all but one to avoid saturating the machine
  cl <- makeCluster(numCores)
  registerDoParallel(cl, c("Sim_study_MS"))
  
  # Define the variables Method_name, Param, Vect_Alpha, Vect_Beta
  # Here you should define the variables used in the code
  #Param = Parameters_Min  # Set of hyperparameters established for the prior distribution and obtained by choosing the bound_var parameter.
  Param = hyperparameters_values
  #n = 0  # Represents the scenario in question; n = 0 corresponds to scenario 1, n = 1 to scenario 2, and n = 2 to scenario 3.
  
  # Initialize the data frames
  Results_Scen_Mu <- data.frame()
  Results_Scen_Sigma <- data.frame()
  Results_Scen_Diff <- data.frame()
  Method_name = c("BM", "BT", "EM1", "ET1", "EM2", "ET2", "EM3", "ET3", "EM4", "ET4")
  
  # Parallelize the for loop
  Results_list <- foreach(n_hiper = 1:length(Method_name), .combine = 'c',.export = c("Sim_study_MS","Gen_Joint_Dist","Gen_FC_X1_X2","Measure_Diagnostic","Mom_Mu_Sigma",
                                                                                      "Vect_Alpha","Vect_Beta"),
                          .packages = c("betafunctions", "coda")) %dopar% {
                            Results_Scen <- Sim_study_MS(N = 10^4, N_FC = 2, prop_prec = 3, a = Param$Va[n_hiper + 10 * n], b = Param$Vb[n_hiper + 10 * n], c = Param$Vc[n_hiper + 10 * n], d = Param$Vd[n_hiper + 10 * n],
                                                      thin1 = 1, X10_given = "random", dig_tol = 15, thin2 = 1, burnin = 2000, n_sample = c(seq(10, 50, 5), seq(100, 200, 50)),
                                                      alpha_real = Vect_Alpha[n + 1], beta_real = Vect_Beta[n + 1], N_Iter_Sim = N_Iter_Sim, sample_size_IS = (10^4-2000)/1)
                            
                            Results_Scen$Result_Mu$Method <- Method_name[n_hiper]
                            Results_Scen$Result_Sigma$Method <- Method_name[n_hiper]
                            Results_Scen$Descriptive_Sample_Prior$Differences$Method <- Method_name[n_hiper]
                            
                            list(
                              Mu = Results_Scen$Result_Mu,
                              Sigma = Results_Scen$Result_Sigma,
                              Diff = Results_Scen$Descriptive_Sample_Prior$Differences
                            )
                          }
  # Stop the cluster
  stopCluster(cl)
  # Combine the results after the foreach loop
  for (result in seq_along(Results_list)) {
    Results_Scen_Mu <- rbind(Results_Scen_Mu, Results_list[result]$Mu)
    Results_Scen_Sigma <- rbind(Results_Scen_Sigma, Results_list[result]$Sigma)
    Results_Scen_Diff <- rbind(Results_Scen_Diff, Results_list[result]$Diff)
  }
  
  
  return(list(Results_Scen_Mu=Results_Scen_Mu, Results_Scen_Sigma=Results_Scen_Sigma,Results_Scen_Diff=Results_Scen_Diff))
}

Vect_Alpha = c(0.5, 16, 3)
Vect_Beta = c(0.5, 4, 12)
Vect_mu = Vect_Alpha / (Vect_Alpha + Vect_Beta)  # Mean of the beta distribution of X.
Vect_var = Vect_Alpha * Vect_Beta / ((Vect_Alpha + Vect_Beta)^2 * (Vect_Alpha + Vect_Beta + 1))  # Variance of the beta distribution of X.


Parameters_Scenario_MS1=Parall_sim_study(n=0,hyperparameters_values=Param,N_Iter_Sim=1000)
Parameters_Scenario_MS2=Parall_sim_study(n=1,hyperparameters_values=Param,N_Iter_Sim=1000)
Parameters_Scenario_MS3=Parall_sim_study(n=2,hyperparameters_values=Param,N_Iter_Sim=1000)



Graphs_Comparison_Hyper_MS=function(n,results_sim_study){
  # Posterior estimates obtained for the Mu Moment of the Beta distribution of the variable X.
  Comparison_Hyper(results_sim_study$Results_Scen_Mu, value_real = Vect_mu[n + 1], lim_x = c(10, 25, 50, 100, 150, 200),
                   title_text = " ", y_Text = "Mu")
  # MinAlphaSigMu01SigV01Scen-3
  
  # Posterior estimates obtained for the Beta parameter of the Beta distribution of the variable X.
  Comparison_Hyper(results_sim_study$Results_Scen_Sigma, value_real = Vect_var[n + 1], lim_x = c(10, 25, 50, 100, 150, 200),
                   title_text = " ", y_Text = "Sigma")
  # MinBetaSigMu01SigV01Scen-3
}


Graphs_Comparison_Hyper_MS(n=0,results_sim_study = Parameters_Scenario_MS1)

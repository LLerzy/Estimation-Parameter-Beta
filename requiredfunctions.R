###########################################################################
#                                                                         #
# Description:                                                            #
# This file contains a set of functions and routines used for generating  #
# random samples from a joint distribution for the shape parameters of    #
# the beta distribution, and for evaluating the convergence of the        #
# generated samples. Additionally, it includes functions for simulating,  #
# monitoring, and comparing posterior estimates obtained using different  #
# hyperparameter sets and sample sizes.                                   #
#                                                                         #
# Author: Llerzy Torres Ome                                               #
# Creation Date: July 23, 2024                                            #
#                                                                         #
# Functions included:                                                     #
# 1. Prior: Defines the proposed prior probability density function.      #
# 2. FC_X1_Given_v: Full conditional distribution of X1 given X2 = v.     #
# 3. Graph_Fc_X1: Plots the full conditional distribution for given       #
#    parameter values.                                                    #
# 4. Gen_FC_X1_X2: Metropolis-Hastings algorithm using random walks for   #
#    generating samples from the full conditional distribution.           #
# 5. Mon_Measure: Monitors acceptance rates and ESS for different values  #
#    of v and precision.                                                  #
# 6. Mon_R_Hat: Monitors the Gelman-Rubin diagnostic (R-hat) for          #
#    different values of v and precision.                                 #
# 7. Graphs: Plots histograms, density, trace, and convergence control    #
#    using the average.                                                   #
# 8. Gen_Joint_Dist: Gibbs sampling for generating joint distributions    #
#    of X1 and X2.                                                        #
# 9. Mtovar_vs2: Generalizes Tovar's method for obtaining hyperparameter  #
#    values.                                                              #
# 10. Mom_Prior_Dist: Calculates joint moments of order l for the proposed  #
#     prior distribution.                                                 #
# 11. Measure_Diagnostic: Compares analytical and numerical results for   #
#     given data samples.                                                 #
# 12. Measure_Analy: Computes analytical results for the proposed prior   #
#     distribution.                                                       #
# 13. Hyperparameters: Obtains hyperparameters using empirical Bayes and  #
#     subjective approaches.                                              #
# 14. Est_Post: Posterior estimation for alpha and beta parameters of the #
#     beta distribution using importance sampling.                        #
# 15. Sim_study: Conducts simulation studies to compare posterior       #
#     estimates using different hyperparameters and sample sizes.         #
# 16. Individual_Graphs: Creates individual graphs to monitor posterior   #
#     estimates using bias and MSE as indicators.                         #
# 17. Comparison_Hyper: Compares joint functions for different hyperparam- #
#     eter sets.                                                          #
#                                                                         #
# Notes:                                                                  #
# 1. It's recommended to review and adapt each function according to the  #
#    specific needs of each analysis.                                     #
# 2. Ensure you understand each function before using it to guarantee     #
#    accurate results and avoid potential errors.                         #
# 3. For any questions or suggestions, contact                            #
#    llerzy.torres@correounivalle.edu.co                                  #
###########################################################################


##########################################################
##########################################################
# Necessary library
##########################################################
##########################################################
library(ggplot2)
library(gridExtra)
library(tidyr) # provides tools for tidying up data and especially useful for its pipe operator (%>%), which streamlines data manipulation and transformation.
library(plotly) # used for interactive visualization of contour plots
library(coda)
library(foreach)
library(doParallel)
library(betafunctions)
library(openxlsx)
library(xtable)
##########################################################

##########################################################
## Proposed Prior Probability Density Function
##########################################################
# alph and bet are the random variables.
# a, b, c, and d are the parameter values of the distribution.
# Note: The relationship between exp and log is used to avoid situations that generate NaN or Inf.
# The prior function returns the prior probability density for the Beta distribution parameters.
Prior=function(alph, bet, a, b, c, d) {
  # We calculate the prior density using the logarithmic expression to avoid numerical issues.
  # The original formula is:
  # 1 / (beta(a, b) * beta(c, d)) * alph^(a-1) * bet^(b-1) * (alph + bet)^(d-a-b) * (alph + bet + 1)^(-c-d)
  # We use log and exp to improve numerical stability:
  return(exp((a-1) * log(alph) + (b-1) * log(bet) + (d-a-b) * log(alph + bet) - (c+d) * log(alph + bet + 1)))
}

# Note:
# Although "prior" is the prior distribution for alpha and beta of the Beta distribution of a variable X,
# we will refer to them as Y1 and Y2, respectively.
# Additionally, X1 and X2 will represent the mean and variance associated with the Beta distribution of X.
##########################################################
##########################################################
# Metropolis-Hastings Method using Random Walks
# for the conditional distribution of X1 given X2
##########################################################
##########################################################
# Full Conditional (FC) of X1 given X2 equals v
# "(a,b,c,d)" is the vector of parameters
# "X1" is within (0,1) and "X2=v" is less than X1(1-X1)
FC_X1_Given_v=function(X1, a, b, c, d, v) {
  # Calculate the full conditional density of X1 given X2 = v
  # The original formula is:
  # result1 = X1^(a-c-d) * (1-X1)^(b-c-d) * (X1*(1-X1)-v)^(d-1)
  # We use log and exp to improve numerical stability:
  return(exp((a-c-d) * log(X1) + (b-c-d) * log(1-X1) + (d-1) * log(X1*(1-X1)-v)))
}

#####
# Plot of the Full Conditional (FC) for given values of a, b, c, d, and three values of v.
#####
# "v1", "v2", and "v3" are given values of the variance.
# "v1name", "v2name", and "v3name" are the names of each plot.
# "(ae, be, ce, de)" is the vector of given parameter values.
Graph_Fc_X1=function(v1, v1name, v2, v2name, v3, v3name, ae, be, ce, de) {
  ggplot() + xlim(c(0, 1)) +
    # Plot FC for v1
    geom_function(fun=function(X1) mapply(FC_X1_Given_v, X1, a=ae, b=be, c=ce, d=de, v=v1), lwd=1,
                  linetype=1, aes(col=v1name)) +
    # Plot FC for v2
    geom_function(fun=function(X1) mapply(FC_X1_Given_v, X1, a=ae, b=be, c=ce, d=de, v=v2), lwd=1,
                  linetype=1, aes(col=v2name)) +
    # Plot FC for v3
    geom_function(fun=function(X1) mapply(FC_X1_Given_v, X1, a=ae, b=be, c=ce, d=de, v=v3), lwd=1,
                  linetype=1, aes(col=v3name)) +
    labs(title="Full Conditional of X1 given X2",
         caption=substitute(
           list("Plot of", f(X1/v)==X1^(a-c-d) * (1-X1)^(b-c-d) * (X1*(1-X1)-v)^(d-1),
                "with", a==ae, b==be, c==ce, d==de), list(ae=ae, be=be, ce=ce, de=de))) +
    xlab(expression(X1)) + ylab(expression(f(X1/X2==v))) +
    scale_colour_manual(values = c("red", "black", "purple"), name="X2=v")
}
#####

#####
# Metropolis-Hastings using Random Walks Algorithm 
#####
# "N" is the sample size to be generated.
# "prop_prec" is the precision set for the algorithm. 
# "a", "b", "c", and "d" are given values for the parameters.
# "v" is the given value for X2.
# "option" allows you to select the entire sample ("all") or just the last value generated ("end").
# "thin" is the thinning interval for MCMC. Every "thin" generated samples, one is stored to reduce autocorrelation.
# "burnin" is the number of iterations to discard.
# The seed type can be specified with "X10_given" to be "random" or "fixed".
# "target_acceptance" is the acceptable tolerance rate for acceptance.
# "dig_tol" is the number of decimal places for -X10^2 + X10 - v, and -yt^2 + yt - v to be different from zero. 
# This criterion is important in the numerical method to avoid numerical problems.
Gen_FC_X1_X2 <- function(N, prop_prec, a, b, c, d, v, option = "end", thin = 1, burnin = 0, 
                         X10_given = "random", target_acceptance = 0.3, dig_tol = 15) {
  X1_lower = 0.5 - 0.5 * sqrt(1 - 4 * v)
  X1_upper = 0.5 + 0.5 * sqrt(1 - 4 * v)
  
  # Initialization
  if (X10_given == "random") {
    # Repeated sampling until the condition is met
    while (TRUE) {
      X10 = rBeta.4P(n = 1, l = X1_lower, u = X1_upper, alpha = a, beta = b)
      if (round(-X10^2 + X10 - v, dig_tol) != 0) {
        break  # Exit loop if condition is met
      }
    }
  } else {
    X10 = X10_given
  }
  
  # Variable initialization
  chain = numeric(N)
  chain[1] = X10
  acc_rate = 0
  burnin_accepted = 0
  post_burnin_accepted = 0
  proposals = 0
  alpha = 0
  
  # Main algorithm
  for (k in 2:N) {
    prop_alpha = (chain[k - 1] - X1_lower) / (X1_upper - X1_lower) * prop_prec
    prop_beta = (X1_upper - chain[k - 1]) / (X1_upper - X1_lower) * prop_prec
    
    yt = rBeta.4P(n = 1, l = X1_lower, u = X1_upper, alpha = prop_alpha, beta = prop_beta)
    if (round(-yt^2 + yt - v, dig_tol) != 0) {
      prop_alpha_proposal = (yt - X1_lower) / (X1_upper - X1_lower) * prop_prec
      prop_beta_proposal = (X1_upper - yt) / (X1_upper - X1_lower) * prop_prec
      alpha[k] = exp((a - c - d) * log(yt / chain[k - 1]) + (b - c - d) * log((1 - yt) / (1 - chain[k - 1])) +
                       (d - 1) * log(yt * (1 - yt) - v) - (d - 1) * log(chain[k - 1] * (1 - chain[k - 1]) - v) +
                       log(dBeta.4P(chain[k - 1], l = X1_lower, u = X1_upper, alpha = prop_alpha_proposal, beta = prop_beta_proposal)) -
                       log(dBeta.4P(yt, l = X1_lower, u = X1_upper, alpha = prop_alpha, beta = prop_beta)))
      if (alpha[k] == Inf) {
        alpha[k] = 1
      }
    } else {
      alpha[k] = 0
    }
    
    # Check if alpha[k] is non-numeric
    if (is.nan(alpha[k]) || is.infinite(alpha[k])) {
      stop(paste("Non-numeric alpha detected at iteration", k, 
                 "with proposal", yt, 
                 "and previous chain value", chain[k - 1], "and value of v", v, "The value alpha is", alpha[k]))
    }
    
    if (runif(1) < alpha[k]) {
      chain[k] = yt
      acc_rate = acc_rate + 1
      if (k <= burnin) {
        burnin_accepted <- burnin_accepted + 1
      } else {
        post_burnin_accepted <- post_burnin_accepted + 1
      }
    } else {
      chain[k] = chain[k - 1]
    }
    
    # Adaptive adjustment of the precision parameter during burn-in
    if (k <= burnin && k %% 100 == 0 && target_acceptance!=0) {
      acceptance_rate <- burnin_accepted / 100
      if (acceptance_rate < target_acceptance) {
        prop_prec <- prop_prec * 0.95
      } else if (acceptance_rate > target_acceptance) {
        prop_prec <- prop_prec * 1.05
      }
      burnin_accepted <- 0
    }
    proposals = proposals + 1
  }
  
  # Burn-in and thinning 
  final_chain = chain[(burnin + 1):N]
  thinned_chain = final_chain[seq(1, length(final_chain), by = thin)]
  
  # Acceptance rate
  acc_rate_pos_burnin = (post_burnin_accepted) / (proposals - burnin)
  acc_rate = acc_rate / proposals
  
  # Output options
  if (option == "end") {
    return(list(thinned_chain = tail(thinned_chain, 1), acc_rate = acc_rate, precision = prop_prec,
                acc_rate_pos_burnin = acc_rate_pos_burnin, proposals = proposals))
  } else if (option == "all") {
    return(list(thinned_chain = thinned_chain, acc_rate = acc_rate, DomInt = c(X1_lower, X1_upper),
                ProbAccept = alpha, precision = prop_prec, acc_rate_pos_burnin = acc_rate_pos_burnin, proposals = proposals))
  }
}

#####
# Function to monitor the  Effective Sample Size and acceptance rate for different values of v and precision (prop_prec) provided
#####
# "N" is the sample size to be generated.
# "prop_prec_values" is the list of values that precision can take.
# "a", "b", "c", and "d" are given values for the parameters.
# "v_values" is the list of values that variance can take.
# "thin" and "burnin" are parameters for the function Gen_FC_X1_X2.
# This function constructs two plots comparing the behavior of the Effective Sample Size (ESS)
# and the acceptance rate for different values of precision.
Mon_Measure = function(N, prop_prec_values, a, b, c, d, v_values, thin = 1, burnin = 1) {
  # Get the total number of cores
  num_cores <- detectCores()
  
  # Use half of the available cores
  cl <- makeCluster(num_cores %/% 2)
  
  # Register the cluster for use with foreach
  registerDoParallel(cl)
  
  # Initialize an empty data.frame
  df <- data.frame()
  
  # Use foreach to iterate in parallel
  results_list <- foreach(v = v_values, .combine = 'rbind', .export = c('Gen_FC_X1_X2', 'FC_X1_Given_v'), 
                          .packages = c('coda', 'betafunctions')) %dopar% {
                            tmp_df <- data.frame()
                            for (prop_prec in prop_prec_values) {
                              results <- Gen_FC_X1_X2(N, prop_prec, a, b, c, d, v, option = "all", thin, burnin)
                              result_AR <- results$acc_rate
                              mcmc_obj <- mcmc(results$thinned_chain)
                              result_ESS <- effectiveSize(mcmc_obj)
                              tmp_df <- rbind(tmp_df, data.frame(v = v, Precision = prop_prec, Accept_Rate = result_AR, ESS = result_ESS))
                            }
                            tmp_df
                          }
  
  # Stop the cluster
  stopCluster(cl)
  #print((N - burnin) / thin)
  
  df <- rbind(df, results_list)
  
  # ESS plot for v and precision values
  measure_quantile = quantile(df$ESS, probs = c(0.25, 0.5, 0.75))
  
  plot_ESS = ggplot(df, aes(x = v, y = Precision, z = ESS)) +
    geom_tile(aes(fill = ESS)) +
    scale_fill_gradientn(colors = c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$ESS), measure_quantile[[3]])),
                         name = "Effective Size",
                         breaks = measure_quantile, # Ensures that the minimum and maximum values are displayed in the legend
                         labels = sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) + # Increases the number of values displayed on the x-axis
    scale_y_continuous(breaks = seq(min(prop_prec_values), max(prop_prec_values), by=1)) + # Ensures that all integers are displayed on the y-axis
    labs(title = " ",
         x = "v",
         y = "Precision",
         fill = "Effective Size")
  
  # Acceptance rate plot for v and precision values
  measure_quantile = quantile(df$Accept_Rate, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  plot_AR = ggplot(df, aes(x = v, y = Precision, z = Accept_Rate)) +
    geom_tile(aes(fill = Accept_Rate)) +
    scale_fill_gradientn(colors = c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$Accept_Rate), measure_quantile[[5]])),
                         name = "Acceptance Rate",
                         breaks = measure_quantile, # Ensures that the minimum and maximum values are displayed in the legend
                         labels = sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) + # Increases the number of values displayed on the x-axis
    scale_y_continuous(breaks = seq(min(prop_prec_values), max(prop_prec_values), by=1)) + # Ensures that all integers are displayed on the y-axis
    labs(title = " ",
         x = "v",
         y = "Precision",
         fill = "Acceptance Rate")
  
  grid.arrange(plot_ESS, plot_AR, nrow = 1, ncol = 2, layout_matrix = rbind(c(1, 2)))
}

#####
# Function to monitor the Gelman-Rubin diagnostic (R-hat) for different values of v and precision (prop_prec) provided
#####
# "N" is the sample size to be generated.
# "prop_prec_values" is the list of values that precision can take.
# "a", "b", "c", and "d" are given values for the parameters.
# "v_values" is the list of values that variance can take.
# "thin" and "burnin" are parameters for the function Gen_FC_X1_X2.
# This function constructs a plot comparing the behavior of the Gelman-Rubin diagnostic (R-hat)
# for different values of precision.
Mon_R_Hat = function(N, prop_prec_values, a, b, c, d, v_values, thin = 1, burnin = 1) {
  # Get the total number of cores
  num_cores <- detectCores()
  
  # Use half of the available cores
  cl <- makeCluster(num_cores %/% 2)
  
  # Register the cluster for use with foreach
  registerDoParallel(cl)
  
  # Initialize an empty data.frame
  df <- data.frame()
  
  # Use foreach to iterate in parallel
  results_list <- foreach(v = v_values, .combine = 'rbind', .export = c('Gen_FC_X1_X2'), 
                          .packages = c('coda', 'betafunctions')) %dopar% {
                            tmp_df <- data.frame()
                            for (prop_prec in prop_prec_values) {
                              sample1 = Gen_FC_X1_X2(N, prop_prec, a, b, c, d, v, option = "all", thin=thin , X10_given = "random", burnin)
                              sample2 = Gen_FC_X1_X2(N, prop_prec, a, b, c, d, v, option = "all", thin=thin , X10_given = "random", burnin)
                              sample3 = Gen_FC_X1_X2(N, prop_prec, a, b, c, d, v, option = "all", thin=thin , X10_given = "random", burnin)
                              Gelm_Rud = gelman.diag(list(mcmc(sample1$thinned_chain), mcmc(sample2$thinned_chain),
                                                          mcmc(sample3$thinned_chain)))$psrf[1]
                              tmp_df <- rbind(tmp_df, data.frame(v = v, Precision = prop_prec, Gelman_Rubin = Gelm_Rud))
                            }
                            tmp_df
                          }
  
  # Stop the cluster
  stopCluster(cl)
  
  df <- rbind(df, results_list)
  
  # R-hat plot for v and precision values
  measure_quantile = quantile(df$Gelman_Rubin, probs = c(0.2, 0.94, 0.96, 0.98, 1))
  
  ggplot(df, aes(x = v, y = Precision, z = Gelman_Rubin)) +
    geom_tile(aes(fill = Gelman_Rubin)) +
    scale_fill_gradientn(colors = c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$Gelman_Rubin), measure_quantile[[5]])),
                         name = "R-hat",
                         breaks = measure_quantile, # Ensures that the minimum and maximum values are displayed in the legend
                         labels = sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) + # Increases the number of values displayed on the x-axis
    scale_y_continuous(breaks = seq(min(prop_prec_values), max(prop_prec_values), by = 1)) + # Ensures that all integers are displayed on the y-axis
    labs(title = " ",
         x = "v",
         y = "Precision",
         fill = "R-hat")
  
  #return(list(plot_H))
}

#####
# Function that plots the histogram, density, trace, and convergence control using the average.
#####
# "nameaxisy" is the name of the vertical axis for trace and convergence monitoring.
# "width" is the width for the confidence intervals of convergence monitoring.
# "lscatt" is an increment to the minimum value generated. It allows plotting the line at an "lscatt" distance from the trace to enhance visualization.
# "uscatt" is an increment to the maximum value generated. It allows plotting the line at an "uscatt" distance from the trace to enhance visualization.
Graphs = function(dataset, nameaxisy, width = 10, lscatt = 0.05, uscatt = 0.05) {
  # Histogram with density
  l = length(dataset[, 1])
  hist = ggplot(dataset, aes(x = dataset[, 1])) + 
    geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "white") +
    geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
    labs(title = "Histogram and Density") + ylab("Density") +
    xlab(if(nameaxisy == "X2") { expression(X[2]^(t)) } 
         else if(nameaxisy == "X1") { expression(X[1]^(t)) }
         else if(nameaxisy == "Y1") { expression(Y[1]^(t)) }
         else if(nameaxisy == "Y2") { expression(Y[2]^(t)) }
         else { substitute(va, list(va = as.name(nameaxisy))) }) +
    theme(plot.title = element_text(size = 11))
  
  # Trace plot with maximum and minimum
  trace = ggplot(dataset, aes(x = 1:l, y = dataset[, 1])) +
    geom_line() + xlab("t") +
    ylab(if(nameaxisy == "X2") { expression(X[2]^(t)) } 
         else if(nameaxisy == "X1") { expression(X[1]^(t)) }
         else if(nameaxisy == "Y1") { expression(Y[1]^(t)) }
         else if(nameaxisy == "Y2") { expression(Y[2]^(t)) }
         else { substitute(va, list(va = as.name(nameaxisy))) }) +
    ylim(c(min(dataset) - lscatt, max(dataset) + uscatt)) +
    geom_hline(aes(yintercept = min(dataset[, 1])), colour = "red", linetype = 2) +
    geom_text(aes(l - l / 10, min(dataset[, 1]), label = round(min(dataset[, 1]), 3), vjust = 2), colour = "red") +
    geom_hline(aes(yintercept = max(dataset[, 1])), colour = "red", linetype = 2) +
    geom_text(aes(l - l / 10, max(dataset[, 1]), label = round(max(dataset[, 1]), 3), vjust = -1), colour = "red") +
    labs(title = substitute(list("Trace of the random sample of size", n), list(n = l))) +
    theme(plot.title = element_text(size = 11))
  
  # Acf plot
  alfa = 0.05
  lim = qnorm((1 - alfa / 2)) / sqrt(l)
  acf_values = acf(dataset, plot = FALSE)
  acf_data = data.frame(Lag = acf_values$lag[-1],  # Remove the first lag value (always 0)
                        ACF = acf_values$acf[-1])  # Remove the first ACF value (always 1)
  acfplot = ggplot(acf_data, aes(x = Lag, y = ACF)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = c(lim, -lim), linetype = "dashed") +
    labs(title = "Autocorrelation Function",
         x = "Lag",
         y = "ACF") +
    theme_minimal() + theme(plot.title = element_text(size = 11))  
  
  # Convergence control using averaging
  dataset$estintden = cumsum(dataset[, 1]) / (1:l)
  dataset$esterrden = sqrt(cumsum((dataset[, 1] - dataset$estintden)^2)) / (1:l)
  
  mean_X1_X2 = ggplot(dataset, aes(x = 1:l, y = estintden)) + geom_line() +
    geom_line(aes(x = 1:l, y = estintden - 1.95 * esterrden, colour = "Upper")) +
    geom_line(aes(x = 1:l, y = estintden + 1.95 * esterrden, colour = "Lower")) +
    ylim(mean(dataset$estintden) + width * c(-dataset$esterrden[l], dataset$esterrden[l])) +
    ylab(if(nameaxisy == "X2") { expression(X[2]^(t)) } 
         else if(nameaxisy == "X1") { expression(X[1]^(t)) }
         else if(nameaxisy == "Y1") { expression(Y[1]^(t)) }
         else if(nameaxisy == "Y2") { expression(Y[2]^(t)) }
         else { substitute(va, list(va = as.name(nameaxisy))) }) +
    xlab("t") + geom_hline(yintercept = mean(dataset[, 1]), colour = "red", linetype = 2) +
    geom_text(aes(l - l / 10, mean(dataset[, 1]), label = round(mean(dataset[, 1]), 3), vjust = -2), colour = "red") +
    labs(title = "Convergence Control using Averaging", color = "Bounds") + 
    scale_shape_discrete(name = " ") + theme(plot.title = element_text(size = 11), legend.position = "top", legend.box = "horizontal")
  
  # Plots of histogram, trace, convergence control, and acf.
  grid.arrange(hist, trace, mean_X1_X2, acfplot, 
                      ncol = 2, nrow = 2, widths = c(4, 4), heights = c(2, 2), layout_matrix = rbind(c(1, 2), c(3, 4)))
}

##########################################################
##########################################################
# Gibbs Sampling
##########################################################
##########################################################
# "N1": Gibbs Sampling sample size
# "N2": Random walks sample size for full conditional
# "a", "b", "c", and "d" are given values for parameters
# "thin" is the thinning interval for Random Walks. Every "thin" generated sample is stored to reduce autocorrelation.
# The seed type can be specified with "X10_given" as "random" or "fixed".
# "lower_epsilon": lower limit for the conditional distribution of the variance given a value for the mean.
# "dig_tol" in Gen_FC_X1_X2: number of decimal places for -X10^2 + X10 - v, and -yt^2 + yt - v to be different from zero.
Gen_Joint_Dist = function(N1, N2, prop_prec, a, b, c, d, thin = 1, X10_given = "random", lower_epsilon = 0, dig_tol = 15) {
  SampleGen = matrix(data = NA, nrow = N1, ncol = 2, dimnames = list(NULL, c("X2", "X1")))
  SampleGen = as.data.frame(SampleGen)
  SampleGen$X1[1] = rbeta(1, a, b)
  SampleGen$X2[1] = rBeta.4P(1, l = 0, u = SampleGen$X1[1] * (1 - SampleGen$X1[1]), alpha = c, beta = d)
  for (t in 2:N1) {
    # Generate X1[t] using the Metropolis-Hastings algorithm
    SampleGen$X1[t] = Gen_FC_X1_X2(N2, prop_prec, a, b, c, d, SampleGen$X2[t-1], option = "end", thin, burnin = 0, X10_given, dig_tol = dig_tol)$thinned_chain
    # Generate X2[t] using the conditional distribution
    SampleGen$X2[t] = rBeta.4P(1, l = lower_epsilon, u = SampleGen$X1[t] * (1 - SampleGen$X1[t]), alpha = c, beta = d)
  }
  return(SampleGen)
}

# Example
# trial0 = Gen_Joint_Dist(N1 = 10, N2 = 5, prop_prec = 3, a = 3, b = 2.5, c = 4, d = 6, thin = 1, X10_given = "random", lower_epsilon = 0)

##########################################################
##########################################################
# Generalization of Tovar's method to obtain hyperparameter values
##########################################################
##########################################################
# "q1" and "q2" are values obtained from a person considered an expert on the topic of interest
# "low" and "upp" are the values where the parameter of interest remains
# "alp" is the confidence level that the expert has that the interval (low, upp) contains the true value of the parameter.
Mtovar_vs2 = function(q1, q2, low, upp, alp) {
  tht0 = (q1 + q2) / 2  # Mean of the expert's interval
  w = (tht0 - low) / (upp - tht0)  # Weighting factor based on the expert's interval
  sig = sqrt(alp) * (q1 - tht0)  # Adjusted standard deviation based on the confidence level
  b = ((upp - low)^2 * w - ((w + 1)^2 * sig^2)) / ((w + 1)^3 * sig^2)  # Hyperparameter b calculation
  a = w * b  # Hyperparameter a calculation
  return(list(a = a, b = b, c = tht0))  # Return the hyperparameters and the mean
}

##########################################################
##########################################################
# Joint moments of order l=l1+l2 for proposed prior distribution
##########################################################
##########################################################
# l1 is the marginal order for alpha
# l2 is the marginal order for beta
# a, b, c, and d are hyperparameter values.
Mom_Prior_Dist = function(l1, l2, a, b, c, d) {
  beta(c - l1 - l2, l1 + l2 + d) * beta(l1 + a, l2 + b)
}

##########################################################
##########################################################
# Comparison of analytic and numeric results
##########################################################
##########################################################
# data1 and data2 are datasets generated by the Gibbs sampling method mentioned earlier.
# thin and burnin are parameters applied to select the data elements to be used for determining numerical measures.
# digits: number of decimal places for numerical measures.
# a, b, c, and d are hyperparameter values of the proposed distribution.
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
    
    return(list(Numerical = Numerical_results))
    
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

##########################################################
##########################################################
# Analytical results
##########################################################
##########################################################
# "a", "b", "c", and "d" are hyperparameter values.
# "digits" is the number of decimal places for the analytical measures.
Measure_Analy = function(a, b, c, d, digits) {
  K = Mom_Prior_Dist(0, 0, a, b, c, d)
  Analytic_results = round(data.frame(
    "Mean.1" = Mom_Prior_Dist(1, 0, a, b, c, d) / K,
    "Var.1" = Mom_Prior_Dist(2, 0, a, b, c, d) / K - (Mom_Prior_Dist(1, 0, a, b, c, d) / K)^2,
    "Mean.2" = Mom_Prior_Dist(0, 1, a, b, c, d) / K,
    "Var.2" = Mom_Prior_Dist(0, 2, a, b, c, d) / K - (Mom_Prior_Dist(0, 1, a, b, c, d) / K)^2,
    "Cov" = Mom_Prior_Dist(1, 1, a, b, c, d) / K - (Mom_Prior_Dist(1, 0, a, b, c, d) / K) * Mom_Prior_Dist(0, 1, a, b, c, d) / K,
    "K" = K
  ), digits)
  
  return(Analytic_results)
}

#####################
## Obtaining hyperparameters from two approaches.
## The first approach is empirical Bayes: it uses the Bootstrap quantile interval.
## The second approach is subjective: it uses quantile intervals from an expert's opinion.
#####################

# ssample: original sample 
# r_boostrap: number of resamples.
# q_boostrap: Bootstrap quantiles
# option_mu: method to obtain hyperparameters a and b for the mean mu, can be "moments" or "tovar".
# sig_mu: significance level for Tovar's method for obtaining hyperparameters for the mean mu.
# bound_var: method to define the upper limit for the variance, can be "min", "mean", "max".
# sig_var: significance level for Tovar's method for obtaining hyperparameters for the variance.
# digits: number of decimal places for the Bootstrap quantile interval for the mean and variance.
# graphs_boot: logical indicator to generate a histogram, T or F.
# Q_E_mu and Q_E_cv: quantiles for the mean and variance obtained from the expert.

Hyperparameters = function(ssample, r_boostrap = 100, q_boostrap = c(0.025, 0.975), option_mu = "moments", 
                           sig_mu = 0.05, bound_var = "max", sig_var = 0.05, digits = 4, 
                           graphs_boot = F, Q_E_mu = 0, Q_E_cv = 0) {
  if (r_boostrap != 0) {
    n_sample = length(ssample)
    boot = matrix(sample(ssample, size = r_boostrap * n_sample, replace = T), nrow = n_sample, ncol = r_boostrap)
    boots_mean = round(apply(boot, 2, mean), digits)
    boots_sd = round(apply(boot, 2, sd), digits)
    boots_cv = round(boots_sd / boots_mean, digits)
    
    # Select the quantiles associated with q_boostrap
    # For the mean
    quantile_mu = round(quantile(boots_mean, probs = q_boostrap), digits)
    # For the CV
    quantile_cv = quantile(boots_cv, probs = q_boostrap)
  } else if (r_boostrap == 0) {
    quantile_mu = Q_E_mu
    quantile_cv = Q_E_cv
  }
  
  # For the variance
  quantile_var = round((quantile_cv * mean(quantile_mu))^2, digits)
  
  #############################
  # Interval for the mean
  #############################
  if (option_mu == "moments") {
    portion = (mean(quantile_mu) * (1 - mean(quantile_mu)) / ((quantile_mu[[2]] - quantile_mu[[1]]) / 4)^2 - 1)
    hiper_mean = data.frame("a" = mean(quantile_mu) * portion, "b" = (1 - mean(quantile_mu)) * portion)
  } else if (option_mu == "tovar") {
    hiper_mean = Mtovar_vs2(quantile_mu[[1]], quantile_mu[[2]], 0, 1, sig_mu)
  }
  
  #############################
  # Interval for the variance
  #############################
  bound_var_value = if (bound_var == "min") { min(quantile_mu) * (1 - min(quantile_mu)) }
  else if (bound_var == "mean") { mean(quantile_mu) * (1 - mean(quantile_mu)) }
  else if (bound_var == "max") { max(quantile_mu) * (1 - max(quantile_mu)) }
  
  hiper_var = Mtovar_vs2(quantile_var[[1]], quantile_var[[2]], 0, bound_var_value, sig_var)
  
  #############################
  # Histogram and Bootstrap density
  #############################
  
  if (graphs_boot == T) {
    hist_orig = ggplot(as.data.frame(ssample), aes(x = ssample)) + 
      geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "white") +
      geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
      labs(title = "Original Sample") + ylab("Density") +
      xlab(substitute(va, list(va = "X")))
    if (r_boostrap != 0) {
      hist_boot_mean = ggplot(as.data.frame(boots_mean), aes(x = boots_mean)) + 
        geom_histogram(aes(y = ..density..), colour = 1, fill = "white") +
        geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
        labs(title = "Bootstrap for the Mean") + ylab("Density") +
        xlab(substitute(va, list(va = "Mean of X")))
      
      hist_boot_cv = ggplot(as.data.frame(boots_cv), aes(x = boots_cv)) + 
        geom_histogram(aes(y = ..density..), colour = 1, fill = "white") +
        geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
        labs(title = "Bootstrap for the CV") + ylab("Density") +
        xlab(substitute(va, list(va = "CV of X")))
      
      grid.arrange(hist_orig, hist_boot_cv, hist_boot_mean, 
                   ncol = 2, nrow = 2, widths = c(2, 2), heights = c(2, 2), layout_matrix = rbind(c(1, 1), c(2, 3)))
    } else if (r_boostrap == 0) {
      hist_orig
    }
  }
  
  return(list(Q_mean = quantile_mu, Q_cv = quantile_cv, Q_var = quantile_var, hiper_mean = hiper_mean, hiper_var = hiper_var,
              Var_Upper_Bound = bound_var_value))
}

#####################
## Posterior estimation for the parameters alpha and beta of the Beta distribution
## using importance sampling.
#####################
# ssample: set of ordered pairs representing the observed sample.
# N: sample size generated by the Gibbs sampling method.
# N_FC: sample size for the Metropolis-Hastings random walk method for the conditional distribution X2 given X1.
# Precision: precision implemented by the instrumental distribution in the Metropolis-Hastings random walk method for the conditional distribution X2 given X1.
# a, b, c, and d are the values of the hyperparameters generated by the function of the same name.
# thin1: parameter for the Gen_Joint_Dist function.
# thin2: parameter used to set the step size stored from the chain of size N generated.
# burnin: number of discarded samples for posterior estimates.
# dig_tol: number of decimal places to avoid numerical issues in the sample generated by the conditional distribution involved in Gen_Joint_Dist.

Est_Post = function(ssample, N, N_FC, Precision, a, b, c, d, thin1, thin2, burnin, dig_tol) {
  Example_Joint_Dist = Gen_Joint_Dist(N1 = N, N2 = N_FC, prop_prec = Precision, a, b, c, d, thin = thin1, X10_given = "random", dig_tol = dig_tol)
  piece = (Example_Joint_Dist$X1[seq((burnin + 1), N, by = thin2)] * (1 - Example_Joint_Dist$X1[seq((burnin + 1), N, by = thin2)]) /
             Example_Joint_Dist$X2[seq((burnin + 1), N, by = thin2)] - 1)
  sample_alpha = Example_Joint_Dist$X1[seq((burnin + 1), N, by = thin2)] * piece
  sample_beta = (1 - Example_Joint_Dist$X1[seq((burnin + 1), N, by = thin2)]) * piece
  x0 = prod(ssample)
  y0 = prod(1 - ssample)
  
  zeta=exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)-(length(ssample))*log(beta(sample_alpha,sample_beta)))
  
  marginalike=sum(zeta)/length(sample_beta)
  
  meanpalpha=sum(zeta*exp(log(sample_alpha)-log(marginalike)))
  varpalpha=sum(zeta*exp(2*log(sample_alpha)-log(marginalike)))-meanpalpha^2
  
  meanpbeta=sum(zeta*exp(log(sample_beta)-log(marginalike))) 
  varpbeta=sum(zeta*exp(2*log(sample_beta)-log(marginalike)))-meanpbeta^2
  
  covalphabeta=sum(zeta*exp(log(sample_alpha)+log(sample_beta)-log(marginalike))) - meanpalpha*meanpbeta
  
  Descriptivo = Measure_Diagnostic(data1 = Example_Joint_Dist$X1, data2 = Example_Joint_Dist$X2, var = "transform",
                                   digits = 4, burnin, thin = thin2, a, b, c, d)
  return(list(EstPost=data.frame(Prior_Lik=marginalike,P_A=meanpalpha,P_B=meanpbeta,VP_A=varpalpha,VP_B=varpbeta,CP_AB=covalphabeta),
              Descriptivo=Descriptivo,
              SA=sample_alpha,SB=sample_beta,SX1=Example_Joint_Dist$X1,SX2=Example_Joint_Dist$X2))
}



#####################
# Simulation study to compare and monitor the behavior of posterior estimates
# produced by different sets of hyperparameters (obtained from the empirical Bayes approach and the subjective approach),
# and observed samples of different sizes.
#####################
# N: sample size generated by the Gibbs sampling method.
# N_FC: sample size for the Metropolis-Hastings random walk method for the conditional distribution X2 given X1.
# Precision: precision implemented by the instrumental distribution in the Metropolis-Hastings random walk method for the conditional distribution X2 given X1.
# a, b, c, and d are the values of the hyperparameters generated by the function of the same name.
# thin1: parameter for the Gen_Joint_Dist function.
# The seed type can be specified with "X10_given" as "random" or "fixed".
# dig_tol: number of decimal places to avoid numerical issues in the sample generated by the conditional distribution involved in Gen_Joint_Dist.
# thin2: parameter used to set the step size stored from the chain of size N generated.
# burnin: number of discarded samples for posterior estimates.
# n_sample: list of sample sizes for samples generated from the Beta distribution using the true parameter values, alpha_real and beta_real.
# alpha_real and beta_real are the true values used to set the simulation scenario.
# N_Iter_Sim: Number of times the estimation is repeated for each sample size.

Sim_study = function(N, N_FC, prop_prec, a, b, c, d, thin1, X10_given = "random", dig_tol = 15, thin2, burnin, n_sample, alpha_real, beta_real, N_Iter_Sim,
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
  # Generate prior sample
  Sample_Prior_Hip = Gen_Joint_Dist(N1 = N, N2 = N_FC, prop_prec = 3, a, b, c, d, thin = thin1, X10_given = "random", dig_tol = dig_tol)
  piece = (Sample_Prior_Hip$X1[seq((burnin + 1), N, by = thin2)] * (1 - Sample_Prior_Hip$X1[seq((burnin + 1), N, by = thin2)]) /
             Sample_Prior_Hip$X2[seq((burnin + 1), N, by = thin2)] - 1)
  sample_alpha = Sample_Prior_Hip$X1[seq((burnin + 1), N, by = thin2)] * piece
  sample_beta = (1 - Sample_Prior_Hip$X1[seq((burnin + 1), N, by = thin2)]) * piece
  
  # Descriptive statistics for the prior sample
  Descriptive_Sample_Prior = Measure_Diagnostic(data1 = Sample_Prior_Hip$X1, data2 = Sample_Prior_Hip$X2, var = "transform",
                                                digits = 4, burnin, thin = thin2, a, b, c, d)
  
  Iter_Alpha = NA
  Iter_Beta = NA
  VIter_Alpha = matrix(data = NA, nrow = 0, ncol = 12, dimnames = list(NULL, c("Min", "Q2.5", "Q50", "Q97.5", "Max", "Mean", "Var", "Bias", "MSE", "SampleSize","Coverage","Length")))
  VIter_Beta = VIter_Alpha
  RC_A=matrix(data = NA, nrow = 0, ncol = 2, dimnames = list(NULL, c("Q2.5", "Q97.5")))
  RC_B=RC_A
  for (j in 1:length(n_sample)) {
    for (i in 1:N_Iter_Sim) {
      sampe_prueba = rbeta(n_sample[j], alpha_real, beta_real)
      x0 = prod(sampe_prueba)
      y0 = prod(1 - sampe_prueba)
      
      zeta=exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)-(n_sample[j])*log(beta(sample_alpha,sample_beta)))
      marginalike=sum(zeta)
      Iter_Alpha[i] = sum(zeta*exp(log(sample_alpha)-log(marginalike)))
      Iter_Beta[i] = sum(zeta*exp(log(sample_beta)-log(marginalike))) 

      # Sampling importance resampling
      selected_elements=sample(1:length(sample_alpha),size = sample_size_IS,replace = T,prob = zeta/marginalike)
      selected_alpha=sample_alpha[selected_elements]
      selected_beta=sample_beta[selected_elements]
      
      # Regions of credibility
      RC_A = rbind(RC_A, quantile(selected_alpha, probs = c(0.025, 0.975), na.rm = T))
      RC_B = rbind(RC_B, quantile(selected_beta, probs = c(0.025, 0.975), na.rm = T))
      
    }
    VIter_Alpha = rbind(VIter_Alpha, c(quantile(Iter_Alpha, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = T), mean(Iter_Alpha), var(Iter_Alpha), (alpha_real - mean(Iter_Alpha)),
                                       (var(Iter_Alpha) + (alpha_real - mean(Iter_Alpha))^2), n_sample[j],covprob(RC_A,alpha_real),meanlenght(RC_A)))
    VIter_Beta = rbind(VIter_Beta, c(quantile(Iter_Beta, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = T), mean(Iter_Beta), var(Iter_Beta), (beta_real - mean(Iter_Beta)),
                                     (var(Iter_Beta) + (beta_real - mean(Iter_Beta))^2), n_sample[j], covprob(RC_B,beta_real), meanlenght(RC_B)))
  }
  
  return(list(Result_Alpha = as.data.frame(VIter_Alpha), Result_Beta = as.data.frame(VIter_Beta), Descriptive_Sample_Prior = Descriptive_Sample_Prior))
}




############
# Graphs for the results obtained from the simulation study.
# These graphs allow monitoring the behavior of the estimates,
# using bias and MSE as indicators.
############
# Individual Graph
# Data: Data contains a data frame with the necessary variables, including SampleSize (x-axis).
# lim_x: limits for the x-axis.
# value_real: true value of the parameter in the simulation study.
# y_text: name of the y-axis.
# title_text: title of the graph.
Individual_Graphs = function(Data, lim_x, value_real, y_text, title_text) {
  # Graph of the mean with bands constructed with the quantiles of the iterations for each sample size.
  Mean_Graph = ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = Q2.5, color = "Q2.5")) +
    geom_line(aes(y = Mean, color = "Mean")) +
    geom_line(aes(y = Q97.5, color = "Q97.5")) +
    labs(title = title_text,
         x = "Sample Size",
         y = y_text) + 
    scale_color_manual(name = " ",
                       values = c("Q2.5" = "blue", "Mean" = "green", "Q97.5" = "red")) +
    scale_x_continuous(breaks = lim_x) +  # Modification to show values at specified intervals
    geom_hline(yintercept = value_real, linetype = "dashed", color = "red") +
    theme_minimal()
  
  # Graph of the MSE
  Mse_Graph = ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = MSE, color = "MSE")) +
    labs(title = " ",
         x = "Sample Size",
         y = "MSE") +
    scale_color_manual(name = " ",
                       values = c("MSE" = "blue")) +
    scale_x_continuous(breaks = lim_x) +  # Modification to show values at specified intervals
    theme_minimal()
  
  # Graph of the Bias
  Bias_Graph = ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = Bias, color = "Bias")) +
    labs(title = " ",
         x = "Sample Size",
         y = "Bias") +
    scale_color_manual(name = " ",
                       values = c("Bias" = "black")) +
    scale_x_continuous(breaks = lim_x) +  # Modification to show values at specified intervals
    theme_minimal()
  
  grid.arrange(Mean_Graph, Mse_Graph, Bias_Graph, nrow = 2, ncol = 2, layout_matrix = rbind(c(1, 1), c(2, 3)))
}




####
# Comparison of joint functions
####
Comparison_Hyper = function(data, value_real, lim_x, title_text, y_Text) {
  # Mean
  Comparacion_Mean = ggplot(data, aes(x = SampleSize, y = Mean, color = Method, shape = Method, linetype = Method)) +
    geom_line() +
    geom_point(size = 2) +
    labs(title = title_text,
         x = " ",
         y = y_Text) +
    scale_x_continuous(breaks = lim_x) +
    geom_hline(yintercept = value_real, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(legend.position = "none") +  # Remove legend
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Different shapes for points
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # Bias
  Comparacion_Bias = ggplot(data, aes(x = SampleSize, y = Bias, color = Method, shape = Method, linetype = Method)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = " ",
         y = "Bias") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() +
    theme(legend.position = "none") +  # Remove legend
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Different shapes for points
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # MSE
  Comparacion_Mse = ggplot(data, aes(x = SampleSize, y = MSE, color = Method, shape = Method, linetype = Method)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = " ",
         y = "MSE") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() + 
    theme(legend.position = "none") +  # Remove legend
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Different shapes for points
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # Coverage
  Comparacion_Coverage = ggplot(data, aes(x = SampleSize, y = Coverage, color = Method, shape = Method, linetype = Method)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = " ",
         y = "Coverage") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() + 
    theme(legend.position = "none") +  # Remove legend
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Different shapes for points
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # Length
  Comparacion_Length = ggplot(data, aes(x = SampleSize, y = Length, color = Method, shape = Method, linetype = Method)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = "Sample Size",
         y = "Length") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() + 
    theme(legend.position = "none") +  # Remove legend
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Different shapes for points
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # Function to extract the legend from a plot
  g_legend <- function(a.gplot) {
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Extract the legend from the mean plot
  legend <- g_legend(Comparacion_Mean + theme(legend.position = "right", legend.title = element_blank()))
  
  # Combine the plots without legend with the legend on the right
  grid.arrange(
    arrangeGrob(Comparacion_Mean, Comparacion_Bias, Comparacion_Mse, Comparacion_Coverage, Comparacion_Length, ncol = 1, heights = c(9, 9, 9, 9, 9)),
    legend, 
    ncol = 2, 
    widths = c(7, 1)
  )
}


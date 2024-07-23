###########################################################################
# File: distribution.R                                                    #
#                                                                         #
# Description:                                                            #
# This file contains a set of functions and routines used for generating  #
# random samples from a joint distribution for the shape parameters of    #
# the beta distribution, and for evaluating the convergence of the        #
# generated samples.                                                      #
#                                                                         #
# Author: Llerzy Torres Ome                                               #
# Creation Date: July 23, 2024                                            #
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
#library(xlsx)
library(tidyr) # provides tools for tidying up data and especially useful for its pipe operator (%>%), which streamlines data manipulation and transformation.
library(plotly) # used for interactive visualization of contour plots
library(coda)
library(foreach)
library(doParallel)
library(betafunctions)
#library(rio)
library(openxlsx)
library(xtable)
##########################################################

##########################################################
## Función de densidad de probabilidad propuesta
##########################################################
prior=function(alph,bet,a,b,c,d){
  #1/(beta(a,b)*beta(c,d))*alph^(a-1)*bet^(b-1)*(alph+bet)^(d-a-b)*(alph+bet+1)^(-c-d)
  return( exp((a-1)*log(alph)+(b-1)*log(bet)+(d-a-b)*log(alph+bet)-(c+d)*log(alph+bet+1)))
}


##########################################################
##########################################################
# Metropolis Hasting Method using Random Walks
# for conditional distribution X1 given X2 squared
##########################################################
##########################################################
# Full Conditional (FC) of X1 given X2 squared (v)
# "(a,b,c,d)" is the vector of parameters
# "X1" is within (0,1) and "X2=v" is less than X1(1-X1)
FC_X1_given_v=function(X1,a,b,c,d,v){
  #result1=X1^(a-c-d)*(1-X1)^(b-c-d)*(X1*(1-X1)-v)^(d-1)
  exp( (a-c-d)*log(X1)+(b-c-d)*log(1-X1)+(d-1)*log(X1*(1-X1)-v))
  #return(list(result1=result1,result2=result2))
  #return(result2)
} 

#####
# Plot of FC for given values of a, b, c, d, and three values of v.
#####
# "v1", "v2" and "v3" are given values of the variance
# "v1name", "v2name", "s3name" are the names of each plot.
# "(ae,be,ce,de)" is the vector of given parameters values.
Graph_Fc_X1=function(v1,v1name,v2,v2name,v3,v3name,ae,be,ce,de){
  ggplot()+xlim(c(0,1))+
    geom_function(fun=function(X1) mapply(FC_X1_given_v,X1,a=ae,b=be,c=ce,d=de,v=v1),lwd=1,
                  linetype=1,aes(col=v1name))+
    geom_function(fun=function(X1) mapply(FC_X1_given_v,X1,a=ae,b=be,c=ce,d=de,v=v2),lwd=1,
                  linetype=1,aes(col=v2name))+
    geom_function(fun=function(X1) mapply(FC_X1_given_v,X1,a=ae,b=be,c=ce,d=de,v=v3),lwd=1,
                  linetype=1,aes(col=v3name))+
    labs(title="Condicional completa de X1 dado X2",
         caption = substitute(
           list("Gráfica de",f(X1/v)==X1^(a-c-d)*(1-X1)^(b-c-d)*(X1*(1-X1)-v)^(d-1),
                "con", a==ae,b==be,c==ce,d==de),list(ae=ae,be=be,ce=ce,de=de)))+
    xlab(expression(X1))+ylab(expression(f(X1/X2==v))) + 
    scale_colour_manual(values = c("red","black","purple"),name="X2=v")
}
#####

#####
# Metropolis Hasting using Random Walks Algorithm 
#####
# "N" is the  sample size to be generated
# "prop_prec" is the precision set for the algorithm. 
# "a","b","c" and "d" are given values for parameters
# "v" is the given value for X2
# "option" allows you to select the entire sample ("all") or just the last value generated ("end").
# "thin" is the thinning interval for MCMC. Every "thin" generated samples, one is stored to reduce autocorrelation.
# "burn in" is the number of iterations to discard.
# the seed type can be specified with "X10_given" to be "random" or "fixed".

Gen_FC_X1_X2 <- function(N, prop_prec, a, b, c, d, v, option = "end", thin = 1, burnin = 0, 
                         X10_given = "random",target_acceptance=0.3,dig_tol=15) {
  X1_lower=0.5-0.5*sqrt(1-4*v)
  X1_upper=0.5+0.5*sqrt(1-4*v)
  # Initialization
  if (X10_given == "random"){
    # Repeated sampling until the condition is met
    while(TRUE) {
      X10 = rBeta.4P(n = 1,l = X1_lower,u = X1_upper, alpha = a, beta = b) 
      if (round(-X10^2 + X10 - v,dig_tol)!=0) {
        break  # Exit loop if condition is met
      }
    }
  } else {
    X10 = X10_given
  }
  
  # Variable initialization
  chain = numeric(N)
  chain[1] = X10
  acc_rate = 0; burnin_accepted=0; post_burnin_accepted=0;
  proposals = 0
  alpha=0
  # Main algorithm
  for (k in 2:N){
    prop_alpha = (chain[k-1]-X1_lower)/(X1_upper-X1_lower) * prop_prec
    prop_beta = (X1_upper-chain[k-1])/(X1_upper-X1_lower) * prop_prec
    
    # Repeated sampling until the condition is met
    #while(TRUE) {
    #  yt = rBeta.4P(n = 1, l = X1_lower, u = X1_upper, alpha = prop_alpha, beta = prop_beta)
    #  if (round(-yt^2 + yt - v,15)!=0) {
    #    break  # Exit loop if condition is met
    #  }
    #  alpha[k]=0
    #  chain[k] = chain[k-1]
    #  proposals = proposals + 1
    #  k=k+1
    #}
    #alpha = exp(log(FC_X1_given_v(yt, a, b, c, d, v)) - log(FC_X1_given_v(chain[k-1], a, b, c, d, v)))
    #alpha[k]=exp((a-c-d)*log(yt/chain[k-1])+(b-c-d)*log((1-yt)/(1-chain[k-1]))+
    #(d-1)*log((yt*(1-yt)-v)/(chain[k-1]*(1-chain[k-1])-v)))
    yt = rBeta.4P(n = 1, l = X1_lower, u = X1_upper, alpha = prop_alpha, beta = prop_beta)
    if(round(-yt^2 + yt - v,dig_tol)!=0){
      prop_alpha_proposal =(yt-X1_lower)/(X1_upper-X1_lower) * prop_prec
      prop_beta_proposal = (X1_upper-yt)/(X1_upper-X1_lower) * prop_prec
      alpha[k]=exp((a-c-d)*log(yt/chain[k-1])+(b-c-d)*log((1-yt)/(1-chain[k-1]))+
                     (d-1)*log((yt*(1-yt)-v))-(d-1)*log((chain[k-1]*(1-chain[k-1])-v))+
                     log(dBeta.4P(chain[k-1],l = X1_lower,u = X1_upper,alpha = prop_alpha_proposal,beta = prop_beta_proposal))
                         -log(dBeta.4P(yt,l = X1_lower,u = X1_upper,alpha = prop_alpha,beta = prop_beta)))
      if(alpha[k]==Inf){
        alpha[k]=1
      }
    }else {
      alpha[k]=0
    }
    
    # Check if alpha[k] is non-numeric
    if (is.nan(alpha[k]) || is.infinite(alpha[k])) {
      stop(paste("Non-numeric alpha detected at iteration", k, 
                 "with proposal", yt, 
                 "and previous chain value", chain[k - 1],"and value of v",v,"The value alpha is", alpha[k]))
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
      chain[k] = chain[k-1]
    }
    
    # Ajuste adaptativo del parámetro de precisión durante el burn-in
    if (k <= burnin && k %% 100 == 0) {
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
  acc_rate_pos_burnin=(post_burnin_accepted)/(proposals-burnin)
  acc_rate = acc_rate / proposals
  # Output options
  if (option == "end") {
    return(list(thinned_chain = tail(thinned_chain, 1), acc_rate = acc_rate,precision=prop_prec,
                acc_rate_pos_burnin=acc_rate_pos_burnin,proposals=proposals))
  } else if (option == "all") {
    return(list(thinned_chain = thinned_chain, acc_rate = acc_rate, DomInt=c(X1_lower,X1_upper),
                ProbAccept=alpha,precision=prop_prec,acc_rate_pos_burnin=acc_rate_pos_burnin,proposals=proposals))
  }
}

#####
# Function to monitor the acceptance rate for different values of v and precision (prop_prec) provided
#####
# "N" is the  sample size to be generated
# "prop_prec_values" is the list of values that precision can take.
# "a","b","c" and "d" are given values for parameters
# "v_values" is the list of values that variance can take.

Mon_measure = function(N, prop_prec_values, a, b, c, d, v_values,thin = 1,burnin=1) {
  # Obtener el número total de núcleos
  num_cores <- detectCores()
  
  # Usar la mitad de los núcleos disponibles
  cl <- makeCluster(num_cores %/% 2)
  
  # Registrar el clúster para su uso con foreach
  registerDoParallel(cl)
  
  # Inicializar un data.frame vacío
  df <- data.frame()
  
  # Utilizar foreach para iterar en paralelo
  results_list <- foreach(v = v_values, .combine = 'rbind', .export = c('Gen_FC_X1_X2', 'FC_X1_given_v'), 
                          .packages = c('coda','betafunctions')) %dopar% {
                            tmp_df <- data.frame()
                            for (prop_prec in prop_prec_values) {
                              results <- Gen_FC_X1_X2(N, prop_prec, a, b, c, d, v, option = "all", thin,burnin)
                              result_AR <- results$acc_rate
                              mcmc_obj <- mcmc(results$thinned_chain)
                              result_ESS <- effectiveSize(mcmc_obj)
                              tmp_df <- rbind(tmp_df, data.frame(v = v, Precision = prop_prec, Accept_Rate = result_AR, ESS = result_ESS))
                              #print(prop_prec)
                            }
                            tmp_df
                            #print(v_values)
                          }
  
  # Detener el clúster
  stopCluster(cl)
  print((N-burnin)/thin)
  
  df <- rbind(df, results_list)
  
  # ESS plot for v and precision values.
  measure_quantile = quantile(df$ESS,probs = c(0.25,0.5,0.75))
  
  plot_ESS=ggplot(df, aes(x=v, y=Precision, z=ESS)) +
    geom_tile(aes(fill=ESS)) +
    scale_fill_gradientn(colors=c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$ESS), measure_quantile[[3]])),
                         name= "Tamaño Efectivo",
                         breaks=measure_quantile, # This ensures that the minimum and maximum values are displayed in the legend
                         labels=sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks=seq(0, 0.25, by=0.05)) + # This increases the number of values displayed on the x-axis
    scale_y_continuous(breaks=seq(min(prop_prec_values), max(prop_prec_values), by=1)) + # This ensures that all integers are displayed on the y-axis
    labs(title=" ",
         x="v",
         y="Precisión",
         fill= "Tamaño Efectivo")
  
  # Acceptance rate plot for v and precision values.
  measure_quantile = quantile(df$Accept_Rate,probs = c(0,0.25,0.5,0.75,1))
  
  plot_AR=ggplot(df, aes(x=v, y=Precision, z=Accept_Rate)) +
    geom_tile(aes(fill=Accept_Rate)) +
    scale_fill_gradientn(colors=c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$Accept_Rate), measure_quantile[[5]])),
                         name= "Tasa de Aceptación",
                         breaks=measure_quantile, # This ensures that the minimum and maximum values are displayed in the legend
                         labels=sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks=seq(0, 0.25, by=0.05)) + # This increases the number of values displayed on the x-axis
    scale_y_continuous(breaks=seq(min(prop_prec_values), max(prop_prec_values), by=1)) + # This ensures that all integers are displayed on the y-axis
    labs(title=" ",
         x="v",
         y="Precisión",
         fill= "Tasa de Aceptación")
  
  return(list(plot_ESS,plot_AR))
}





Mon_R_Hat=function(N, prop_prec_values, a, b, c, d, v_values,thin = 1,burnin=1) {
  # Obtener el número total de núcleos
  num_cores <- detectCores()
  
  # Usar la mitad de los núcleos disponibles
  cl <- makeCluster(num_cores %/% 2)
  
  # Registrar el clúster para su uso con foreach
  registerDoParallel(cl)
  
  # Inicializar un data.frame vacío
  df <- data.frame()
  
  # Utilizar foreach para iterar en paralelo
  results_list <- foreach(v = v_values, .combine = 'rbind', .export = c('Gen_FC_X1_X2', 'FC_X1_given_v'), 
                          .packages = c('coda','betafunctions')) %dopar% {
                            tmp_df <- data.frame()
                            for (prop_prec in prop_prec_values) {
                              sample1=Gen_FC_X1_X2(N,prop_prec,a,b,c,d,v,option="all",thin = 1,X10_given = "random",burnin)
                              sample2=Gen_FC_X1_X2(N,prop_prec,a,b,c,d,v,option="all",thin = 1,X10_given = "random",burnin)
                              sample3=Gen_FC_X1_X2(N,prop_prec,a,b,c,d,v,option="all",thin = 1,X10_given = "random",burnin)
                              Gelm_Rud=gelman.diag(list(mcmc(sample1$thinned_chain),mcmc(sample2$thinned_chain),
                                                        mcmc(sample3$thinned_chain)))$psrf[1]
                              tmp_df <- rbind(tmp_df, data.frame(v = v, Precision = prop_prec, Gelman_Rudin = Gelm_Rud))
                              #print(prop_prec)
                            }
                            tmp_df
                            #print(v_values)
                          }
  
  # Detener el clúster
  stopCluster(cl)
  
  df <- rbind(df, results_list)
  
  # RHat plot for v and precision values.
  
  measure_quantile = quantile(df$Gelman_Rudin,probs = c(0.2,0.94,0.96,0.98,1))
  
  plot_H=ggplot(df, aes(x=v, y=Precision, z=Gelman_Rudin)) +
    geom_tile(aes(fill=Gelman_Rudin)) +
    scale_fill_gradientn(colors=c("red", "white", "blue"), 
                         values = scales::rescale(c(measure_quantile[[1]], mean(df$Gelman_Rudin), measure_quantile[[5]])),
                         name= "R-hat",
                         breaks=measure_quantile, # This ensures that the minimum and maximum values are displayed in the legend
                         labels=sprintf("%.2f", measure_quantile)) + # Formats the values to 2 decimal places
    scale_x_continuous(breaks=seq(0, 0.25, by=0.05)) + # This increases the number of values displayed on the x-axis
    scale_y_continuous(breaks=seq(min(prop_prec_values), max(prop_prec_values), by=1)) + # This ensures that all integers are displayed on the y-axis
    labs(title=" ",
         x="v",
         y="Precisión",
         fill= "R-hat")
  
  return(list(plot_H))
}









#####
# Function that plots the histogram, density, trace and convergence control using the average.
#####
# "nameaxisy" is the name of the vertical axis for trace and convergence monitoring.
# "witdh" is the width for the confidence intervals of convergence monitoring.
# "lscatt" is an increment to the minimum value generated. It allows plotting the line at an "lscatt" distance from the trace to enhance visualization.
# "uscatt" is an increment to the maximum value generated. It allows plotting the line at an "uscatt" distance from the trace to enhance visualization.
Graphs=function(dataset,nameaxisy,width=10,lscatt=0.05,uscatt=0.05){
  # Histogram with density
  l=length(dataset[,1])
  hist=ggplot(dataset, aes(x = dataset[,1])) + 
    geom_histogram(aes(y =after_stat(density)),colour = 1, fill = "white") +
    geom_density(lwd = 1.2,linetype = 2,colour = 2,fill=4,alpha=0.25)+
    labs(title = "Histograma y Densidad")+ylab("Densidad")+
    xlab(if(nameaxisy=="X2"){expression(X[2]^(t))}else if(nameaxisy=="X1"){expression(X[1]^(t))}
         else if(nameaxisy=="Y1"){expression(Y[1]^(t))}else if(nameaxisy=="Y2"){expression(Y[2]^(t))}else{
           substitute(va,list(va=as.name(nameaxisy)))})+
    theme(plot.title = element_text(size=11))
  
  # Trace plot with maximum and minimum
  trace=ggplot(dataset, aes(x=1:l,y=dataset[,1])) +
    geom_line()+xlab("t")+
    ylab(if(nameaxisy=="X2"){expression(X[2]^(t))}else if(nameaxisy=="X1"){expression(X[1]^(t))}
         else if(nameaxisy=="Y1"){expression(Y[1]^(t))}else if(nameaxisy=="Y2"){expression(Y[2]^(t))}else{
           substitute(va,list(va=as.name(nameaxisy)))})+
    ylim(c(min(dataset)-lscatt,max(dataset)+uscatt))+
    geom_hline(aes(yintercept = min(dataset[,1])),colour="red",linetype=2)+
    geom_text(aes(l-l/10,min(dataset[,1]),label = round(min(dataset[,1]),3), vjust = 2),colour="red")+
    geom_hline(aes(yintercept = max(dataset[,1])),colour="red",linetype=2)+
    geom_text(aes(l-l/10,max(dataset[,1]),label = round(max(dataset[,1]),3), vjust = -1),colour="red")+
    labs(title = substitute(list("Traza de la muestra aleatoria de tamaño",n),list(n=l)))+
    theme(plot.title = element_text(size=11))
  
  # Acf plot
  alfa=0.05
  lim= qnorm((1 - alfa / 2)) / sqrt(l)
  acf_values=acf(dataset,plot = FALSE)
  acf_data <- data.frame(Lag = acf_values$lag[-1],  # Eliminar el primer valor de lag (siempre 0)
                         ACF = acf_values$acf[-1])  # Eliminar el primer valor de ACF (siempre 1)
  acfplot=ggplot(acf_data, aes(x = Lag, y = ACF)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = c(lim, -lim), linetype = "dashed") +
    labs(title = "Función de Autocorrelación",
         x = "Lag",
         y = "ACF") +
    theme_minimal()+ theme(plot.title = element_text(size=11))  
  
  # Convergence control using averaging
  dataset$estintden=cumsum(dataset[,1])/(1:l)
  dataset$esterrden=sqrt(cumsum((dataset[,1]-dataset$estintden)^2))/(1:l)
  
  mean_X1_X2=ggplot(dataset,aes(x=1:l,y=estintden))+geom_line()+
    geom_line(aes(x=1:l,y=estintden-1.95*esterrden,colour="Upper"))+
    geom_line(aes(x=1:l,y=estintden+1.95*esterrden,colour="Lower"))+
    ylim(mean(dataset$estintden)+width*c(-dataset$esterrden[l],dataset$esterrden[l]))+
    ylab(if(nameaxisy=="X2"){expression(X[2]^(t))}else if(nameaxisy=="X1"){expression(X[1]^(t))}
         else if(nameaxisy=="Y1"){expression(Y[1]^(t))}else if(nameaxisy=="Y2"){expression(Y[2]^(t))}else{
           substitute(va,list(va=as.name(nameaxisy)))})+
    xlab("t")+geom_hline(yintercept = mean(dataset[,1]),colour="red",linetype=2)+
    geom_text(aes(l-l/10,mean(dataset[,1]),label = round(mean(dataset[,1]),3), vjust = -3),colour="red")+
    labs(title = "Control de convergencia usando promedio",color="Bounds")+ 
    scale_shape_discrete(name = " ")+ theme(plot.title = element_text(size=11),legend.position="top", legend.box="horizontal")
  
  # Plots of histogram, trace, convergence control and acf.
  allgraphs=grid.arrange(hist, trace, mean_X1_X2,acfplot, 
                         ncol=2, nrow=2, widths=c(4, 4), heights=c(2,2),layout_matrix=rbind(c(1,2),c(3,4)))
  return(allgraphs)
}

##########################################################
##########################################################
# Gibbs Sampling
##########################################################
##########################################################
# "N1": Gibbs Sampling sample size
# "N2": Random walks sample size for full conditional
# "a","b","c" and "d" are given values for parameters
# "thin" is the thinning interval for Random Walks. Every "thin" generated samples, one is stored to reduce autocorrelation.
# The seed type can be specified with "mu0_given", "random" or "fixed".
Gen_Joint_Dist=function(N1,N2,prop_prec,a,b,c,d,thin=1,X10_given="random",lower_epsilon=0,dig_tol=15){
  SampleGen=matrix(data = NA,nrow = N1,ncol = 2,dimnames = list(list(),list("X2","X1")))
  SampleGen=as.data.frame(SampleGen)
  SampleGen$X1[1]=rbeta(1,a,b)
  SampleGen$X2[1]=rBeta.4P(1,l = 0,u = SampleGen$X1[1]*(1-SampleGen$X1[1]),alpha = c,beta = d)
  for (t in 2:N1) {
    #SampleGen$X2[t]=rBeta.4P(1,l = 0,u = SampleGen$X1[t-1]*(1-SampleGen$X1[t-1]),alpha = c,beta = d)
    #result = Gen_FC_X1_X2(N2, prop_prec, a, b, c, d, SampleGen$X2[t], option = "end", thin, burnin = 0, X10_given)
    #SampleGen$X1[t] = result$thinned_chain
    
    SampleGen$X1[t] = Gen_FC_X1_X2(N2, prop_prec, a, b, c, d, SampleGen$X2[t-1], option = "end", thin, burnin = 0, X10_given,dig_tol = dig_tol)$thinned_chain
    SampleGen$X2[t]=rBeta.4P(1,l = lower_epsilon,u = SampleGen$X1[t]*(1-SampleGen$X1[t]),alpha = c,beta = d)
  }
  return(SampleGen)
}

#Example
#prueba0=Gen_Joint_Dist(N1=10,N2=5,prop_prec=3,a=3,b=2.5,c=4,d=6,thin=1,X10_given="random",lower_epsilon=0)

##########################################################
##########################################################
# Generalization of Tovar's method to obtain hyperparameter values
##########################################################
##########################################################
# "x1" and "x2" are values obtained from a person expert considered an expert on topic of interest
# "a" and "b" is the values where the parameter of interest remains
# "alpha" is the complement of the expert's belief that the interval (x1,x2) contains the true parameter value.
mtovar_vs2=function(x1,x2,low,upp,alp){
  tht0=(x1+x2)/2
  w=(tht0-low)/(upp-tht0)
  sig=sqrt(alp)*(x1-tht0)
  b= ((upp-low)^2*w-((w+1)^2*sig^2))/((w+1)^3*sig^2)
  a=w*b
  return(list(a=a,b=b,c=tht0))
}

##########################################################
##########################################################
# Joint moments of order m for proposed a priori distribution.
##########################################################
##########################################################
mompriordist=function(l1,l2,a,b,c,d){
  beta(c-l1-l2,l1+l2+d)*beta(l1+a,l2+b)
}




##########################################################
##########################################################
# Comparation of analytic and numeric results
##########################################################
##########################################################


Measure_Diagnostic=function(data1,data2,var="original",burnin,thin,digits=5,a,b,c,d){
  
  N=length(data1)
  data1 <- data1[seq((burnin+1), N, by=thin)]
  data2 <- data2[seq((burnin+1), N, by=thin)]
  
  if(var=="original"){
    new_names=c("Mean_X1","Var_X1","ESS_X1","Mean_X2","Var_X2","ESS_X2","Cov","Length")
    
    # Numerical results
    Numerical_results=round(data.frame( "mean1"=mean(data1),"var1"=var(data1),"ESS1"=effectiveSize(mcmc(data1))[[1]],
                                        "mean2"=mean(data2),"var2"=var(data2),"ESS2"=effectiveSize(mcmc(data2))[[1]],
                                        "cov12"=cov(data1,data2),"lenght"=length(data1)
    ),digits)
    
    names(Numerical_results)=new_names
    
    return(list(Numerical=Numerical_results))
    
  } else if(var=="transform") {
    piece=(data1*(1-data1)/data2-1)
    new_data1=data1*piece
    new_data2=(1-data1)*piece
    new_names=c("Mean_Y1","Var_Y1","ESS_Y1","Mean_Y2","Var_Y2","ESS_Y2","Cov","Length")
    
    # Numerical results
    Numerical_results=round(data.frame( "mean1"=mean(new_data1),"var1"=var(new_data1),"ESS1"=effectiveSize(mcmc(new_data1))[[1]],
                                        "mean2"=mean(new_data2),"var2"=var(new_data2),"ESS2"=effectiveSize(mcmc(new_data2))[[1]],
                                        "cov12"=cov(new_data1,new_data2),"lenght"=length(new_data1)
    ),digits)
    names(Numerical_results)=new_names
    
    # Analytically results
    K=mompriordist(0,0,a,b,c,d)
    Analytic_results=round(data.frame("Mean.1"=mompriordist(1,0,a,b,c,d)/K,
                                      "Var.1"=mompriordist(2,0,a,b,c,d)/K-(mompriordist(1,0,a,b,c,d)/K)^2,
                                      "ESS.1"=length(new_data1),
                                      "Mean.2"=mompriordist(0,1,a,b,c,d)/K,
                                      "Var.2"=mompriordist(0,2,a,b,c,d)/K-(mompriordist(0,1,a,b,c,d)/K)^2,
                                      "ESS.2"=length(new_data1),
                                      "Cov"=mompriordist(1,1,a,b,c,d)/K-(mompriordist(1,0,a,b,c,d)/K)*mompriordist(0,1,a,b,c,d)/K,
                                      "length"=length(new_data1)),digits)
    names(Analytic_results)=new_names
    
    # Differences between analytical and numerical results.
    Differences=round(Analytic_results-Numerical_results,digits)
    return(list(Numerical=Numerical_results,Analytical=Analytic_results,Differences=Differences))
  }
  
  #Heidel1=heidel.diag(mcmc(data1))
  #Heidel2=heidel.diag(mcmc(data2))
  #acf1=acf(data1);
  #acf2=acf(data2);
  #return(list(Diagnostic,Heidel1=Heidel1,Heidel2=Heidel2,acf1,acf2))
  #return(list(Numerical=Numerical_results,Analytical=Analytic_results,Differences=Differences))
}



# Analytically results
Measure_Analy=function(a,b,c,d,digits){
  
  K=mompriordist(0,0,a,b,c,d)
  Analytic_results=round(data.frame("Mean.1"=mompriordist(1,0,a,b,c,d)/K,
                                    "Var.1"=mompriordist(2,0,a,b,c,d)/K-(mompriordist(1,0,a,b,c,d)/K)^2,
                                    "Mean.2"=mompriordist(0,1,a,b,c,d)/K,
                                    "Var.2"=mompriordist(0,2,a,b,c,d)/K-(mompriordist(0,1,a,b,c,d)/K)^2,
                                    "Cov"=mompriordist(1,1,a,b,c,d)/K-(mompriordist(1,0,a,b,c,d)/K)*mompriordist(0,1,a,b,c,d)/K,
                                    "K"=K),
                         digits)
  return(Analytic_results)
}

#####################
## Obtención de hiperparámetros utilizando Intervalo Boostrap
#####################

# ssample: muestra origianl 
# r_boostrap: número de remuestras.
# q_boostrap: cuantiles boostrap
# option_mu: método para obtener hiperparámetros a y b de la media mu, puede ser "moments" o "tovar.
# sig_mu: Nivel de significancia para el método de obtención de hiperparámetros de tovar para la media mu.
# bound_var: método con el que se define el límite superior para la varianza, puede ser "min", "mean", "max".
# sig_var: Nivel de significancia para el método de obtención de hiperparámetros de tovar para la media la varianza.
# digits: cantidad de decimales para el intervalo de cuantiles boostrap de la media y la varianza.
# grpahs_boot: Indicador lógico para generar histograma, T o F.
# Q_E_mu y Q_E_cv: cuantiles de la media y la varianza obtenidos del Especialista.


hiperparameters=function(ssample,r_boostrap=100,q_boostrap=c(0.025,0.975),option_mu="moments",sig_mu=0.05,bound_var="max",sig_var=0.05,
                         digits=4,grpahs_boot=F,Q_E_mu=0,Q_E_cv=0){
  if(r_boostrap!=0){
    n_sample=length(ssample)
    boot=matrix(sample(ssample,size = r_boostrap*n_sample, replace = T),nrow = n_sample,ncol=r_boostrap)
    boots_mean=round(apply(boot, 2, mean),digits)
    boots_sd=round(apply(boot, 2, sd),digits)
    boots_cv=round(boots_sd/boots_mean,digits)
    
    #Se escogen los cuantiles asociados q_boostrap
    #Para la media
    quantile_mu=round(quantile(boots_mean,probs = q_boostrap),digits)
    #Para el CV
    quantile_cv=quantile(boots_cv,probs = q_boostrap)
  } else if(r_boostrap==0){
    quantile_mu =Q_E_mu
    quantile_cv =Q_E_cv
  }
  
  #Para la varianza
  quantile_var=round((quantile_cv*mean(quantile_mu))^2,digits)
  #############################
  #Intervalo para la media
  #############################
  if(option_mu=="moments"){
    portion=(mean(quantile_mu)*(1-mean(quantile_mu)) /((quantile_mu[[2]]-quantile_mu[[1]])/4)-1 )
    hiper_mean=data.frame("a"=mean(quantile_mu) * portion,"b"=(1-mean(quantile_mu)) * portion)
  } else if(option_mu=="tovar"){hiper_mean=mtovar_vs2(quantile_mu[[1]],quantile_mu[[2]],0,1,sig_mu)}
  
  #############################
  #Intervalo para la varianza
  #############################
  bound_var_value=if(bound_var=="min"){min(quantile_mu)*(1-min(quantile_mu))}
  else if(bound_var=="mean"){mean(quantile_mu)*(1-mean(quantile_mu))} 
  else if(bound_var=="max"){max(quantile_mu)*(1-max(quantile_mu))}
  
  hiper_var=mtovar_vs2(quantile_var[[1]],quantile_var[[2]],0,bound_var_value,sig_var)
  
  #############################
  # Histograma y densidad boostrap
  #############################
  
  if(grpahs_boot==T){
    hist_orig=ggplot(as.data.frame(ssample), aes(x = ssample)) + 
      geom_histogram(aes(y = after_stat(density)),colour = 1, fill = "white")+
      geom_density(lwd = 1.2,linetype = 2,colour = 2,fill=4,alpha=0.25)+
      labs(title = "Muestra Original")+ylab("Densidad")+
      xlab(substitute(va,list(va="X")))
    if(r_boostrap!=0){
      hist_boot_mean=ggplot(as.data.frame(boots_mean), aes(x = boots_mean)) + 
        geom_histogram(aes(y = ..density..),
                       colour = 1, fill = "white") +
        geom_density(lwd = 1.2,
                     linetype = 2,
                     colour = 2,fill=4,alpha=0.25)+
        labs(title = "Boostrapping para la media")+ylab("Densidad")+
        xlab(substitute(va,list(va="Media de X")))
      
      hist_boot_cv=ggplot(as.data.frame(boots_cv), aes(x = boots_cv)) + 
        geom_histogram(aes(y = ..density..),
                       colour = 1, fill = "white") +
        geom_density(lwd = 1.2,
                     linetype = 2,
                     colour = 2,fill=4,alpha=0.25)+
        labs(title = "Boostrapping para el CV")+ylab("Densidad")+
        xlab(substitute(va,list(va="CV de X")))
      
      grid.arrange(hist_orig, hist_boot_cv,hist_boot_mean, 
                   ncol=2, nrow=2, widths=c(2,2), heights=c(2,2),layout_matrix=rbind(c(1,1),c(2,3)))
    }else if(r_boostrap==0){
      hist_orig
    }
  }
  
  
  return(list(Q_mean=quantile_mu,Q_cv=quantile_cv,Q_var=quantile_var,hiper_mean=hiper_mean,hiper_var=hiper_var,
              Var_Upper_Bound=bound_var_value))
}





#####################
## Estimación posterior utilizando muestreo por importancia.
#####################
Est_Post=function(ssample,N,N_FC,Precision,a,b,c,d,thin1,thin2,burnin,dig_tol){
  Example_Joint_Dist=Gen_Joint_Dist(N1 = N,N2 = N_FC,prop_prec = Precision,a,b,c,d,thin=thin1,X10_given = "random",dig_tol = dig_tol)
  piece=(Example_Joint_Dist$X1[seq((burnin+1), N, by=thin2)]*(1-Example_Joint_Dist$X1[seq((burnin+1), N, by=thin2)])/
           Example_Joint_Dist$X2[seq((burnin+1), N, by=thin2)]-1)
  sample_alpha=Example_Joint_Dist$X1[seq((burnin+1), N, by=thin2)]*piece
  sample_beta=(1-Example_Joint_Dist$X1[seq((burnin+1), N, by=thin2)])*piece
  x0=prod(ssample)
  y0=prod(1-ssample)
  #marginalike=sum(x0^(sample_alpha-1)*y0^(sample_beta-1)*prior(sample_alpha,sample_beta,a,b,c,d)/
  #                  beta(sample_alpha,sample_beta)^(length(ssample)))
  
  marginalike=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+log(prior(sample_alpha,sample_beta,a,b,c,d))-
                        (length(ssample))*log(beta(sample_alpha,sample_beta))))
  
  meanpalph=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+log(sample_alpha)+
                      log(prior(sample_alpha,sample_beta,a,b,c,d))-
                      (length(ssample))*log(beta(sample_alpha,sample_beta))-log(marginalike)))
  
  meanpbet=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+log(sample_beta)+
                     log(prior(sample_alpha,sample_beta,a,b,c,d))-
                     (length(ssample))*log(beta(sample_alpha,sample_beta))-log(marginalike)))
  
  Descriptivo=Measure_Diagnostic(data1 = Example_Joint_Dist$X1,data2 = (Example_Joint_Dist$X2),var ="transform",
                                 digits = 4,burnin,thin = thin2,a,b,c,d)
  return(list(EstPost=data.frame(Prior_Lik=marginalike,P_A=meanpalph,P_B=meanpbet),Descriptivo=Descriptivo,
              SA=sample_alpha,SB=sample_beta,SX1=Example_Joint_Dist$X1,SX2=Example_Joint_Dist$X2))
}















####
# Función que automatiza lo anterior
# Estudio de simulación para comparar y monitorear el comportamiento de las estimaciones posteriores
# producidas por los diferentes conjuntos de hiperparámetros y muestras observadas de distintos tamaño.
####

Estudio_Sim=function(N1,N2,prop_prec,a,b,c,d,thin1,X10_given = "random",dig_tol = 15,thin2,burnin,n_sample,alpha_real,beta_real,N_Iter_Sim){
  
  Sample_Prior_Hip=Gen_Joint_Dist(N1,N2,prop_prec = 3, a, b, c, d,thin=thin1,X10_given = "random",dig_tol=dig_tol)
  piece=(Sample_Prior_Hip$X1[seq((burnin+1),N1, by=thin2)]*(1-Sample_Prior_Hip$X1[seq((burnin+1), N1, by=thin2)])/
           Sample_Prior_Hip$X2[seq((burnin+1),N1, by=thin2)]-1)
  sample_alpha=Sample_Prior_Hip$X1[seq((burnin+1), N1, by=thin2)]*piece
  sample_beta=(1-Sample_Prior_Hip$X1[seq((burnin+1), N1, by=thin2)])*piece
  
  Descriptive_Sample_Prior=Measure_Diagnostic(data1 = Sample_Prior_Hip$X1,data2 = Sample_Prior_Hip$X2,var ="transform",
                                              digits = 4,burnin,thin = thin2,a,b,c,d)
  
  Iter_Alpha=NA;Iter_Beta=NA
  VIter_Alpha=matrix(data = NA,nrow = 0,ncol = 10,dimnames = list(list(),list("Min","Q2.5","Q50","Q97.5","Max","Mean","Var","Bias","MSE","SampleSize")))
  VIter_Beta=VIter_Alpha
  
  for (j in 1:length(n_sample) ) {
    for (i in 1:N_Iter_Sim) {
      sampe_prueba=rbeta(n_sample[j],alpha_real,beta_real)
      x0=prod(sampe_prueba)
      y0=prod(1-sampe_prueba)
      marginalike=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+
                            log(prior(sample_alpha,sample_beta,a,b,c,d))-
                            (length(sampe_prueba))*log(beta(sample_alpha,sample_beta))))
      
      Iter_Alpha[i]=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+log(sample_alpha)+
                              log(prior(sample_alpha,sample_beta,a,b,c,d))-
                              (length(sampe_prueba))*log(beta(sample_alpha,sample_beta))-log(marginalike)))
      
      Iter_Beta[i]=sum(exp((sample_alpha-1)*log(x0)+(sample_beta-1)*log(y0)+log(sample_beta)+
                             log(prior(sample_alpha,sample_beta,a,b,c,d))-
                             (length(sampe_prueba))*log(beta(sample_alpha,sample_beta))-log(marginalike)))
    }
    VIter_Alpha=rbind(VIter_Alpha,c(quantile(Iter_Alpha,probs = c(0,0.025,0.5,0.975,1),na.rm = T),mean(Iter_Alpha),var(Iter_Alpha),( alpha_real-mean(Iter_Alpha) ),
                                    ( var(Iter_Alpha)+(alpha_real-mean(Iter_Alpha))^2),n_sample[j]))
    VIter_Beta=rbind(VIter_Beta,c(quantile(Iter_Beta,probs = c(0,0.025,0.5,0.975,1),na.rm = T),mean(Iter_Beta),var(Iter_Beta),( beta_real-mean(Iter_Beta) ),
                                  ( var(Iter_Beta)+(beta_real-mean(Iter_Beta))^2),n_sample[j]))
  }
  return(list(Result_Alpha=as.data.frame(VIter_Alpha),Result_Beta=as.data.frame(VIter_Beta),Descriptive_Sample_Prior=Descriptive_Sample_Prior))
}


#Gráficas individuales para el monitoreo de la convergencia de las estimaciones.
Individual_Graphs=function(Data,lim_x,value_real,y_text,title_text){
  # Gráfica de la media con bandas construidas con los cuantiles de las iteraciones de cada tamaño de muestra.
  Mean_Graph=ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = Q2.5, color = "Q2.5")) +
    geom_line(aes(y = Mean, color = "Mean")) +
    geom_line(aes(y = Q97.5, color = "Q97.5")) +
    labs(title = title_text,
         x = "Tamaño de Muestra",
         y = y_text) + 
    scale_color_manual(name = " ",
                       values = c("Q2.5" = "blue", "Mean" = "green", "Q97.5" = "red")) +# ylim(c(7,18))+
    scale_x_continuous(breaks = lim_x) +  # Modificación para mostrar valores de 5 en 5
    geom_hline(yintercept=value_real, linetype="dashed", color = "red")+
    theme_minimal()
  
  # Gráfica del MSE
  Mse_Graph=ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = MSE, color = "MSE")) +
    labs(title = " ",
         x = "Tamaño de Muestra",
         y = "MSE") +
    scale_color_manual(name = " ",
                       values = c("MSE" = "blue")) +# ylim(c(7,18))+
    scale_x_continuous(breaks = lim_x) +  # Modificación para mostrar valores de 5 en 5
    theme_minimal()
  
  # Gráfica del Bias
  Bias_Grap=ggplot(Data, aes(x = SampleSize)) +
    geom_line(aes(y = Bias, color = "Bias")) +
    labs(title = " ",
         x = "Tamaño de Muestra",
         y = "Bias") +
    scale_color_manual(name = " ",
                       values = c("Bias" = "black")) +# ylim(c(7,18))+
    scale_x_continuous(breaks = lim_x) +  # Modificación para mostrar valores de 5 en 5
    theme_minimal()
  
  grid.arrange(Mean_Graph,Mse_Graph,Bias_Grap,nrow = 2,ncol=2,layout_matrix = rbind(c(1, 1),c(2, 3)))
}


####
# Gráfica de funciones conjuntas
####
Comparacion_hip=function(data,value_real,lim_x,title_text,y_Text){
  # Mean
  Comparacion_Mean=ggplot(data, aes(x = SampleSize, y = Mean, color = Source, shape = Source, linetype = Source)) +
    geom_line() +
    geom_point(size = 2) +
    labs(title = title_text,
         x = " ",
         y = y_Text) +
    scale_x_continuous(breaks = lim_x) +
    geom_hline(yintercept=value_real, linetype="dashed", color = "red")+
    theme_minimal() +
    theme(legend.position = "none") +  # Quitar leyenda
    #theme(legend.title = element_blank())+
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Diferentes formas para los puntos
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  
  
  # Bias
  Comparacion_Bias=ggplot(data, aes(x = SampleSize, y = Bias, color = Source, shape = Source, linetype = Source)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = " ",
         y = "Bias") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() +
    theme(legend.position = "none") +  # Quitar leyenda
    #theme(legend.title = element_blank())+
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Diferentes formas para los puntos
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # MSE
  Comparacion_Mse=ggplot(data, aes(x = SampleSize, y = MSE, color = Source, shape = Source, linetype = Source)) +
    geom_line() +
    geom_point() +
    labs(title = " ",
         x = "Tamaño de Muestra",
         y = "Mse") +
    scale_x_continuous(breaks = lim_x) +
    theme_minimal() + 
    theme(legend.position = "none") +  # Quitar leyenda
    #theme(legend.title = element_blank())+
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 3, 4, 5, 6)) +  # Diferentes formas para los puntos
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "twodash", "dotted", "solid", "dashed"))
  
  # Función para extraer la leyenda de un gráfico
  g_legend <- function(a.gplot){
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Extraer la leyenda del gráfico de medias
  legend <- g_legend(Comparacion_Mean + theme(legend.position = "right",legend.title = element_blank()))
  
  # Combinar los gráficos sin leyenda con la leyenda a la derecha
  grid.arrange(
    arrangeGrob(Comparacion_Mean, Comparacion_Bias, Comparacion_Mse, ncol = 1,heights = c(6, 6,6)),
    legend, 
    ncol = 2, 
    widths = c(7, 1)
  )
  
  
  #grid.arrange(Comparacion_Mean,Comparacion_Bias,Comparacion_Mse,nrow = 2,ncol=2,layout_matrix = rbind(c(1, 1),c(2, 3)))
}


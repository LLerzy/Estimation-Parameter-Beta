# Posterior Estimation of the Shape Parameters of the Beta Distribution

This repository is a branch of a main repository called [**Estimation-Parameter-Beta**](https://github.com/LLerzy/Estimation-Parameter-Beta).

## Table of Contents

-   [Overview](#overview)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Included Functions](#included-functions)
-   [Branches](#branches)
-   [Contributions](#contributions)

## Overview

This branch contains a set of routines designed to estimate the shape parameters of the Beta distribution for the random variable \(X\) from a Bayesian perspective. The methodology is based on a new bivariate prior distribution for the mean and variance, given by:  

\[
f_{\mu,\sigma^2}(u,v|\phi=(a,b,c,d))=\dfrac{1}{\text{beta}(a,b) \text{beta}(c,d)} \dfrac{v^{c-1} (u (1-u)-v)^{d-1}}{u^{c+d-a} (1-u)^{c+d-b} },
\]

along with a corresponding method for generating random samples, which can be found in the [**Algorithm-Sim-Samples**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples) branch.  

The main objective of this branch is to conduct a simulation study to analyze the behavior of the estimations obtained for these parameters.  

The scripts include both empirical and subjective Bayesian approaches for hyperparameter estimation, utilizing bootstrap intervals as well as expert-defined intervals with different biases and semi-amplitudes.  

All code is written in **R-Project**, leveraging the statistical and computational tools available in this software. The project provides both numerical and theoretical results.  


## Installation

To run the code, you need to have `R` installed with the following packages:

- `ggplot2` 
- `gridExtra` 
- `tidyr` 
- `plotly` 
- `coda` 
- `foreach` 
- `doParallel` 
- `betafunctions` 
- `openxlsx` 
- `xtable`

You can install the required packages using the following command:

``` r
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Usage

To run the designed algorithms, follow these steps:

1. Clone this repository:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/main.git
    ```

2. Open the R script in any compatible Integrated Development Environment (IDE):

    -   For the script that contains the main functions in the shape parameter estimation process, refer to [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/requiredfunctions.R).

3. You can run the script and load all the functions with the following command:

    ``` bash
    Rscript requiredfunctions.R
    ```

## Included Functions

The functions defined in the `requiredfunctions.R` script include:

-   `Prior`: Defines the proposed prior probability density function, referred to as combination 1 in the [New-Biv-Dist](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist) branch.
-   `FC_X1_Given_v`: Full conditional distribution of $X_1$ given $X_2 = v$.
-   `Graph_Fc_X1`: Plots the full conditional distribution for given parameter values.
-   `Gen_FC_X1_X2`: Metropolis-Hastings algorithm with random walk to generate samples from the full conditional distribution.
-   `Mon_Measure`: Monitors acceptance rates and ESS (effective sample size) for different values of $v$ and precision.
-   `Mon_R_Hat`: Monitors the Gelman-Rubin diagnostic ($R$-hat) for different values of $v$ and precision.
-   `Graphs`: Generates histograms, density plots, trace plots, and convergence diagnostics using averages.
-   `Gen_Joint_Dist`: Gibbs sampling to generate samples from the joint distribution of the random vector $(X_1,X_2)$.
-   `Mtovar_vs2`: Generalizes Tovar’s method to obtain hyperparameter values.
-   `Mom_Prior_Dist`: Computes the joint moments of order $l$ for the proposed prior distribution.
-   `Measure_Diagnostic`: Compares analytical and numerical results for user-provided data samples.
-   `Measure_Analy`: Computes the analytical results for the proposed prior distribution.
-   `Hyperparameters`: Obtains hyperparameter values using empirical and subjective Bayesian approaches.
-   `Est_Post`: Posterior estimation of the Beta distribution's shape parameters $\alpha$ and $\beta$ using importance sampling.
-   `Sim_study`: Conducts simulation studies to compare posterior estimates using different hyperparameters and sample sizes.
-   `Individual_Graphs`: Creates individual graphs to monitor posterior estimates using bias and MSE as indicators.
-   `Comparison_Hyper`: Compares joint density functions for different sets of hyperparameters.

### Notes:

-   It is recommended to review and adapt each function according to the specific needs of each analysis.
-   Make sure you understand each function before using it to ensure accurate results and avoid potential errors.

## Branches

The following branches are available in this repository:

-   [**New-Biv-Dist**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist): Presents numerically approximated characteristics for seven of the 25 proposed probability distributions representing the random behavior of the Beta distribution’s shape parameters. Additionally, some numerical and theoretical approximations for combination 1 are compared.
-   [**Algorithm-Sim-Samples**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples): Monitors the convergence of chains generated by an algorithm using MCMC methods to generate random samples from a new proposed probability distribution, referred to as combination 1 in the [New-Biv-Dist](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist) branch.
-   [**Post-Estimate**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Post-Estimate): Conducts a simulation study to estimate the shape parameters of the Beta distribution. This study involves: two questions that can be used in an elicitation process, hyperparameter determination methods, three scenarios for the shape parameters, hypothetical expert information with different levels of bias, and graphs comparing posterior estimators through mean estimation, bias, mean squared error, coverage probability, and average length.

## Contributions

Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.

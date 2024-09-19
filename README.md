# MCMC-Based Algorithm for Simulating Samples from a New Bivariate Distribution

This repository is a branch of a main repository called **[Estimation-Parameter-Beta](https://github.com/LLerzy/Estimation-Parameter-Beta)**.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Files and Documents](#files-and-documents)
- [Contributions](#contributions)

## Overview

This project presents a set of functions implementing MCMC algorithms designed to generate random samples $(x_1,x_2)_{(n)}$ from the random vector \((X_1, X_2)\) with the density function 

$$f(x_1,x_2|\phi)=\dfrac{1}{\text{beta}(a,b)\text{beta}(c,d)}\dfrac{x_2^{c-1}[x_1(1-x_1)-x_2]^{d-1}}{x_1^{c+d-a}(1-x_1)^{c+d-b}}.$$

Subsequently, the generated sample is transformed into the space of the random vector $(Y_1,Y_2)$ with the density function,

$$f_{\phi}(y_1,y_2)= \dfrac{1}{\text{beta}(a, b) \text{beta}(c, d)}\ y_1^{a-1}y_2^{b-1}(y_1+y_2)^{d-(a+b)} (y_1+y_2+1)^{-c-d},\hspace{1cm}y_1,y_2\in\mathbb{R}_+.$$ 

All the code is written in R-Project, leveraging the statistical and computational tools of this software. The project includes both numerical and theoretical results, focusing on the analysis of the generated chains through eight different convergence criteria. These diagnostics provide a comprehensive evaluation of the chain performance, ensuring reliability and accuracy of the estimates. The presented methods are suitable for a wide range of applications involving Bayesian inference and MCMC methods.

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
```R
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Usage
To replicate the analysis and run the algorithms, follow these steps:

1. Clone this repository:
   ```bash
   git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples.git
   ```
2. Open the R script or RMarkdown file in RStudio or any compatible Integrated Development Environment (IDE):
   - For the script containing the main functions, refer to [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/requiredfunctions.R).
   - Alternatively, you can view the rendered Markdown output in [`randomsamples.md`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/randomsamples.md).

3. Run the script to generate the results:
   ```bash
   Rscript -e "rmarkdown::render('randomsamples.Rmd')"
   ```

The code will generate convergence diagnostics such as histograms versus density plots, heatmaps for effective sample size, acceptance rate and R-Hat, cumulative mean, trace plots, autocorrelation, and contour plots along with scatter plots of the generated sample.

## Files and Documents
- [`randomsamples.Rmd`](randomsamples.Rmd): The main RMarkdown file containing the code to generate the results.
- [`randomsamples.md`](randomsamples.md): A Markdown version of the results.
- [`requiredfunctions.R`](requiredfunctions.R): Contains the custom functions used in the algorithm, such as Gibbs Sampling, Metropolis Hastings, and the new bivariate distribution.
- [`randomsamples_files/`](randomsamples_files/figure-gfm): Folder containing the figures generated by the RMarkdown file.

### Key Figures:
- **Effective Sample Size and Acceptance Rate vs. Sample Size and Precision**: See the generated figure [`EffectiveSize-AcceptanceRate.png`](randomsamples_files/figure-gfm/unnamed-chunk-3-1.png) to analyze the algorithm's performance in different configurations of precision and sample size.
- **Contour Curves vs. Scatter Plot of the Random Vector $(Y_1,Y_2)$**: See the graph in [`Curvas-de-Nivel-vs-Diagrama-Puntos.png`](randomsamples_files/figure-gfm/unnamed-chunk-14-1.png) for an evaluation of the chain convergence.

## Contributions
Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.
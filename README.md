# New Bivariate Probability Distributions

This repository is a branch of a main repository called [**Estimation-Parameter-Beta**](https://github.com/LLerzy/Estimation-Parameter-Beta).

## Table of Contents

-   [Overview](#overview)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Files and Documents](#files-and-documents)
-   [Contributions](#contributions)

## Overview

This project presents several bivariate distributions constructed for the random vector $(Y_1,Y_2)$ using a continuous and differentiable transformation,

$$\begin{matrix}
   T^{-1}: &  \mathbb{R}^2_+ & \Longrightarrow &(0,1)\times (0,U(X_1))\\
        &\left(Y_1,Y_2\right)&\longrightarrow & (X_1,X_2):=\left(\dfrac{Y_1}{Y_1+Y_2},\dfrac{Y_1Y_2}{(Y_1+Y_2)^2(Y_1+Y_2+1)}\right)
\end{matrix},$$

along with the joint density function representation of the vector $(X_1,X_2)$,

$$f_{X_1,X_2}(x_1,x_2|\phi=(\phi_1,\phi_2)) = f_{X_1}(x_1|\phi_1)f_{X_2|X_1}(x_2|x_1,\phi_2).$$

Although five distributions for $X_1$ and five for $X_2|X_1$ were considered, for a total of 25 combinations, this branch presents the characteristics of only seven distributions, all with a parameter vector $\phi=(a,b,c,d)$. For each one, the density surface, contour curves, marginal density, marginal cumulative probability distribution, and moments (mean, variance) are constructed. All the code is written in Mathematica, taking advantage of the symbolic manipulation and computational power of this software. The project includes numerical results obtained using the NIntegrate function, as well as some theoretical results.

## Installation

To run the code, you need to have `Mathematica` installed.

## Usage

To replicate the analysis and execute the algorithms, follow these steps:

1.  Clone this repository:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist.git
    ```

2.  Open the Mathematica script in any compatible Integrated Development Environment (IDE):

    -   For the script that contains the characteristics (mentioned in the overview) of the seven selected bivariate distributions for the vector $(Y_1,Y_2)$, refer to [`SelectedBivariateDistributions.nb`](SelectedBivariateDistributions.nb).
    -   For the script that contains the characteristics of the bivariate distribution of the vector $(Y_1,Y_2)$, constructed by considering $X_1 \sim Beta(a,b)$ and $X_2|X_1 \sim Beta(0,U(X_1),c,d)$, refer to [`DistribucionBetaBeta4P.nb`](DistribucionBetaBeta4P.nb).

3.  The script must be run within the IDE.

The code in `DistribucionBetaBeta4P.nb` considers four configurations for the parameter vector $(a,b,c,d)$ of the bivariate distribution constructed using $X_1 \sim Beta(a,b)$ and $X_2|X_1 \sim Beta(0,U(X_1),c,d)$. It presents the density surface, contour curves, marginal cumulative distribution, survival surface, marginal density function, and compares the theoretical moments with the numerical approximations obtained using the NIntegrate function of the `Mathematica` software.

## Files and Documents

-   [`DistBivSelecc.pdf`](DistBivSelecc.pdf): Contains numerical results obtained for each of the seven selected bivariate distributions.
-   [`DistribucionBetaBeta4P.pdf`](DistribucionBetaBeta4P.pdf): Contains numerical and theoretical results obtained for the bivariate distribution of $(Y_1,Y_2)$, constructed by considering $X_1 \sim Beta(a,b)$ and $X_2|X_1 \sim Beta(0,U(X_1),c,d)$.

## Contributions

Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.
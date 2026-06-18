# BASIN: Bayesian Sediment Source Apportionment Tool

**BASIN** is an open-source Bayesian modelling framework designed for robust sediment source apportionment. It addresses inherent methodological limitations of classical sediment fingerprinting, such as tracer non-conservatism and the "equifinality trap", by shifting from the analysis of isolated tracers to the evaluation of structural geochemical covariance.

The tool provides a user-friendly graphical interface (GUI) built with R and Shiny, while relying on the Hamiltonian Monte Carlo (HMC) algorithm via Stan for fast, robust posterior sampling.

## Key Features

* **Intelligent Tracer Pre-filtering:** Includes an SSA-aware 1D range test, the Geochemical Association Preservation Test (GAPT), and a PCA Convex Hull Penalty to mitigate equifinality by preserving structural integrity.
* **Flexible Error Architecture:** Evaluates elemental covariance directly within the Markov chains via *Mixture Covariance* and *Full Bayes* modes, transforming tracer multicollinearity into a valuable discriminatory signal.
* **Dynamic Particle Size Correction:** Evaluates non-linear exponential corrections (Beta-correction) iteratively during MCMC sampling.
* **Robust Bayesian Bias Absorption:** Utilizes a heavy-tailed Student's t-distribution and a random displacement vector to buffer localized hydrodynamic fluctuations, un-sampled "ghost hotspots," and secondary geochemical shifts.
* **Compositional Data Analysis (CoDA):** Optional isometric log-ratio (ILR) transformation module to eliminate spurious correlations inherent to closed geochemical data.
* **Built-in Model Comparison:** Natively integrates Pareto Smoothed Importance Sampling Leave-One-Out cross-validation (PSIS-LOO) and Exact LOO-CV to evaluate out-of-sample predictive accuracy and prevent overfitting.
* **Virtual Mixtures Generation:** Integrated tools to generate stochastic, semi-stochastic, and deterministic virtual mixtures with controlled non-conservatism for rigorous model validation.

## Prerequisites and Installation

BASIN is run locally in your web browser. To deploy the model, you only need an installed R environment. The Stan code of the MCMC sampler is dynamically compiled "on the fly" with each run, ensuring strict reproducibility across platforms.

### 1. Install R
* Download and install [R](https://cran.r-project.org/).
* (Optional) Download and install [RStudio](https://posit.co/download/rstudio-desktop/).

### 2. Install Rtools (Windows Only) / Xcode (Mac)
Because BASIN compiles C++ code dynamically via Stan, you must have a C++ toolchain installed:
* **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (make sure the version matches your R version).
* **Mac:** Install Xcode command line tools by running `xcode-select --install` in your terminal.
* **Linux:** Install `build-essential` via your package manager.

### 3. Dependencies Installation
All required R packages are configured to be **installed automatically** upon the first launch of the application. 

**⚠️ Fallback (Manual Installation):**
If the automatic installation fails (e.g., due to firewall restrictions or lack of administrative privileges), you can install the dependencies manually. Open your command line / terminal, launch R, select a CRAN mirror when prompted, and execute the following commands:

```R
install.packages("pacman")
pacman::p_load(shiny, rstan, tidyverse, caret, reshape2, dplyr, ggplot2, DT, shinythemes, data.table, MASS, zip, compositions, scoringRules)

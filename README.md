# BASIN: Bayesian Sediment Source Apportionment Tool

**BASIN** is an open-source, context-dependent modular Bayesian modelling framework designed for robust sediment source fingerprinting (SSF). Built on the **Stan** probabilistic programming language (utilizing Hamiltonian Monte Carlo and the No-U-Turn Sampler, HMC/NUTS) and an accessible **R/Shiny** graphical user interface, BASIN addresses the classic "equifinality trap", observation-error misweighting, and transport-induced tracer non-conservatism.

Recognizing that no single statistical parameterization is universally optimal across contrasting sampling designs, simplex geometries, and geomorphological regimes, BASIN replaces rigid "black-box" pipelines with **four independently configurable design axes**:
1. **Cross-Sample Information Sharing (Pooling Modes)**
2. **Residual Error-Scale & Inter-Tracer Covariance Architecture**
3. **Source End-Member Representation**
4. **In-Sampler Non-Conservatism & Compositional Buffers**

---

## Key Features

### 1. Geometry-Preserving Tracer Pre-Filtering
Instead of rigid univariate exclusion that dismantles the multidimensional mixing space, BASIN provides an optional, structure-preserving pre-filtering pipeline:
- **Proxy-Aware 1D Range Test:** Evaluates source–mixture range compatibility with optional normalization for hydrodynamic sorting or enrichment proxies (e.g., $\text{SSA}$, $d_{50}$, $\text{TOC}$).
- **Inter-Tracer Association Preservation Test (APT / GAPT):** Retains co-migrating tracers that preserve their structural pairwise relationships between sources and target mixtures.
- **Multivariate PCA Convex-Hull Leverage Filter:** An $n$-dimensional drop-out diagnostic that identifies and removes individual tracers pulling target mixtures outside the source mixing polygon.
- **Kruskal–Wallis & Collinearity Pruning:** Confirms inter-source discriminatory power ($p < 0.05$) while purging near-duplicate signals ($|\rho| > 0.95$).

### 2. Four-Axis Bayesian Unmixing Engine
- **Selectable Residual Error-Scale Specifications:**
  - *Free Regularized Residual Scales (Default):* Assigns each tracer a weakly regularized scale parameter ($\sigma_{\text{mix},t} \sim \text{Half-Student-}t(3, 0, \sigma_{\text{prior}})$), acting as an adaptive Bayesian M-estimator. Decoupling variance normalization from squared proportions ($P_{i,g}^2$) eliminates the mathematical penalty against pure end-members ($P_g = 0$ or $1$) and prevents centroid shrinkage on boundary-rich mixtures, even for isolated single targets ($N_{\text{mix}} = 1$).
  - *Process-Error Only (Source-Variance Propagated):* Deterministically propagates intra-source variances ($\sum_g P_{i,g}^2 \sigma_{\text{src},g,t}^2$) into the mixture likelihood—optimal for interior-dominated, low-source-count ($G = 2$) mixtures with strongly heterogeneous source variances.
- **Heavy-Tailed Multivariate Likelihood & Covariance Modes:** Models target observations using a multivariate Student's $t$-distribution ($\nu = 4$) to downweight erratic biogeochemical outliers, paired with three Cholesky-factorized error structures:
  - *Independent (Diagonal)*
  - *Target Mixture Covariance ($\Sigma_{\text{mix}}$, LKJ-regularized)*
  - *Full Bayes (Simultaneous Source + Mixture Covariance, $\Sigma_{\text{src}} + \Sigma_{\text{mix}}$)*
- **Controlled Gradient of Information Sharing:** Seamlessly switch between *Strictly Isolated Inference* ($N_{\text{mix}} = 1$), *Non-Hierarchical Joint Noise-Learning* (independent sample proportions with shared error parameters), *Hierarchical Partial Pooling* (multi-level Dirichlet hierarchy for stationary ensembles), and *Complete Pooling* (technical replicates), with direct support for prior geographical/connectivity weights ($w_{\text{conn}}$).
- **Mechanistic Non-Conservatism & Compositional Buffers:**
  - *Sorting Correction (Beta-module):* Evaluates tracer-specific non-linear exponential sorting adjustments ($\beta_{\text{sorting},t}$) directly within the MCMC sampler.
  - *Systematic Bias Absorption ($\delta_{\text{bias}}$):* A heavy-tailed displacement vector that isolates coherent ensemble-wide offsets (e.g., heavy-mineral lag or minor unsampled erosion hotspots contributing $\le 5\text{--}10\%$) without skewing known source proportions.
  - *Compositional Data Analysis (CoDA):* Isometric log-ratio (ILR) transformation for strictly conservative closed compositional datasets.

### 3. Integrated Diagnostics, Cross-Validation & Virtual Mixtures
- **Multi-Tracer Mass-Balance Reconstruction:** Evaluates internal mathematical closure (`sMAPE`, `Mean Bias`, and `NRMSE_src` per tracer) to diagnose discrepancies between internal fit and external accuracy.
- **Bayesian Model Comparison:** Natively integrates approximate **PSIS-LOO** (with Pareto-$k$ outlier diagnostics) and **Exact LOO-CV** to rank competing model configurations by Expected Log Predictive Density (`elpd_loo`) and prevent overfitting.
- **Stochastic & Deterministic Virtual Mixtures:** Generates virtual mixtures across 4 regimes (Full Stochastic, Semi-Stochastic Simplex Grid, Semi-Deterministic, Full Deterministic) and reports point accuracy (`MAE`, `RMSE`, `R2`, `NSE`), proper scoring rules (`CRPS`), and boundary-adjusted empirical interval coverage (`W50`, `W95`, `P50/P95 adj`).

---

## Practitioner Deployment Guidelines

| Sampling Design & Process Regime | Recommended Pooling Mode | Error-Scale & Covariance Specification | Non-Conservatism Buffers |
| :--- | :--- | :--- | :--- |
| **1. Single conservative / weakly non-conservative target ($N_{\text{mix}} = 1$)** | No Pooling (Isolated) | Free Regularized Scales + **Full Bayes** Covariance | None (keep parsimonious) |
| **2. Single target with suspected transport shifts ($N_{\text{mix}} = 1$)** | No Pooling (Isolated) | Free Regularized Scales + Independent or MixCov | **Beta-Correction** &/or **Bias Absorption** |
| **3. Stationary ensemble (replicate traps, background core layers)** | **Hierarchical Partial Pooling** | Free Scales (or Process-Error if interior-rich) + MixCov | Optional; verify via PSIS-LOO |
| **4. Multi-event ensemble or core with extreme flood/typhoon pulses** | **No Pooling** (preserves pulse amplitude) | Free Regularized Scales + MixCov or Full Bayes | **Bias Absorption** + Beta (if proxy available) |
| **5. Interior-rich, $G=2$ design with heterogeneous source variances** | No Pooling or Hierarchical | **Process-Error Only** (calibrates 50% CI coverage) | None unless sorting is severe |

---

## Prerequisites and Installation

BASIN runs locally in your standard web browser. With each execution, the Stan C++ code corresponding to the selected module combination is dynamically generated and compiled "on the fly", ensuring cross-platform reproducibility and full mathematical transparency.

### 1. Install R
- Download and install [R](https://cran.r-project.org/) (version $\ge 4.1$ recommended).
- *(Optional)* Download and install [RStudio](https://posit.co/download/rstudio-desktop/).

### 2. Install C++ Toolchain (Rtools / Xcode / build-essential)
Because BASIN compiles Stan models to C++, a system C++ compiler is required:
- **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.
- **macOS:** Install Xcode Command Line Tools by running `xcode-select --install` in the Terminal.
- **Linux:** Install `build-essential` (e.g., `sudo apt-get install build-essential` on Ubuntu/Debian).

### 3. Install R Package Dependencies
All required R packages are configured to be **installed automatically** via `pacman` upon the first launch of the application.

**⚠️ Manual Installation (Fallback):**
If automatic installation is restricted by your network or permissions, open R and run:

```R
install.packages("pacman")
pacman::p_load(
  shiny, rstan, tidyverse, caret, reshape2, dplyr, ggplot2,
  DT, shinythemes, data.table, MASS, zip, compositions,
  scoringRules, loo
)

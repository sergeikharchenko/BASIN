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
- **Proxy-Aware 1D Range Test:** Evaluates source–mixture range compatibility with optional normalization for hydrodynamic sorting or enrichment proxies (e.g., SSA, <i>d</i><sub>50</sub>, TOC).
- **Inter-Tracer Association Preservation Test (APT / GAPT):** Retains co-migrating tracers that preserve their structural pairwise relationships between sources and target mixtures.
- **Multivariate PCA Convex-Hull Leverage Filter:** An <i>n</i>-dimensional drop-out diagnostic that identifies and removes individual tracers pulling target mixtures outside the source mixing polygon.
- **Kruskal–Wallis & Collinearity Pruning:** Confirms inter-source discriminatory power (<i>p</i> &lt; 0.05) while purging near-duplicate signals (|&rho;| &gt; 0.95).

### 2. Four-Axis Bayesian Unmixing Engine
- **Selectable Residual Error-Scale Specifications:**
  - *Free Regularized Residual Scales (Default):* Assigns each tracer a weakly regularized scale parameter (&sigma;<sub>mix,<i>t</i></sub> &sim; Half-Student-<i>t</i>(3, 0, &sigma;<sub>prior</sub>)), acting as an adaptive Bayesian M-estimator. Decoupling variance normalization from squared proportions (<i>P</i><sub><i>i,g</i></sub><sup>2</sup>) eliminates the mathematical penalty against pure end-members (<i>P<sub>g</sub></i> = 0 or 1) and prevents centroid shrinkage on boundary-rich mixtures, even for isolated single targets (<i>N</i><sub>mix</sub> = 1).
  - *Process-Error Only (Source-Variance Propagated):* Deterministically propagates intra-source variances (&sum;<sub><i>g</i></sub> <i>P</i><sub><i>i,g</i></sub><sup>2</sup> &sigma;<sub>src,<i>g,t</i></sub><sup>2</sup>) into the mixture likelihood—optimal for interior-dominated, low-source-count (<i>G</i> = 2) mixtures with strongly heterogeneous source variances.
- **Heavy-Tailed Multivariate Likelihood & Covariance Modes:** Models target observations using a multivariate Student's <i>t</i>-distribution (&nu; = 4) to downweight erratic biogeochemical outliers, paired with three Cholesky-factorized error structures:
  - *Independent (Diagonal)*
  - *Target Mixture Covariance (&Sigma;<sub>mix</sub>, LKJ-regularized)*
  - *Full Bayes (Simultaneous Source + Mixture Covariance, &Sigma;<sub>src</sub> + &Sigma;<sub>mix</sub>)*
- **Controlled Gradient of Information Sharing:** Seamlessly switch between *Strictly Isolated Inference* (<i>N</i><sub>mix</sub> = 1), *Non-Hierarchical Joint Noise-Learning* (independent sample proportions with shared error parameters), *Hierarchical Partial Pooling* (multi-level Dirichlet hierarchy for stationary ensembles), and *Complete Pooling* (technical replicates), with direct support for prior geographical/connectivity weights (<i>w</i><sub>conn</sub>).
- **Mechanistic Non-Conservatism & Compositional Buffers:**
  - *Sorting Correction (Beta-module):* Evaluates tracer-specific non-linear exponential sorting adjustments (&beta;<sub>sorting,<i>t</i></sub>) directly within the MCMC sampler.
  - *Systematic Bias Absorption (&delta;<sub>bias</sub>):* A heavy-tailed displacement vector that isolates coherent ensemble-wide offsets (e.g., heavy-mineral lag or minor unsampled erosion hotspots contributing &le; 5–10%) without skewing known source proportions.
  - *Compositional Data Analysis (CoDA):* Isometric log-ratio (ILR) transformation for strictly conservative closed compositional datasets.

### 3. Integrated Diagnostics, Cross-Validation & Virtual Mixtures
- **Multi-Tracer Mass-Balance Reconstruction:** Evaluates internal mathematical closure (`sMAPE`, `Mean Bias`, and `NRMSE_src` per tracer) to diagnose discrepancies between internal fit and external accuracy.
- **Bayesian Model Comparison:** Natively integrates approximate **PSIS-LOO** (with Pareto-<i>k</i> outlier diagnostics) and **Exact LOO-CV** to rank competing model configurations by Expected Log Predictive Density (`elpd_loo`) and prevent overfitting.
- **Stochastic & Deterministic Virtual Mixtures:** Generates virtual mixtures across 4 regimes (Full Stochastic, Semi-Stochastic Simplex Grid, Semi-Deterministic, Full Deterministic) and reports point accuracy (`MAE`, `RMSE`, `R2`, `NSE`), proper scoring rules (`CRPS`), and boundary-adjusted empirical interval coverage (`W50`, `W95`, `P50/P95 adj`).

---

## Practitioner Deployment Guidelines

| Sampling Design & Process Regime | Recommended Pooling Mode | Error-Scale & Covariance Specification | Non-Conservatism Buffers |
| :--- | :--- | :--- | :--- |
| **1. Single conservative / weakly non-conservative target (<i>N</i><sub>mix</sub> = 1)** | No Pooling (Isolated) | Free Regularized Scales + **Full Bayes** Covariance | None (keep parsimonious) |
| **2. Single target with suspected transport shifts (<i>N</i><sub>mix</sub> = 1)** | No Pooling (Isolated) | Free Regularized Scales + Independent or MixCov | **Beta-Correction** &/or **Bias Absorption** |
| **3. Stationary ensemble (replicate traps, background core layers)** | **Hierarchical Partial Pooling** | Free Scales (or Process-Error if interior-rich) + MixCov | Optional; verify via PSIS-LOO |
| **4. Multi-event ensemble or core with extreme flood/typhoon pulses** | **No Pooling** (preserves pulse amplitude) | Free Regularized Scales + MixCov or Full Bayes | **Bias Absorption** + Beta (if proxy available) |
| **5. Interior-rich, <i>G</i> = 2 design with heterogeneous source variances** | No Pooling or Hierarchical | **Process-Error Only** (calibrates 50% CI coverage) | None unless sorting is severe |

---

## Prerequisites and Installation

BASIN runs locally in your standard web browser. With each execution, the Stan C++ code corresponding to the selected module combination is dynamically generated and compiled "on the fly", ensuring cross-platform reproducibility and full mathematical transparency.

### 1. Install R
- Download and install [R](https://cran.r-project.org/) (version &ge; 4.1 recommended).
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
```

---

## Running BASIN & Input Data Format

1. Clone or download this repository to your local machine.
2. Open R or RStudio, set the working directory to the repository folder, and launch the app:

```R
shiny::runApp("app.R")
```

3. Prepare your input dataset as a standard comma-separated values (`*.csv`) file:
   - **Column 1 (`ID`):** Unique sample identifiers (e.g., `S_01`, `Mix_01`).
   - **Column 2 (`Group`):** Source category names (e.g., `Forest`, `Cropland`, `Subsoil`) and a single label for target mixtures (e.g., `Mix`).
   - **Optional Columns:** A sorting/enrichment proxy column (e.g., `SSA`, `d50`, `TOC`) and/or a grouping covariate column for hierarchical models (e.g., `Season`, `Core_Zone`).
   - **Remaining Numeric Columns:** Candidate tracer values (geochemical elements, colorimetric indices, radionuclides, magnetic parameters, etc.).

Every completed model run is automatically archived in the local `BASIN_Workspaces/` and `BASIN_Saved_Results/` directories (including raw MCMC chains, `.rds` Stanfit objects, posterior summary tables, diagnostic plots, and runtime logs).

---

## Data Availability, Citation & Acknowledgements

- **Benchmark Datasets & Reproducible Scripts:** Complete validation datasets (including the Mano Dam catchment, Fukushima, Japan; the Thio River mining catchment, New Caledonia; computationally spoiled datasets; and the Hayama Lake sediment core archive) alongside all benchmarking scripts are openly archived on Zenodo: **[https://doi.org/10.5281/zenodo.24013349](https://doi.org/10.5281/zenodo.24013349)**.
- **Developer:** Sergey Kharchenko (*Faculty of Geography, Lomonosov Moscow State University* & *Institute of Geography, Russian Academy of Sciences*). Contact: `xar4enkkoff@yandex.ru`.
- **Funding:** This research was supported by the Russian Science Foundation (RSF), grant No. 25-17-00143.

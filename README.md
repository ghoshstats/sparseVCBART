# sparseVCBART
R package for a sparse version of VCBART. Uses Bayesian Additive Regression Trees to estimate covariate effect functions in the linear varying coefficient model, with global–local shrinkage.

## Installation

### System requirements

* **R ≥ 4.1**
* A working C/C++ toolchain for compiling R packages:
  * **macOS:** Xcode Command Line Tools (`xcode-select --install`)
  * **Windows:** Rtools (matching your R version)
  * **Linux:** build essentials (`gcc g++ make`)

### R package dependencies
* `Rcpp`, `RcppArmadillo` (C++ interface & linear algebra)



From a local checkout:

```r
# In the package root directory:
R CMD build .
R CMD INSTALL sparseVCBART_*.tar.gz

```
---

## Reproducing the experiments

All scripts live in the `codes/` directory. A shared helper file, `codes/common_sim.R`, contains common simulation settings and utility functions used by multiple scripts. Wrappers for competitor methods are in `codes/competitor_wrappers.R`.

### Simulation 1

* **VCBART variants:** `codes/sim1_VCBART_variants.R`
* **Other competitors:** `codes/sim1_others.R`

### Simulation 2

* **VCBART variants:** `codes/sim2_VCBART_variants.R`
* **Other competitors:** `codes/sim2_others.R`

### Semi-synthetic and real data analyses
The scripts for the applied data analyses are provided below:

* Semi-synthetic Political Science Example: `codes/PolScidata.R`

* Health and Retirement Study (HRS) Analysis: `codes/hrs-analysis.R`

* Yeast Transcription Factor (TF) Example: `codes/yeast_analysis.R`

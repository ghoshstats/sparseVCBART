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

### Install `sparseVCBART`

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("ghoshstats/sparseVCBART")
```

From a local checkout:

```r
# In the package root:
R CMD build .
R CMD INSTALL sparseVCBART_*.tar.gz
# …or from R:
# install.packages("path/to/sparseVCBART_*.tar.gz", repos = NULL, type = "source")
```

---

## Reproducing the experiments

All scripts live in the `codes/` directory. A shared helper file, `codes/common_sim.R`, contains common simulation settings and utility
functions used by multiple scripts. Wrappers for competitor methods are in `codes/competitor_wrappers.R`.

### Simulation 1 (Section 4.1)

* **VCBART variants:** `codes/sim1_VCBART_variants.R`
* **Other competitors:** `codes/sim1_others.R`

### Simulation 2 (Section 4.1, high-dimensional)

* **VCBART variants:** `codes/sim2_VCBART_variants.R`
* **Other competitors:** `codes/sim2_others.R`

### Semi-synthetic analysis (Section 4.2)
The script in the directory `codes/realdata.R` reproduces the semi-synthetic study. 


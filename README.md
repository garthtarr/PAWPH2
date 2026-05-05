# PAWPH2: Penalized Adaptive Robust Proportional Hazard Model

## Overview

**PAWPH2** is an R package implementing the Penalized Adaptive Weighted Proportional Hazard (PAWPH) model for robust survival analysis. The method simultaneously performs variable selection and outlier detection in Cox proportional hazards regression by augmenting the model with per-observation intercepts that are penalised via an adaptive lasso.

---

## Origins and Source Materials

This package consolidates two previously separate codebases into a single, self-contained R package:

### 1. `prcox.R` (PAWPH2 project)

The core PAWPH methodology was implemented in a standalone R script (`prcox.R`) in the PAWPH2 project. This script defined eight functions:

| Original function | Role |
|---|---|
| `prcox()` | Base model fit via augmented design matrix |
| `cvf.prcox()` | Single cross-validation fold worker |
| `cv.prcox()` | 10-fold stratified cross-validation |
| `cv.prcox.adap()` | Two-step adaptive lasso cross-validation |
| `prcoxreg()` | Main 3-step PAWPH procedure |
| `SE.boot()` | Parallel bootstrap standard errors |
| `getCI()` | Confidence intervals from bootstrap output |
| `lik.surv()` | Survival likelihood on a held-out test set |

### 2. `ncvreg2` (custom fork of Patrick Breheny's ncvreg)

The PAWPH method relies on a modified version of `ncvsurv()` from Patrick Breheny's [ncvreg](https://github.com/pbreheny/ncvreg) package (v3.13.0). The standard `ncvsurv()` does not accept per-observation intercepts; the `ncvreg2` fork added two new arguments:

- `intercept`: a numeric vector of length *n* giving fixed per-observation offsets (γ_i) incorporated into the linear predictor at the start of coordinate descent.
- `permu.weight`: a numeric vector of length *n* used to weight the event indicator.

Rather than maintain a dependency on the local `ncvreg2` package, PAWPH2 copies across only the functions and C code that are actually needed, with full attribution to the original author.

**Functions adapted from ncvreg2 (all internal, not part of the public API):**

| Function | ncvreg2 source file |
|---|---|
| `ncvsurv()` | `R/ncvsurv.R` — modified |
| `coef.ncvreg()` | `R/predict.R` |
| `predict.ncvreg()` | `R/predict.R` |
| `predict.ncvsurv()` | `R/predict-surv.R` |
| `std()` | `R/std.R` |
| `setupLambdaCox()` | `R/setupLambdaCox.R` |
| `convexMin()` | `R/convexMin.R` |
| `lamNames()` | `R/lamNames.R` |

**C source files adapted from ncvreg2:**

| File | C routines |
|---|---|
| `src/cox-dh.c` | `cdfit_cox_dh()` — modified |
| `src/standardize.c` | `standardize()` |
| `src/maxprod.c` | `maxprod()` |
| `src/init.c` | `R_init_PAWPH2()` — modified |

All adapted files carry attribution comments identifying the original source (https://github.com/pbreheny/ncvreg) and are covered under the GPL-3 licence, consistent with the original.

---

## Changes Made to Source Material

### C code

**`cox-dh.c`**

The functional logic is unchanged. Two updates were required to compile with modern toolchains:

1. `Calloc(n, type)` → `R_Calloc(n, type)`: The `Calloc` macro is deprecated in R ≥ 4.4 and not recognised by Apple Clang 21+ under `-std=gnu2x`. Replaced throughout with the current `R_Calloc` form.
2. `Free(ptr)` → `R_Free(ptr)`: Same deprecation; updated throughout `cleanupCox()`.

**`init.c`** (from `ncvreg_init.c`)

- Renamed `R_init_ncvreg` → `R_init_PAWPH2` so the dynamic library entry point matches the package name.
- Removed `extern` declarations and `CallEntries` registrations for C routines not used by PAWPH2: `cdfit_binomial`, `cdfit_gaussian`, `cdfit_poisson`, `rawfit_gaussian`, `mfdr_binomial`, `mfdr_cox`, `mfdr_gaussian`.
- Retained: `cdfit_cox_dh`, `maxprod`, `standardize`, and all shared helper functions (`crossprod`, `wcrossprod`, `wsqsum`, `sqsum`, `sum`, `MCP`, `SCAD`, `lasso`).

`standardize.c` and `maxprod.c` are copied verbatim (no changes needed).

### R code — ncvreg2 internals

- All `.Call("name", ...)` invocations updated to `.Call("name", ..., PACKAGE = "PAWPH2")` to ensure the C symbols are resolved from this package's shared library, not from any other loaded package.
- `coef.ncvreg()`, `predict.ncvreg()`, `predict.ncvsurv()` decorated with `@method` and `@export` roxygen2 tags so they are properly registered as S3 methods in `NAMESPACE`.
- `approx()` and `approxfun()` calls made explicit with `stats::` prefix.
- `convexMin()` trimmed to the Cox family branch only (other families — gaussian, binomial, poisson — are not used in PAWPH2 but the structure is preserved with a comment).

### R code — PAWPH functions

**`SE.boot()`** (from `prcox.R`):

- Removed `source("prcox.R")` from the `foreach` body — now unnecessary because all functions are part of the package.
- Replaced the hardcoded `ncvreg2` in `.packages` with `PAWPH2`.
- Replaced the hardcoded `makeCluster(16)` with `makeCluster(ncores)` where `ncores` defaults to `max(1, parallel::detectCores() - 1)`, making it portable across machines.

**`prcoxreg()`** (from `prcox.R`):

- `all.equal(alpha, 1) == T` → `isTRUE(all.equal(alpha, 1))` to avoid the fragile `== T` comparison.

**`cv.prcox()`**:

- `quantile(x, 0.8)` → `stats::quantile(x, 0.8, na.rm = TRUE)` to use explicit namespace and handle edge cases.

**`lik.surv()`** (internal, `lik-surv.R`):

- Added `val.basehaz <- pmax(val.basehaz, .Machine$double.eps)` after the `gbm::basehaz.gbm()` call for the instantaneous hazard. The GBM smoother can return slightly negative values due to smoothing artefacts; passing these to `log()` produced `NaN` warnings. Clamping to `.Machine$double.eps` before the log is taken suppresses the warning and is physically correct (hazard rates cannot be negative).

---

## Method Description

The PAWPH model extends the standard Cox proportional hazards model to be robust against outlying survival times. It augments the predictor matrix with an *n* × *n* identity matrix, giving each observation its own intercept γ_i. These intercepts are penalised via an adaptive lasso so that only observations whose survival time is genuinely surprising given their covariates receive a non-zero γ_i.

### The three-step procedure (`prcoxreg`)

**Step 1 — Variable selection** (when `alpha < 1`):
An adaptive lasso cross-validation is run over a grid of `alpha` values (the mixing weight between covariate and intercept penalties). This identifies relevant covariates. A two-step adaptive lasso is used: an initial fit provides weights `1/|β̂_j|` and `1/|γ̂_i|` for a second, reweighted fit.

**Step 2 — Outlier detection**:
The event indicator is temporarily set to 1 for all observations (treating censored observations as events). The adaptive lasso is then run on the selected covariates, allowing the per-observation intercepts to absorb any unexplained survival heterogeneity. Observations with large negative γ_i are unexpectedly long survivors (masking outliers); positive γ_i indicates unexpectedly short survival (swamping outliers). Censored observations with `γ_i > −2` have their intercept zeroed (they are unlikely to be genuine outliers).

**Step 3 — Refitting**:
A standard unpenalised Cox regression is fit on the selected covariates with the estimated γ_i values held fixed. This gives the refined estimates `betaHat_re`, unconfounded by the outlying observations.

After Step 3, any intercept with `|γ_i| < 0.5` is set to zero (a thresholding step to avoid spurious outlier flags).

### Cross-validation detail

- Folds are assigned by stratifying on event status to preserve the event rate within each fold.
- The validation metric is the mean negative log-likelihood, with the top 20% of individual losses clipped before averaging. This reduces the influence of observations with extreme likelihood values on lambda selection.
- The baseline hazard for the test-set likelihood is estimated using the GBM smooth estimator (`gbm::basehaz.gbm()`).

---

## Installation

```r
# Install from source
install.packages("~/PAWPH2", repos = NULL, type = "source")

# Or, if using devtools:
devtools::install("~/PAWPH2")
```

**Dependencies:** `survival`, `gbm`, `dplyr`, `doParallel`, `foreach`, `parallel`

---

## Usage

### Simulating example data

The examples below use a simulated dataset with `n = 200` observations, `p = 8` covariates, two truly active predictors (X1 and X2), and five injected outliers with unexpectedly long survival times.

```r
library(PAWPH2)
library(survival)

set.seed(8471)
n <- 200; p <- 8

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# True model: X1 (β = 1.2) and X2 (β = −0.8) are active; rest are noise
beta_true <- c(1.2, -0.8, rep(0, p - 2))
lp        <- X %*% beta_true

time   <- rexp(n, rate = exp(lp))   # exponential survival times
status <- rbinom(n, 1, 0.8)         # ~20% censoring

# Inject 5 outliers: unexpectedly long survivors
outlier_idx         <- sample(n, 5)
time[outlier_idx]   <- time[outlier_idx] * 20

y <- cbind(time = time, status = status)
```

### Basic usage — outlier detection only

Set `alpha = 1` to run outlier detection without variable selection.

```r
res <- prcoxreg(y, X, seed = 5312, alpha = 1)

# Refined coefficient estimates (Step 3)
res$betaHat_re
#>     X1     X2     X3     X4     X5     X6     X7     X8
#>  1.217 -0.893  0.001  0.097  0.125 -0.164  0.128  0.000

# Per-observation outlier intercepts (non-zero = flagged as outlier)
# Large negative values indicate unexpectedly long survivors
res$gammaHat[res$gammaHat != 0]
#>      13      41      43     133     151     168
#>  -1.007  -1.657  -1.061  -2.486  -1.219  -0.569
```

### With variable selection

Supply a vector of `alpha` values to sweep. Values less than 1 allocate some penalty weight to covariate selection; `alpha = 0.5` gives equal weight to both. The optimal `alpha` is chosen by cross-validation.

```r
res <- prcoxreg(y, X, seed = 5312, alpha = seq(0.1, 1, by = 0.1))

# Cross-validation selected alpha
res$opt.alpha
#> [1] 0.6

# Step-1 variable selection estimates (sparse — only X1 and X2 selected)
res$betaHat
#>    X1     X2     X3     X4     X5     X6     X7     X8
#>  1.509 -1.062  0.000  0.000  0.000  0.000  0.000  0.000

# Step-3 refined estimates — primary estimates for inference
res$betaHat_re
#>    X1     X2     X3     X4     X5     X6     X7     X8
#>  1.166 -0.899  0.000  0.000  0.000  0.000  0.000  0.000
```

### With bootstrap confidence intervals

Set `cond.CI = TRUE` to compute bootstrap standard errors and confidence intervals for the refined estimates. This is computationally intensive; adjust `B` (number of bootstrap replicates) and `ncores` as needed.

```r
res <- prcoxreg(y, X, seed = 5312, alpha = 1, cond.CI = TRUE, B = 200)

# CI matrix rows: SE, lower_n, upper_n, lower_q, upper_q, z.stat, pvalue
res$CI
```

### Lower-level functions

```r
# Cross-validation only (returns full CV path)
cv_res <- cv.prcox(y, X, seed = 5312)

# Adaptive lasso CV over an alpha grid
adap_res <- cv.prcox.adap(y, X, seed = 5312, alpha = c(0.5, 0.75, 1))

# Base model fit (returns full lambda path)
fit <- prcox(y, X)

# Bootstrap SEs separately, then extract a CI for the first coefficient
se_res <- SE.boot(y, X, B = 100, seed = 5312, alpha = 1)
getCI(se_res$betaHat[, 1], b = se_res$betaHat_re[1])
```

---

## Package structure

```
PAWPH2/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── PAWPH2-package.R      # package-level imports
│   ├── prcox.R               # prcox(), cvf.prcox()
│   ├── cv-prcox.R            # cv.prcox(), cv.prcox.adap()
│   ├── prcoxreg.R            # prcoxreg() — main entry point
│   ├── bootstrap.R           # SE.boot(), getCI()
│   ├── lik-surv.R            # lik.surv() — internal
│   ├── ncvsurv.R             # ncvsurv() — adapted from ncvreg2
│   ├── ncvreg-internals.R    # std(), setupLambdaCox(), convexMin(), lamNames()
│   ├── predict-ncvreg.R      # coef.ncvreg(), predict.ncvreg() S3 methods
│   └── predict-ncvsurv.R     # predict.ncvsurv() S3 method
└── src/
    ├── cox-dh.c              # coordinate descent for Cox (modified)
    ├── standardize.c         # design matrix standardisation
    ├── maxprod.c             # lambda_max computation
    ├── init.c                # C routine registration (modified)
    └── Makevars
```

---

## Attribution

The C source files `cox-dh.c`, `standardize.c`, and `maxprod.c`, along with the R functions `ncvsurv()`, `std()`, `setupLambdaCox()`, `convexMin()`, `lamNames()`, `coef.ncvreg()`, `predict.ncvreg()`, and `predict.ncvsurv()` are adapted from **ncvreg v3.13.0** by Patrick Breheny (University of Iowa), licensed under GPL-3.

> Breheny P and Huang J (2011). Coordinate descent algorithms for nonconvex penalized regression. *Annals of Applied Statistics*, 5(1), 232–253. https://doi.org/10.1214/10-AOAS388

Original ncvreg source: https://github.com/pbreheny/ncvreg

---

## Licence

GPL-3

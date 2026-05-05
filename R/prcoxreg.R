#' Fit the PAWPH robust Cox regression model
#'
#' The main user-facing function for the Penalized Adaptive Weighted
#' Proportional Hazard (PAWPH) model. Implements a three-step procedure for
#' simultaneous variable selection and outlier detection in survival data:
#'
#' **Step 1 — Variable selection** (when `alpha < 1`): An adaptive lasso
#' cross-validation ([cv.prcox.adap()]) is run with the supplied `alpha` vector
#' to identify relevant covariates.
#'
#' **Step 2 — Outlier detection**: The event indicator is temporarily set to 1
#' for all observations (treating censored observations as events) and an
#' adaptive lasso is run on the selected covariates (or all covariates if none
#' were selected in Step 1). The resulting per-observation intercepts
#' (gamma_i) flag outliers: a large negative gamma_i indicates an observation
#' whose survival time is unexpectedly long given its covariates.
#'
#' **Step 3 — Refitting**: A standard unpenalised Cox regression is fit on the
#' selected covariates with the estimated gamma_i intercepts held fixed. This
#' gives refined coefficient estimates `betaHat_re` unconfounded by outliers.
#'
#' After Step 2, gamma_i values for censored observations that are not strongly
#' outlying (`gamma_i > -2`) are set to zero, and after Step 3 all small
#' intercepts (`|gamma_i| < 0.5`) are zeroed.
#'
#' @param y Survival outcome matrix (n x 2): time and event indicator. Typically
#'   created with [survival::Surv()].
#' @param X Predictor matrix (n x p).
#' @param seed Integer seed for reproducible cross-validation fold assignment.
#' @param alpha Numeric vector of penalty mixing weights to sweep over via
#'   [cv.prcox.adap()]. Values in (0, 1): proportion of penalty allocated to
#'   per-observation intercepts. Use `alpha = 1` (default) to skip variable
#'   selection and detect outliers only.
#' @param cond.CI Logical. If `TRUE`, compute bootstrap confidence intervals for
#'   `betaHat_re` via [SE.boot()]. Default `FALSE`.
#' @param pout.max Maximum proportion of observations that may be flagged as
#'   outliers. Default `0.2` (20%). Set `NULL` to disable the constraint.
#' @param permu.weight Numeric vector of length n giving per-observation weights
#'   applied to the event indicator in [ncvsurv()]. Default is `rep(1, n)`.
#' @param ... Additional arguments passed to [cv.prcox.adap()].
#'
#' @return A list (the result of [cv.prcox.adap()]) with the following
#'   additional or modified elements:
#'   \describe{
#'     \item{`betaHat`}{Coefficient estimates from Step 1 variable selection
#'       (length p). Zero for unselected variables.}
#'     \item{`betaHat_re`}{Refined coefficient estimates from Step 3 (length p).
#'       These are the primary estimates for inference.}
#'     \item{`gammaHat`}{Per-observation intercept estimates (length n). Values
#'       far from zero indicate outliers; negative values correspond to
#'       unexpectedly long survival (masking outliers), positive values to
#'       unexpectedly short survival (swamping outliers).}
#'     \item{`opt.alpha`}{Selected alpha value.}
#'     \item{`cve`}{Cross-validation error path.}
#'     \item{`CI`}{(If `cond.CI = TRUE`) Matrix (7 x p) of confidence interval
#'       summaries with rows: SE, lower_n, upper_n, lower_q, upper_q, z.stat,
#'       pvalue.}
#'     \item{`se.res`}{(If `cond.CI = TRUE`) Full bootstrap result from
#'       [SE.boot()].}
#'   }
#'
#' @seealso [cv.prcox()], [cv.prcox.adap()], [prcox()], [SE.boot()], [getCI()]
#'
#' @references
#'   Breheny P and Huang J (2011). Coordinate descent algorithms for nonconvex
#'   penalized regression. *Annals of Applied Statistics*, 5(1), 232–253.
#'   \doi{10.1214/10-AOAS388}
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(MASS)
#'
#' set.seed(42)
#' n <- 300; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("X", 1:p)
#' # True coefficients: X1 and X2 are active, rest zero
#' beta_true <- c(1, -0.8, rep(0, p - 2))
#' lp <- X %*% beta_true
#' time <- rexp(n, rate = exp(lp))
#' status <- rbinom(n, 1, 0.8)
#' y <- cbind(time, status)
#'
#' res <- prcoxreg(y, X, seed = 100, alpha = seq(0.1, 1, by = 0.1))
#' res$betaHat_re   # refined estimates
#' res$gammaHat     # outlier scores
#' }
#'
#' @export
prcoxreg <- function(y, X,
                     seed,
                     alpha = 1,
                     cond.CI = FALSE,
                     pout.max = 0.2,
                     permu.weight = rep(1, n),
                     ...) {
  p <- ncol(X)
  n <- nrow(X)
  delta <- y[, 2]

  is.vs    <- !(isTRUE(all.equal(alpha, 1)))
  is.cens  <- sum(delta) != n
  prop.cens <- sum(1 - delta) / n

  # Allow extra room for the outlier proportion constraint when censoring is present
  pout.max2 <- NULL
  if (!is.null(pout.max)) {
    pout.max2 <- pout.max + prop.cens
  }

  if (is.vs) {
    # Step 1: variable selection
    res.vs <- cv.prcox.adap(y, X, seed, alpha = alpha, ...)
    betaHat.vs <- res.vs$betaHat

    # Step 2: outlier detection on selected variables (treat all obs as events)
    y[, 2] <- rep(1, n)
    if (sum(betaHat.vs != 0) > 0) {
      res <- cv.prcox.adap(y, as.matrix(X[, which(betaHat.vs != 0)]),
                           seed, alpha = 1, pout.max = pout.max2)
      res$betaHat  <- betaHat.vs
      res$opt.alpha <- res.vs$opt.alpha
    } else {
      res <- cv.prcox.adap(y, X, seed, alpha = alpha, pout.max = pout.max2, ...)
    }

  } else {
    # No variable selection: go directly to outlier detection
    y[, 2] <- rep(1, n)
    res <- cv.prcox.adap(y, X, seed, alpha = alpha, pout.max = pout.max2,
                         permu.weight = permu.weight)
  }

  # Zero out intercepts for censored observations that are not strongly outlying
  ind.cens <- which(delta == 0)
  res$gammaHat[ind.cens] <- ifelse(res$gammaHat[ind.cens] > -2, 0, res$gammaHat[ind.cens])

  # Restore original event indicator
  y[, 2] <- delta

  # Step 3: refit on selected variables with fixed intercepts
  if (sum(res$betaHat != 0) > 0) {
    fit <- ncvsurv(
      as.matrix(X[, which(res$betaHat != 0)]),
      y,
      penalty    = "lasso",
      intercept  = res$gammaHat,
      lambda     = 0,
      returnX    = FALSE,
      permu.weight = permu.weight
    )
    betaHat_re <- rep(0, p)
    betaHat_re[which(res$betaHat != 0)] <- as.numeric(fit$beta)
    res$betaHat_re <- betaHat_re
  } else {
    res$betaHat_re <- res$betaHat
  }
  names(res$betaHat_re) <- names(res$betaHat)

  # Apply threshold: small intercepts are treated as zero
  res$gammaHat <- ifelse(abs(res$gammaHat) < 0.5, 0, res$gammaHat)

  # Conditional bootstrap CIs
  if (cond.CI) {
    se.res <- SE.boot(y, X,
                      pout.max = pout.max,
                      alpha    = ifelse(is.vs, res$opt.alpha, 1),
                      seed     = seed,
                      ...)
    res$CI <- sapply(1:ncol(X), function(i) {
      getCI(se.res$betaHat[, i], b = res$betaHat_re[i])
    })
    rownames(res$CI) <- c("SE", "lower_n", "upper_n", "lower_q", "upper_q", "z.stat", "pvalue")
    colnames(res$CI) <- names(res$betaHat)
    res$se.res <- se.res
  }

  res
}

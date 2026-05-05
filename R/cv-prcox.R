#' Cross-validation for the PAWPH base model
#'
#' Performs k-fold cross-validation for [prcox()] to select the regularisation
#' parameter lambda. Event status is stratified across folds so that each fold
#' has approximately the same proportion of events. The validation metric is the
#' mean negative log-likelihood, clipped at the 80th percentile to reduce the
#' influence of extreme values.
#'
#' @param y Survival outcome matrix (n x 2): time and event indicator.
#' @param X Predictor matrix (n x p).
#' @param nfolds Number of CV folds. Default is `10`.
#' @param seed Optional integer seed for reproducible fold assignment.
#' @param penalty.factor Per-variable penalty multipliers (length `p + n`).
#'   Defaults to `c(rep(0, p), rep(1, n))`.
#' @param pout.max Maximum allowable proportion of observations flagged as
#'   outliers. Lambda values that exceed this proportion have their CV error set
#'   to `Inf` and are effectively excluded from selection. Set to `NULL` to
#'   disable (default `NULL`).
#' @param dfmax Maximum number of active variables in the model. Passed to
#'   [prcox()].
#' @param ... Additional arguments passed to [prcox()].
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`betaHat`}{Estimated coefficients for the `p` covariates at the
#'       optimal lambda.}
#'     \item{`gammaHat`}{Estimated per-observation intercepts (length n) at the
#'       optimal lambda.}
#'     \item{`fit`}{The full [prcox()] fit object (all lambda values).}
#'     \item{`cve`}{Vector of cross-validation errors (mean clipped neg
#'       log-likelihood) for each lambda.}
#'     \item{`cve.min`}{Minimum CV error.}
#'     \item{`opt.lam`}{Optimal lambda value.}
#'     \item{`lik.val`}{Matrix of per-observation validation likelihoods.}
#'     \item{`lambda`}{Lambda sequence used.}
#'   }
#'
#' @seealso [prcox()], [cv.prcox.adap()], [prcoxreg()]
#'
#' @export
cv.prcox <- function(y, X,
                     nfolds = 10,
                     seed,
                     penalty.factor = c(rep(0, p), rep(1, n)),
                     pout.max = NULL,
                     dfmax = 100,
                     ...) {
  n <- nrow(X)
  p <- ncol(X)
  fit <- prcox(y, X, penalty.factor = penalty.factor, ...)

  # Set up stratified folds (maintains event rate within each fold)
  fail <- y[, 2]
  if (!missing(seed)) set.seed(seed)
  ind1 <- which(fail == 1)
  ind0 <- which(fail == 0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1 == 0] <- nfolds
  fold0[fold0 == 0] <- nfolds
  fold <- integer(n)
  fold[fail == 1] <- sample(fold1)
  fold[fail == 0] <- sample(fold0)

  # Compute validation likelihoods across all folds
  lik.val <- matrix(Inf, nrow = n, ncol = length(fit$lambda))
  colnames(lik.val) <- fit$lambda
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$penalty.factor <- penalty.factor
  for (i in 1:nfolds) {
    res <- cvf.prcox(i, X, y, fold, cv.args)
    lik.val[fold == i, colnames(res$lik)] <- res$lik
  }

  # Clipped mean: exclude top 20% to reduce influence of extreme values
  cve <- apply(lik.val, 2, function(x) {
    xx <- x[x < stats::quantile(x, 0.8, na.rm = TRUE)]
    mean(xx, na.rm = TRUE)
  })

  # Enforce outlier proportion constraint if requested
  if (!is.null(pout.max)) {
    prop.out <- apply(fit$beta[-(1:p), ], 2, function(b) sum(b != 0)) / n
    lam.ind <- which(prop.out <= pout.max)
    cve[-lam.ind] <- Inf
  }

  # Exclude lambda values that produced any NA coefficients
  cve[which(apply(fit$beta, 2, anyNA))] <- Inf

  lam <- fit$lambda[which.min(cve)]
  betaHat <- coef.ncvreg(fit, lambda = lam)
  gammaHat <- betaHat[(p + 1):(n + p)]
  betaHat <- betaHat[1:p]

  list(
    betaHat  = betaHat,
    gammaHat = gammaHat,
    fit      = fit,
    cve      = cve,
    cve.min  = min(cve, na.rm = TRUE),
    opt.lam  = lam,
    lik.val  = lik.val,
    lambda   = fit$lambda
  )
}


#' Adaptive lasso cross-validation for the PAWPH model
#'
#' Implements a two-step adaptive lasso procedure for the PAWPH model over a
#' grid of `alpha` values (mixing weight between covariate and intercept
#' penalties). For each `alpha`:
#'
#' 1. An initial [cv.prcox()] fit is computed.
#' 2. Adaptive lasso weights are derived from the initial estimates:
#'    `1/|beta_j|` for covariates, `1/|gamma_i|` for intercepts.
#' 3. A second [cv.prcox()] is fit with these weights.
#'
#' The `alpha` value minimising the second-step CV error is selected.
#'
#' @param y Survival outcome matrix (n x 2): time and event indicator.
#' @param X Predictor matrix (n x p).
#' @param seed Optional integer seed for reproducible fold assignment.
#' @param alpha Numeric vector of penalty mixing weights (between 0 and 1) to
#'   sweep over. Default is `1` (intercepts only penalised). Values close to 0
#'   favour covariate selection; values close to 1 favour intercept selection.
#' @param penalty.factor Optional additional per-variable multipliers (length
#'   `p + n`). Defaults to all-ones.
#' @param lambda.min Smallest lambda as a fraction of lambda.max. Default `0.05`.
#' @param pout.max Maximum proportion of observations that may be flagged as
#'   outliers. Passed to [cv.prcox()].
#' @param ... Additional arguments passed to [cv.prcox()].
#'
#' @return The [cv.prcox()] result list from the best `alpha`, with an
#'   additional element `opt.alpha` giving the selected `alpha` value.
#'
#' @seealso [cv.prcox()], [prcoxreg()]
#'
#' @export
cv.prcox.adap <- function(y, X,
                           seed,
                           alpha = 1,
                           penalty.factor = rep(1, p + n),
                           lambda.min = 0.05,
                           pout.max = NULL,
                           ...) {
  p <- ncol(X)
  n <- nrow(X)
  res_list <- array(list(), length(alpha))

  for (i in seq_along(alpha)) {
    # Step 1: initial fit with uniform weights scaled by alpha
    res0 <- cv.prcox(y, X,
                     seed = seed,
                     penalty.factor = c(rep(1 - alpha[i], p), rep(alpha[i], n)) * penalty.factor,
                     lambda.min = lambda.min,
                     pout.max = pout.max,
                     ...)
    gammaHat <- res0$gammaHat
    betaHat  <- res0$betaHat

    # Step 2: adaptive lasso refit with inverse-magnitude weights
    res_list[[i]] <- cv.prcox(y, X,
                               seed = seed,
                               penalty.factor = c(
                                 (1 - alpha[i]) / pmax(abs(betaHat), 1e-5),
                                 alpha[i] / pmax(abs(gammaHat), 1e-5)
                               ) * penalty.factor,
                               lambda.min = lambda.min,
                               pout.max = pout.max,
                               ...)
  }

  cve <- sapply(res_list, getElement, "cve.min")
  i.alpha <- which.min(cve)
  res <- res_list[[i.alpha]]
  res$opt.alpha <- alpha[i.alpha]
  res
}

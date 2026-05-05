#' Fit a PAWPH base model
#'
#' Augments the predictor matrix `X` with an n x n identity matrix (one column
#' per observation), then fits a penalised Cox regression path via [ncvsurv()].
#' The augmented columns allow the model to learn per-observation intercepts
#' (gamma_i), which are the basis of PAWPH outlier detection.
#'
#' The covariate columns in `X` are left unpenalised by default
#' (`penalty.factor = 0`), while the augmented intercept columns are penalised
#' (`penalty.factor = 1`). This can be overridden via the `penalty.factor`
#' argument.
#'
#' @param y Survival outcome matrix with two columns: time and event indicator.
#'   Typically created with [survival::Surv()].
#' @param X Predictor matrix (n x p).
#' @param penalty.factor Numeric vector of length `p + n` giving per-variable
#'   penalty multipliers. Defaults to `c(rep(0, p), rep(1, n))` — covariates
#'   unpenalised, intercepts penalised.
#' @param lambda.min Smallest lambda value as a fraction of lambda.max.
#'   Default is `0.05`.
#' @param lambda Optional user-supplied lambda sequence. If missing, an
#'   automatic sequence is computed.
#' @param ... Additional arguments passed to [ncvsurv()].
#'
#' @return An object of class `c("ncvsurv", "ncvreg")` as returned by
#'   [ncvsurv()], fit on the augmented design matrix `[X | I_n]`.
#'
#' @seealso [prcoxreg()], [cv.prcox()]
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- cbind(time = rexp(n), status = rbinom(n, 1, 0.8))
#' fit <- prcox(y, X)
#' }
#'
#' @export
prcox <- function(y, X,
                  penalty.factor = c(rep(0, p), rep(1, n)),
                  lambda.min = 0.05,
                  lambda,
                  ...) {
  p <- ncol(X)
  n <- nrow(X)
  X.aug <- cbind(X, diag(n))
  if (missing(lambda)) {
    ncvsurv(X.aug, y, penalty = "lasso", lambda.min = lambda.min,
            penalty.factor = penalty.factor, returnX = FALSE, ...)
  } else {
    ncvsurv(X.aug, y, penalty = "lasso", penalty.factor = penalty.factor,
            lambda = lambda, returnX = FALSE, ...)
  }
}


# Internal worker: fit one CV fold and return test-set likelihoods.
cvf.prcox <- function(i, XX, y, fold, cv.args) {
  p <- ncol(XX)
  cv.args$X <- XX[fold != i, , drop=FALSE]
  cv.args$y <- y[fold != i, ]
  cv.args$permu.weight <- cv.args$permu.weight[fold != i]
  cv.args$penalty.factor <- cv.args$penalty.factor[c(1:p, p + which(fold != i))]
  fit.i <- do.call("prcox", cv.args)
  prop.out <- apply(fit.i$beta, 2, function(beta) {
    sum(beta[-(1:p)] != 0) / length(cv.args$y)
  })
  X.aug <- cbind(cv.args$X, diag(nrow(cv.args$X)))
  lp.train <- X.aug %*% fit.i$beta
  test.X <- XX[fold == i, , drop=FALSE]
  test.y <- y[fold == i, ]
  lp.test <- test.X %*% fit.i$beta[1:p, ]
  lik <- lik.surv(
    t.train     = cv.args$y[, 1],
    delta.train = cv.args$y[, 2],
    lp.train    = lp.train,
    t.test      = test.y[, 1],
    delta.test  = test.y[, 2],
    lp.test     = lp.test
  )
  colnames(lik) <- fit.i$lambda
  list(lik = lik, prop.out = prop.out)
}

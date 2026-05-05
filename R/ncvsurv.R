# Adapted from ncvreg v3.13.0 by Patrick Breheny
# Original package: https://github.com/pbreheny/ncvreg
# Licensed under GPL-3.
#
# Modification: Added `intercept` (per-observation intercept vector) and
# `permu.weight` (permutation weights) arguments to enable the PAWPH outlier
# detection mechanism. These are passed through to the underlying C routine
# cdfit_cox_dh().

#' Fit a penalised Cox regression model
#'
#' Fits a regularisation path for Cox proportional hazards regression using
#' lasso, MCP, or SCAD penalties. This is an internal function adapted from
#' `ncvreg::ncvsurv()` with the addition of per-observation `intercept` and
#' `permu.weight` arguments required by the PAWPH method.
#'
#' @param X Predictor matrix (n x p). Columns corresponding to augmented
#'   identity matrix (per-observation intercepts) are added by [prcox()].
#' @param y Survival outcome matrix with two columns: time and event indicator.
#' @param penalty Penalty type: `"MCP"`, `"SCAD"`, or `"lasso"`.
#' @param gamma Concavity parameter for MCP/SCAD.
#' @param alpha Per-variable penalty mixing parameters (vector of length p).
#' @param lambda.min Smallest lambda as a fraction of lambda.max.
#' @param nlambda Number of lambda values in the path.
#' @param lambda Optional user-supplied lambda sequence.
#' @param eps Convergence tolerance.
#' @param max.iter Maximum coordinate descent iterations.
#' @param convex Check for local convexity? Default `TRUE`.
#' @param dfmax Maximum number of active variables.
#' @param penalty.factor Per-variable penalty multipliers. Set to 0 to leave
#'   a variable unpenalised.
#' @param warn Warn if algorithm does not converge?
#' @param returnX Return the standardised design matrix in the output?
#' @param intercept Numeric vector of length n giving per-observation intercepts
#'   (gamma_i values). Defaults to `rep(0, n)`. This is the key modification
#'   enabling PAWPH outlier detection.
#' @param permu.weight Numeric vector of length n giving per-observation weights
#'   applied to the event indicator. Defaults to `rep(1, n)`.
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class `c("ncvsurv", "ncvreg")`.
#'
#' @references Breheny P and Huang J (2011). Coordinate descent algorithms for
#'   nonconvex penalized regression. *Annals of Applied Statistics*, 5(1),
#'   232–253. \doi{10.1214/10-AOAS388}
#'
#' @export
ncvsurv <- function(X, y,
                    penalty = c("MCP", "SCAD", "lasso"),
                    gamma = switch(penalty, SCAD = 3.7, 3),
                    alpha = rep(1, p),
                    lambda.min = ifelse(n > p, .001, .05),
                    nlambda = 100,
                    lambda,
                    eps = 1e-4,
                    max.iter = 10000,
                    convex = TRUE,
                    dfmax = p,
                    penalty.factor = rep(1, ncol(X)),
                    warn = TRUE,
                    returnX,
                    intercept = rep(0, n),
                    permu.weight = rep(1, n),
                    ...) {
  p <- ncol(X)
  n <- nrow(X)

  # Coercion
  penalty <- match.arg(penalty)
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X) == "integer") storage.mode(X) <- "double"
  if (!inherits(y, "matrix")) {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be a matrix or able to be coerced to a matrix", call.=FALSE)
    if (ncol(y) != 2) stop("y must have two columns for survival data: time-on-study and a censoring indicator", call.=FALSE)
  }
  if (typeof(y) == "integer") storage.mode(y) <- "double"
  if (typeof(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"

  # Error checking
  if (gamma <= 1 & penalty == "MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty == "SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (max(alpha) <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor) != ncol(X)) stop("penalty.factor does not match up with X", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvsurv", call.=FALSE)

  # Set up XX, yy, lambda
  tOrder <- order(y[, 1])
  yy <- as.double(y[tOrder, 1])
  intercept <- intercept[tOrder]
  Delta <- y[tOrder, 2]
  n <- length(yy)
  XX <- std(X[tOrder, , drop=FALSE])
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  # Fit via C routine (cdfit_cox_dh with intercept support)
  res <- .Call("cdfit_cox_dh", XX, Delta * permu.weight, penalty, lambda,
               eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               as.double(alpha), as.integer(dfmax),
               as.integer(user.lambda | any(penalty.factor == 0)),
               as.integer(warn),
               as.double(intercept),
               PACKAGE = "PAWPH2")
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1 * res[[2]]
  iter <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)

  # Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (warn & sum(iter) == max.iter) warning("Algorithm failed to converge for some values of lambda")

  # Local convexity check
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda * (1 - max(alpha)), "cox", penalty.factor, Delta = Delta) else NULL

  # Unstandardize
  beta <- matrix(0, nrow = ncol(X), ncol = length(lambda))
  bb <- b / attr(XX, "scale")[ns]
  beta[ns, ] <- bb
  offset <- -crossprod(attr(XX, "center")[ns], bb)

  # Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  dimnames(beta) <- list(varnames, lamNames(lambda))

  # Output
  val <- structure(
    list(
      beta         = beta,
      iter         = iter,
      lambda       = lambda,
      penalty      = penalty,
      gamma        = gamma,
      alpha        = alpha,
      convex.min   = convex.min,
      loss         = loss,
      penalty.factor = penalty.factor,
      n            = n,
      time         = yy,
      fail         = Delta,
      order        = tOrder
    ),
    class = c("ncvsurv", "ncvreg")
  )
  val$Eta <- sweep(Eta, 2, offset, "-")

  if (missing(returnX)) {
    if (utils::object.size(XX) > 1e8) {
      warning("Due to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE.")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) val$X <- XX
  val
}

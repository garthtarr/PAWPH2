# Adapted from ncvreg v3.13.0 by Patrick Breheny
# Original package: https://github.com/pbreheny/ncvreg
# Licensed under GPL-3. Copied with minor formatting changes only.
#
# Functions: std(), setupLambdaCox(), convexMin(), lamNames()
# These are internal helper functions used by ncvsurv().

#' @keywords internal
std <- function(X) {
  if (typeof(X) == 'integer') storage.mode(X) <- 'double'
  if (!inherits(X, "matrix")) {
    if (is.numeric(X)) {
      X <- matrix(as.double(X), ncol=1)
    } else {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
  }
  STD <- .Call("standardize", X, PACKAGE = "PAWPH2")
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}

#' @keywords internal
setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)

  # Fit to unpenalized covariates
  ind <- which(penalty.factor != 0)
  if (length(ind) != p) {
    nullFit <- survival::coxph(survival::Surv(y, Delta) ~ X[, -ind, drop=FALSE])
    eta <- nullFit$linear.predictors
    rsk <- rev(cumsum(rev(exp(eta))))
    s <- Delta - exp(eta) * cumsum(Delta / rsk)
  } else {
    w <- 1 / (n - (1:n) + 1)
    s <- Delta - cumsum(Delta * w)
  }

  # Determine lambda.max
  zmax <- .Call("maxprod", X, s, ind, penalty.factor, PACKAGE = "PAWPH2") / n
  lambda.max <- zmax / max(alpha)

  if (lambda.min == 0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001 * lambda.max), len = nlambda - 1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), len = nlambda))
  }
  lambda
}

#' @keywords internal
convexMin <- function(b, X, penalty, gamma, l2, family, penalty.factor, a, Delta=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)

  if (penalty == "MCP") {
    k <- 1 / gamma
  } else if (penalty == "SCAD") {
    k <- 1 / (gamma - 1)
  } else if (penalty == "lasso") {
    return(NULL)
  }
  if (l == 0) return(NULL)

  val <- NULL
  for (i in 1:l) {
    A1 <- if (i == 1) rep(1, p) else b[, i] == 0
    if (i == l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[, i + 1] == 0
      U <- A1 & A2
      L2 <- l2[i + 1]
    }
    if (sum(!U) == 0) next
    Xu <- X[, !U]
    p.. <- k * (penalty.factor[!U] != 0) - L2 * penalty.factor[!U]
    if (family == "cox") {
      eta <- if (i == l) X %*% b[, i] else X %*% b[, i + 1]
      haz <- drop(exp(eta))
      rsk <- rev(cumsum(rev(haz)))
      h <- haz * cumsum(Delta / rsk)
      xwxn <- crossprod(sqrt(h) * Xu) / n
      eigen.min <- min(eigen(xwxn - diag(diag(xwxn) * p.., nrow(xwxn), ncol(xwxn)))$values)
    } else {
      # Other families not needed in PAWPH2 but preserved for completeness
      next
    }
    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}

#' @keywords internal
lamNames <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d, 4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}

# Adapted from ncvreg v3.13.0 by Patrick Breheny
# Original package: https://github.com/pbreheny/ncvreg
# Licensed under GPL-3. Copied with minor formatting changes only.
#
# Functions: coef.ncvreg(), predict.ncvreg()
# These S3 methods support coefficient extraction from fitted ncvsurv objects.

#' @method coef ncvreg
#' @export
#' @keywords internal
coef.ncvreg <- function(object, lambda, which = 1:length(object$lambda), drop = TRUE, ...) {
  if (!missing(lambda)) {
    if (max(lambda) > max(object$lambda) | min(lambda) < min(object$lambda)) {
      stop('Supplied lambda value(s) are outside the range of the model fit.', call.=FALSE)
    }
    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta <- (1 - w) * object$beta[, l, drop=FALSE] + w * object$beta[, r, drop=FALSE]
    colnames(beta) <- lamNames(lambda)
  } else {
    beta <- object$beta[, which, drop=FALSE]
  }
  if (drop) return(drop(beta)) else return(beta)
}

#' @method predict ncvreg
#' @export
#' @keywords internal
predict.ncvreg <- function(object, X,
                           type = c("link", "response", "class", "coefficients", "vars", "nvars"),
                           lambda, which = 1:length(object$lambda), ...) {
  type <- match.arg(type)
  beta <- coef.ncvreg(object, lambda = lambda, which = which, drop = FALSE)
  if (type == "coefficients") return(beta)
  if (!inherits(object, 'ncvsurv')) {
    alpha <- beta[1, ]
    beta <- beta[-1, , drop=FALSE]
  }
  if (type == "nvars") return(apply(beta != 0, 2, sum))
  if (type == "vars") return(drop(apply(beta != 0, 2, FUN = which)))
  eta <- sweep(X %*% beta, 2, alpha, "+")
  if (type == "link") return(drop(eta))
}

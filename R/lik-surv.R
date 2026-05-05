# Internal helper: compute survival likelihood on a test set.

#' Survival likelihood on a held-out test set
#'
#' For each lambda value, computes the negative log-likelihood of test
#' observations given a model fitted on training data. The baseline hazard is
#' estimated using the GBM smoother via [gbm::basehaz.gbm()].
#'
#' @param t.train Numeric vector of training event/censoring times.
#' @param delta.train Integer vector of training event indicators (1 = event).
#' @param lp.train Matrix (n_train x n_lambda) of linear predictors from the
#'   training fit, one column per lambda value.
#' @param t.test Numeric vector of test event/censoring times.
#' @param delta.test Integer vector of test event indicators.
#' @param lp.test Matrix (n_test x n_lambda) of linear predictors for test
#'   observations.
#'
#' @return A matrix of negative log-likelihoods with rows corresponding to
#'   test observations and columns corresponding to lambda values.
#'
#' @keywords internal
lik.surv <- function(t.train, delta.train, lp.train, t.test, delta.test, lp.test) {
  lik.mat <- -log(sapply(1:ncol(lp.train), function(j) {
    if (sum(!is.na(lp.train[, j])) == 0) {
      rep(0, length(t.test))
    } else {
      cum.basehaz <- gbm::basehaz.gbm(t.train, delta.train, f.x = lp.train[, j],
                                       t.eval = t.test, cumulative = TRUE)
      val.basehaz <- gbm::basehaz.gbm(t.train, delta.train, f.x = lp.train[, j],
                                       t.eval = t.test, smooth = TRUE, cumulative = FALSE)
      # GBM smoother can produce slightly negative values; clamp to avoid log(x <= 0)
      val.basehaz <- pmax(val.basehaz, .Machine$double.eps)
      St.base <- exp(-cum.basehaz)
      St <- St.base^(exp(lp.test[, j]))
      dplyr::if_else(delta.test == 1, St * val.basehaz * exp(lp.test[, j]), St)
    }
  }))
  lik.mat[lik.mat == Inf] <- NA
  lik.mat
}

# Suppress R CMD check NOTE for foreach loop variable 'ind'
utils::globalVariables("ind")

#' Bootstrap standard errors for PAWPH coefficient estimates
#'
#' Computes bootstrap standard errors for the refined coefficient estimates
#' (`betaHat_re`) from [prcoxreg()] by refitting the model on `B` bootstrap
#' samples in parallel.
#'
#' @param orig_y Survival outcome matrix (n x 2) from the original dataset.
#' @param orig_X Predictor matrix (n x p) from the original dataset.
#' @param B Number of bootstrap replicates. Default `100`.
#' @param seed Integer seed for reproducible bootstrap sampling. Default `1234`.
#' @param ncores Number of parallel workers. Default is
#'   `max(1, parallel::detectCores() - 1)`. Set to `1` to run sequentially.
#' @param ... Additional arguments passed to [prcoxreg()] for each bootstrap
#'   replicate (e.g., `alpha`, `pout.max`).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`betaHat`}{Matrix (B x p) of bootstrap `betaHat_re` estimates.}
#'     \item{`res`}{List of length B containing the full [prcoxreg()] result
#'       for each bootstrap replicate.}
#'   }
#'
#' @seealso [prcoxreg()], [getCI()]
#'
#' @export
SE.boot <- function(orig_y, orig_X, B = 100, seed = 1234,
                    ncores = max(1, parallel::detectCores() - 1),
                    ...) {
  set.seed(seed)
  n <- nrow(orig_X)
  p <- ncol(orig_X)
  boot.args <- list(...)

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  res <- foreach::foreach(
    ind = replicate(B, sample(n, n, replace = TRUE), simplify = FALSE),
    .packages = c("PAWPH2", "gbm", "dplyr", "survival")
  ) %dopar% {
    args <- boot.args
    if (!is.null(args$penalty.factor)) {
      args$penalty.factor <- args$penalty.factor[c(1:p, p + ind)]
    }
    args$X    <- orig_X[ind, ]
    args$y    <- orig_y[ind, ]
    args$seed <- seed
    do.call("prcoxreg", args)
  }

  betaHat <- t(sapply(res, getElement, "betaHat_re"))
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  list(betaHat = betaHat, res = res)
}


#' Confidence intervals from bootstrap samples
#'
#' Computes normal-based and quantile-based confidence intervals, a z-statistic,
#' and a two-sided p-value from a vector of bootstrap coefficient estimates.
#'
#' @param betas Numeric vector of length B: bootstrap estimates of a single
#'   coefficient (one column of the `betaHat` matrix returned by [SE.boot()]).
#' @param b Numeric scalar: the point estimate around which CIs are centered
#'   (typically `betaHat_re[j]` from [prcoxreg()]).
#' @param alpha Significance level. Default `0.05` gives 95% CIs.
#'
#' @return A named numeric vector with elements:
#'   \describe{
#'     \item{`SE`}{Bootstrap standard error.}
#'     \item{`lower_n`}{Normal-based lower CI bound.}
#'     \item{`upper_n`}{Normal-based upper CI bound.}
#'     \item{`lower_q`}{Quantile-based lower CI bound.}
#'     \item{`upper_q`}{Quantile-based upper CI bound.}
#'     \item{`z.stat`}{z-statistic (`b / SE`).}
#'     \item{`pvalue`}{Two-sided p-value under the standard normal.}
#'   }
#'
#' @seealso [SE.boot()], [prcoxreg()]
#'
#' @export
getCI <- function(betas, b, alpha = 0.05) {
  se <- stats::sd(betas)
  c(
    SE      = se,
    lower_n = b - stats::qnorm(1 - alpha / 2) * se,
    upper_n = b + stats::qnorm(1 - alpha / 2) * se,
    lower_q = stats::quantile(betas, probs = alpha / 2),
    upper_q = stats::quantile(betas, probs = 1 - alpha / 2),
    z.stat  = b / se,
    pvalue  = stats::pnorm(abs(b / se), lower.tail = FALSE) * 2
  )
}

#' Draw one sample from a truncated multivariate normal (TMVN) distribution
#' using sequential nearest neighbor (SNN) method
#'
#' @param cens_lb lower bound vector for TMVN of length n
#' @param cens_ub upper bound vector for TMVN of length n
#' @param m positive integer for the number of nearest neighbors used
#' @param covmat n-by-n dense covariance matrix, either `covmat` or `locs`,
#' `cov_name`, and `cov_parms` need to be provided
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for the covariance function from the `GpGp` package
#' @param NN n X m matrix for nearest neighbors. i-th row is the nearest neighbor indices of y_i. `NN[i, 1]` should be `i`
#' @param ordering `0` for do not reorder, `1` for variance descending order, `2` for maximin ordering
#' @param seed set seed for reproducibility
#' @return a vector of length n representing the underlying GP responses
#' @export
#' @examples
#' library(nntmvn)
#' library(TruncatedNormal)
#' set.seed(123)
#' x <- matrix(seq(from = 0, to = 1, length.out = 51), ncol = 1)
#' cov_name <- "matern15_isotropic"
#' cov_parm <- c(1.0, 0.1, 0.001) #'' variance, range, nugget
#' cov_func <- getFromNamespace(cov_name, "GpGp")
#' covmat <- cov_func(cov_parm, x)
#' lb <- rep(-Inf, nrow(x))
#' ub <- rep(-1, nrow(x))
#' m <- 30
#' samp_SNN <- matrix(NA, 3, nrow(x))
#' for (i in 1:3) {
#'   samp_SNN[i, ] <- nntmvn::rtmvn(lb, ub, m = m, covmat = covmat, locs = x, ordering = 0)
#' }
#' samp_TN <- TruncatedNormal::rtmvnorm(3, rep(0, nrow(x)), covmat, lb, ub)
#' qqplot(samp_SNN, samp_TN, xlim = range(samp_SNN, samp_TN), ylim = range(samp_SNN, samp_TN))
#' abline(a = 0, b = 1, lty = "dashed", col = "red")
#'
rtmvn <- function(cens_lb, cens_ub, m = 30, covmat = NULL,
                  locs = NULL, cov_name = NULL, cov_parm = NULL, NN = NULL,
                  ordering = 0, seed = NULL) {
  if (any(cens_lb > cens_ub)) {
    stop("There exists entry in cens_lb greater than its counterpart in cens_ub\n")
  }
  mask_cens <- rep(TRUE, length(cens_lb))
  y <- rep(NA, length(cens_lb))
  ind_equal <- which(cens_lb == cens_ub)
  mask_cens[ind_equal] <- FALSE
  y[ind_equal] <- cens_lb[ind_equal]
  return(rptmvn(
    y, cens_lb, cens_ub, mask_cens, m, covmat,
    locs, cov_name, cov_parm, NN, ordering, seed
  ))
}

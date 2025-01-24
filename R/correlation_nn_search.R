#' Find ordered nearest neighbors based on correlation, assuming the
#' absolute value of the correlation is monotonically decreasing with distance.
#' Returns an n X (m + 1) matrix, each row indicating the m + 1 nearest
#' neighbors including itself.
#'
#' @param covmat the covariance matrix
#' @param m the number of nearest neighbors
#' @return an n X (m + 1) matrix
#' @examples
#' library(RANN)
#' library(nntmvn)
#' set.seed(123)
#' d <- 3
#' n <- 100
#' locs <- matrix(runif(d * n), n, d)
#' covparms <- c(2, 0.01, 0)
#' covmat <- GpGp::matern15_isotropic(covparms, locs)
#' m <- 10
#' NNarray_test <- RANN::nn2(locs, k = m + 1)[[1]]
#' NNarray <- nntmvn::corr_nn(covmat, m)
#' cat("Number of mismatch is", sum(NNarray != NNarray_test, na.rm = TRUE))
#'
#' @export
corr_nn <- function(covmat, m) {
  inv_sd <- 1 / sqrt(diag(covmat))
  cormat <- outer(inv_sd, inv_sd) * covmat
  n <- nrow(covmat)
  if (m > n - 1) {
    warning("m cannot be greater than n - 1 when calling find_nn_corr\n")
    m <- n - 1
  }
  NN <- find_nn_corr_internal(cormat, m) + 1
  if (any(NN[, 1] != 1:n)) {
    warning("there exists correlation 1 at machine precision
            when calling find_nn_corr\n")
  }
  return(NN)
}

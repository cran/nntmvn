#' Simulate the underlying GP responses for censored responses using
#'  nearest neighbors
#'
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom utils getFromNamespace
#' @importFrom RANN nn2
#' @import GpGp
#' @param y uncensored responses of length n, where n is the number of all responses
#' @param cens_lb lower bound vector for TMVN of length n
#' @param cens_ub upper bound vector for TMVN of length n
#' @param mask_cens mask for censored responses (also locations) of length n
#' @param m positive integer for the number of nearest neighbors used
#' @param covmat n-by-n dense covariance matrix, either `covmat` or `locs`,
#' `cov_name`, and `cov_parms` need to be provided
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for the covariance function from the `GpGp` package
#' @param NN n X m matrix for nearest neighbors. i-th row is the nearest neighbor indices of y_i. `NN[i, 1]` should be `i`
#' @param seed set seed for reproducibility
#' @return a vector of length n representing the underlying GP responses
#' @export
#' @examples
#' library(GpGp)
#' library(RANN)
#' library(nntmvn)
#' set.seed(123)
#' x <- matrix(seq(from = 0, to = 1, length.out = 51), ncol = 1)
#' cov_name <- "matern15_isotropic"
#' cov_parm <- c(1.0, 0.1, 0.001) #' variance, range, nugget
#' cov_func <- getFromNamespace(cov_name, "GpGp")
#' covmat <- cov_func(cov_parm, x)
#' y <- t(chol(covmat)) %*% rnorm(length(x))
#' mask <- y < 0.3
#' y_cens <- y
#' y_cens[mask] <- NA
#' lb <- rep(-Inf, 100)
#' ub <- rep(0.3, 100)
#' m <- 10
#' y_samp_mtd1 <- rtmvn_snn(y_cens, lb, ub, mask, m = m, locs = x, 
#'                          cov_name = cov_name, cov_parm = cov_parm, seed = 123)
#' y_samp_mtd2 <- rtmvn_snn(y_cens, lb, ub, mask, m = m, covmat = covmat, 
#'                          seed = 123)
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y[mask], col = "blue")
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y_samp_mtd1[mask], col = "red")
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y_samp_mtd2[mask], col = "brown")
#'
rtmvn_snn <- function(y, cens_lb, cens_ub, mask_cens, m = 30, covmat = NULL,
                      locs = NULL, cov_name = NULL, cov_parm = NULL, NN = NULL,
                      seed = NULL) {
  if (is.null(covmat)) {
    if (is.null(locs) || is.null(cov_name) || is.null(cov_parm)) {
      stop("locs, cov_name, cov_parm cannot be NULL when covmat is NULL\n")
    }
    covfunc <- getFromNamespace(cov_name, "GpGp")
  } else {
    covfunc <- NULL
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # find NN
  if (is.null(NN)) {
    if (!is.null(locs)) {
      NN <- RANN::nn2(locs, k = m + 1)[[1]]
    } else {
      NN <- corr_nn(covmat, m)
    }
  }

  ind_cens <- which(mask_cens)
  for (i in ind_cens) {
    NN_row <- NN[i, ]
    if (is.null(covfunc)) {
      covmat_sub <- covmat[NN_row, NN_row]
    } else {
      covmat_sub <- covfunc(cov_parm, locs[NN_row, , drop = FALSE])
    }
    mask_cens_sub <- mask_cens[NN_row]
    y_sub <- y[NN_row]
    cens_lb_sub <- cens_lb[NN_row]
    cens_ub_sub <- cens_ub[NN_row]
    n_cens_sub <- sum(mask_cens_sub)
    if (n_cens_sub == length(NN_row)) {
      cond_covmat_sub_cens <- covmat_sub
      cond_mean_sub_cens <- rep(0, n_cens_sub)
    } else {
      tmp_mat <- covmat_sub[mask_cens_sub, !mask_cens_sub] %*%
        solve(covmat_sub[!mask_cens_sub, !mask_cens_sub])
      cond_covmat_sub_cens <- covmat_sub[mask_cens_sub, mask_cens_sub] -
        tmp_mat %*% covmat_sub[!mask_cens_sub, mask_cens_sub]
      cond_covmat_sub_cens[lower.tri(cond_covmat_sub_cens)] <-
        t(cond_covmat_sub_cens)[lower.tri(cond_covmat_sub_cens)]
      cond_mean_sub_cens <- as.vector(tmp_mat %*% y_sub[!mask_cens_sub])
    }
    samp_cens_sub <- t(TruncatedNormal::rtmvnorm(
      1, cond_mean_sub_cens,
      cond_covmat_sub_cens,
      cens_lb_sub[mask_cens_sub], cens_ub_sub[mask_cens_sub]
    ))
    y[i] <- samp_cens_sub[1]
    mask_cens[i] <- FALSE
  }
  y
}

#' Simulate the underlying GP responses for censored responses using
#'  nearest neighbors
#'
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom utils getFromNamespace
#' @import GpGp
#' @param y uncensored responses of length n, where n is the number of all responses
#' @param cens_lb lower bound vector for TMVN of length n
#' @param cens_ub upper bound vector for TMVN of length n
#' @param mask_cens mask for censored responses (also locations) of length n
#' @param NN n X m matrix for nearest neighbors. i-th row is the nearest neighbor indices of y_i. `NN[i, 1]` should be `i`
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for the covariance function from the `GpGp` package
#' @param covmat (optional) n-by-n dense covariance matrix, not needed if `locs`, 
#' `cov_name`, and `cov_parms` are provided
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
#' cov_parm <- c(1.0, 0.1, 0.001) # variance, range, nugget
#' cov_func <- getFromNamespace(cov_name, "GpGp")
#' covmat <- cov_func(cov_parm, x)
#' y <- t(chol(covmat)) %*% rnorm(length(x))
#' mask <- y < 0.3
#' y_cens <- y
#' y_cens[mask] <- NA
#' lb <- rep(-Inf, 100)
#' ub <- rep(0.3, 100)
#' m <- 10
#' NN <- RANN::nn2(x, k = m + 1)[[1]]
#' y_samp <- rtmvn_snn(y_cens, lb, ub, mask, NN, x, cov_name, cov_parm)
#' 
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y[mask], col = "blue")
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y_samp[mask], col = "red")
#' 
rtmvn_snn <- function(y, cens_lb, cens_ub, mask_cens, NN, locs, cov_name,
                               cov_parm, covmat = NULL,
                               seed = NULL) {
  if (is.null(covmat)) {
    covfunc <- getFromNamespace(cov_name, "GpGp")
  } else {
    covfunc <- NULL
  }
  if (!is.null(seed)) {
    set.seed(seed)
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
      cond_covmat_sub_cens[lower.tri(cond_covmat_sub_cens)] = 
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

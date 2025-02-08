#' Draw one sample of the underlying GP responses for a partially censored Gaussian
#' process using sequential nearest neighbor (SNN) method
#'
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom utils getFromNamespace
#' @importFrom RANN nn2
#' @importFrom stats dnorm
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
#' @param ordering `0` for do not reorder, `1` for variance descending order, `2` for maximin ordering
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
#' y_samp_mtd1 <- rptmvn(y_cens, lb, ub, mask,
#'   m = m, locs = x,
#'   cov_name = cov_name, cov_parm = cov_parm, seed = 123
#' )
#' y_samp_mtd2 <- rptmvn(y_cens, lb, ub, mask,
#'   m = m, covmat = covmat,
#'   seed = 123
#' )
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y[mask], col = "blue")
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y_samp_mtd1[mask], col = "red")
#' plot(x, y_cens, ylim = range(y))
#' points(x[mask, ], y_samp_mtd2[mask], col = "brown")
#'
rptmvn <- function(y, cens_lb, cens_ub, mask_cens, m = 30, covmat = NULL,
                   locs = NULL, cov_name = NULL, cov_parm = NULL, NN = NULL,
                   ordering = 0, seed = NULL) {
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
  if (ordering == 0) {

  } else if (ordering == 1) {
    if (!is.null(covmat)) {
      scaled_ub <- cens_ub[mask_cens] / sqrt(diag(covmat)[mask_cens])
      scaled_lb <- cens_lb[mask_cens] / sqrt(diag(covmat)[mask_cens])
    } else {
      scaled_ub <- cens_ub[mask_cens]
      scaled_lb <- cens_lb[mask_cens]
    }
    order_new <- 1:length(y)
    marginal_pr <- exp(TruncatedNormal::lnNpr(scaled_lb, scaled_ub))
    if (any(marginal_pr < 1e-20)) {
      warning(paste(
        "Reordering based on marginal variance cannot be",
        "performed due to marginal probability being too small\n"
      ))
    } else {
      marginal_var <- 1 - (x_times_dnorm(scaled_ub) -
        x_times_dnorm(scaled_lb)) / marginal_pr -
        ((dnorm(scaled_ub) - dnorm(scaled_lb)) / marginal_pr)^2
      order_new[mask_cens] <- order_new[mask_cens][order(marginal_var,
        decreasing = TRUE
      )]
      y <- y[order_new]
      cens_lb <- cens_lb[order_new]
      cens_ub <- cens_ub[order_new]
      # mask_cens <- mask_cens[order_new]
      if (!is.null(locs)) {
        locs <- locs[order_new, , drop = FALSE]
      }
      if (!is.null(covmat)) {
        covmat <- covmat[order_new, order_new, drop = FALSE]
      }
      if (!is.null(NN)) {
        warning(paste("When ordering is", ordering, "the input NN is ignored\n"))
        NN <- NULL
      }
    }
  } else if (ordering == 2) {
    if (is.null(locs)) {
      stop("locs must be provided if ordering = 2 for maximin ordering\n")
    }
    order_new <- GpGp::order_maxmin(locs)
    y <- y[order_new]
    cens_lb <- cens_lb[order_new]
    cens_ub <- cens_ub[order_new]
    mask_cens <- mask_cens[order_new]
    if (!is.null(locs)) {
      locs <- locs[order_new, , drop = FALSE]
    }
    if (!is.null(covmat)) {
      covmat <- covmat[order_new, order_new, drop = FALSE]
    }
    if (!is.null(NN)) {
      warning(paste("When ordering is", ordering, "the input NN is ignored\n"))
      NN <- NULL
    }
  } else {
    stop("Undefined ordering. Allowed input for ordering is 0 or 1")
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
  if (exists("order_new")) {
    order_new_rev <- c(1:length(y))
    order_new_rev[order_new] <- c(1:length(y))
    return(y[order_new_rev])
  } else {
    return(y)
  }
}

x_times_dnorm <- function(x) {
  y <- x * dnorm(x)
  y[is.nan(y)] <- 0
  y
}

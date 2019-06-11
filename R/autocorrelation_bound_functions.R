# diagonal_autocov_approx_0 computes the intergral of y^hat_0(t,t) with respect to \mu(dt).
# Input: f_data = the functional data matrix with functions in columns
# Output: scalar value of the integral of y^hat_0(t,t) with respect to \mu(dt).
#
# roxygen comments:
#' Compute the diagonal covariance
#'
#' \code{diagonal_autocov_approx_0} Computes the diagonal covariance of the given functional data.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @return A numeric value; integral approximation of the diagonal covariance of the functional data.
diagonal_autocov_approx_0 <- function(f_data) {
  J <- NROW(f_data)
  gamma_hat_0 <- autocov_approx_h(f_data, 0)
  sum(diag(gamma_hat_0)) / J
}

# autocorrelation_coff_h computes the approximate functional autocorrelation coefficient
#   rho^hat_h at lag h, defined in (17)
# Input: f_data = the functional data matrix with functions in columns
#        lag = lag for which to compute the coefficient
# Output: scalar value of the approximate functional autocorrelation coefficient at lag h.
#
# roxygen comments:
#' Computes the approximate functional autocorrelation coefficient at a given lag.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#' @return numeric value; the approximate functional autocorrelation coefficient at lag h.
autocorrelation_coeff_h <- function(f_data, lag) {
  N <- NCOL(f_data)
  num <- sqrt(t_statistic_Q(f_data, lag))
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  coefficient <- num / denom
  coefficient
}

# B_h_bound returns an approximate asymptotic upper 1-alpha confidence bound for the functional
#   autocorrelation coefficient at lag h under the assumption that f_data forms a weak white
#   noise.
# Input: f_data = the functional data matrix with functions in columns
#        lag = the lag for which to ccmpute the bound
#        alpha = significance level of the bound
#        M = optional argument specifying the sampling size in the related Monte Carlo method
#        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
#                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# Output: scalar value of the 1-alpha confidence bound for the functional autocorrelation
#         coefficient at lag h under a weak white noise assumption.
#
# roxygen comments:
#' Compute weak white noise confidence bound for autocorrelation coefficient.
#'
#' \code{B_h_bound} Computes an approximate asymptotic upper 1-alpha confidence bound for the functional
#' autocorrelation coefficient at lag h under a weak white noise assumption.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#' @param alpha the significance level to be used in the hypothesis test
#' @param M Number of samples to take when applying a Monte-Carlo approximation
#' @param low_disc Boolean value indicating whether or not to use low-discrepancy sampling in the Monte
#' Carlo method. Note, low-discrepancy sampling will yield deterministic results.
#' @return numeric value; the 1-alpha confidence bound for the functional autocorrelation
#' coefficient at lag h under a weak white noise assumption.
B_h_bound <- function(f_data, lag, alpha=0.05, M=NULL, low_disc=FALSE) {
  N <- NCOL(f_data)
  quantile = Q_WS_quantile(f_data, lag, alpha=alpha, M=M, low_disc=low_disc)$quantile
  num <- sqrt(quantile)
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  bound <- num / denom
  bound
}

# B_h_bound returns an approximate asymptotic upper 1-alpha confidence bound for the functional
#   autocorrelation coefficient at lag h under the assumption that f_data forms a strong
#   white noise.
# Input: f_data = the functional data matrix with functions in columns
#        alpha = significance level of the bound
# Output: scalar value of the 1-alpha confidence bound for the functional autocorrelation
#         coefficient at lag h under a strong white noise assumption.
#
# roxygen comments:
#' Compute strong white noise confidence bound for autocorrelation coefficient.
#'
#' \code{B_iid_bound} Computes an approximate asymptotic upper 1-alpha confidence bound for the functional
#' autocorrelation coefficient at lag h under the assumption that f_data forms a strong white noise
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param alpha the significance level to be used in the hypothesis test
#' @return Numeric value; the 1-alpha confidence bound for the functional autocorrelation coefficient
#' at lag h under a strong white noise assumption.
#' @rdname B_iid_bound
B_iid_bound <- function(f_data, alpha=0.05) {
  N <- NCOL(f_data)
  quantile_iid = Q_WS_quantile_iid(f_data, alpha=alpha)$quantile
  num <- sqrt(quantile_iid)
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  bound <- num / denom
  bound
}

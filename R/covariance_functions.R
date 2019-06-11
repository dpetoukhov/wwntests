# autocov_approx_h computes the approximate autocovariance gamma^hat_h(t,s) for a given lag (h)
#   for every (t,s) in U_J X U_J.
# Input: f_data = the functional data matrix with observed functions in columns
#        lag = the fixed lag for which to compute gamma^hat_h(t,s)
# Output: a 2-D array encoding the values of gamma^hat_h(t,s) for every (t,s) in U_J X U_J.
#
# roxygen comments:
#' Compute the approximate autocovariance at specified lag
#'
#' \code{autocov_approx_h} Computes the approximate autocovariance for a given lag h of the functional
#' data
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#' @return A 2-dimensional array encoding the autocovariance matrix for a given lag h.
autocov_approx_h <- function(f_data, lag) {
  N = NCOL(f_data)
  c_f_data <- center(f_data)
  gamma_hat_sum <- 0
  for (i in 1:(N-lag)) {
    gamma_hat_sum <- gamma_hat_sum + c_f_data[,i] %o% c_f_data[,i+lag]
  }
  gamma_hat <- gamma_hat_sum / N
  gamma_hat
}

# covariance_i_j returns the approximate covariance c^hat_i_j(t,s,u,v), encoding all values of
#   (t,s,u,v) in U_J X U_J X U_J X U_J, i,j in 1:T. Both T and J are inferred from f_data.
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing the covariance for
# Output: a 4-D array encoding c^hat_i_j(t,s,u,v) evaluated at all (t,s,u,v) in U_JxU_JxU_JxU_J.
#
# roxygen comments:
#' Compute the approximate covariance tensor for lag windows defined by i,j
#'
#' \code{covariance_i_j} Computes the approximate covariance tensor of the functional data for lag
#' windows defined by i,j.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @return A 4-dimensional array, encoding the covariance tensor of the functional data for lag
#' windows defined by i,j.
covariance_i_j <- function(f_data, i, j) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)
  sum <- array(0, c(J, J, J, J))
  for (k in (1+max(i,j)):N) {
    sum <- sum + c_f_data[,k-i] %o% c_f_data[,k] %o% c_f_data[,k-j] %o% c_f_data[,k]
  }
  cov <- sum / N
  cov
}

# covariance_i_j_vec is a vectorized version of the function covariance_i_j
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing the covariance for
# Output: a 4-D array encoding c^hat_i_j(t,s,u,v) evaluated at all (t,s,u,v) in U_JxU_JxU_JxU_J.
#
# roxygen comments:
#' Compute the approximate covariance tensor for lag windows defined by i,j
#'
#' \code{covariance_i_j_vec} Computes the approximate covariance tensor of the functional data for lag
#' windows defined by i,j; a vectorized version of covariance_i_j.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @return A 4-dimensional array, encoding the covariance tensor of the functional data for lag
#' windows defined by i,j.
covariance_i_j_vec <- function(f_data, i, j) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)
  sum_parts <- as.list((1+max(i,j)):N)
  sum_parts <- lapply(sum_parts,
                      function(k) c_f_data[,k-i] %o% c_f_data[,k] %o%
                        c_f_data[,k-j] %o% c_f_data[,k])
  cov <- Reduce('+', sum_parts)
  cov / N
}

# diagonal_covariance_i_j returns the approximate covariance c^hat_i_i(t,s,t,s), encoding all
#   values of t,s in U_J X U_J, i in 1:T.
# Input: f_data = the functional data matrix with functions in columns
#        i = the index i in 1:T that we are computing the covariance diagonal for
# Output: a 2-D array encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in U_JxU_J.
#
# roxygen comments:
#' Compute the approximate diagonal covariance matrix for lag windows defined by i
#'
#' \code{diagonal_covariance_i} Computes the approximate diagonal covariance matrix of the functional
#' data for lag windows defined by i.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i the index in 1:T that we are computing the covariance for
#' @return A 2-dimensional array, encoding the covariance matrix of the functional data for lag
#' windows defined by i.
diagonal_covariance_i <- function(f_data, i) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)
  sum1 <- array(0, c(J, J))
  for (k in (1+i):N) {
    sum1 <- sum1 + ((c_f_data[,k-i])^2 %o% (c_f_data[,k])^2)
  }
  cov <- (1 / N) * sum1
  cov
}

# scalar_covariance_i_j returns the approximate covariance c^hat_i_j(t,s,u,v) evaluated at a
#   given t,s,u,v in U_J X U_J X U_J X U_J (for use in MCint method).
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing the covariance for
#        times = a 4-element vector representing the values (t,s,u,v)
# Output: scalar value of the computed covariance c^hat_i_j(t,s,u,v).
#
# roxygen comments:
#' Compute the approximate covariance at a point for lag windows defined by i,j
#'
#' \code{scalar_covariance_i_j} Computes the approximate covariance at a point of the functional data
#' for lag windows defined by i,j; a scalarized version of covariance_i_j that takes point estimates.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @param times A vector with 4 columns containing indices specifying which subset of f_data to consider
#' @return A numeric value; the covariance of the functional data at a point for lag
#' windows defined by i,j.
scalar_covariance_i_j <- function(f_data, i, j, times) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum1 <- 0
  for (k in (1+max(i,j)):N) {
    sum1 <- sum1 + c_f_data[times[1],k-i] * c_f_data[times[2],k] * c_f_data[times[3],k-j]  *
      c_f_data[times[4],k]
  }
  cov <- (1/N) * sum1
  cov
}

# scalar_covariance_i_j_vec is a vectorized version of the function scalar_covariance_i_j
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing the covariance for
#        times = a 4-element vector representing the values (t,s,u,v)
# Output: scalar value of the computed covariance c^hat_i_j(t,s,u,v).
#
# roxygen comments:
#' Compute the approximate covariance at a point for lag windows defined by i,j
#'
#' \code{scalar_covariance_i_j_vec} Computes the approximate covariance at a point of the functional data
#' for lag windows defined by i,j; a vectorized version of scalar_covariance_i_j.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @param times A vector with 4 columns containing indices specifying which subset of f_data to consider
#' @return A numeric value; the covariance of the functional data at a point for lag
#' windows defined by i,j.
scalar_covariance_i_j_vec <- function(f_data, i, j, times) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum_parts <- list((1+max(i,j)):N)
  sum_parts <- lapply(sum_parts,
                      function(k) c_f_data[times[1],k-i] * c_f_data[times[2],k] *
                        c_f_data[times[3],k-j]  * c_f_data[times[4],k])
  cov <- (1/N) * Reduce('+', sum_parts)
  cov
}

# iid_covariance returns one of the two independent sum terms in the approximate covariance
#   c^*_0(t,s,u,v) definition, encoding all values of (t,s) in U_J X U_J, i,j in 1:T.
# Input: f_data = the functional data matrix with functions in columns
# Output: returns a 2-D tensor of c^*(t,s), one of the two independent sums in the computation
#         of c^*(t,s,u,v).
#
# roxygen comments:
#' Compute part of the covariance under a strong white noise assumption
#'
#' \code{iid_covariance} A helper function used to compute one of the two independent sum terms in the
#' computation of the approximate covariance of the functional data under a strong white noise assumption.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @return A 2-dimensional matrix containing one of the two independent sums in the computation of the
#' covariance.
iid_covariance <- function(f_data) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum1 <- 0
  for (i in 1:N) {
    sum1 <- sum1 + c_f_data[,i] %o% c_f_data[,i]
  }
  sum1 / N
}

#
# roxygen comments:
#' Compute part of the covariance under a strong white noise assumption
#'
#' \code{iid_covariance_vec} A helper function used to compute one of the two independent sum terms in the
#' computation of the approximate covariance of the functional data under a strong white noise assumption;
#' a vectorized version of iid_covariance.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @return A 2-dimensional matrix containing one of the two independent sums in the computation of the
#' covariance.
iid_covariance_vec <- function(f_data) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum_parts <- 1:N
  sum_parts <- lapply(sum_parts,
                      function(i) c_f_data[,i] %o% c_f_data[,i])
  cov <- (1 / N) * Reduce('+', sum_parts)
  cov
}

# covariance_diag_store returns a list storage of the approximate covariances c^hat_i_i(t,s,t,s),
#   for all i in 1:K, for each encoding all values of t,s in U_J X U_J.
# Input: f_data = the functional data matrix with functions in columns
#        K = the maximum value of i in the range 1:K for which to compute c^hat_i_i(t,s,t,s)
# Output: a list containing K 2-D arrays encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in
#         U_JxU_J, for i in 1:K
#
# roxygen comments:
#' List storage of diagonal covariances.
#'
#' \code{covariance_diag_store} Creates a list storage of approximate diagonal covariances computed
#' by the function diagonal_covariance_i
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param K the range of lags 1:K to use
#' @return A list containing K 2-dimensional arrays containing the diagonal covariance matrices of the
#' functional data, for lags h in the range 1:K.
covariance_diag_store <- function(f_data, K) {
  cov_i_store <- list()
  for (j in 1:K) {
    cov_i_store[[j]] <- diagonal_covariance_i(f_data, j)
  }
  cov_i_store
}

# true_eta_approx is a non-stochaistic approximation of eta_i_j (see (15)) using a Riemann sum.
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
# Output: scalar value of eta^hat_i_j computed using a simple Riemann sum.
true_eta_approx_i_j <- function(f_data, i, j) {
  J <- NROW(f_data)
  cov_tensor <- covariance_i_j(f_data, i, j)
  2 * sum(cov_tensor^2) / (J^4)
}

# MCint_eta_approx_i_j computes an approximation of eta_i_j (defined under (15)) using the second
#   Monte Carlo integration method "MCint" defined on page 8.
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
#        M = number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
#        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
#                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# Output: scalar value of eta^_hat_i_j computed using the MCint method.
MCint_eta_approx_i_j <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
  J <- NROW(f_data)
  T <- NCOL(f_data)
  if (is.null(M)) {
    M = floor((max(150 - T, 0) + max(100-J,0) + (J / sqrt(2))))
  }
  if (low_disc == TRUE) {
    if (requireNamespace('fOptions', quietly = TRUE)) {
      rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
      rand_samp_mat[which(rand_samp_mat == 0)] = 1
    } else {
      stop("Please install the 'fOptions' package for low discrepancy sampling.")
    }
  } else {
    rand_samp_mat <- matrix(nrow=M, ncol=4)
    for (k in 1:4) {
      rand_samp <- floor(J * runif(M, 0, 1))
      rand_samp[which(rand_samp == 0)] = 1
      rand_samp_mat[,k] <- rand_samp
    }
  }
  eta_hat_i_j_sum <- 0
  for (k in 1:M) {
    cov <- scalar_covariance_i_j(f_data, i, j, rand_samp_mat[k,])
    eta_hat_i_j_sum <- eta_hat_i_j_sum + (cov^2)
  }
  eta_hat_i_j <- (2/M) * eta_hat_i_j_sum
  eta_hat_i_j
}

# MCint_eta_approx_i_j_vec is a vectorized version of MCint_eta_approx_i_j.
# Input: f_data = the functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
#        M = number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
#        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
#                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# Output: scalar value of eta^_hat_i_j computed using the MCint method.
MCint_eta_approx_i_j_vec <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  M = floor((max(150 - N, 0) + max(100-J,0) + (J / sqrt(2))))
  if (low_disc == TRUE) {
    if (requireNamespace('fOptions', quietly = TRUE)) {
      rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
      rand_samp_mat[which(rand_samp_mat == 0)] = 1
    } else {
      stop("Please install the 'fOptions' package for low discrepancy sampling.")
    }
  } else {
    rand_samp_mat <- matrix(nrow=M, ncol=4)
    for (k in 1:4) {
      rand_samp <- floor(J * runif(M, 0, 1))
      rand_samp[which(rand_samp == 0)] = 1
      rand_samp_mat[,k] <- rand_samp
    }
  }
  eta_parts <- as.list(1:M)
  eta_parts <- lapply(eta_parts, function(k) scalar_covariance_i_j(f_data, i, j,
                                                          rand_samp_mat[k,]) ^ 2)
  eta_hat_i_j <- (2 / M) * Reduce('+', eta_parts)
  eta_hat_i_j
}


# mean_hat_V_K computes the approximation of the mean defined in (15) which is used in the Welch-
#   Satterthwaite approximation as mean of the chi-squared random variable approximating V_K.
# Input: f_data = the functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for for the test statistic V_K
# Output: scalar approximation of the mean of the test statistic V_K.
mean_hat_V_K <- function(f_data, K) {
  J <- NROW(f_data)
  sum1 <- 0
  store <- covariance_diag_store(f_data, K)
  for (i in 1:K) {
    sum1 <- sum1 + sum(store[[i]])
  }
  mu_hat_V_K <- (1 / (J^2)) * sum1
  mu_hat_V_K
}

# mean_hat_V_K_iid computes the approximation of the mean defined in (15) which is used in the
#   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#   strong white noise.
# Input: f_data = the functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for for the test statistic V_K
# Output: scalar approximation of the mean of the test statistic V_K under a strong white noise
#         assumption.
mean_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)
  mu_hat_Q_h <- K * ((sum(diag(cov)) / J)^2)
  mu_hat_Q_h
}

# mean_hat_Q_h computes the approximation of the mean defined in (15) which is used in the Welch-
#   Satterthwaite approximation as mean of the chi-squared random variable approximating Q_h.
# Input: f_data = the functional data matrix with functions in columns
#        lag = specifies the lag use in the test statistic Q_h (lag = h in paper)
# Output: scalar approximation of the mean of the test statistic Q_h.
mean_hat_Q_h <- function(f_data, lag) {
  J <- NROW(f_data)
  cov <- diagonal_covariance_i(f_data, lag)
  mu_hat_Q_h <- (1 / (J^2)) * sum(cov)
  mu_hat_Q_h
}

# mean_hat_Q_h_iid computes the approximation of the mean defined in (15) which is used in the
#   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#   strong white noise.
# Input: f_data = the functional data matrix with functions in columns
# Output: scalar approximation of the mean of the test statistic Q_h under a strong white noise
#         assumption.
mean_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)
  mu_hat_Q_h_iid <- (sum(diag(cov)) / J)^2
  mu_hat_Q_h_iid
}

# variance_hat_V_K computes the approximation of the variance defined in (15) which is used in
#   the Welch- Satterthwaite approximation as the variance of the chi-squared random variable
#   approximating V_K.
# Input: f_data = the functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
#        M = optional argument specifying the sampling size in the related Monte Carlo method
#        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
#                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# Output: scalar approximation of the variance of the test statistic V_K.
variance_hat_V_K <- function(f_data, K, M=NULL, low_disc=FALSE) {
  N <- NCOL(f_data)
  sum1 <- 0
  for (i in 1:K) {
    sum1 <- sum1 + MCint_eta_approx_i_j(f_data, i, i, M=M, low_disc=low_disc)
  }
  bandwidth <- ceiling(0.25 * (N ^ (1/3)))
  if (K > 1) {
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        if (abs(i-j) > bandwidth) { # empirically, past a lag of 15, error is less than 1%
          next
        }
        sum1 <- sum1 + (2 * MCint_eta_approx_i_j(f_data, i, j, M=M, low_disc=low_disc))
      }
    }
  }
  variance_V_K <- sum1
  variance_V_K
}


# variance_hat_V_K_iid computes the approximation of the variance defined in (15) which is used
#   in the Welch- Satterthwaite approximation under the assumption that the functional data
#   follows a strong white noise.
# Input: f_data = the functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
# Output: scalar approximation of the variance of the test statistic V_K
variance_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)
  variance_V_K_iid <- K * 2 * ( sum(cov_iid^2) / (J^2) )^2
  variance_V_K_iid
}

# variance_hat_Q_h computes the approximation of the variance defined in (15) which is used in
#   the Welch- Satterthwaite approximation as variance of the chi-squared random variable
#   approximating Q_h.
# Input: f_data = the functional data matrix with functions in columns
#        lag = specifies the lag use in the test statistic Q_h (lag = h in paper)
#        M = optional argument specifying the sampling size in the related Monte Carlo method
#        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
#                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# Output: scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h <- function(f_data, lag, M=NULL, low_disc=FALSE) {
  variance_Q_h <- MCint_eta_approx_i_j(f_data, lag, lag, M=M, low_disc=low_disc)
  variance_Q_h
}

# variance_hat_Q_h_iid computes the approximation of the variance defined in (15) which is used
#   in the Welch- Satterthwaite approximation under the assumption that the functional data
#   follows a strong white noise.
# Input: f_data = the functional data matrix with functions in columns
#        lag = specifies the lag use in the test statistic Q_h (lag = h in paper)
# Output: scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)
  variance_Q_h_iid <-2 * ( sum(cov_iid^2) / (J^2) )^2
  variance_Q_h_iid
}

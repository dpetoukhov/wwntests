# t_statistic_Q computes the test statistic Q_{T,h} = T*||y^hat_h||^2 for fixed h and for T
#   inferred from the functional data f_data that is passed.
# Input: f_data = the functional data matrix with observed functions in columns
#        lag = the fixed time lag used in the computation of the statistic
# Output: scalar value of the statistic Q_{T,h} to test the hypothesis H_{0,h} : y_h(t,s) = 0.
t_statistic_Q <- function(f_data, lag) {
  N = NCOL(f_data)
  J = NROW(f_data)
  gamma_hat <- autocov_approx_h(f_data, lag)
  Q_T_h <- N * sum(gamma_hat^2) / (J^2)
  Q_T_h
}

# t_statistic_V computes the statistic V_{T,K} = T*sum_h(||y^hat_h||^2) or h in 1:K and for T
#   inferred from the functional data f_data that is passed to the function.
# Input: f_data = the functional data with functions in columns
#        K = the max value in the range of time lags (1:K) used
# Output: scalar value of the statistic V_{T,K} to test the hypothesis
#         H'_{0,K} : for all h in 1:K y_h(t,s) = 0.
t_statistic_V <- function(f_data, K) {
  V_T_K <- 0
  for (h in 1:K) {
    V_T_K <- V_T_K + t_statistic_Q(f_data, h)
  }
  V_T_K
}

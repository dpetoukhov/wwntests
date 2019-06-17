# imhof_test returns the the SVD of the tensor c^hat_i_j(t,s,u,v) and the p-value computing
#   the probability that the observed value of the statistic Q_h is larger than the 1-alpha
#   quantile of the quadratic form in normal variables described in (15)
# Input: f_data = the functional data matrix with functions in columns
#        lag = the lag for which to compute the imhof test
# Output: a list containing the SVD of tensor c^hat_i_j(t,s,u,v) and the p-value computing the
#   probability that the observed value of the statistic Q_h is larger than the 1-alpha quantile
#   of the quadratic form in normal variables.
imhof_test <- function(f_data, lag) {
  if (!requireNamespace('tensorA')) {
    stop("Please install the 'tensorA' package to perform the imhof test.")
  }
  if (!requireNamespace('CompQuadForm')) {
    stop("Please install the 'CompQuadForm' package to perform the imhof test.")
  }
  if ((lag < 1) | (lag %% 1 != 0)) {
    stop("The 'lag' parameter must a positive integer.")
  }
  N = NCOL(f_data)
  J = NROW(f_data)
  t_statistic_val = t_statistic_Q(f_data, lag)
  c_f_data <- center(f_data)
  tensor <- array(0, c(J, J, J, J))
  sum1 <- 0
  for (k in 1:(N-lag)) {
    tensor <- tensor + c_f_data[,k] %o% c_f_data[,k+lag] %o% c_f_data[,k] %o% c_f_data[,k+lag]
  }
  tensor <- tensor / N
  temp_tensor <- as.numeric(tensor)
  tensor_numeric <- tensorA::to.tensor(temp_tensor, c(J,J,J,J))
  names(tensor_numeric) = c("a", "b", "c", "d")
  SVD <- tensorA::svd.tensor(tensor_numeric, i=c("a", "b"))
  eigenvalues <- as.numeric(SVD$d / (J^2))
  pval_imhof <- CompQuadForm::imhof(t_statistic_val, lambda = eigenvalues)$Qq
  list(statistic = t_statistic_val, p_value = pval_imhof)
}

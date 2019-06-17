# spectral_t_statistic computes the spectral density operator based test statistic of the functional
#   data f_data.
# Input: f_data = The functional data matrix with observed functions in columns
#        kernel = The kernel function to use. The currently supported kernels are 'Bartlett' and 'Parzen'.
#                 The default kernel is 'Bartlett'.
#        bandwidth = specifies the bandwidth to use. Currently admitted arguments are positive
#                    integers, 'static' which computes the bandwith p via p = n^(1/(2q+1)) where
#                    n is the sample size and q is the kernel order, or 'adaptive' which uses a
#                    bandwith selection method that is based on the functional data.
# Output: scalar value of the spectral density based test statistic.
spectral_t_statistic <- function(f_data, kernel = 'Bartlett', bandwidth = 'adaptive') {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  kernel_string <- kernel
  if (kernel == 'Bartlett') {
    kernel <- bartlett_kernel
    kernel_order <- 1
  } else if (kernel == 'Parzen') {
    kernel <- parzen_kernel
    kernel_order <- 2
  #} else if (kernel == 'Daniell') {
    #kernel <- daniell_kernel
    #kernel_order <- 2
  } else {
    stop("This kernel is not supported. Please see the documentation for supported kernel functions.")
  }
  if (bandwidth == 'static') {
    bandwidth <- N^(1 / (2 * kernel_order + 1))
  } else if (bandwidth == 'adaptive') {
    bandwidth <- max(2, adaptive_bandwidth(f_data, kernel_string))
  } else if (!is.numeric(bandwidth)) {
    stop("Please see the documentation for valid bandwith arguments.")
  }
  data_inner_prod <- crossprod(f_data) / J
  C_hat_HS_norm <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS_norm[j+1] <- N^(-2) * sum(data_inner_prod[(j+1):N, (j+1):N] * data_inner_prod[1:(N-j), 1:(N-j)])
  }
  kernel_vals <- sapply(1:(N-1) / bandwidth, kernel)
  spectral_distance_Q_sq <- 2 * sum((kernel_vals^2) * C_hat_HS_norm[-1])
  C_n_k <- sum( (1 - 1:(N-1)/N) * kernel_vals^2 )
  D_n_k <- sum( (1 - 1:(N-2)/N) * (1 - 2:(N-1)/N) * kernel_vals[-(N-1)]^4 )
  sigma_squared_hat <- sum(diag(data_inner_prod)) / N

  t_stat_term <- sigma_squared_hat^(-2) * C_hat_HS_norm[1] * sqrt(2 * D_n_k) # for convenience
  untrans_num <- 2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq

  ### TODO: add case for when H is not R. This is denoted by const in original codebase
  beta <- 1 - (2/3) * sum(kernel_vals^2) * sum(kernel_vals^6) / (sum(kernel_vals^4)^2)
  t_stat <- ((2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq)^beta -
    (C_n_k^beta + 2^(-1) * beta * (beta - 1) * C_n_k^(beta-2) * t_stat_term^2)) / (beta * C_n_k^(beta-1) * t_stat_term)
  list(stat = t_stat, band = bandwidth)
}


# adaptive_bandwidth computes the "optimal" bandwidth using a bandwidth selection method based on the
#   spectral density operator which adapts to the functional data.
# Input: f_data = the functional data matrix with observed functions in columns
#        kernel = the kernel function to use. The currently supported kernels are 'Bartlett' and 'Parzen'.
#                 The default kernel is 'Bartlett'.
# Output: a scalar value of the "optimal" data-adapted bandwidth.
adaptive_bandwidth <- function(f_data, kernel) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  if (kernel == 'Bartlett') {
    kernel <- bartlett_kernel
    order <- 1
    xi <- 1
    kern_int <- 2 / 3
  } else if (kernel == 'Parzen') {
    kernel <- parzen_kernel
    order <- 2
    xi <- 6
    kern_int <- 151 / 280
  #} else if (kernel == 'Daniell') {
    #kernel <- daniell_kernel
    #order <- 2
    #xi <- (pi^2) / 6
    #kern_int <- 1
  } else {
    stop('Please see the documentation for supported kernels.')
  }
  data_inner_prod <- crossprod(f_data) / (N * J)
  C_hat_HS <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS[j+1]<-sum(data_inner_prod[(j+1):N,(j+1):N] * data_inner_prod[1:(N-j),1:(N-j)])
  }
  initial_band_q <- 4 * N^(1 / (2*order +1))
  k_n_j_q <- kernel(1:(N-1) / initial_band_q)
  initial_band_0 <- initial_band_q / 4
  k_n_j_0 <- kernel(1:(N-1) / initial_band_0)
  Q_hat_sq <- 2 * sum(k_n_j_0^2 * C_hat_HS[-1])
  Q_hat_sq <- Q_hat_sq + sum(C_hat_HS[0])
  Term2 <- Q_hat_sq
  Term1 <- 2 * sum(k_n_j_q^2 * ((1:(N-1))^(2*order)) * C_hat_HS[-1])
  C_hat_TR <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_TR[j+1] <- sum(data_inner_prod[1:(N-j), (j+1):N])
  }
  trace <- 2 * sum(k_n_j_0^2 * C_hat_TR[-1])
  Term3 <- trace + sum(C_hat_TR[1])
  band_constant <- (2 * order * xi^2 * Term1 / (kern_int * (Term2 + Term3)))^(1/(2*order + 1))
  band_constant * N^(1 / (2*order + 1))
}


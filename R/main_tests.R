#' Single-Lag Hypothesis Test
#'
#' \code{single_lag_test} Computes the single-lag hypothesis test at a single user-specified lag.
#'
#' @param f_data The functional data matrix with observed functions in the columns
#' @param lag Positive integer value. The lag to use to compute the single lag test statistic.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param iid A Boolean value, FALSE by default. If given TRUE, the hypothesis test will use a strong-white
#' noise assumption (instead of a weak-white noise assumption).
#' @param M Positive integer value. Number of Monte-Carlo simulations for the Welch-Satterthwaite approximation.
#' @param low_disc A Boolean value, FALSE by default. If given TRUE, uses low-discrepancy sampling in the
#' Monte-Carlo method. Note, low-discrepancy sampling will yield deterministic results.
#' Requires the 'fOptions' package.
#' @param bootstrap A Boolean value, FALSE by default If given TRUE, the hypothesis test is done by
#' approximating the limiting distribution of the test statistic via a block bootstrap process.
#' @param block_size A positive Integer value, with the default value being computed via the adaptive
#' bandwidth selection method in the "spectral" test. Determines the block size (of each block in each
#' bootstrap sample) if the test is being bootstrapped.
#' @param straps A positive Integer, with a default value of 300. Determines the number of bootstrap samples
#' to take if the test is being bootstrapped. Only used if 'bootstrap' == TRUE.
#' @param moving A Boolean value, FALSE by default If given TRUE, the performed block bootstrap will be moving
#' rather than stationary.
#' @param suppress_raw_output Boolean value, FALSE by default. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @details The "single-lag" portmanteau test is based on the sample autocovariance function computed from the
#' functional data. This test assesses the significance of lagged autocovariance operators at a single,
#' user-specified lag h. More specifically, it tests the null hypothesis that the lag-h autocovariance
#' operator is equal to 0. This test is designed for stationary functional time-series, and is valid under
#' conditional heteroscedasticity conditions.
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#'
#' @references
#' [1] Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' f <- far_1_S(150, 50, S = 0.75)
#' single_lag_test(f, lag = 1)
#' single_lag_test(f, lag = 2, M=100)
#'
#' @import stats
#' @export
single_lag_test <- function(f_data, lag=1, alpha=0.05, iid=FALSE,
                          M=NULL, low_disc=FALSE, bootstrap=FALSE,
                          block_size='adaptive', straps=300, moving = FALSE,
                          suppress_raw_output=FALSE, suppress_print_output=FALSE) {
  if (bootstrap == TRUE & (iid == TRUE | low_disc == TRUE)) {
    stop("Bootstrapping this test only requires the lag parameter
         (and optionally, a significance level).")
  }
  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. Atleast one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }
  if (bootstrap == TRUE) {
    results <- Q_WS_hyp_test(f_data, lag, alpha = alpha, bootstrap = TRUE, block_size = block_size,
                             moving = moving, straps = straps)
    if (suppress_print_output == FALSE) {
      if (moving == TRUE) {
        title_print <- sprintf(" Moving Block Bootstrapped Single-Lag Test\n\n")
      } else if (moving == FALSE) {
        title_print <- sprintf("Block Bootstrapped Single-Lag Test\n\n")
      }
      test_type <- 'the series is a weak white noise\n'
      null_print <- sprintf("null hypothesis: %s", test_type)
      p_val_print <- sprintf("p-value = %f\n", results$p_value)
      samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      lag_print <- sprintf("lag = %d\n", lag)
      boot_num <- sprintf("number of bootstrap samples = %d\n", straps)
      block_sze <- sprintf("block size = %d\n\n\n", results$block_size)
      message(c(title_print, null_print, p_val_print, samp_print,
            lag_print, boot_num, block_sze))
    }
    if (suppress_raw_output == FALSE) {
      results[-4]
    }
  } else if (iid  == FALSE) {
    results <- Q_WS_hyp_test(f_data, lag, alpha=alpha, M=M, low_disc=low_disc)
    if (suppress_print_output == FALSE) {
      title_print <- sprintf("Single-Lag Test\n\n")
      null_print <- sprintf("null hypothesis: the series is uncorrelated at lag %d\n", lag)
      p_val_print <- sprintf("p-value = %f\n", results$p_value)
      samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      lag_print <- sprintf("lag = %d\n\n\n", lag)
      message(c(title_print, null_print, p_val_print, samp_print,
            lag_print))
    }
    if (suppress_raw_output == FALSE) {
      results
    }
  } else if (iid == TRUE) {
    results <- Q_WS_hyp_test(f_data, iid = TRUE, lag = lag, alpha=alpha)
    if (suppress_print_output == FALSE) {
      title_print <- sprintf("Single-Lag Test (iid assumption)\n\n")
      test_type <- 'the series is a strong white noise\n'
      null_print <- sprintf("null hypothesis: %s", test_type)
      p_val_print <- sprintf("p-value = %f\n", results$p_value)
      samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      lag_print <- sprintf("lag = %d\n\n\n", lag)
      message(c(title_print, null_print, p_val_print, samp_print,
            lag_print))
    }
    if (suppress_raw_output == FALSE) {
      results
    }
  }
}


#' Multi-Lag Hypothesis Test
#'
#' \code{multi_lag_test} Computes the multi-lag hypothesis test over a range of user-specified lags.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag Positive integer value. The lag to use to compute the single lag test statistic
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param iid A Boolean value, FALSE by default. If given TRUE, the hypothesis test will use a strong-white
#' noise assumption (instead of a weak-white noise assumption).
#' @param M Positive integer value. Number of Monte-Carlo simulation for Welch-Satterthwaite approximation.
#' @param low_disc A Boolean value, FALSE by default. If given TRUE, uses low-discrepancy sampling in the
#' Monte-Carlo method. Note, low-discrepancy sampling will yield deterministic results.
#' Requires the 'fOptions' package.
#' @param suppress_raw_output Boolean value, FALSE by default. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @details The "multi-lag" portmanteau test is also based on the sample autocovariance function computed from the
#' functional data. This test assesses the cumulative significance of lagged autocovariance operators, up to a
#' user-selected maximum lag K. More specifically, it tests the null hypothesis that the first K lag-h autocovariance
#' operators (h going from 1 to K) is equal to 0. This test is designed for stationary functional time-series, and
#' is valid under conditional heteroscedasticity conditions.
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#'
#' @references
#' [1] Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' b <- brown_motion(150, 50)
#' multi_lag_test(b, lag = 5)
#' multi_lag_test(b, lag = 10, M = 50)
#'
#' @import stats
#' @export
multi_lag_test <- function(f_data, lag = 20, M=NULL, low_disc=FALSE, iid=FALSE,
                           alpha=0.05, suppress_raw_output=FALSE,
                           suppress_print_output=FALSE) {
  K <- lag
  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. Atleast one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }
  if (iid == FALSE) {
    results <- V_WS_quantile(f_data, K, alpha=alpha, M=M, low_disc=low_disc)
    if (suppress_print_output == FALSE) {
      title_print <- sprintf("Multi-Lag Test\n\n")
      test_type <- 'the series is a weak white noise\n'
      null_print <- sprintf("null hypothesis: %s", test_type)
      p_val_print <- sprintf("p-value = %f\n", results$p_value)
      samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      lag_print <- sprintf("maximum lag = %d\n", K)
      mc_print <- sprintf("number of monte-carlo simulations = %d\n\n\n", M)
      message(c(title_print, null_print, p_val_print, samp_print,
            lag_print, mc_print))
    }
    if (suppress_raw_output == FALSE) {
      results
    }
  } else {
    results <- V_WS_quantile_iid(f_data, K, alpha=alpha)
    if (suppress_print_output == FALSE) {
      title_print <- sprintf("Multi-Lag Test (iid assumption)\n\n")
      test_type <- 'the series is a strong white noise\n'
      null_print <- sprintf("null hypothesis: %s", test_type)
      p_val_print <- sprintf("p-value = %f\n", results$p_value)
      samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      lag_print <- sprintf("maximum lag = %d\n\n\n", K)
      message(c(title_print, null_print, p_val_print, samp_print,
            lag_print))
    }
    if (suppress_raw_output == FALSE) {
      results
    }
  }
}


#' Spectral Density Test
#'
#' \code{spectral_test} Computes the spectral hypothesis test under a user-specified kernel function and
#' bandwidth; automatic bandwidth selection methods are provided.
#'
#' @param f_data The functional data matrix with observed functions in the columns
#' @param kernel A String specifying the kernel function to use. The currently supported kernels are the
#' 'Bartlett' and  'Parzen' kernels. The default kernel is 'Bartlett'.
#' @param bandwidth A String or positive Integer value which specifies the bandwidth to use. Currently admitted
#' string handles are 'static' which computes the bandwidth p via p = n^(1/(2q+1)) where n is the sample size
#' and q is the kernel order, or 'adaptive' which uses a bandwidth selection method that is based on the
#' functional data.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used for the test.
#' The significance level is 0.05 by default. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param suppress_raw_output Boolean value, FALSE by default. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @description The "spectral" portmanteau test is based on the spectral density operator. It essentially measures
#' the proximity of a functional time series to a white noise - the constant spectral density operator of an
#' uncorrelated series. Unlike the "single-lag" and "multi-lag" tests, this test is not for general white noise
#' series, and may not hold under functional conditionally heteroscedastic assumptions.
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#'
#' @references
#' [1] Characiejus V., & Rice G. (2019). A general white noise test based on kernel lag-window estimates of the
#' spectral density operator. Econometrics and Statistics, submitted.
#'
#' [2] Chen W.W. & Deo R.S. (2004). Power transformations to induce normality and their applications.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66, 117–130.
#'
#' @examples
#' b <- brown_motion(100, 50)
#' spectral_test(b)
#' spectral_test(b, kernel = 'Parzen', bandwidth = 'adaptive')
#' spectral_test(b, kernel = 'Bartlett', bandwidth = 2)
#'
#' @export
spectral_test <- function(f_data, kernel = 'Bartlett', bandwidth = 'adaptive', alpha = 0.05,
                          suppress_raw_output=FALSE, suppress_print_output=FALSE) {
  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. Atleast one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }
  quantile <- qnorm(1 - alpha)
  statistic <- spectral_t_statistic(f_data, kernel = kernel, bandwidth = bandwidth)
  band <- statistic$band
  statistic <- statistic$stat
  p_val <- 1 - pnorm(statistic)
  results <- list(statistic = statistic, quantile = quantile, p_value = p_val, band = band)
  if (suppress_print_output == FALSE) {
    title_print <- sprintf("Spectral Test\n\n")
    test_type <- 'the series is iid\n'
    null_print <- sprintf("null hypothesis: %s", test_type)
    p_val_print <- sprintf("p-value = %f\n", results$p_value)
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    kern_print <- sprintf("kernel function = %s\n", kernel)
    band_print <- sprintf("bandwidth = %f\n", results$band)
    if (is.numeric(bandwidth)) {
      band_sel <- sprintf("bandwidth selection = %d\n\n\n", bandwidth)
    } else {
      band_sel <- sprintf("bandwidth selection = %s\n\n\n", bandwidth)
    }
    message(c(title_print, null_print, p_val_print, samp_print, kern_print,
              band_print, band_sel))
  }
  if (suppress_raw_output == FALSE) {
    results[-4]
  }
}


#' Independence Test
#'
#' \code{independence_test} Computes the independence test with a user-specified number of principal components
#' and range of lags.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param components A positive Integer specifying the number of principal components to project the data on;
#' ranked in order of importance (importance is determined by the proportion of the variance that is explained
#' by the individual principal component.)
#' @param lag A positive Integer value, specifying the maximum lag to include - this can be seen as the bandwidth
#' or lag-window.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param suppress_raw_output Boolean value, FALSE by default. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @details The "independence" portmanteau test is a test of independence and identical distribution based on a
#' dimensionality reduction by projecting the data onto the most important functional principal components.
#' It is based on the resulting lagged cross-variances. This test is not for general white noise series, and
#' may not hold under functional conditionally heteroscedastic assumptions. Please consult the vignette for a
#' deeper exposition, and consult the reference for a complete treatment.
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#' @references
#' [1] Gabrys R., & Kokoszka P. (2007). Portmanteau Test of Independence for Functional Observations.
#' Journal of the American Statistical Association, 102:480, 1338-1348, DOI: 10.1198/016214507000001111.
#'
#' @examples
#' b <- brown_motion(250, 100)
#' independence_test(b, components = 3, lag = 5)
#'
#' @importFrom rainbow fts
#' @importFrom ftsa ftsm
#'
#' @export
independence_test <- function(f_data, components, lag, alpha = 0.05,
                              suppress_raw_output=FALSE, suppress_print_output=FALSE) {
  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. Atleast one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }
  if ((components < 1) | (components %% 1 != 0)) {
    stop("The 'components parameter must be a positive integer.")
  }
  if ((lag < 1) | (lag %% 1 != 0)) {
    stop("The 'components lag must be a positive integer.")
  }
  N <- NCOL(f_data)
  J <- NROW(f_data)
  f_data <- center(f_data)
  suppressWarnings(pc_decomp <- ftsa::ftsm(rainbow::fts(1:J, f_data), order = components, mean = FALSE))
  scores <- pc_decomp$coeff
  C_0 <- crossprod(scores) / N
  c_h <- array(0, dim=c(components,components,lag))
  for (h in 1:lag) {
    for (k in 1:components) {
      for (l in 1:components) {
        score_uni <- 0
        for (t in 1:(N-h)) {
          score_uni <- score_uni + (scores[t,k] * scores[t+h,l])
        }
        c_h[k,l,h] <- score_uni / N
      }
    }
  }
  r_f_h <- r_b_h <- array(0, dim=c(components,components,lag))
  summand <- vector('numeric', lag)
  for (h in 1:lag) {
    r_f_h[,,h] <- solve(C_0) %*% c_h[,,h]
    r_b_h[,,h] <- c_h[,,h] %*% solve(C_0)
    summand[h] <- sum(r_f_h[,,h] * r_b_h[,,h])
  }
  Q_n <- N * sum(summand)
  p_val <- as.numeric(1 - pchisq(Q_n, df = components^2 * lag))
  quantile <- as.numeric(qchisq(1 - alpha, df = components^2 * lag))
  results <- list(statistic = Q_n, quantile = quantile, p_value = p_val)
  if (suppress_print_output == FALSE) {
    title_print <- sprintf("Independence Test\n\n")
    test_type <- 'the series is iid\n'
    null_print <- sprintf("null hypothesis: %s", test_type)
    p_val_print <- sprintf("p-value = %f\n", results$p_value)
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    comp_print <- sprintf('number of principal components = %d\n', components)
    lag_print <- sprintf("maximum lag = %d\n\n\n", lag)
    message(c(title_print, null_print, p_val_print, comp_print,
          lag_print))
  }
  if (suppress_raw_output == FALSE) {
    results
  }
}

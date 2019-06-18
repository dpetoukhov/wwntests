#' Compute Functional Hypothesis Tests
#'
#' \code{fport_test} Computes a variety of functional portmanteau hypothesis tests. All hypothesis tests in this
#' package are accessible through this function.
#'
#' @param f_data The functional data matrix with observed functions in the columns.
#' @param test A String specifying the hypothesis test. Currently available tests are referred to by their
#' string handles: "single-lag", "multi-lag", "spectral", "independence", and "imhof". Please see the Details
#' section of the documentation, or the vignette, for a short overview of the available tests. For a more
#' complete treatment of these hypothesis tests, please consult the references.
#' @param lag A positive integer value. Only used for the "single-lag", "multi-lag", "independence", and "imhof" tests.
#' This parameter specifies the single lag, or maximum lag, to be used by the specified test.
#' @param iid Only used for the "single-lag" and "multi-lag" tests. A Boolean value, FALSE by default. If given TRUE,
#' the hypothesis test will use a strong-white noise assumption (instead of a weak-white noise assumption).
#' @param M Only used for the "single-lag" and "multi-lag" tests. A positive Integer. Determines the number of
#' Monte-Carlo simulations employed in the Welch-Satterthwaite approximation of the limiting distribution of the
#' test statistic.
#' @param low_disc Only used for the "single-lag" and "multi-lag" tests. A Boolean value, FALSE by default.
#' If given TRUE, uses low-discrepancy sampling in the Monte-Carlo method. Note, low-discrepancy sampling will
#' yield deterministic results. Requires the 'fOptions' package.
#' @param kernel Only used for the "spectral" test. A String, 'Bartlett' by default. Specifies the kernel to be
#' used in the "spectral" test. Currently supported kernels are the 'Bartlett' and 'Parzen' kernels.
#' @param bandwidth Only used for the "spectral" test. Either a String or a positive Integer value, 'adaptive' by
#' default. Determines the bandwidth (or lag-window) to be used for the test. Given the string handle 'adaptive',
#' the bandwidth is computed via a bandwidth selection method which aims to minimize the integrated normed
#' error of the spectral density operator. If the given string handle is 'static', the bandwidth is computed
#' to be n^(1/(2q + 1)), where n is the sample size and q is the kernel order. If a positive integer is
#' given, that will be the bandwidth that is used.
#' @param components Only used for the "independence" test. A positive Integer value. Determines the number of
#' functional principal components to use (ranked by their importance).
#' @param bootstrap Only used for the "single-lag" test. A Boolean value, FALSE by default. If given TRUE, the
#' hypothesis test is evaluated by approximating the limiting distribution of the test statistic via a block
#' bootstrapping process.
#' @param block_size Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A positive Integer
#' value, with the default value being computed via the adaptive bandwidth selection method in the "spectral" test.
#' Determines the block size (of each block in each bootstrap sample) if the test is being bootstrapped.
#' @param straps Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A positive Integer with
#' a default value of 300. Determines the number of bootstrap samples to take if the test is being bootstrapped.
#' @param moving Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A Boolean value, FALSE
#' by default If given TRUE, the performed block bootstrap will be moving rather than stationary.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param complete_test A Boolean value, FALSE by default. If TRUE, the function requires no other parameters
#' other than f_data, and will return a table with a single column containing p-values from an array of tests
#' contained in the rows.
#' @param suppress_raw_output A Boolean value, FALSE by default. If given TRUE, the function will not return a
#' list containing the p-value, quantile and statistic, and instead only prints output to the console.
#' @param suppress_print_output A Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @details The "single-lag" portmanteau test is based on the sample autocovariance function computed from the
#' functional data. This test assesses the significance of lagged autocovariance operators at a single, user-specified
#' lag h. More specifically, it tests the null hypothesis that the lag-h autocovariance operator is equal to 0.
#' This test is designed for stationary functional time-series, and is valid under conditional heteroscedasticity
#' conditions.
#' The required parameter for this test are 'lag', which determines the lag at which the test is evaluated. If this
#' parameter is left blank, it will take a default of 1.
#' The optional parameters for this test are 'iid', 'M', 'low_disc', 'bootstrap', 'block_size', 'straps', 'moving',
#' and 'alpha'.
#'
#' The "multi-lag" portmanteau test is also based on the sample autocovariance function computed from the functional
#' data. This test assesses the cumulative significance of lagged autocovariance operators, up to a user-selected
#' maximum lag K. More specifically, it tests the null hypothesis that the first K lag-h autocovariance operators
#' (h going from 1 to K) is equal to 0. This test is designed for stationary functional time-series, and is valid
#' under conditional heteroscedasticity conditions.
#' The required parameter for this test is 'lag', which determines the maximum lag at which the test is evaluated.
#' If this parameter is left blank, it will take a default of 20.
#' The optional parameters for this test are 'iid', 'M', 'low_disc', 'bootstrap', 'block_size', 'straps', 'moving',
#' and 'alpha'.
#'
#' The "spectral" portmanteau test is based on the spectral density operator. It essentially measures the proximity of a
#' functional time series to a white noise - the constant spectral density operator of an uncorrelated series.
#' Unlike the "single-lag" and "multi-lag" tests, this test is not for general white noise series, and may not hold
#' under functional conditionally heteroscedastic assumptions.
#' The optional parameters for this test are 'kernel', 'bandwidth', and 'alpha'.
#'
#' The "independence" portmanteau test is a test of independence and identical distribution based on a dimensionality
#' reduction by projecting the data onto the most important functional principal components. It is based on the
#' resulting lagged cross-variances. This test is not for general white noise series, and may not hold under
#' functional conditionally heteroscedastic assumptions.
#' The required parameters for this test are 'lag' and 'components'. The 'lag' parameter determines the maximum lag at
#' which the test is evaluated. The 'components' parameter determines the number of the most important principal
#' components to use (importance is determined by the proportion of the variance that is explained by the
#' individual principal component.)
#'
#' The "imhof" portmanteau test is an analogue of the "single-lag" test. While the "single-lag" test computes the
#' limiting distribution of the test statistic via a Welch-Satterthwaite approximation, the "imhof" test directly
#' computes the coefficients of the quadratic form in Normal variables which the test statistic converges too as
#' the sample size goes to infinity. We warn the user that this test is extremely computationally expensive, and
#' is only recommended for small datasets as a means of cross-verification against the single-lag test.
#' The required parameter for this test is 'lag', which determines the lag at which the test is evaluated.
#' The "imhof" test requires the "tensorA" and "CompQuadForm" packages. Note also that the imhof test does not
#' return a statistic, and thus returns a list with only 2 elements if suppress_raw_output = FALSE.
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE. If 'complete-test' = TRUE, will return a 1-column table instead containing
#' the p-values for a variety of tests, which are given short descriptions in the index of the table.
#'
#' @references
#' [1] Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' [2] Characiejus V., & Rice G. (2019). A general white noise test based on kernel lag-window estimates of the
#' spectral density operator. Econometrics and Statistics, submitted.
#'
#' [3] Gabrys R., & Kokoszka P. (2007). Portmanteau Test of Independence for Functional Observations.
#' Journal of the American Statistical Association, 102:480, 1338-1348, DOI: 10.1198/016214507000001111.
#'
#' [4] Zhang X. (2016). White noise testing and model diagnostic checking for functional time series.
#' Journal of Econometrics, 194, 76-95.
#'
#' [5] Chen W.W. & Deo R.S. (2004). Power transformations to induce normality and their applications.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66, 117–130.
#'
#' @examples
#' b <- brown_motion(250, 50)
#' fport_test(b, test = 'single-lag', lag = 10)
#' fport_test(b, test = 'multi-lag', lag = 10, alpha = 0.01)
#' fport_test(b, test = 'single-lag', lag = 1, M = 250, low_disc = TRUE)
#' fport_test(b, test = 'spectral', kernel = 'Bartlett', bandwidth = 'static', alpha = 0.05)
#' fport_test(b, test = 'spectral', alpha = 0.1, kernel = 'Parzen', bandwidth = 'adaptive')
#' fport_test(b, test = 'independence', components = 3, lag = 3)
#'
#' @export
#' @import stats
fport_test <- function(f_data, test = 'multi-lag', lag=NULL, iid=FALSE, M=NULL,
                       low_disc=FALSE, kernel = "Bartlett", bandwidth = "adaptive",
                       components = 3, bootstrap=FALSE, block_size = "adaptive", moving=FALSE,
                       straps = 300, alpha=0.05, complete_test=FALSE,
                       suppress_raw_output = FALSE, suppress_print_output = FALSE) {
  tests = c('single-lag', 'multi-lag', 'spectral', 'independence', 'imhof')
  if (test == 'multi-lag' & is.null(lag) & complete_test==FALSE) {
    warning("You did not specify a maximum lag for the multi-lag test. We use a default of lag = 20")
    lag = 20
  }
  if (test == 'single-lag' & is.null(lag) & complete_test==FALSE) {
    warning("You did not specify a maximum lag for the single-lag test. We use a default of lag = 1")
    lag = 1
  }
  if (!(test %in% tests)) {
    stop("Please see the documentation for available tests.")
  }
  if (!is.matrix(f_data)) {
    stop("Invalid arguments, functional data f_data must be passed in matrix form.")
  }
  if (!is.null(lag)) {
    if (!all.equal(lag, as.integer(lag)) | lag <= 0) {
      stop("Invalid arguments, lag must be a positive integer for the single-lag and multi-lag tests.")
    }
  }
  if (alpha < 0 | alpha > 1) {
    stop("Invalid arguments, the significance level alpha must be between 0 and 1.")
  }
  if (!is.logical(iid) | !is.logical(low_disc)) {
    stop("Invalid arguments, the iid and low_disc parameters must be logical values.")
  }
  if (!is.null(M)) {
    if (!all.equal(M, as.integer(M)) | M < 0) {
      stop("Invalid arguments, M must be a positive integer or NULL.")
    }
  }
  iid_error = base::simpleError("When iid = true, this function does not use Monte Carlo methods,
and thus also does not support low-discrepancy sequence sampling or parallelization. Please change the parameters.")
  if ((iid == TRUE) & ((low_disc == TRUE) | (!is.null(M)))) {
    stop(iid_error)
  }
  if (complete_test == TRUE) {
    m <- as.table(matrix(0, ncol = 1, 10))
    colnames(m) <- c('p_value')
    rownames(m) <- c('single-lag, lag = 1', 'single-lag, lag = 2',
                     'single-lag, lag = 3', 'multi-lag, lag = 5',
                     'multi-lag, lag = 10', 'multi-lag, lag = 20',
                     'spectral, static bandwidth',
                     'spectral, adaptive bandwidth',
                     'independence, 3 components, lag = 3',
                     'independence, 16 components, lag = 10')
    m[1] <- fport_test(f_data, test = 'single-lag', lag = 1)$p_value
    m[2] <- fport_test(f_data, test = 'single-lag', lag = 2)$p_value
    m[3] <- fport_test(f_data, test = 'single-lag', lag = 3)$p_value

    m[4] <- fport_test(f_data, test = 'multi-lag', lag = 5)$p_value
    m[5] <- fport_test(f_data, test = 'multi-lag', lag = 10)$p_value
    m[6] <- fport_test(f_data, test = 'multi-lag', lag = 20)$p_value

    m[7] <- fport_test(f_data, test = 'spectral',
                       bandwidth = 'static')$p_value
    m[8] <- fport_test(f_data, test = 'spectral',
                       bandwidth = 'adaptive')$p_value

    m[9] <- fport_test(f_data, test = 'independence',
                       components = 3, lag = 3)$p_value

    m[10] <- fport_test(f_data, test = 'independence',
                        components = 16, lag = 10)$p_value
    m
  } else if (test == 'multi-lag') {
    multi_lag_test(f_data, lag, M=M, low_disc = low_disc, iid=iid,
                   suppress_raw_output = suppress_raw_output,
                   suppress_print_output = suppress_print_output)
  } else if (test == 'single-lag') {
    single_lag_test(f_data, lag, alpha=alpha, iid=iid,
                    M=M, low_disc=low_disc, bootstrap=bootstrap,
                    block_size=block_size, straps=straps, moving = moving,
                    suppress_raw_output = suppress_raw_output,
                    suppress_print_output = suppress_print_output)
  } else if (test == 'spectral') {
    spectral_test(f_data, kernel = kernel, bandwidth = bandwidth, alpha = alpha,
                  suppress_raw_output=suppress_raw_output,
                  suppress_print_output = suppress_print_output)
  } else if (test == 'independence') {
    independence_test(f_data, components = components, lag = lag,
                      suppress_raw_output = suppress_raw_output,
                      suppress_print_output = suppress_print_output)
  } else if (test == 'imhof') {
    input <- readline("We warn the user that the imhof test is extremely computationally expensive. \n
Press [enter] if you would like to continue.")
    if (input != '') {
      stop("User cancelled the test.")
    }
    results <- imhof_test(f_data, lag)
    if (suppress_print_output == FALSE) {
      title_print <- sprintf("Imhof Test\n\n")
      test_type <- 'the series is a weak white noise\n'
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


#' Plot Confidence Bounds of Estimated Functional Autocorrelation Coefficients
#'
#' \code{autocorrelation_coeff_plot} Computes the 1-alpha upper confidence bounds for the functional
#' autocorrelation coefficients at lags h = 1:K under both weak white noise (WWN) and strong white
#' noise (SWN) assumptions. It plots the coefficients as well as the bounds for all lags h = 1:K.
#' Note, the SWN bound is constant, while the WWN is dependent on the lag.
#'
#' @param f_data The functional data matrix with observed functions in the columns.
#' @param K A positive Integer value. The maximum lag for which to compute the single-lag test (tests
#' will be computed for lags h in 1:K).
#' @param alpha A numeric value between 0 and 1 specifying the significance level to be used in the single-lag
#' test. The default value is 0.05.
#' @param M A positive Integer value. Determines the number of Monte-Carlo simulations employed in the
#' Welch-Satterthwaite approximation of the limiting distribution of the test statistics, for each test.
#' @param low_disc A Boolean value, FALSE by default. If given TRUE, uses low-discrepancy sampling in the
#' Monte-Carlo method. Note, low-discrepancy sampling will yield deterministic results.
#' Requires the 'fOptions' package.
#' @details This function computes and plots autocorrelation coefficients at lag h, for h in 1:K. It also
#' computes an estimated asymptotic 1 - alpha confidence bound, under the assumption that the series
#' forms a weak white noise. Additionally, it computes a similar (constant) bound under the assumption the
#' series form a strong white noise. Please see the vignette or the references for a more complete treatment.
#' @return Plot of the estimated autocorrelation coefficients for lags h in 1:K with the weak
#' white noise 1-alpha upper  confidence bound for each lag, as well as the constant strong white
#' noise 1-alpha confidence bound.
#'
#' @references
#' [1] Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' b <- brown_motion(75, 40)
#' autocorrelation_coeff_plot(b)
#' autocorrelation_coeff_plot(b, M = 200, low_disc = TRUE)
#'
#' @export
#' @import sde
#' @importFrom graphics legend lines par plot
autocorrelation_coeff_plot <- function(f_data, K=20, alpha=0.05, M=NULL, low_disc=FALSE) {
  if ((K < 1) | (K %% 1 != 0)) {
    stop("The parameter 'K' must be a positive integer.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }
  J = NROW(f_data)
  coefficients = array(0, K)
  B_h_bounds = array(0,K)
  B_iid_bounds = array(0,K)
  lags = 1:K
  for (h in lags){
    coefficients[h] <- autocorrelation_coeff_h(f_data, h)
    B_h_bounds[h] <- B_h_bound(f_data, h, M=M, low_disc=low_disc)
  }
  par(mfrow=c(1,1))
  plot(lags, coefficients, ylim=c(0,2 * max(coefficients)), type='h', xlab='Lag',
         ylab='Autocorrelation Coefficient', main = 'Autocorrelation Bounds')
  lines(B_h_bounds, col='blue', lty='dotted')
  lines(rep(B_iid_bound(f_data), K), col='red', lty='solid')
  legend('topleft',
         legend=c('Estimated Autocorrelation Coefficients', 'WWN Bound', 'SWN Bound'),
         col=c('black', 'blue', 'red'), lty=c('solid', 'dotted', 'solid'), cex=0.75)
}

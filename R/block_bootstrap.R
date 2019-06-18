#' \code{block_bootstrap} Performs a block bootstrap on the functional data f_data with block size b.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param b the block size (of each block in each bootstrap sample)
#' @param B the number of bootstraps samples
#' @param moving boolean value specifying whether the block bootstrap should be moving or not. A moving black
#' bootstrap samples individual functional observations and adds on the consequent block, rather than sampling
#' blocks of the data.
#' @return Returns a list of B elements, each element being a block bootstrap sample in the same format
#' as the original functional data f_data.
#'
#' @export
#'
block_bootsrap <- function(f_data, b, B = 300, moving = FALSE) {
  N <- NCOL(f_data)
  if (b > N) {
    stop("Please select a block size that is less than or equal to the sample size of
         the functional data. It is best to select a block size that evenly divides the
         sample size.")
  } else if (b < 1) {
    stop("The block size must be a positive integer.")
  } else if (B < 1) {
    stop("The number of bootstrap samples must be a positive integer.")
  }
  blocks <- list()
  M <- floor(N / b)
  for (s in 1:M) {
    blocks[[s]] <- (b*(s - 1) + 1):(b*s)
  }
  bootstrap_samples <- list()
  for (j in 1:B) {
    if (moving == FALSE) {
      samples <- sample(1:M, M, replace = TRUE)
      bootstrapped_data <- f_data[,blocks[[samples[1]]]]
      for (i in samples[-1]) {
        bootstrapped_data <- cbind(bootstrapped_data, f_data[,blocks[[samples[i]]]])
      }
    } else if (moving == TRUE) {
      samples <- sample(1:(N - b), M, replace = TRUE)
      bootstrapped_data <- f_data[, samples[1]:(samples[1] + b)]
      for (i in 2:M) {
        bootstrapped_data <- cbind(bootstrapped_data, f_data[,samples[i]:(samples[i] + b)])
      }
    }
    bootstrap_samples[[j]] <- bootstrapped_data
  }
  bootstrap_samples
}

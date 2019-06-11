# center centers the functional data f_data by substracting the row means from the data.
# Input: f_data = the functional data matrix with observed functions in columns
# Output: a matrix containing the centered functional data
# roxygen comments:
#' Center functional data
#'
#' \code{center} Centers the given functional data
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @return A matrix of the same form as f_data containing the centered functional data.
center <- function(f_data) {
  c_f_data <- f_data - rowMeans(f_data)
  c_f_data
}


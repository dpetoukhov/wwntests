#' Bartlett Kernel Function
#'
#' \code{bartlett_kernel} Computes the Bartlett kernel function at a given point value.
#' @param x the point value at which the kernel function is evaluated
#' @return A scalar value; the value of the Bartlett kernel function at the point value x.
bartlett_kernel <- function(x) {
  len <- length(x)
  for (i in 1:len) {
    if (abs(x[i]) <= 1) {
      x[i] <- 1 - abs(x[i])
    } else {
      x[i] <- 0
    }
  }
  x
}


#' Parzen Kernel Function
#'
#' \code{parzen_kernel} Computes the Parzen kernel function at a given point value.
#' @param x the point value at which the kernel function is evaluated
#' @return A scalar value; the value of the Parzen kernel function at the point value x.
parzen_kernel <- function(x) {
  len <- length(x)
  for (i in 1:len) {
    if (abs(x[i]) <= 1) {
      if (abs(x[i]) <= 0.5) {
        x[i] <- 1 - 6 * x[i]^2 + 6 * abs(x[i])^3
      } else {
        x[i] <- 2 * (1 - abs(x[i]))^3
      }
    } else {
      x[i] <- 0
    }
  }
  x
}


#' Daniell Kernel Function
#'
#' \code{daniell_kernel} Computes the Daniell kernel function at a given point value.
#' @param x the point value at which the kernel function is evaluated
#' @return A scalar value; the value of the Daniell kernel function at the point value x.
daniell_kernel <- function(x) {
  len <- length(x)
  for (i in 1:len) {
    if (x[i] == 0) {
      x[i] <- 0
    } else {
      x[i] <- sin(pi * x[i]) / (pi * x[i])
    }
  }
  x
}

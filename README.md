
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wwntests

The goal of the wwntests package is to provide an array of functional
hypothesis tests tests and related visualizations. Currently supported
hypothesis tests include weak and strong white noise tests, which are
built for general white noise and thus hold under conditional
heteroscedasticity assumptions, as well as tests of independence and
identical distribution.

## Installation

Using the *devtools* package:

``` r
install.packages("devtools")
library(devtools)
install_github("dpetoukhov/wwntests")
```

## CRAN

You can install the released version of wwntests from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("wwntests")
```

## Examples

We provide some basic examples which shows you how to apply the package
to functional data. We start with a sequence of Brownian motion curves
and a sequence of FAR(1, 0.75)-IID (Functional Auto Regressive Process
with IID errors) curves, both of which are available in *wwntests* via
the *brown\_motion* and *far\_1\_S* functions. For both functions, *N*
determines the number of samples and *J* determines the number of times
each individual curve is sampled. The FAR(1,S) process is not a weak
white noise (nor a strong white noise) when S \> 0.

``` r
set.seed(1997)
b <- brown_motion(N = 150, J = 100)
f <- far_1_S(N = 150, J = 100, S = 0.75)
```

We might wish to test whether or not our functional data is a weak-white
noise. The ‘single-lag’ and ‘multi-lag’ tests test this null hypothesis.
For the sake of brevity, we’ll look just at the ‘single-lag’ test,
please refer to the vignette for a more in depth look into both tests.

Note, all tests can be accessed through the *fport\_test* by passing
their string handles (see documentation) to the argument *test*. If you
are interested in a particular test, stand-alone functions are also
available for each test (see documentation).

To be more specific, the ‘single-lag’ test at lag h tests the null
hypothesis that the lag-h autocovariance function of the data is equal
to 0. We apply the single-lag test at lag =
1.

``` r
# Note, we pass TRUE to the suppress_raw_output to suppress the list output of the function.
fport_test(f_data = b, test = 'single-lag', lag = 1, suppress_raw_output = TRUE)
#>  Single-Lag Test
#> 
#>  null hypothesis: the series is a weak white noise
#>  p-value = 0.656092
#>  sample size = 150
#>  lag = 1
fport_test(f_data = f, test = 'single-lag', lag = 1, suppress_raw_output = TRUE)
#>  Single-Lag Test
#> 
#>  null hypothesis: the series is a weak white noise
#>  p-value = 0.000000
#>  sample size = 150
#>  lag = 1
```

Our p-values are as expected, given that a Brownian motion is indeed a
weak white noise while the FAR data is not.

Now suppose we want to test whether or not the functional data is
independent and identically distributed. The ‘spectral’ and
‘independence’ tests test this as a null hypothesis. For brevity, we
will give an example using the ‘spectral’ test, which is based on the
spectral density operator of the functional data, and essentially
measures the proximity of a functional time series to a white noise -
the constant spectral density operator of an uncorrelated
series.

``` r
fport_test(f_data = b, test = 'spectral', kernel = 'Bartlett', bandwidth = 'adaptive', 
           suppress_raw_output = TRUE)
#>  Spectral Test
#> 
#>  null hypothesis: the series is iid
#>  p-value = 0.962314
#>  sample size = 150
#>  kernel function = Bartlett
#>  bandwidth = 11.347351
#>  bandwidth selection = adaptive
fport_test(f_data = f, test = 'spectral', kernel = 'Bartlett', bandwidth = 'adaptive', 
           suppress_raw_output = TRUE)
#>  Spectral Test
#> 
#>  null hypothesis: the series is iid
#>  p-value = 0.000000
#>  sample size = 150
#>  kernel function = Bartlett
#>  bandwidth = 6.951734
#>  bandwidth selection = adaptive
```

Again, our p-values are as expected. The ‘spectral’ and ‘independence’
test both provide the advantage that they are very quick, much less
computationally burdensome than the ‘single-lag’ and (especially) the
‘multi-lag’ test. However, the ‘spectral’ and ‘independence’ test have
the drawback that they are not designed for general white noise and may
not hold under functional conditionally heteroscedastic assumptions. To
see this in action, we will consider an fGARCH(1,1) (Functional
Generalized Autoregressive Conditional Heteroscedasticity) process. This
process is a weak-white noise, but not a strong white noise, and as a
result, we may expect the ‘spectral-test’ to falsely assert that the
fGARCH(1,1) process is independent and identically distributed.

``` r
fgarch <- fgarch_1_1(N = 150, J = 100)
fport_test(f_data = fgarch, test = 'spectral', kernel = 'Bartlett', bandwidth = 'adaptive',
           suppress_raw_output = TRUE)
#>  Spectral Test
#> 
#>  null hypothesis: the series is iid
#>  p-value = 0.671692
#>  sample size = 150
#>  kernel function = Bartlett
#>  bandwidth = 5.513880
#>  bandwidth selection = adaptive
```

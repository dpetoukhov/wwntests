context('Testing the Hypothesis Test Function')

set.seed(42)
b1 <- brown_motion(100, 50)
f1 <- far_1_S(100, 50, 0.7)

test_that("fport_test Producing Valid Output for All Valid Combinations of Parameters", {
  # Testing single-lag test
  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', lag = 12), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', lag = 12), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', lag = 2, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', lag = 2, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', lag = 1, M=250), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', lag = 1, M=250), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', lag = 2, low_disc=TRUE, M=100), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', lag = 2, low_disc=TRUE, M=100), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', lag = 20, alpha=0.1, M=500, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', lag = 20, alpha=0.1, M=500, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', iid=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', iid=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', bootstrap = TRUE, lag = 1), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', bootstrap = TRUE, lag = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', bootstrap = TRUE, lag = 1, moving = TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', bootstrap = TRUE, lag = 5, moving = TRUE, straps = 150), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', bootstrap = TRUE, lag = 1, block_size = 1), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', bootstrap = TRUE, lag = 5, block_size = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'single-lag', bootstrap = TRUE, lag = 1, straps = 150, block_size = 1), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'single-lag', bootstrap = TRUE, lag = 5, straps = 100, block_size = 10), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))




  # Testing multi-lag test
  expect_success(expect_is(res <- fport_test(b1, test = 'multi-lag', lag = 10), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'multi-lag', lag = 10), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'multi-lag', lag = 20, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'multi-lag', lag = 20, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'multi-lag', lag = 8, M=250), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'multi-lag', lag = 8, M=250), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'multi-lag', lag = 5, alpha=0.1, M=100, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'multi-lag', lag = 5, alpha=0.1, M=100, low_disc=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'multi-lag', lag = 10, iid=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'multi-lag', lag = 10, iid=TRUE), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))


  # Testing spectral test
  expect_success(expect_is(res <- fport_test(b1, test = 'spectral'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', bandwidth = 2), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Parzen', bandwidth = 'static'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Parzen', bandwidth = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Daniell', bandwidth = 'static'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Daniell', bandwidth = 10), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Daniell', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Parzen', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'spectral', kernel = 'Bartlett', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', bandwidth = 2), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Parzen', bandwidth = 'static'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Parzen', bandwidth = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Daniell', bandwidth = 'static'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Daniell', bandwidth = 10), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Daniell', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Parzen', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'spectral', kernel = 'Bartlett', bandwidth = 'adaptive'), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))


  # Testing Independence Test
  expect_success(expect_is(res <- fport_test(b1, test = 'independence', lag = 3), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'independence', components = 5, lag = 20), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'independence', lag = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(b1, test = 'independence', components = 2, lag = 2), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'independence', lag = 3), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'independence', components = 5, lag = 1), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'independence', lag = 5), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))

  expect_success(expect_is(res <- fport_test(f1, test = 'independence', components = 2, lag = 2), "list"))
  expect_success(expect_length(res, 3))
  expect_success(expect_equal(sum(is.na(res)), 0))
})


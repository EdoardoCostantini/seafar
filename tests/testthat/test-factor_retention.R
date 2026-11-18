testthat::skip_if_not_installed("lavaan")

data(HolzingerSwineford1939, package = 'lavaan')
HS_1939_std <- as.matrix(scale(HolzingerSwineford1939[, paste0('x', 1:9)]))
stopifnot(is.matrix(HS_1939_std))

test_that("Parallel Analysis is equivalent to nFactors::parallel for N > J case", {
  set.seed(10)
  parallel_analysis_output <- parallel_analysis(
    N = nrow(HS_1939_std),
    J = ncol(HS_1939_std),
    iterations = 1e3,
    quant = 0.05
  )
  set.seed(10)
  nFactors_pa <- nFactors::parallel(
    subject = nrow(HS_1939_std),
    var = ncol(HS_1939_std),
    rep = 1e3,
    quantile = 0.05
  )

  expect_true(
    cor(
      parallel_analysis_output$summary_statistics_eigenvalues$mean,
      nFactors_pa$eigen$mevpea
    ) > 0.99
  )
})

test_that("Parallel Analysis works for N < J case", {
  set.seed(10)
  parallel_analysis_output <- parallel_analysis(
    N = nrow(autism),
    J = ncol(autism),
    iterations = 2,
    quantile = 0.05
  )

  expect_type(parallel_analysis_output, "list")
  expect_true(all(
    c("sampled_eigenvalues", "summary_statistics_eigenvalues") %in%
      names(parallel_analysis_output)
  ))
})

test_that("Factor retention works on N > J", {
  set.seed(10)
  Q <- factor_retention(data = HS_1939_std)

  expect_s3_class(Q, "nfactors")
  expect_type(Q, "list")
  expect_true(all(
    c("kaiser", "acceleration_factor", "optimal_coordinates", "parallel") %in%
      names(Q)
  ))
  expect_true(all(sapply(Q, is.numeric)))
})

test_that("Factor retention works on N < J", {
  set.seed(10)
  Q <- factor_retention(data = autism, iterations = 2)

  expect_s3_class(Q, "nfactors")
  expect_type(Q, "list")
  expect_true(all(
    c("kaiser", "acceleration_factor", "optimal_coordinates", "parallel") %in%
      names(Q)
  ))
  expect_true(all(sapply(Q, is.numeric)))
})

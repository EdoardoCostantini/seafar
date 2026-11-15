testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

J <- ncol(hs)
N <- nrow(hs)
Q <- 3

test_that("seafar_lasso works", {
  fit <- seafar_lasso(hs, nfactors = Q, lambda = 100, INIT = 'svd')

  expect_s3_class(fit, "lasso")
  expect_named(fit, c("scores", "loadings", "PVE", "Residual", "converged", "center", "scale"))
  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(dim(fit$scores), c(N, Q))
  expect_true(fit$converged)
  expect_true(any(fit$loadings == 0))
})

test_that("seafar_lasso_multistart works", {
  fit <- seafar_lasso_multistart(hs, nfactors = Q, lambda = 100, nstarts = 5)

  expect_s3_class(fit, "lasso_multistart")
  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(fit$lambda, 100)
  expect_length(fit$all_losses, 5)

  # unshrink refits the zero pattern of the lasso solution without penalty
  expect_equal(fit$loadings == 0, fit$loadings_lasso == 0)
})

test_that("summary prints the selected model", {
  fit <- seafar_lasso_multistart(hs, nfactors = Q, lambda = 100, nstarts = 5)

  expect_output(summary(fit), "The LASSO penalty is: 100.000")
  expect_output(summary(fit), sprintf("The PEV is: %.3f", fit$PVE))
})

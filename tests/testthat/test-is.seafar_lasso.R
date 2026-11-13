testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

J <- ncol(hs)
Q <- 3

test_that("is.seafar_lasso with warm starts works", {
  res <- is.seafar_lasso(hs, nfactors = Q)

  expect_s3_class(res, "is_seafar_lasso")
  expect_named(res, c("sellambda", "selmodel", "pve", "pvepca"))
  expect_equal(sum(res$selmodel$loadings != 0), 12)
})

test_that("is.seafar_lasso without warm starts works", {
  res <- is.seafar_lasso(hs, nfactors = Q, warmstart = FALSE, INIT = 'svd', lambda.length = 5)
  expect_s3_class(res, "is_seafar_lasso")
  expect_named(res, c("sellambda", "selmodel", "pve", "pvepca"))
})

test_that("summary prints the selection results", {
  res <- is.seafar_lasso(hs, nfactors = Q)

  expect_output(summary(res), sprintf("%.3f", res$sellambda), fixed = TRUE)
  expect_output(summary(res), "selected by the Index of Sparseness")
  expect_output(summary(res$selmodel), "The LASSO penalty is")
})

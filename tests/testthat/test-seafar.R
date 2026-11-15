testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

J <- ncol(hs)
N <- nrow(hs)
Q <- 3
C <- 9

expect_valid_seafar <- function(fit, C_total, orthogonal = FALSE) {
  expect_s3_class(fit, "seafar")
  expect_type(fit, "list")
  expect_named(fit, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(dim(fit$scores), c(N, Q))
  expect_equal(sum(fit$loadings != 0), C_total)

  expect_true(all(fit$PVE <= 1))
  expect_length(fit$Residual, 1)
  expect_true(fit$Residual >= 0)
}

test_that("seafar with a single start works for orthogonal factors", {
  fit <- seafar(
    hs,
    nfactors = Q,
    C = C,
    orthogonal = TRUE,
    INIT = "svd",
    nstarts = 1
    )

  expect_s3_class(fit, "seafar")
  expect_type(fit, "list")
  expect_named(fit, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(dim(fit$scores), c(N, Q))
  expect_equal(sum(fit$loadings != 0), C)
})

test_that("seafar with multiple starts works for orthogonal factors", {
  fit <- seafar(hs, nfactors = Q, C = C, orthogonal = TRUE, nstarts = 10)
  expect_s3_class(fit, "seafar")
  expect_type(fit, "list")
  expect_named(fit, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(dim(fit$scores), c(N, Q))
  expect_equal(sum(fit$loadings != 0), C)
})

test_that("seafar with multiple starts works for correlated factors", {
  fit <- seafar(hs, nfactors = Q, C = C, nstarts = 10)

  expect_s3_class(fit, "seafar")
  expect_type(fit, "list")
  expect_named(fit, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit$loadings), c(J, Q))
  expect_equal(dim(fit$scores), c(N, Q))
  expect_equal(sum(fit$loadings != 0), C)
})

test_that("seafar accepts a vector of nonzero loadings per factor", {
  C_vec <- c(3, 3, 3)
  fit_ortho <- seafar(hs, nfactors = Q, C = C_vec, orthogonal = TRUE, nstarts = 5)
  fit_corr <- seafar(hs, nfactors = Q, C = C_vec, nstarts = 5)

  expect_s3_class(fit_ortho, "seafar")
  expect_type(fit_ortho, "list")
  expect_named(fit_ortho, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit_ortho$loadings), c(J, Q))
  expect_equal(dim(fit_ortho$scores), c(N, Q))

  expect_s3_class(fit_corr, "seafar")
  expect_type(fit_corr, "list")
  expect_named(fit_corr, c("scores", "loadings", "PVE", "Residual", "center", "scale"))

  expect_equal(dim(fit_corr$loadings), c(J, Q))
  expect_equal(dim(fit_corr$scores), c(N, Q))

  expect_equal(unname(colSums(fit_ortho$loadings != 0)), C_vec)
  expect_equal(unname(colSums(fit_corr$loadings != 0)), C_vec)
})

test_that("seafar errors when every start fails", {
  hs_const <- cbind(hs, x10 = 1)
  expect_error(
    seafar(hs_const, nfactors = Q, C = C, nstarts = 3),
    "all 3 starts failed.*zero variance"
  )
})

test_that("INIT defaults to semirational and rejects unknown methods", {

  expect_s3_class(seafar(hs, Q, C), "seafar")

  expect_error(seafar(hs, Q, C, INIT = "rnadom"), "INIT must be one of")
})

testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

stopifnot(is.matrix(hs))

J <- ncol(hs)
N <- nrow(hs)
Q <- 3
candidates <- seq(J * Q - 1, 3 * Q) # all cardinalities tested when card.length is large

test_that("seafar.IS with warm starts works with orthogonal factors", {
  res <- seafar.IS(data = hs, nfactors = Q, orthogonal = TRUE)
  expect_equal(res$selcard, 9)

  expect_true(all(c(
    "cardinality", "value", "pve", "pvepca", "propzero", "smallestP",
    "selcard", "center", "scale", "selmodel"
  ) %in% names(res)))
})

test_that("seafar.IS with warm starts works with correlated factors", {
  res <- seafar.IS(data = hs, nfactors = Q)
  expect_equal(res$selcard, 9)

  expect_true(all(c(
    "cardinality", "value", "pve", "pvepca", "propzero", "smallestP",
    "selcard", "center", "scale", "selmodel"
  ) %in% names(res)))
})

test_that("seafar.IS without warm start works with orthogonal factors", {
  res <- seafar.IS(data = hs, nfactors = Q, warmstart = FALSE, orthogonal = TRUE,
                   INIT = 'svd', THR = .3)

  expect_all_equal(res$selcard, 11)

  expect_true(all(c(
    "selcard", "selmodel", "pve", "pvepca"
  ) %in% names(res)))
})

test_that("summary works on seafar.IS output", {
  res_wst <- seafar.IS(hs, nfactors = Q, warmstart = TRUE)
  expect_output(summary(res_wst), paste("selected by the Index of Sparseness is:", res_wst$selcard))

  res_cold <- seafar.IS(hs, nfactors = Q, warmstart = FALSE, nstarts = 5)
  expect_output(summary(res_cold), paste("selected by the Index of Sparseness is:", res_cold$selcard))
  expect_output(summary(res_cold), "PEV of the unconstrained model")
})

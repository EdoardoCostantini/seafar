testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

J <- ncol(hs)
N <- nrow(hs)
Q <- 3
C <- 9

test_that("IS with warm starts works", {
  res <- is.seafar_wst(
    data = hs,
    nfactors = Q
  )

  expect_s3_class(res, "ISwarmstart")
  expect_type(res, "list")
  expect_true(all(c(
    "cardinality", "value", "pve", "pvepca", "propzero", "smallestP",
    "selcard", "center", "scale"
  ) %in% names(res)))

  expect_equal(res$selcard, 9)
})

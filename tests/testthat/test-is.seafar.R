testthat::skip_if_not_installed("lavaan")

data("HolzingerSwineford1939", package = "lavaan")
hs <- as.matrix(HolzingerSwineford1939[, paste0("x", 1:9)])

J <- ncol(hs)
N <- nrow(hs)
Q <- 3
C <- 9

test_that("is.seafar_original() returns the IS components", {
  IS <- is.seafar_original(
    data = hs,
    nfactors = Q,
    C = 9,
    INIT = "svd",
    orthogonal = TRUE,
    nstarts = 1
  )

  expect_type(IS, "list")
  expect_named(
    IS,
    c("value", "vaf", "propzero", "smallestP", "maxsdP", "center", "scale")
  )

  expect_true(is.numeric(IS$value) && length(IS$value) == 1)
  expect_gt(IS$value, 0)
  expect_gt(IS$vaf, 0)
  expect_lte(IS$vaf, 1)

})

test_that("is.seafar() returns a single total cardinality", {

  selcard <- is.seafar(
    data = hs,
    nfactors = Q,
    INIT = "svd",
    orthogonal = TRUE,
    nstarts = 1
  )

  expect_s3_class(selcard, "is_seafar")
  expect_true(is.numeric(selcard))
  expect_length(selcard, 1)

  expect_equal(unname(unclass(selcard)), 25)
})


test_that("summary.is_seafar() prints the selected cardinality", {
  selcard <- is.seafar(
    data = hs,
    nfactors = Q,
    INIT = "svd",
    orthogonal = TRUE,
    nstarts = 1,
    THR = 0.2
  )

  expect_output(
    summary(selcard),
    "The cardinality selected by the Index of Sparseness is: 12"
  )
})

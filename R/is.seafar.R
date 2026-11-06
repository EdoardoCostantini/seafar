#' Model selection for SEAFAR using the Index of Sparseness.
#'
#' @param data A data frame or matrix (NxJ).
#' @param nfactors Number of factors.
#' @param C Number of nonzero loadings.
#' @param maxiter Maximum number of iterations for the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix.
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param nstarts Number of starts.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance
#'   after centering. Centering is applied regardless.
#'
#' @importFrom stats sd
#'
#' @returns A list with
#' \item{value}{IS value.}
#' \item{vaf}{Variance accounted for.}
#' \item{propzero}{Proportion of zero loadings.}
#' \item{smallestP}{Smallest nonzero loading.}
#' \item{maxsdP}{Maximum standard deviation of loading.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' IS <- is.seafar_original(
#'   data = X,
#'   nfactors = 5,
#'   C = 240,
#'   INIT = "random",
#'   orthogonal = TRUE,
#'   nstarts = 50
#' )
#' }
#'
is.seafar_original <- function(data,
                               nfactors,
                               C,
                               maxiter = 50,
                               eps = 10^-4,
                               INIT = "semirational",
                               orthogonal = FALSE,
                               nstarts = 50,
                               standardize = TRUE) {
  INIT <- check_init(INIT)
  N <- dim(data)[1]
  J <- dim(data)[2]

  # 0. Center (always) and optionally scale to unit variance (N denominator), so that
  #    vzero and va are computed on the same data as the seafar fits
  data <- scale(data, center = TRUE, scale = FALSE)
  xcenter <- attr(data, "scaled:center")
  xscale <- rep(1, J)
  if (standardize) {
    xscale <- sqrt(colSums(data^2) / N)
    if (any(xscale < sqrt(.Machine$double.eps))) {
      stop("data contains item(s) with zero variance; remove them or use standardize = FALSE")
    }
    data <- scale(data, center = FALSE, scale = xscale)
  }

  seafar_result <- seafar_multistart(
    data = data,
    nfactors = nfactors,
    C = C,
    maxiter = maxiter,
    eps = eps,
    INIT = INIT,
    orthogonal = orthogonal,
    nstarts = nstarts,
    show_progress = FALSE,
    standardize = FALSE # data already preprocessed above
  )
  Hmat <- seafar_result$scores
  Pmat <- seafar_result$loadings

  vzero <- sum(data^2)
  d <- svd(data)$d
  va <- sum((d[1:nfactors])^2)
  vs <- sum((Hmat %*% t(Pmat))^2)
  nrzeqcoef <- sum(round(Pmat, 3) == 0)
  nrcoef <- J * nfactors

  IS <- list()
  IS$value <- va * vs / vzero^2 * nrzeqcoef / nrcoef
  IS$vaf <- vs / vzero
  IS$propzero <- nrzeqcoef / nrcoef

  IS$smallestP <- ifelse(sum(rowSums(Pmat != 0)) < sum(C),
                         0, min(abs(Pmat[Pmat != 0]))
  )
  IS$maxsdP <- max(apply(Pmat, 2, sd))
  IS$center <- xcenter
  IS$scale <- xscale

  return(IS)
}


#' Model selection for SEAFAR using modified Index of Sparseness
#' (Van Deun et al., 2026).
#'
#' @param data A data frame or matrix (NxJ).
#' @param nfactors Number of factors.
#' @param maxiter Maximum number of iterations for the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix.
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param nstarts Number of starts.
#' @param TOL Number of decimals used to decide whether two PEV values are the same.
#' @param THR Thresholding selection for the smallest nonzero loading.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance
#'   after centering. Centering is applied regardless.
#' @param card.length Maximum number of candidate cardinality values to test, default is 100.
#'
#' @returns The selected cardinality value.
#' @export
#'
#' @examples
#' \dontrun{
#' selcard <- is.seafar(
#'   data = X,
#'   nfactors = 5,
#'   INIT = "random",
#'   orthogonal = TRUE,
#'   nstarts = 50,
#'   TOL = 2,
#'   THR = 0.15
#' )
#' }
is.seafar <- function(data,
                      nfactors,
                      maxiter = 50,
                      eps = 10^-4,
                      INIT = "semirational",
                      orthogonal = FALSE,
                      nstarts,
                      TOL,
                      THR,
                      standardize = TRUE) {
                      nstarts = 50,
                      TOL = 3,
                      THR = 0.03,
                      standardize = TRUE,
                      card.length = 100) {
  INIT <- check_init(INIT)
  J <- dim(data)[2]

  col_names <- c("K", "IS", "PEV", "Prop0", "MinNonZeroL", "MaxSDL")
  rows <- list()

  # --- cardinalities (total number of nonzero loadings) ---
  if (nfactors * (J - 3) > card.length) {
    cardvec <- round(seq(3 * nfactors, J * nfactors - 1, length.out = card.length))
  } else {
    cardvec <- seq(3 * nfactors, J * nfactors - 1, by = 1)
  }

  for (l in seq_along(cardvec)) {
    a <- tryCatch(
      is.seafar_original(
        data = data,
        nfactors = nfactors,
        C = cardvec[l],
        maxiter = maxiter,
        eps = eps,
        INIT = INIT,
        orthogonal = orthogonal,
        nstarts = nstarts,
        standardize = standardize
      ),
      error = function(e) NULL
    )

    if (is.null(a)) next

    rows[[length(rows) + 1]] <- c(
      cardvec[l], a$value, a$vaf, a$propzero,
      a$smallestP, a$maxsdP
    )
  }

  # --- combine successful rows only ---
  if (length(rows) > 0) {
    avec <- do.call(rbind, rows)
    avec <- as.matrix(avec)
    colnames(avec) <- col_names
    avec <- round(avec, 3)
  } else {
    stop("is.seafar_original() failed for every cardinality; run it once to see the error")
  }

  round(avec, 4)

  KL <- dim(avec)[1]

  placeholder <- seq(1:KL)
  indexnonzeroL <- placeholder[avec[, "MinNonZeroL"] < THR][1] - 1 # THR holding selection
  if (is.na(indexnonzeroL)) {
    indexnonzeroL <- KL # no smallest loading below THR: keep the largest cardinality
  } else if (indexnonzeroL == 0) {
    indexnonzeroL <- NA # below THR already at the smallest cardinality: fall back to the maximum IS
  }
  selcardinality <- avec[indexnonzeroL, 1]

  # TOL = 2
  samepve <- rep(TRUE, KL)
  if (!is.na(indexnonzeroL)) {
    samepve <- round(avec[, 3], TOL) == round(avec[indexnonzeroL, 3], TOL)
  }

  index <- indexnonzeroL
  if (is.na(indexnonzeroL) || max(round(avec[samepve, 2], 3)) > round(avec[indexnonzeroL, 2], 3)) {
    index <- placeholder[samepve & avec[, 2] == max(avec[samepve, 2])][1]
    selcardinality <- avec[index, 1]
  }
  attr(selcardinality, "class") <- "is_seafar"

  return(selcardinality)
}

#' Display a summary of the results of \code{is.seafar()}.
#'
#' @param object Object of class inheriting from 'is_seafar'.
#' @param ... Argument to be passed to or from other methods.
#'
#' @returns Summary of the results.
#' @export
#'
#' @examples
#' \dontrun{
#' summary(IS_result)
#' }
summary.is_seafar <- function(object, ...) {
  cat(sprintf(
    "The cardinality selected by the Index of Sparseness is: %s\n",
    object
  ))
}

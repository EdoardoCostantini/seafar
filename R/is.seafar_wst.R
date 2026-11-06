#' Model selection for SEAFAR using Index of Sparseness with warm start
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param maxiter Maximum number of iterations for the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix, default is \code{"rational"}.
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param card.length Number of candidate cardinality values to test, default is 100.
#' @param show_progress Option of print progress bar, default is FALSE.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns A list with
#' \item{cardinality}{Candidate cardinalities (number of nonzero loadings).}
#' \item{value}{IS value.}
#' \item{pve}{Proportion of variance explained for each cardinality.}
#' \item{pvepca}{Proportion of variance explained by the unconstrained model (first nfactors principal components).}
#' \item{propzero}{Proportion of zero loadings.}
#' \item{smallestP}{Smallest nonzero loading.}
#' \item{selcard}{Cardinality with the highest IS value.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' IS_result <- is.seafar_wst(
#'   data = X,
#'   nfactors = 5,
#'   INIT = "random",
#'   orthogonal = TRUE,
#'   card.length = 100
#' )
#' }
is.seafar_wst <- function(data,
                          nfactors,
                          eps = 1e-4,
                          maxiter = 50,
                          INIT = "rational",
                          orthogonal = FALSE,
                          card.length = 100,
                          show_progress = FALSE,
                          standardize = TRUE) {
  INIT <- check_init(INIT)
  IS <- list()
  n <- dim(data)[1]
  J <- dim(data)[2]

  # 0. Center (always) and optionally scale to unit variance (N denominator), once for the whole path,
  #    so that pvepca and the initial loadings are computed on the same data as the seafar() fits
  data <- scale(data, center = TRUE, scale = FALSE)
  xcenter <- attr(data, "scaled:center")
  xscale <- rep(1, J)
  if (standardize) {
    xscale <- sqrt(colSums(data^2) / n)
    if (any(xscale < sqrt(.Machine$double.eps))) {
      stop("data contains item(s) with zero variance; remove them or use standardize = FALSE")
    }
    data <- scale(data, center = FALSE, scale = xscale)
  }
  nrcoef <- J * nfactors

  # obtain initial loading matrix
  loadings <- seafar_init(data, nfactors, C = nrcoef, INIT = INIT)
  svd1 <- svd(data, 0, 0)
  pvepca <- sum(svd1$d[1:nfactors]^2) / sum(svd1$d^2) # PVE of the unconstrained model, independent of INIT and scaling
  PVE <- c()
  propzero <- c()
  min_nzero_loading <- c()
  isvalue <- c()
  if (nfactors * (J - 3) > card.length) {
    regpath <- round(seq(((J * nfactors) - 1), nfactors * 3, length = card.length))
  } else {
    regpath <- seq(((J * nfactors) - 1), nfactors * 3)
  }

  if (show_progress == TRUE) {
    pb <- txtProgressBar(min = 0, max = length(regpath), style = 3)
  }


  fit_fun <- if (orthogonal == TRUE) seafar_orthogonal else seafar_general

  for (i in 1:length(regpath)) {
    card <- regpath[i]
    seafa_result <- fit_fun(data, nfactors, C = card, eps = eps, maxiter = maxiter, initloadings = loadings, INIT = INIT, standardize = FALSE) # data already preprocessed above
    loadings <- seafa_result$loadings
    # scores <- seafa_result$scores
    pveseafa <- seafa_result$PVE[length(seafa_result$PVE)]
    PVE <- c(PVE, pveseafa)
    propzero <- c(propzero, 1 - card / nrcoef)
    is <- pvepca * pveseafa * (1 - card / nrcoef)
    isvalue <- c(isvalue, is)
    minl <- ifelse(sum(colSums(loadings != 0)) > card, 0,
                   min(abs(loadings[loadings != 0]))
    )
    min_nzero_loading <- c(min_nzero_loading, minl)

    if (show_progress == TRUE) {
      setTxtProgressBar(pb, i)
    }
  }

  if (show_progress == TRUE) {
    close(pb)
  }

  cardinality <- regpath

  IS$cardinality <- cardinality
  IS$value <- isvalue
  IS$pve <- PVE
  IS$pvepca <- pvepca
  IS$propzero <- propzero
  IS$smallestP <- min_nzero_loading
  IS$selcard <- cardinality[which.max(isvalue)]
  IS$center <- xcenter
  IS$scale <- xscale
  attr(IS, "class") <- "ISwarmstart"

  return(IS)
}

#' Display a summary of the results of \code{is.seafar_wst()}.
#'
#' @param object Object of class inheriting from 'ISwarmstart'.
#' @param ... Argument to be passed to or from other methods.
#'
#' @returns Summary of the results.
#' @export
#'
#' @examples
#' \dontrun{
#' summary(IS_result)
#' }
summary.ISwarmstart <- function(object, ...) {
  cat(sprintf(
    "The cardinality selected by the Index of Sparseness is: %s\n",
    object$cardinality[which.max(object$value)]
  ))

  cat(sprintf(
    "The PEV with this cardinality is: %.3f\n",
    object$pve[which.max(object$value)]
  ))

  cat(sprintf(
    "The PEV of the unconstrained model (no zero loadings) is: %.3f\n",
    object$pvepca
  ))
}

#' Model selection for SEAFAR using Index of Sparseness
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param warmstart Logical. If \code{TRUE} (default), use warm starts.
#' @param maxiter Maximum number of iterations for the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix. The default is \code{"rational"} if \code{warmstart = TRUE}
#'   and \code{"semirational"} if \code{warmstart = FALSE}.
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param nstarts Number of starts for fitting the selected model and, if \code{warmstart = FALSE}, for model selection.
#' @param TOL Number of decimals used to decide whether two PEV values are the same. Only used if \code{warmstart = FALSE}.
#' @param THR Thresholding selection for the smallest nonzero loading. Only used if \code{warmstart = FALSE}.
#' @param card.length Number of candidate cardinality values to test, default is 100.
#' @param show_progress Option of print progress bar, default is FALSE. Only used if \code{warmstart = TRUE}.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance
#'   after centering. Centering is applied regardless.
#'
#' @returns A list with
#' \item{selcard}{The selected cardinality.}
#' \item{selmodel}{The output of \code{seafar()} fitted with the selected cardinality.}
#' If \code{warmstart = TRUE}, the list also contains the output of \code{is.seafar_wst()}.
#' If \code{warmstart = FALSE}, the list also contains
#' \item{pve}{Proportion of variance explained by the selected model.}
#' \item{pvepca}{Proportion of variance explained by the unconstrained model (first nfactors principal components).}
#' @export
#'
#' @examples
#' \dontrun{
#' IS_result <- seafar.IS(
#'   data = X,
#'   nfactors = 5,
#'   warmstart = TRUE,
#'   card.length = 100
#' )
#' }
seafar.IS <- function(data,
                      nfactors,
                      warmstart = TRUE,
                      maxiter = 50,
                      eps = 10^-4,
                      INIT = if (warmstart) "rational" else "semirational",
                      orthogonal = FALSE,
                      nstarts = 50,
                      TOL = 3,
                      THR = 0.03,
                      card.length = 100,
                      show_progress = FALSE,
                      standardize = TRUE) {
  if (warmstart == TRUE) {
    result <- is.seafar_wst(
      data = data,
      nfactors = nfactors,
      eps = eps,
      maxiter = maxiter,
      INIT = INIT,
      orthogonal = orthogonal,
      card.length = card.length,
      show_progress = show_progress,
      standardize = standardize
    )
  } else {
    selcard <- is.seafar(
      data = data,
      nfactors = nfactors,
      maxiter = maxiter,
      eps = eps,
      INIT = INIT,
      orthogonal = orthogonal,
      nstarts = nstarts,
      TOL = TOL,
      THR = THR,
      standardize = standardize,
      card.length = card.length
    )
    result <- list(selcard = as.numeric(selcard))
    attr(result, "class") <- "is_seafar"
  }
  result$selmodel <- seafar(
    data = data,
    nfactors = nfactors,
    C = result$selcard,
    maxiter = maxiter,
    eps = eps,
    INIT = INIT,
    orthogonal = orthogonal,
    nstarts = nstarts,
    standardize = standardize
  )
  if (warmstart == FALSE) {
    selmodel <- result$selmodel
    svd1 <- svd(scale(data, center = selmodel$center, scale = selmodel$scale), 0, 0)
    result$pve <- selmodel$PVE[length(selmodel$PVE)]
    result$pvepca <- sum(svd1$d[1:nfactors]^2) / sum(svd1$d^2)
  }
  return(result)
}

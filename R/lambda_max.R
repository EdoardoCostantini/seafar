#' Function to find the maximum LASSO penalty parameter for a given data set
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param maxiter Maximum iterations for SEAFA with LASSO.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param nstarts Number of starts
#' @param min_items Number of minimum items per factor.
#' @param TOL Tolerance for the binary search convergence criterion.
#' @param INIT Method to initialize loading matrix.
#' @param max_iter Maximum iterations for the binary search.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns The maximum value for LASSO penalty parameter.
#' @export
#'
#' @examples
#'\dontrun{
#' max_lambda <- lambda_max(
#'   data = scale(big5),
#'   nfactors = 5,
#'   INIT = 'semirational')
#'}
lambda_max <- function(data,
                       nfactors,
                       maxiter = 50,
                       eps = 10^-4,
                       nstarts = 1,
                       min_items = 3,
                       TOL = 10^-4,
                       INIT = 'semirational',
                       max_iter = 25,
                       standardize = TRUE){
  INIT <- check_init(INIT)
  N <- dim(data)[1]
  J <- dim(data)[2]

  # 0. Center (always) and optionally scale to unit variance (N denominator)
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

  svd1 <- svd(data, nfactors, nfactors)
  scores_0 <- sqrt(N) * svd1$u

  # rough upper bound based on SVD
  C  <- abs(crossprod(scores_0, data))
  lambda_high <- 2 * min(apply(C, 1L, function(v) sort(v, TRUE)[min_items]))


  # check to make sure the upper bound produces error
  repeat {
    loadings <- try(fit_once(data = data,
                             nfactors = nfactors,
                             lambda = lambda_high,
                             maxiter = maxiter,
                             eps = eps,
                             nstarts = nstarts,
                             INIT = INIT,
                             standardize = FALSE), # data already preprocessed above
                    silent = TRUE)
    ok <- !(inherits(loadings, "try-error") ||
              any(!enough_items(loadings, min_items)))
    if (!ok) {
      break
    } else {
      lambda_high <- lambda_high * 2
    }
  }

  lambda_low <- 0

  ## 2. binary search
  for (k in seq_len(max_iter)) {
    lambda_mid <- 0.5 * (lambda_high + lambda_low)

    loadings <- try(fit_once(data = data,
                             nfactors = nfactors,
                             lambda = lambda_mid,
                             maxiter = maxiter,
                             eps = eps,
                             nstarts = nstarts,
                             INIT = INIT,
                             standardize = FALSE), # data already preprocessed above
                    silent = TRUE)
    ok <- !(inherits(loadings, "try-error") ||
              any(!enough_items(loadings, min_items)))

    if (ok) {
      lambda_low <- lambda_mid
    } else {
      lambda_high <- lambda_mid
    }

    if ((lambda_high - lambda_low) < TOL * (1 + lambda_low)) {
      break
    }
  }

  lambda_max <- lambda_low
  return(lambda_max)
}

#' Function to check there are at least certain number of items per factor
#'
#' @param loadings Loading matrix resulting from SEAFA.
#' @param min_items Minimum number of items per factors.
#' @param eps Numerical tolerance used to treat very small loadings as zero.
#'
#' @returns A logical vector indicating whether each factor has at least `min_items` nonzero loadings.
#' @export
#'
#' @examples
#'\dontrun{
#' check <- enough_items(
#'   loadings)
#'}
enough_items <- function(loadings,
                         min_items = 3,
                         eps = 1e-8){
  colSums(abs(loadings) > eps) >= min_items
}

#' Function to fit SEAFA with LASSO penalty once and get loading matrix
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param lambda LASSO penalty parameter.
#' @param maxiter Maximum iterations for SEAFA with LASSO.
#' @param eps Convergence criterion based on difference in loss between iterations.
#' @param INIT Method to initialize loading matrix.
#' @param nstarts Number of starts.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns Estimated loading matrix.
#' @export
#'
#' @examples
#'\dontrun{
#' P <- seafar_lasso_multistart(
#'   scale(big5),
#'   nfactors = 5,
#'   lambda = 10,
#'   INIT = 'semirational',
#'   nstarts = 1)
#'}
fit_once <- function(data,
                     nfactors,
                     lambda,
                     maxiter = 50,
                     eps = 10^-4,
                     INIT = 'semirational',
                     nstarts = 1,
                     standardize = TRUE){
  INIT <- check_init(INIT)
  fit <- seafar_lasso_multistart(data = data,
                                 nfactors = nfactors,
                                 maxiter = maxiter,
                                 eps = eps,
                                 nstarts = nstarts,
                                 lambda  = lambda,
                                 INIT = INIT,
                                 standardize = standardize,
                                 unshrink = FALSE)
  loadings <- fit$loadings
  return(loadings)
}

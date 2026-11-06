#' Function to calculate the IS for each LASSO penalty parameter.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors.
#' @param lambda LASSO penalty tuning parameter.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loadings.
#' @param nstarts Number of starts.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns A list with
#' \item{value}{IS value.}
#' \item{vaf}{Variance accounted for.}
#' \item{propzero}{Proportion of zero loadings.}
#' \item{smallestP}{Smallest nonzero loading.}
#' \item{maxsdP}{Maximum standard deviation of loading.}
#' @export
#'
#' @examples
#' \dontrun{
#' IS <- IS_lasso(
#'   data = X,
#'   nfactors = 5,
#'   lambda = 100,
#'   INIT = "semirational",
#'   nstarts = 50
#' )
#' }
IS_lasso <- function(data,
                     nfactors,
                     lambda,
                     maxiter = 50,
                     eps = 10^-4,
                     INIT = "semirational",
                     nstarts = 50,
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

  VarSelect0 <- seafar_lasso_multistart(
    data = data,
    nfactors = nfactors,
    lambda = 0,
    maxiter = maxiter,
    eps = eps,
    INIT = INIT,
    nstarts = nstarts,
    standardize = FALSE, # data already preprocessed above
    unshrink = FALSE
  )
  P_hat0 <- VarSelect0$loadings
  T_hat0 <- VarSelect0$scores

  V_oo <- sum(data^2)
  V_s <- sum((T_hat0 %*% t(P_hat0))^2)

  VarSelect <- seafar_lasso_multistart(
    data = data,
    nfactors = nfactors,
    lambda = lambda,
    maxiter = maxiter,
    eps = eps,
    INIT = INIT,
    nstarts = nstarts,
    standardize = FALSE, # data already preprocessed above
    unshrink = FALSE
  )
  P_hat <- VarSelect$loadings
  T_hat <- VarSelect$scores

  card <- sum(P_hat != 0)

  V_a <- sum((T_hat %*% t(P_hat))^2)
  IS <- list()
  IS$value <- (V_a * V_s / V_oo^2) * (sum(round(P_hat, 3) == 0) / (J * nfactors))
  IS$vaf <- V_a / V_oo
  IS$propzero <- sum(round(P_hat, 3) == 0) / (J * nfactors)

  IS$smallestP <- ifelse(sum(rowSums(P_hat != 0)) < sum(card),
                         0, min(abs(P_hat[P_hat != 0]))
  )
  IS$maxsdP <- max(apply(P_hat, 2, sd))
  return(IS)
}


#' Model selection procedure with IS to choose the LASSO penalty paramter.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors.
#' @param max_lambda Maximum LASSO penalty parameter.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix.
#' @param nstarts Number of starts.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns The selected LASSO penalty parameter.
#' @export
#'
#' @examples
#'\dontrun{
#' lambda <- is.seafar_lasso(
#'   data = X,
#'   nfactors = 5,
#'   warmstart = FALSE,
#'   nstarts = 50)
#'}
is.seafar_lasso <- function(data,
                            nfactors,
                            warmstart = TRUE,
                            maxiter = 50,
                            eps = 10^-4,
                            INIT = if (warmstart) "rational" else "semirational",
                            nstarts = 50,
                            lambda.length = 50,
                            standardize = TRUE) {
  INIT <- check_init(INIT)
  N <- dim(data)[2]
  I <- dim(data)[1]

  lasso_vec <- seq(0.1*N, max_lambda, length.out = 50)
  K <- length(lasso_vec)

  avec <- matrix(nrow = K, ncol = 6)
  colnames(avec) <- c( "lambda", "IS", "PEV", "Prop0", "MinNonZeroL", "MaxSDL")
  for (k in 1:K) {
    # retry loop
    success <- FALSE
    ntry <- 0
    while (!success) {
      ntry <- ntry + 1
      res <- tryCatch(
        IS_lasso(
          data = data,
          nfactors = nfactors,
          lambda = lasso_vec[k],
          maxiter = maxiter,
          eps = eps,
          INIT = INIT,
          nstarts = nstarts,
          standardize = standardize
        ),
        error = function(e) {
          message(sprintf("k = %d: IS_lasso failed (lambda = %g): %s, retrying...",
                          k, lasso_vec[k], conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(res)) {
        a       <- res
        success <- TRUE
      } else if (ntry >= 10) {
        stop(sprintf("IS_lasso failed %d times for lambda = %g; see the messages above", ntry, lasso_vec[k]))
      }
      # else it loops and tries again
    }
    avec[k, ] <- c(
      lasso_vec[k], a$value, a$vaf, a$propzero,
      a$smallestP, a$maxsdP
    )
  }

  sel_lasso <- lasso_vec[which.max(avec[,2])]

  return(sel_lasso)
}

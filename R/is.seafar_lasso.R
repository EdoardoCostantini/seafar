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

  return(is_lasso_stats(P_hat, T_hat, V_s, V_oo))
}

# IS and related statistics of one solution, given the explained variance of the
# lambda = 0 solution (V_s) and the total variance (V_oo) of the preprocessed data
is_lasso_stats <- function(P_hat, T_hat, V_s, V_oo) {
  card <- sum(P_hat != 0)

  V_a <- sum((T_hat %*% t(P_hat))^2)
  IS <- list()
  IS$value <- (V_a * V_s / V_oo^2) * (sum(round(P_hat, 3) == 0) / length(P_hat))
  IS$vaf <- V_a / V_oo
  IS$propzero <- sum(round(P_hat, 3) == 0) / length(P_hat)

  IS$smallestP <- ifelse(sum(rowSums(P_hat != 0)) < sum(card),
                         0, min(abs(P_hat[P_hat != 0]))
  )
  IS$maxsdP <- max(apply(P_hat, 2, sd))
  return(IS)
}


#' Model selection procedure with IS to choose the LASSO penalty paramter.
#'
#' The grid of candidate penalties runs up to the largest value that still leaves
#' at least 3 nonzero loadings on every factor, which is found
#' internally with \code{\link{lambda_max}}.
#'
#' With \code{warmstart = TRUE} (the default) the candidate penalties are evaluated once each, in
#' increasing order, using the loadings of the previous penalty as starting values.
#' This is much faster, but the results depend more on the initial solution. The
#' path stops at the first penalty that leaves fewer than 3 nonzero loadings on a factor.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors.
#' @param warmstart Logical. If \code{TRUE} (default), use warm starts.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loading matrix. The default is \code{"rational"} if \code{warmstart = TRUE}
#'   and \code{"semirational"} if \code{warmstart = FALSE}.
#' @param nstarts Number of starts for fitting the selected model and, if \code{warmstart = FALSE}, for model selection.
#' @param lambda.length Number of candidate lambda values to evaluate.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns A list with
#' \item{sellambda}{The selected LASSO penalty parameter.}
#' \item{selmodel}{The output of \code{seafar_lasso_multistart()} fitted with the selected penalty.}
#' \item{pve}{Proportion of variance explained by the selected model.}
#' \item{pvepca}{Proportion of variance explained by the unconstrained model (first nfactors principal components).}
#' @export
#'
#' @examples
#'\dontrun{
#' IS_result <- is.seafar_lasso(
#'   data = X,
#'   nfactors = 5,
#'   nstarts = 50)
#' summary(IS_result)
#' summary(IS_result$selmodel)
#'
#' # slower, using multiple starts for every candidate penalty
#' IS_result <- is.seafar_lasso(
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

  max_lambda <- lambda_max(
    data = data,
    nfactors = nfactors,
    maxiter = maxiter,
    eps = eps,
    INIT = INIT,
    standardize = standardize
  )

  if (max_lambda <= 0) {
    stop("no positive penalty leaves at least 3 nonzero loadings on every factor; try fewer factors")
  }

  # With few observations and many items max_lambda can fall below the usual lower bound
  lambda_min <- if (max_lambda > 0.1*N) 0.1*N else 0.1*max_lambda
  lasso_vec <- seq(lambda_min, max_lambda, length.out = lambda.length)
  K <- length(lasso_vec)

  if (warmstart) {
    avec <- is_lasso_path(
      data = data,
      nfactors = nfactors,
      lasso_vec = lasso_vec,
      maxiter = maxiter,
      eps = eps,
      INIT = INIT,
      standardize = standardize
    )
  } else {
    avec <- lasso_is_table(K)
    for (k in 1:K) {
      a <- retry_lasso(
        function() {
          IS_lasso(
            data = data,
            nfactors = nfactors,
            lambda = lasso_vec[k],
            maxiter = maxiter,
            eps = eps,
            INIT = INIT,
            nstarts = nstarts,
            standardize = standardize
          )
        },
        what = "IS_lasso",
        lambda = lasso_vec[k]
      )
      avec[k, ] <- c(
        lasso_vec[k], a$value, a$vaf, a$propzero,
        a$smallestP, a$maxsdP
      )
    }
  }

  sel_lasso <- lasso_vec[which.max(avec[,2])]

  selmodel <- retry_lasso(
    function() {
      seafar_lasso_multistart(
        data = data,
        nfactors = nfactors,
        lambda = sel_lasso,
        maxiter = maxiter,
        eps = eps,
        INIT = INIT,
        nstarts = nstarts,
        standardize = standardize
      )
    },
    what = "seafar_lasso_multistart",
    lambda = sel_lasso
  )

  svd1 <- svd(scale(data, center = selmodel$center, scale = selmodel$scale), 0, 0)

  result <- list(
    sellambda = sel_lasso,
    selmodel = selmodel,
    pve = selmodel$PVE,
    pvepca = sum(svd1$d[1:nfactors]^2) / sum(svd1$d^2)
  )
  attr(result, "class") <- "is_seafar_lasso"

  return(result)
}

lasso_is_table <- function(K) {
  avec <- matrix(nrow = K, ncol = 6)
  colnames(avec) <- c( "lambda", "IS", "PEV", "Prop0", "MinNonZeroL", "MaxSDL")
  return(avec)
}

# Evaluates the IS along an increasing sequence of penalties, each fit starting from the
# loadings of the previous one. Penalties after the path stops are left as NA.
is_lasso_path <- function(data, nfactors, lasso_vec, maxiter, eps, INIT, standardize) {
  N <- dim(data)[1]
  J <- dim(data)[2]

  # 0. Center (always) and optionally scale to unit variance (N denominator), once for the whole path
  data <- scale(data, center = TRUE, scale = FALSE)
  xscale <- rep(1, J)
  if (standardize) {
    xscale <- sqrt(colSums(data^2) / N)
    if (any(xscale < sqrt(.Machine$double.eps))) {
      stop("data contains item(s) with zero variance; remove them or use standardize = FALSE")
    }
    data <- scale(data, center = FALSE, scale = xscale)
  }

  fit0 <- retry_lasso(
    function() {
      seafar_lasso(
        data = data,
        nfactors = nfactors,
        lambda = 0,
        maxiter = maxiter,
        eps = eps,
        INIT = INIT,
        standardize = FALSE # data already preprocessed above
      )
    },
    what = "seafar_lasso",
    lambda = 0
  )

  V_oo <- sum(data^2)
  V_s <- sum((fit0$scores %*% t(fit0$loadings))^2)
  loadings <- fit0$loadings

  avec <- lasso_is_table(length(lasso_vec))
  for (k in seq_along(lasso_vec)) {
    fit <- seafar_lasso(
      data = data,
      nfactors = nfactors,
      lambda = lasso_vec[k],
      maxiter = maxiter,
      eps = eps,
      INIT = INIT,
      initloadings = loadings,
      standardize = FALSE # data already preprocessed above
    )
    loadings <- fit$loadings

    # a factor that lost its loadings cannot recover along the path
    if (any(!enough_items(loadings, min_items = 3))) {
      break
    }

    a <- is_lasso_stats(loadings, fit$scores, V_s, V_oo)
    avec[k, ] <- c(
      lasso_vec[k], a$value, a$vaf, a$propzero,
      a$smallestP, a$maxsdP
    )
  }

  if (is.na(avec[1, "IS"])) {
    stop("no candidate penalty leaves at least 3 nonzero loadings on every factor; try fewer factors or warmstart = FALSE")
  }

  return(avec)
}

# Fits can fail for an unlucky random start, so retry before giving up
retry_lasso <- function(fit, what, lambda, times = 10) {
  for (ntry in seq_len(times)) {
    res <- tryCatch(
      fit(),
      error = function(e) {
        message(sprintf("%s failed (lambda = %g): %s, retrying...",
                        what, lambda, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) {
      return(res)
    }
  }
  stop(sprintf("%s failed %d times for lambda = %g; see the messages above", what, times, lambda))
}

#' Display a summary of the results of \code{is.seafar_lasso()}.
#'
#' Use \code{summary(object$selmodel)} to display the loading matrix of the selected model.
#'
#' @param object Object of class inheriting from 'is_seafar_lasso'.
#' @param ... Argument to be passed to or from other methods.
#'
#' @returns Summary of the results.
#' @export
#'
#' @examples
#' \dontrun{
#' summary(IS_result)
#' }
summary.is_seafar_lasso <- function(object, ...) {
  cat(sprintf(
    "The LASSO penalty selected by the Index of Sparseness is: %.3f\n",
    object$sellambda
  ))

  cat(sprintf(
    "The number of nonzero loadings with this penalty is: %s\n",
    sum(round(object$selmodel$loadings, 3) != 0)
  ))

  cat(sprintf(
    "The PEV with this penalty is: %.3f\n",
    object$pve
  ))

  cat(sprintf(
    "The PEV of the unconstrained model (no zero loadings) is: %.3f\n",
    object$pvepca
  ))
}

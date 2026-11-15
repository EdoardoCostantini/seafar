#' Sparse exploratory approximate factor analysis using LASSO penalty
#'
#' Function for exploratory approximate factor analysis with a LASSO penalty on the loadings,
#' resulting in a sparse measurement model. The model is fitted from multiple starting values
#' and the solution with the lowest loss is returned.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors.
#' @param lambda LASSO penalty tuning parameter.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loadings.
#' @param initloadings Initial loading matrix if available.
#' @param nstarts Number of starts, default is 50. Only one start is run if \code{initloadings} is provided
#'   or \code{INIT = "svd"}, because every start would be identical.
#' @param show_progress Option of print progress bar, default is FALSE.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#' @param unshrink Logical. If \code{TRUE} (default), the zero pattern of the best solution is refitted with
#'   \code{lambda = 0} (warm-started from that solution) to remove the shrinkage of the nonzero loadings.
#'
#' @returns A list with
#' \item{loadings}{The best estimated loading matrix (unshrunk if \code{unshrink = TRUE}).}
#' \item{scores}{The best estimated factor score matrix (from the unshrunk refit if \code{unshrink = TRUE}).}
#' \item{Lossvec}{A list of vectors of loss values of each successful start (lasso fits).}
#' \item{Loss}{Loss value of the best start; the residual sum of squares of the refit if \code{unshrink = TRUE}.}
#' \item{PVE}{Proportion of variance explained by the returned loadings and scores (no penalty term included).}
#' \item{lambda}{The LASSO penalty tuning parameter used.}
#' \item{all_losses}{Penalized loss of every successful start (lasso fits).}
#' \item{n_best}{Number of starts that reached the lowest lasso loss (within 0.01).}
#' \item{n_distinct}{Number of distinct lasso loss values across starts.}
#' \item{converged}{Whether the returned solution converged (the refit if \code{unshrink = TRUE}).}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#' \item{loadings_lasso, scores_lasso, Loss_lasso}{If \code{unshrink = TRUE}: the original lasso solution that
#'   defined the zero pattern, and its penalized loss.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- seafar_lasso(
#'   data = X,
#'   nfactors = 3,
#'   lambda = 10,
#'   INIT = "random",
#'   nstarts = 50
#' )
#' summary(fit)
#' }
seafar_lasso <- function(data,
                         nfactors,
                         lambda,
                         maxiter = 50,
                         eps = 10^-4,
                         INIT = "semirational",
                         initloadings = NULL,
                         nstarts = 50,
                         show_progress = FALSE,
                         standardize = TRUE,
                         unshrink = TRUE) {
  INIT <- check_init(INIT)
  if (!is.null(initloadings) || INIT == "svd") {
    nstarts <- 1 # every start would be identical
  }

  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  converged <- array()
  if (show_progress == TRUE) {
    pb <- txtProgressBar(min = 0, max = nstarts, style = 3)
  }

  valid_index <- 0 # counts only successful runs
  last_error <- NULL

  for (n in 1:nstarts) {
    # catch when a start throws an error to skip that start
    result <- tryCatch(
      {
        seafar_lasso_single(
          data = data,
          nfactors = nfactors,
          lambda = lambda,
          maxiter = maxiter,
          eps = eps,
          INIT = INIT,
          initloadings = initloadings,
          standardize = standardize
        )
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        return(NULL)
      }
    )

    # Skip failed runs
    if (is.null(result)) {
      if (show_progress) setTxtProgressBar(pb, n)
      next
    }

    # Store successful run
    valid_index <- valid_index + 1
    Pout3d[[valid_index]] <- result$loadings
    Tout3d[[valid_index]] <- result$scores
    LOSS[valid_index] <- result$Residual
    LOSSvec[[valid_index]] <- 1 - result$PVE
    converged[valid_index] <- result$converged
    xcenter <- result$center # same for every start; taken from a successful run
    xscale <- result$scale

    if (show_progress) setTxtProgressBar(pb, n)
  }

  if (show_progress == TRUE) close(pb)

  if (valid_index == 0) {
    stop("all ", nstarts, " starts failed; the last error was: ", last_error)
  }
  if (valid_index < nstarts) {
    warning(nstarts - valid_index, " of ", nstarts, " starts failed and were skipped")
  }

  # check how many times the minimum loss was achieved
  best <- min(LOSS)
  tol <- 1e-2
  n_best <- sum(abs(LOSS - best) < tol)
  n_distinct <- length(unique(round(LOSS, 2)))

  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k) > 1) {
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }

  return_varselect <- list()
  return_varselect$loadings <- Pout3d[[k]]
  return_varselect$scores <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  return_varselect$Loss <- LOSS[k]
  return_varselect$all_losses <- LOSS
  return_varselect$n_best <- n_best
  return_varselect$n_distinct <- n_distinct
  return_varselect$converged <- converged[k]
  return_varselect$lambda <- lambda
  return_varselect$center <- xcenter
  return_varselect$scale <- xscale

  # undo shrinkage: refit the zero pattern of the best start with lambda = 0, warm-started from it
  if (unshrink) {
    refit <- seafar_lasso_single(data = data,
                                 nfactors = nfactors,
                                 lambda = 0,
                                 maxiter = maxiter,
                                 eps = eps,
                                 initloadings = Pout3d[[k]],
                                 standardize = standardize,
                                 pattern = Pout3d[[k]])

    return_varselect$loadings <- refit$loadings
    return_varselect$scores <- refit$scores
    return_varselect$Loss <- refit$Residual
    return_varselect$loadings_lasso <- Pout3d[[k]] # keep the shrunk solution the pattern came from
    return_varselect$scores_lasso <- Tout3d[[k]]
    return_varselect$Loss_lasso <- LOSS[k]
    return_varselect$converged <- refit$converged
  }

  # computed directly because the PVE of a lasso fit has the penalty included in its loss
  data_scaled <- scale(data, center = return_varselect$center, scale = return_varselect$scale)
  return_varselect$PVE <- 1 - ssres(data_scaled, return_varselect$scores, return_varselect$loadings) / sum(data_scaled^2)

  attr(return_varselect, "class") <- "seafar_lasso"

  return(return_varselect)
}

#' Single start of sparse exploratory approximate factor analysis using LASSO penalty
#'
#' @param data A data matrix containing the dataset.
#' @param nfactors The number of factors.
#' @param lambda LASSO penalty tuning parameter.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param INIT Method to initialize loadings.
#' @param initloadings Optional user-supplied initial loading matrix.
#' @param standardize Logical. If \code{TRUE} (default) variables are scaled to unit variance after centering.
#'    Centering is applied regardless.
#' @param pattern Optional JxQ matrix; loadings where \code{pattern} is 0 (or \code{FALSE}) are fixed at zero.
#'    With \code{lambda = 0} this refits a given zero pattern without shrinkage (used by \code{unshrink} in
#'    \code{seafar_lasso()}).
#'
#' @returns
#' \item{scores}{The estimated factor score matrix.}
#' \item{loadings}{The estimated loading matrix.}
#' \item{PVE}{A vector of 1 - penalized loss in each iteration of the AO procedure.}
#' \item{Residual}{The penalized loss of the solution.}
#' \item{converged}{Whether the AO procedure converged before \code{maxiter}.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @examples
#' \dontrun{
#' fit <- seafar_lasso_single(ddata = X, nfactors = 3, lambda = 10, INIT = "svd")
#' }
seafar_lasso_single <- function(data,
                                nfactors,
                                lambda,
                                maxiter = 50,
                                eps = 10^-4,
                                INIT = "semirational",
                                initloadings = NULL,
                                standardize = TRUE,
                                pattern = NULL) {
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

  converged <- FALSE
  ssx <- sum(data^2)
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc # Used to compute convergence criterium
  # 1. Initialize loading matrix
  if (is.null(initloadings)) {
    loadings <- seafar_init_l1(data, nfactors, INIT)
  } else {
    loadings <- initloadings
  }

  svd1 <- svd(data, nfactors, nfactors)
  scores <- sqrt(N) * svd1$u
  diffT <- 0
  diffP <- 0
  while (convAO == 0) {
    iter0 <- 1
    Losst <- 1
    Lossvec0 <- c()
    stopcritT0 <- 0
    # Lossvec1 <- 1
    # 2.1. Update factor scores
    while (stopcritT0 == 0) {
      Lossu1old <- ssres(data, scores, loadings) / ssx
      for (q in 1:nfactors) {
        E <- data - scores %*% t(loadings)
        Er <- E + scores[, q] %*% t(loadings[, q])
        num <- Er %*% loadings[, q]
        if (sum(num^2) > 0) {
          scores[, q] <- sqrt(N) * num / sqrt(sum(num^2))
        }
      }
      # t(scores)%*%scores
      # Calculate loss
      Lossu0 <- LOSS(data, scores, loadings, lambda) / ssx
      Lossvec0 <- c(Lossvec0, Lossu0)
      # check convergence
      if (iter0 > maxiter) {
        stopcritT0 <- 1
      }
      if (abs(Losst - Lossu0) < eps) {
        stopcritT0 <- 1
      }
      iter0 <- iter0 + 1
      Losst <- Lossu0
    }
    # Loss <- ssres(DATA, scores, loadings)/ssx
    Lossu <- LOSS(data, scores, loadings, lambda) / ssx

    # 2. Update loadings

    for (q in 1:nfactors) {
      E <- data - scores %*% t(loadings)
      Er <- E + scores[, q] %*% t(loadings[, q])
      crosstEr <- t(Er) %*% scores[, q]
      loadings[, q] <- sign(crosstEr) * apply(cbind(abs(crosstEr) - lambda / 2, 0), 1, max) / N
    }
    if (!is.null(pattern)) {
      loadings[pattern == 0] <- 0 # fixed zeros
    }

    # Calculate loss
    Lossu <- LOSS(data, scores, loadings, lambda) / ssx
    Lossvec <- c(Lossvec, Lossu)
    if (iter > maxiter) {
      convAO <- 1
    }
    # if (Lossc-Lossu < -1e-12) {
    #  warning('Increase in Loss')
    #  break
    # }
    if (abs(Lossc - Lossu) < eps) {
      convAO <- 1
    }
    iter <- iter + 1
    Lossc <- Lossu
  }
  if (iter < maxiter) {
    converged <- TRUE
  }

  result <- list("scores" = scores, "loadings" = loadings, "PVE" = 1 - Lossvec, "Residual" = Lossu * ssx, "converged" = converged, "center" = xcenter, "scale" = xscale)
  return(result)
}

#' Display a summary of the results of \code{seafar_lasso()}.
#'
#' @param object Object of class inheriting from 'seafar_lasso'.
#' @param disp The default is \code{"loadings"} which returns only the estimated loadings matrix.
#'             Otherwise if \code{"full"}, function returns both the estimated loadings and factor scores matrices.
#' @param nrow Number of rows of the loading matrix to be displayed, default is 10.
#' @param ...  Argument to be passed to or from other methods.
#' @export
#' @examples
#' \dontrun{
#' summary(fit)
#' }
summary.seafar_lasso <- function(object, disp = "loadings", nrow = 10, ...) {
  if (missing(disp)) {
    disp <- "loadings"
  }

  if (disp == "loadings") {
    cat(sprintf(
      "\nThe LASSO penalty is: %.3f\n",
      object$lambda
    ))
    cat(sprintf(
      "The number of nonzero loadings is: %s\n",
      sum(round(object$loadings, 3) != 0)
    ))
    cat(sprintf(
      "The PEV is: %.3f\n",
      object$PVE
    ))
    cat(sprintf("\nThe estimated loadings matrix is \n"))
    print(utils::head(round(object$loadings, 3), nrow))
  } else if (disp == "full") {
    cat(sprintf(
      "\nThe LASSO penalty is: %.3f\n",
      object$lambda
    ))
    cat(sprintf(
      "The number of nonzero loadings is: %s\n",
      sum(round(object$loadings, 3) != 0)
    ))
    cat(sprintf(
      "The PEV is: %.3f\n",
      object$PVE
    ))

    cat(sprintf("\nThe estimated loadings matrix is \n"))
    print(utils::head(round(object$loadings, 3), nrow))

    cat(sprintf("\nThe estimated factor scores matrix is \n"))
    print(object$scores)
  }
}

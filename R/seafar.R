#' Sparse exploratory approximate factor analysis
#'
#' Function for exploratory approximate factor analysis resulting in a sparse
#' measurement model
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors
#' @param C The number of nonzero loadings
#' @param eps Convergence criterion based on difference in loss between iterates
#' @param maxiter Maximum number of iterations of the AO procedure
#' @param initloadings Initial loading matrix if available.
#' @param INIT Method to initialize loadings
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#' @return Factor loading and factor score matrices, plus the \code{center} and \code{scale} used to transform the data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' big5_seafa <- seafar(USArrests, 2, 4, INIT = "svd", orthogonal = TRUE)
#' }
seafar <- function(data,
                   nfactors,
                   C,
                   eps = 1e-4,
                   maxiter = 50,
                   initloadings = NULL,
                   INIT,
                   INIT = "semirational",
                   orthogonal = FALSE,
                   nstarts = 50,
                   show_progress = FALSE,
                   standardize = TRUE) {
  if (orthogonal == TRUE) {
    result <- seafar_orthogonal(
      data = data,
      nfactors = nfactors,
      C = C,
      eps = eps,
      maxiter = maxiter,
      initloadings = initloadings,
      INIT = INIT,
      standardize = standardize
    )
  } else {
    result <- seafar_general(
      data = data,
      nfactors = nfactors,
      C = C,
      eps = eps,
      maxiter = maxiter,
      initloadings = initloadings,
      INIT = INIT,
      standardize = standardize
  INIT <- check_init(INIT)
  if (!is.null(initloadings) || INIT == "svd") {
    nstarts <- 1 # every start would be identical
  }
    )
  }
  # 3. Return output
  attr(result, "class") <- "seafar"
  return(result)
}

#' Function for exploratory approximate factor analysis resulting in a sparse
#' measurement model for orthogonal factors.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param C Number of nonzero loadings.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param initloadings Initial loading matrix if available.
#' @param INIT Method to initialize loadings.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns
#' \item{loadings}{The best estimated loading matrix.}
#' \item{scores}{The best estimated factor score matrix.}
#' \item{PVE}{A vector of PVE in each iteration of the AO procedure.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' big5_seafa <- seafar(ocean, 5, 240, INIT = "svd")
#' }
seafar_orthogonal <- function(data,
                              nfactors,
                              C,
                              eps = 1e-4,
                              maxiter = 50,
                              initloadings = NULL,
                              INIT = "semirational",
                              standardize = TRUE) {
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
  ssx <- sum(data^2)

  stopcrit <- 0
  Lossvec <- c()
  iter <- 1

  verbose <- FALSE ##added to make future debugging easier

  tdata <- t(data)
  # 1. Initialize loading matrix
  if (is.null(initloadings)) {
    loadings <- seafar_init(data, nfactors, C, INIT)
  } else {
    loadings <- initloadings
  }
  Lossc <- 1

  if (length(C) > 1) {
    C_c <- J - C

    while (stopcrit == 0) {

      # 2.1 Update factor scores
      scores <- orthprocr(data, loadings)
      A <- tdata %*% scores / N
      # check monotonicity AO: can be removed later on
      Lossu_H <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update H: Diff loss ', Lossc-Lossu_H)
      }

      # 2.2 Update factor loadings
      if (sum(C_c) == 0) {
        loadings <- A
      } else {
        for (q in 1:nfactors) {
          ind <- sort(abs(A[, q]), index.return = TRUE)
          A[ind$ix[seq_len(C_c[q])], q] <- 0
          loadings <- A
        }
      }

      # 2.3 Check stopping criteria
      # Calculate loss
      Lossu <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update L: Diff loss ', Lossu_H-Lossu)
      }
      Lossvec <- c(Lossvec, Lossu)
      if (iter > maxiter) {
        stopcrit <- 1
      }

      if (Lossc - Lossu < eps) {
        stopcrit <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  } else {
    C_c <- J * nfactors - C
    # 2. Alternating optimization scheme
    while (stopcrit == 0) {
      # 2.1 Update factor scores
      scores <- orthprocr(data, loadings)
      A <- tdata %*% scores / N
      # check monotonicity AO: can be removed later on
      Lossu <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update H: loss ', Lossu)
      }
      if (verbose){
        message('Iter ', iter, ' Update H: Diff loss ', Lossc-Lossu)
      }

      # 2.2 Update factor loadings
      if (C_c == 0) {
        loadings <- A
      } else {
        loadings <- hardthr(A, C_c)
      }

      # 2.3 Check stopping criteria
      # Calculate loss
      Lossu <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update L: Diff loss ', Lossc-Lossu)
      }
      Lossvec <- c(Lossvec, Lossu)
      if (iter > maxiter) {
        stopcrit <- 1
      }

      if (Lossc - Lossu < eps) {
        stopcrit <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  }

  # 3. Return output
  result <- list("scores" = scores, "loadings" = loadings, "PVE" = 1 - Lossvec, "Residual" = Lossu * ssx, "center" = xcenter, "scale" = xscale)
}

#' Function for exploratory approximate factor analysis resulting in a sparse
#' measurement model for non-orthogonal factors.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors Number of factors.
#' @param C Number of nonzero loadings.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param initloadings Initial loading matrix if available.
#' @param INIT Method to initialize loadings.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#'
#' @returns
#' \item{loadings}{The best estimated loading matrix.}
#' \item{scores}{The best estimated factor score matrix.}
#' \item{PVE}{A vector of PVE in each iteration of the AO procedure.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' big5_seafa <- seafar(ocean, 5, 240, INIT = "svd")
#' }
seafar_general <- function(data,
                           nfactors,
                           C,
                           eps = 1e-4,
                           maxiter = 50,
                           initloadings = NULL,
                           INIT = "semirational",
                           standardize = TRUE) {
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

  N <- dim(data)[1]
  J <- dim(data)[2]
  ssx <- sum(data^2)
  svd1 <- svd(data, nfactors, nfactors)
  scores <- sqrt(N)*svd1$u

  # 1. Initialize loading matrix
  if (is.null(initloadings)) {
    loadings <- seafar_init(data, nfactors, C, INIT)
  } else {
    loadings <- initloadings
  }

  verbose <- FALSE
  stopcrit <- 0
  Lossvec <- c()
  Lossc <- Inf # never stop at iteration 1 (random starts can have loss > 1; warm starts can be less sparse than C)
  iter <- 1

  # 2. Alternating optimization scheme
  tdata <- t(data)
  if (length(C) > 1) {
    C_c <- J - C
    while (stopcrit == 0) {

      # 2.1 Update factor scores
      scores <- oblprocr(data, scores, loadings, maxiter, eps)
      Lossu_H <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update H: loss ', Lossu_H)
      }
      if (verbose){
        message('Iter ', iter, ' Update H: Diff loss ', Lossc-Lossu_H)
      }
      #Lossc <- Lossu
      # if (Lossc-Lossu < -1e-12) {
      #  warning('Increase in Loss: update scores')
      #  break
      # }

      # 2.2 Update factor loadings
      TT <- t(scores) %*% scores
      svdTT <- svd(TT)
      alpha <- svdTT$d[1]
      A <- loadings - (loadings %*% TT - tdata %*% scores) / (alpha)
      if (sum(C_c) == 0) {
        loadings <- A
      } else {
        for (q in 1:nfactors) {
          ind <- sort(abs(A[, q]), index.return = TRUE)
          A[ind$ix[seq_len(C_c[q])], q] <- 0
          loadings <- A
        }
      }

      # 2.3 Check stopping criteria
      # Calculate loss
      Lossu <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update L: Diff loss ', Lossu_H-Lossu)
      }
      Lossvec <- c(Lossvec, Lossu)
      if (iter > maxiter) {
        stopcrit <- 1
      }
      # if (Lossc-Lossu < -1e-12) {
      #  warning('Increase in Loss')
      #  break
      # }
      if (Lossc - Lossu < eps) {
        stopcrit <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  } else {
    C_c <- J * nfactors - C
    while (stopcrit == 0) {

      # 2.1. Update factor scores
      scores <- oblprocr(data, scores, loadings, maxiter, eps)
      Lossu_H <- ssres(data, scores, loadings) / ssx

      if (verbose){
        message('Iter ', iter, ' Update H: Diff loss ', Lossc-Lossu_H)
      }

      # 2.2 Update factor loadings
      TT <- t(scores) %*% scores
      svdTT <- svd(TT)
      alpha <- svdTT$d[1]
      A <- loadings - (loadings %*% TT - tdata %*% scores) / (alpha)
      if (C_c == 0) {
        loadings <- A
      } else {
        loadings <- hardthr(A, C_c)
      }

      # 2.3 Check stopping criteria
      # Calculate loss
      Lossu <- ssres(data, scores, loadings) / ssx
      if (verbose){
        message('Iter ', iter, ' Update L: Diff loss ', Lossu_H-Lossu)
      }
      Lossvec <- c(Lossvec, Lossu)
      if (iter > maxiter) {
        stopcrit <- 1
      }
      if (Lossc - Lossu < -1e-12) {
        warning("Increase in Loss")
        # break
      }
      if (Lossc - Lossu < eps) {
        stopcrit <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  }

  # 3. Return output
  result <- list("scores" = scores, "loadings" = loadings, "PVE" = 1 - Lossvec, "Residual" = Lossu * ssx, "center" = xcenter, "scale" = xscale)
}

#' Display a summary of the results of \code{seafar()}.
#'
#' @param object Object of class inheriting from 'seafar'.
#' @param disp The default is \code{"loadings"} which returns only the estimated loadings matrix.
#'             Otherwise if \code{"full"}, function returns both the estimated loadings and factor scores matrices.
#' @param ...  Argument to be passed to or from other methods.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- seafar(X)
#' summary(results, disp = "full")
#' }
summary.seafar <- function(object, disp = "loadings", ...) {
  if (is.null(disp)) {
    disp <- "loadings"
  }

  if (disp == "loadings") {
    cat(sprintf(
      "\nThe number of nonzero loadings is: %s\n",
      sum(object$loadings != 0)
    ))

    cat(sprintf("\nThe estimated loadings matrix is \n"))
    print(utils::head(round(object$loadings, 3), 10))
  } else {
    cat(sprintf(
      "\nThe number of nonzero loadings is: %s\n",
      sum(object$loadings != 0)
    ))

    cat(sprintf("\nThe estimated loadings matrix is \n"))
    print(utils::head(round(object$loadings, 3), 10))

    cat(sprintf("\nThe estimated factor scores matrix is \n"))
    print(object$scores)
  }
}

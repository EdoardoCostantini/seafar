#' Sparse exploratory approximate factor analysis
#'
#' Function for exploratory approximate factor analysis resulting in a sparse
#' measurement model. The model is fitted from multiple starting values and the
#' solution with the lowest loss is returned.
#'
#' @param data A data frame or matrix (NxJ). Data are always centered internally.
#' @param nfactors The number of factors
#' @param C The number of nonzero loadings
#' @param eps Convergence criterion based on difference in loss between iterates
#' @param maxiter Maximum number of iterations of the AO procedure
#' @param initloadings Initial loading matrix if available.
#' @param INIT Method to initialize loadings
#' @param orthogonal Orthogonal or non-orthogonal factors, default is FALSE.
#' @param nstarts Number of starts, default is 50. Only one start is run if \code{initloadings} is provided
#'   or \code{INIT = "svd"}, because every start would be identical.
#' @param show_progress Option of print progress bar, default is FALSE.
#' @param standardize Logical. If \code{TRUE} (default), items are scaled to unit variance (N denominator)
#'   after centering. Centering is applied regardless.
#' @returns A list with
#' \item{scores}{The best estimated factor score matrix.}
#' \item{loadings}{The best estimated loading matrix.}
#' \item{PVE}{A vector of PVE in each iteration of the AO procedure of the best start.}
#' \item{Residual}{The residual sum of squares of the best start.}
#' \item{center}{Item means used to center the data.}
#' \item{scale}{Item standard deviations used to scale the data (1s if \code{standardize = FALSE}).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' big5_seafa <- seafar(USArrests, 2, 4, INIT = "random", orthogonal = TRUE, nstarts = 50)
#' }
seafar <- function(data,
                   nfactors,
                   C,
                   eps = 1e-4,
                   maxiter = 50,
                   initloadings = NULL,
                   INIT = "semirational",
                   orthogonal = FALSE,
                   nstarts = 50,
                   show_progress = FALSE,
                   standardize = TRUE) {
  INIT <- check_init(INIT)
  if (!is.null(initloadings) || INIT == "svd") {
    nstarts <- 1 # every start would be identical
  }
  fit_fun <- if (orthogonal == TRUE) seafar_orthogonal else seafar_general

  Pout3d <- list()
  Hout3d <- list()
  LOSS <- array()
  PVE <- list()
  if (show_progress == TRUE) {
    pb <- txtProgressBar(min = 0, max = nstarts, style = 3)
  }

  valid_index <- 0 # counts only successful runs
  last_error <- NULL

  for (n in 1:nstarts) {
    # catch when a start throws an error to skip that start
    result <- tryCatch(
      {
        fit_fun(
          data = data,
          nfactors = nfactors,
          C = C,
          eps = eps,
          maxiter = maxiter,
          initloadings = initloadings,
          INIT = INIT,
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
    Hout3d[[valid_index]] <- result$scores
    LOSS[valid_index] <- result$Residual
    PVE[[valid_index]] <- result$PVE
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

  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k) > 1) {
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }

  result <- list(
    "scores" = Hout3d[[k]],
    "loadings" = Pout3d[[k]],
    "PVE" = PVE[[k]],
    "Residual" = LOSS[k],
    "center" = xcenter,
    "scale" = xscale
  )
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

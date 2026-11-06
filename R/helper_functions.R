#' Function for hard thresholding.
#'
#' @param loadings Loadings matrix JxR.
#' @param C_c Number of zero loadings.
#'
#' @returns Loadings with exactly C_c nonzero loadings.
#' @export
#' @author Katrijn Van Deun
#'
#' @examples
#' \dontrun{
#' hardthr(loadings, C_c)
#' }
hardthr <- function(loadings, C_c) {
  ind <- sort(abs(loadings), index.return = TRUE) # without replacement
  loadings[ind$ix[1:C_c]] <- 0
  # loadings <- as(loadings, "sparseMatrix")

  loadings
}


#' Loss function.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param scores Factor scores matrix NxQ.
#' @param loadings Loadings matrix JxQ.
#'
#' @returns Sum of squared residuals.
#' @export
#' @author Katrijn Van Deun.
#'
#' @examples
#' \dontrun{
#' ssres(X, scores, loadings)
#' }
ssres <- function(data, scores, loadings) {
  XHAT <- scores %*% t(loadings)
  res <- sum((XHAT - data)^2)
  return(res)
}


#' Function to calculate the solution to the orthogonal procrustes problem.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param loadings Loadings matrix JxQ.
#'
#' @returns Orthogonal factor scores matrix NxJ.
#' @export
#' @author Katrijn Van Deun.
#'
#' @examples
#' \dontrun{
#' orthprocr(X, loadings)
#' }
orthprocr <- function(data, loadings) {
  N <- dim(data)[1]
  Q <- dim(loadings)[2]
  XP <- data %*% loadings / (sqrt(N))
  svdXP <- svd(XP, Q, Q)
  scores <- sqrt(N) * svdXP$u %*% t(svdXP$v)

  scores
}

#' Function to calculate the solution to the olique procrustes problem.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param loadings Loadings matrix JxQ.
#' @param maxiter Maximum number of iterations.
#' @param eps Tolerance for difference in loss.
#'
#' @returns Orthogonal factor scores matrix NxJ.
#' @export
#' @author Katrijn Van Deun.
#'
#' @examples
#' \dontrun{
#' orthprocr(X, loadings)
#' }
oblprocr <- function(data, scores, loadings, maxiter, eps){

  verbose <- FALSE

  N <- dim(data)[1]
  Q <- dim(loadings)[2]
  ssx <- sum(data^2)

  stopcritT0 <- 0
  iter0 <- 1
  Losst_c <- ssres(data, scores, loadings) / ssx

  XL <- data %*% loadings  #cache
  LL <- t(loadings) %*% loadings #cache

  if (verbose){
    message('Initial loss H ', Losst_c)
  }

  while (stopcritT0 == 0) {
    scoresLL <- scores %*% LL
    for (q in 1:Q) {
      oldscores <- scores[, q]
      num <- XL[, q] - scoresLL[, q] + scores[, q] * LL[q, q]
      if (sum(num^2) > 0) {
        scores[, q] <- sqrt(N) * num / sqrt(sum(num^2))
      }
      delta <- scores[, q] - oldscores
      scoresLL <- scoresLL + delta %*% LL[q, , drop = FALSE]
    }
    # t(scores)%*%scores
    # Calculate loss
    Losst_u <- ssres(data, scores, loadings) / ssx
    if (verbose){
      message('IterH ', iter0, ' Update H: Diff loss ', Losst_c-Losst_u)
    }

    # check convergence
    if (iter0 > maxiter) {
      stopcritT0 <- 1
    }
    if (abs(Losst_c - Losst_u) < eps) {
      stopcritT0 <- 1
    }

    iter0 <- iter0 + 1
    Losst_c <- Losst_u
  }
  scores
}

#' Function to calculate the solution to the olique procrustes problem.
#' This version is slightly less efficient.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param loadings Loadings matrix JxQ.
#' @param maxiter Maximum number of iterations.
#' @param eps Tolerance for difference in loss.
#'
#' @returns Orthogonal factor scores matrix NxJ.
#' @export
#'
#' @examples
#' \dontrun{
#' orthprocr(X, loadings)
#' }
oblprocr_old <- function(data, scores, loadings, maxiter, eps){
  N <- dim(data)[1]
  Q <- dim(loadings)[2]
  ssx <- sum(rowSums(data^2))
  iter0 <- 1
  Losst <- 1
  Lossvec0 <- c()
  stopcritT0 <- 0
  while (stopcritT0 == 0) {
    Lossu1old <- ssres(data, scores, loadings) / ssx
    for (q in 1:Q) {
      E <- data - scores %*% t(loadings)
      Er <- E + scores[, q] %*% t(loadings[, q])
      num <- Er %*% loadings[, q]
      scores[, q] <- sqrt(N) * num / sqrt(sum(num^2))
    }
    # t(scores)%*%scores
    # Calculate loss
    Lossu0 <- ssres(data, scores, loadings) / ssx
    message('Iter ', iter0, ' Update H: Diff loss ', Lossu1old-Lossu0)
    Lossu1old <- Lossu0
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
  scores
}


#' Function for obtaining an initial loading matrix.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param nfactors The number of factors Q.
#' @param C The number of zero loadings.
#' @param INIT Methods to initialize loadings c("svd", "random", "rational", "semirational").
#' Rational start uses the result of Ding and He (2004) to partition the
#' variables into nfactors clusters.
#' Mixed start combines random (20%) and SVD (80%).
#'
#' @importFrom stats kmeans rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @return A matrix of initial loadings.
#' @examples
#' \dontrun{
#' initialP <- seafar_init(ocean, 5, 240, INIT = "svd")
#' }
#' @author Tra Le and Katrijn Van Deun

seafar_init <- function(data,
                        nfactors,
                        C,
                        INIT) {
  INIT <- check_init(INIT)
  n <- dim(data)[1]
  J <- dim(data)[2]
  svd1 <- svd(data, nfactors, nfactors)
  P1 <- matrix(rnorm(J * nfactors), ncol = nfactors, nrow = J)
  P2 <- svd1$v %*% diag(svd1$d[1:nfactors]) / sqrt(n)
  if (INIT == "svd") {
    P <- P2
  } else if (INIT == "random") {
    P <- P1
  } else if (INIT == "rational") {
    P <- P2
    # 1.2 Clustering of pca loadings
    c <- kmeans(P, nfactors, nstart = 5)
    TARGET <- matrix(0, nrow = J, ncol = nfactors)
    for (q in 1:nfactors) {
      TARGET[c$cluster == q, q] <- 1
    }
    # 1.3 Rotate pca solution to cluster structure
    WEIGHTS <- TARGET
    B <- pstr(P, TARGET, WEIGHTS, 50, 1e-4)
    P <- P %*% B$Bmatrix
  } else {
    P <- P2 * 0.8 + P1 * 0.2
  }
  if (length(C) > 1) {
    C_c <- J - C
    for (q in 1:nfactors) {
      ind <- sort(abs(P[, q]), index.return = TRUE)
      P[ind$ix[seq_len(C_c[q])], q] <- 0
    }
  } else {
    C_c <- J * nfactors - C
    if (C_c == 0) {
      P <- P
    } else {
      ind <- sort(abs(P), index.return = TRUE)
      P[ind$ix[1:C_c]] <- 0
    }
  }
  loadings <- P
  return(loadings)
}

#' Partially specified target rotation.
#'
#'
#' @param loadings Matrix that needs to be rotated.
#' @param target Target towards which to rotated.
#' @param weights elementwise weights for rotation.
#' @param maxiter Maximum number of iterations of the AO procedure.
#' @param eps Convergence criterion based on difference in loss between iterates.
#' @return Rotation matrix and loss of the pstr criterion.
#' @examples
#' \dontrun{
#' ocean_pstr <- pstr(pcaloadings, target, weights)
#' }
pstr <- function(loadings, target, weights, maxiter = 50, eps = 1e-4) {
  n <- dim(loadings)[1]
  m <- dim(loadings)[2]

  L <- array()
  Bmat <- list()
  REFL <- reflexmat(m) # requires reflexmat.R

  for (i in 1:dim(REFL)[1]) {
    k <- which(REFL[i, ] == -1)
    Binit <- diag(m)
    Binit[, k] <- -1 * Binit[, k]

    B1 <- t(loadings) %*% loadings # ask katrijn
    alpha <- max(eigen(B1)$values)
    iter <- 1

    stop <- 0

    Bcurrent <- Binit
    Lossc <- pstrLoss(Binit, loadings, target, weights)

    while (stop == 0) {
      Pw <- weights * target + loadings %*% Bcurrent - weights * (loadings %*% Bcurrent)
      A <- -2 * t(Pw) %*% loadings
      Fmat <- A + 2 * t(Bcurrent) %*% t(B1) - 2 * alpha * t(Bcurrent)
      F_svd <- svd(-Fmat)
      B <- F_svd$v %*% t(F_svd$u)

      if (iter == maxiter) {
        stop <- 1
      }

      Loss <- pstrLoss(B, loadings, target, weights)
      Diff <- Lossc - Loss

      if (abs(Diff) < eps) {
        stop <- 1
      }

      iter <- iter + 1
      Lossc <- Loss
      Bcurrent <- B
    }

    L[i] <- Lossc
    Bmat[[i]] <- Bcurrent
  }

  k <- which(L == min(L))
  Loss <- L[k[1]]
  B <- Bmat[[k[1]]]


  results <- list()
  results$Bmatrix <- B
  results$Loss <- Loss
  return(results)
}


#' Function to construct matrix of reflections.
#'
#' @param m Number of factors.
#'
#' @returns Matrix of reflections.
#' @export
#'
#' @author Katrijn Van Deun and Zhengguo Gu.
#' @examples
#' \dontrun{
#' matrix <- reflexmat(5)
#' }
reflexmat <- function(m) {
  mat <- rep(1, m)

  for (i in 1:(m - 1)) {
    B <- utils::combn(1:m, i)

    for (j in 1:dim(B)[2]) {
      v <- rep(1, m)
      v[t(B[, j])] <- -1

      mat <- rbind(mat, v)
    }
  }

  return(mat)
}


#' Function to calculate the objective function associated with pstr.
#'
#' @param B Rotation matrix.
#' @param Tmat Original loading matrix.
#' @param target Target matrix for rotation.
#' @param weights Weight.
#'
#' @returns Loss value.
#' @export
#'
#' @examples
#' \dontrun{
#' loss <- pstrLoss(B, Tmat, target, weights)
#' }
pstrLoss <- function(B, Tmat, target, weights) {
  DEV <- Tmat %*% B - target
  wDEV <- weights * DEV
  Loss <- sum(wDEV^2)

  return(Loss)
}

#' Title Loss function for LASSO
#'
#' @param data Data matrix of standardized items NxJ.
#' @param scores Factor scores matrix NxQ.
#' @param loadings Loadings matrix JxQ.
#' @param lambda LASSO penalty tuning parameter.
#'
#' @returns Loss value when using LASSO penalty.
#' @export
#'
#' @examples
#' \dontrun{
#' LOSS(X, scores, loadings, lambda)
#' }
LOSS <- function(data, scores, loadings, lambda) {
  XHAT <- scores %*% t(loadings)
  res <- sum(rowSums((XHAT - data)^2))
  penalty <- sum(abs(loadings))
  loss <- res + lambda * penalty
  return(loss)
}

#' Function for obtaining an initial loading matrix for LASSO penalty.
#'
#' @param data Data matrix of standardized items NxJ.
#' @param nfactors Number of factors.
#' @param INIT Methods to initialize loadings c("svd", "random", "rational", "semirational").
#'
#' @return A matrix of initial loadings.
#' @examples
#' \dontrun{
#' initialP_lasso <- seafar_init_l1(ocean, 5, INIT = "svd")
#' }
seafar_init_l1 <- function(data, nfactors, INIT) {
  INIT <- check_init(INIT)
  n <- dim(data)[1]
  J <- dim(data)[2]
  if (INIT == "rational") {
    # no zero loadings to start from: the penalty decides which ones become zero
    return(seafar_init(data, nfactors, C = J * nfactors, INIT = INIT))
  }
  svd1 <- svd(data, nfactors, nfactors)
  P1 <- matrix(rnorm(J * nfactors), ncol = nfactors, nrow = J)
  P2 <- svd1$v %*% diag(svd1$d[1:nfactors]) / sqrt(n)
  if (INIT == "svd") {
    P <- P2
  } else if (INIT == "random") {
    P <- P1
  } else {
    P <- P2 * 0.8 + P1 * 0.2
  }

  loadings <- P
  return(loadings)
}


#' Function to validate INIT argument
#'
#' @param INIT Method to initialize loading matrix.
#'
#' @returns The correct INIT value.
#'
#' @examples
#' \dontrun{
#' check_init(INIT = 'svd')
#' }
check_init <- function(INIT) {
  methods <- c("semirational", "svd", "random", "rational")
  if (!is.character(INIT) || length(INIT) != 1 || !INIT %in% methods) {
    stop("INIT must be one of: ", paste0('"', methods, '"', collapse = ", "), call. = FALSE)
  }
  INIT
}

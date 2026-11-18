#' Parallel analysis based on SVD decomposition
#'
#' This function samples an \code{N} x \code{J} data matrix \code{iterations}
#' times and extracts its eigenvalues. 
#' The eigenvalues are obtained via the SVD decomposition of the data matrix 
#' rather than eigen-decomposition of the correlation or covariance matrix of 
#' the data matrix. 
#' These eigenvalues are equivalent to those obtained with eigen-decomposition 
#' of the correlation matrix when N > J, but it allows the program to handle 
#' N < J cases.
#' When N < J, the \code{base::svd()} function used to perform the SVD decomposition
#' obtains as many eigenvalues as the matrix rank allows.
#' The function returns the eigenvalues and their summary statistics. 
#' As the number of iterations grows, the mean, standard deviations and values 
#' of the distribution of eigenvalues.
#'
#' @param N The integer number of rows for the data matrix
#' @param J The integer number of independent columns for the data matrix
#' @param iterations The integer number of data matrices with independent columns that will be generated
#' @param quantile The proportion for the quantile of interest. This should be a single value in the interval \[0, 1\].
#'
#' @return
#' \item{sampled_eigenvalues}{ An \code{iterations} x \code{J} matrix of sampled
#' eigenvalues from the eigenvalues distribution }
#' \item{summary_statistics_eigenvalues}{ A \code{data.frame} of summary statistics of
#' the sampled eigenvalues.}
#' \item{summary_statistics_eigenvalues$mean}{ Mean of the sampled eigenvalues }
#' \item{summary_statistics_eigenvalues$sd}{ Standard deviation of the sampled
#' eigenvalues }
#' \item{summary_statistics_eigenvalues$quantile}{ The eigenvalue that leaves the \code{quantile} proportion of the sampled eigenvalues to its left. }
#' @export
#'
#' @examples
#' parallel_analysis(N = 100, J = 5, iterations = 5, quantile = 0.05)
parallel_analysis <- function(
  N = 10,
  J = 5,
  iterations = 100,
  quantile = 0.05
) {
  matrix_rank <- min(N, J)
  # Because the base::svd functions caps the number of eigenvalues to the rank
  # you will only be able to store at most `rank` eigenvalues
  sampled_eigenvalues <- base::matrix(data = NA, nrow = iterations, ncol = matrix_rank)
  for (repetition in 1:iterations) {
    sampled_independent_variables <- base::matrix(
      data = stats::rnorm(N * J),
      nrow = N,
      ncol = J
    )
    uSVt <- base::svd(sampled_independent_variables)
    # Divide by N-1 to scale the eigenvalues the same way as if they came from
    # the correlation matrix
    sampled_eigenvalues[repetition, ] <- uSVt$d^2 / (N - 1)
  }

  sampled_eigenvalues_quantiles <- base::apply(
    X = sampled_eigenvalues,
    MARGIN = 2,
    FUN = stats::quantile,
    probs = quantile
  )
  sampled_eigenvalues_mean <- base::colMeans(sampled_eigenvalues)
  sampled_eigenvalues_sd <- base::apply(
    X = sampled_eigenvalues,
    MARGIN = 2,
    FUN = sd
  )
  summary_statistics_eigenvalues <- base::data.frame(
    mean = sampled_eigenvalues_mean,
    sd = sampled_eigenvalues_sd,
    quantile = sampled_eigenvalues_quantiles
  )
  output <- base::list(
    sampled_eigenvalues = sampled_eigenvalues,
    summary_statistics_eigenvalues = summary_statistics_eigenvalues
  )

  return(output)
}

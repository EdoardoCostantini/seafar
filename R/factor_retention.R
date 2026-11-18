#' Function to retain number of factors/components using various methods.
#'
#' The function provides the number of principal components of an input data matrix to 
#' retain based on four different criteria. The use of the SVD decomposition to 
#' obtain eigenvalues makes it possible to process cases where \code{J} > \code{N}.
#' 
#' @param data An \code{N} x \code{J} matrix or data frame of observed variables.
#' @param iterations The integer number of data matrices with independent columns 
#' that will be generated for the parallel analysis criteria
#' @param quantile The proportion for the quantile of interest for the parallel 
#' analysis criteria. This should be a single value in the interval \[0, 1\].
#' 
#' @returns A list of containing the number of principal components of 
#' \code{data} retained according to four different criteria:
#' \item{kaiser}{Kaiser rule}
#' \item{parallel}{Horn's parallel analysis}
#' \item{acceleration_factor}{Acceleration factor}
#' \item{optimal_coordinates}{Optimal coordinates}
#' For more details on the criteria, please see the documentation of function 
#' \code{\link[nFactors]{nScree}}, which is the engine used provide the numbers of 
#' components to retain.
#' 
#' @seealso \code{\link[nFactors]{nScree}}, \code{\link[nFactors]{parallel}}, \code{\link{parallel_analysis}}
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data(HolzingerSwineford1939, package = 'lavaan')
#' HS_1939_std <- as.matrix(scale(HolzingerSwineford1939[, paste0('x', 1:9)]))
#' npcs <- factor_retention(HS_1939_std, iterations = 10, quantile = 0.05)
#' }
factor_retention <- function(data, iterations = 100, quantile = 0.05) {
  data <- scale(data, center = TRUE, scale = TRUE)

  uSVt <- svd(data)
  # Divide by N-1 to scale the eigenvalues the same way as if they came from
  # the eigen-decomposition of the correlation/covariance matrix
  eigenvalues <- uSVt$d^2 / (nrow(data) - 1)
  parallel_analysis_output <- parallel_analysis(
    N = nrow(data),
    J = ncol(data),
    iterations = iterations,
    quantile = quantile
  )
  non_graphical_solutions <- nFactors::nScree(
    x = eigenvalues,
    aparallel = parallel_analysis_output$summary_statistics_eigenvalues$quantile,
    model = "components"
  )

  number_factors <- list(
    kaiser = non_graphical_solutions$Components$nkaiser,
    parallel = non_graphical_solutions$Components$nparallel,
    acceleration_factor = non_graphical_solutions$Components$naf,
    optimal_coordinates = non_graphical_solutions$Components$noc
  )

  attr(number_factors, "class") <- "nfactors"

  return(number_factors)
}

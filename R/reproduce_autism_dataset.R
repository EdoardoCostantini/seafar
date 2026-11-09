#' Reproduce the autism dataset from GEO accession GSE7329
#'
#' Downloads the original Series Matrix file from NCBI GEO, reproduces the
#' existing preprocessing steps.
#'
#' Running `autism <- reproduce_autism_dataset()` is equivalent to `data("autism")`
#'
#' For more information, consult `data/README.md`.
#'
#' @param source_url Source URL of the original Series Matrix file on the NCBI GEO ftp server.
#' @importFrom utils read.csv
#'
#' @return Invisibly returns the processed matrix `autism`.
#' @examples
#' autism <- reproduce_autism_dataset()
#' @export
reproduce_autism_dataset <- function(
  source_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz",
  envir = parent.frame()
) {
  autism <- .read_autism_matrix_from_source(source_url)
  autism <- .extract_autism_labels_and_values(autism)
  autism <- .remove_autism_mislabeled_individuals(autism)
  autism <- .transpose_and_standardize_autism_matrix(autism)
  autism <- .remove_autism_features_with_missing_values(autism)

  # Return the dataset without printing it
  invisible(autism)
}

.read_autism_matrix_from_source <- function(source_url) {
  con <- gzcon(url(source_url, open = "rb"), text = TRUE)
  on.exit(close(con), add = TRUE)
  read.csv(con, header = FALSE, sep = "\t")
}

.extract_autism_labels_and_values <- function(full_data) {
  col_names <- c()
  for (i in 736:751) {
    col_names <- c(col_names, sapply(full_data[i, ], as.character))
  }

  character <- c()
  for (j in 287:302) {
    character <- c(character, sapply(full_data[j, ], as.character))
  }
  character <- character[-c(1, 32)]
  character <- character[-c(12, 13, 27)]
  character_new <- c(rep("dup15", 7), rep("FMR1", 6), rep("control", 14))

  matrix1 <- full_data[752:703647, ]
  matrix2 <- matrix(rep(0, prod(dim(matrix1))), ncol = dim(matrix1)[2])
  for (i in 1:nrow(matrix1)) {
    matrix2[i, ] <- as.numeric(sapply(matrix1[i, ], as.character))
  }

  autism <- matrix(data = t(matrix2), byrow = TRUE, ncol = length(col_names))
  colnames(autism) <- col_names
  IndexCol <- autism[, 1]
  autism <- autism[, -c(1, dim(autism)[2])]

  autism
}

.remove_autism_mislabeled_individuals <- function(autism) {
  autism[, -c(12, 13, 27)]
}

.transpose_and_standardize_autism_matrix <- function(autism) {
  autism <- t(autism)
  autism <- scale(autism, center = TRUE, scale = TRUE)

  autism
}

.remove_autism_features_with_missing_values <- function(autism) {
  NA_autism <- is.na(autism)
  NA_index <- colSums(NA_autism) != 0
  autism <- autism[, -which(NA_index != 0)]

  autism <- scale(autism, center = TRUE, scale = TRUE)

  autism
}

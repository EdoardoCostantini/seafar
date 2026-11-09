#' Autism Genetic Dataset
#'
#' The dataset contains 43893 genetic marker variables for 27 subjects, in which
#' 14 individuals are in the control group, six are in the FMR1-FM autism group,
#' and seven are in the dup15q autism group.
#'
#' @details
#' The original data was collected by Nishimura et. al. (2007) and is available from [NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7329). The version used in this package was obtained by removing three individuals who were mislabeled and standardizing each column of the data to have a mean of zero and a variance of one (see Guerra-Urzola et. al., 2023).
#' In particular, these processing steps were followed:
#'
#' 1.  download the original compressed Series Matrix file from NCBI GEO,
#' 2.  extract sample labels and expression values,
#' 3.  remove mislabeled individuals,
#' 4.  transpose and standardizes the expression matrix,
#' 5.  remove features with missing values after scaling, and
#' 6.  save the resulting processed matrix for further use.
#'
#' [reproduce_autism_dataset()] implements these steps and can be used to recreate the preprocessed dataset from its source.
#'
#' @section Terms of use:
#' The original dataset is hosted by NCBI on their public ftp infrastructure.
#' Their general terms of use state that data may be [used and distributed without restrictions](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
#' On their ftp server itself, they state that all data there is "[public, non-sensitive, unrestricted scientific data intended for sharing among scientific communities](https://ftp.ncbi.nlm.nih.gov/README.ftp)".
#' Please adhere to these terms when sharing `seafar`, the preprocessed datasets or modified versions of it.
#'
#' @format A matrix with 27 rows and 438893 variables
#' @source [NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7329)
#' @references
#' 
#' Nishimura Y, Martin CL, Vazquez-Lopez A, Spence SJ, Alvarez-Retuerto AI, Sigman M, Steindler C, Pellegrini S, Schanen NC, Warren ST, Geschwind DH. Genome-wide expression profiling of lymphoblastoid cell lines distinguishes different forms of autism and reveals shared pathways. Hum Mol Genet. 2007 Jul 15;16(14):1682-98. doi: 10.1093/hmg/ddm116. Epub 2007 May 21. PMID: 17519220.
#' 
#' Guerra-Urzola, R., Schipper, N. C. de, Tonne, A., Sijtsma, K., Vera, J. C., & Van Deun, K. (2023). Sparsifying the least-squares approach to PCA: Comparison of lasso and cardinality constraint. Advances in Data Analysis and Classification, 17(1), 269–286. https://doi.org/10.1007/s11634-022-00499-2
#' 
#' @examples
#' head(autism)
"autism"

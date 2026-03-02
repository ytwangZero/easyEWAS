#' Example Sample Metadata
#'
#' Example phenotype and covariate data used by `easyEWAS`.
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{SampleName}{Sample identifier.}
#'   \item{var}{Example exposure variable.}
#'   \item{cov1}{Example covariate 1.}
#'   \item{cov2}{Example covariate 2.}
#'   \item{batch}{Example batch variable.}
#' }
#' @source Simulated example data bundled with the package.
"sampledata"

#' Example Methylation Matrix
#'
#' Example CpG methylation matrix used by `easyEWAS`.
#'
#' @format A data frame with 1232 rows and 101 columns containing CpG probe
#' identifiers and sample-level methylation values.
#' @source Simulated example data bundled with the package.
"methydata"

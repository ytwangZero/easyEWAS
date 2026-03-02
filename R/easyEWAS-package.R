#' easyEWAS package
#'
#' @keywords internal
#' @importFrom R6 R6Class
#' @importFrom readxl read_xlsx
#' @importFrom stats IQR as.formula model.matrix p.adjust reformulate sd
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics barplot
#' @importFrom utils data globalVariables
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("chr", "genename", "methydata", "probe", "sampledata"))
}

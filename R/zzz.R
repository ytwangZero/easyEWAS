.onAttach <- function(libname, pkgname) {

  options(
    download.file.method.GEOquery = "auto",
    GEOquery.inmemory.gpl = FALSE
  )

  packageStartupMessage("Welcome! Thank you for using the R package easyEWAS! Detailed guidance can be found at https://easyewas.github.io/")
}

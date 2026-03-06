.onAttach <- function(libname, pkgname) {
  if (!interactive()) {
    return(invisible(NULL))
  }

  version <- as.character(utils::packageVersion(pkgname))
  packageStartupMessage(
    sprintf(
      "easyEWAS %s loaded. Tutorial: https://easyewas-tutorial.github.io/",
      version
    )
  )
}

.onAttach <- function(libname, pkgname) {

  version <- utils::packageVersion(pkgname)

  msg <- paste0(
    "Welcome to easyEWAS ", version, "!\n",
    strrep("=", 24 + nchar(as.character(version))), "\n",
    "Flexible and user-friendly tools for Epigenome-Wide Association Studies.\n\n",
    "Authors: Yuting Wang & Xu Gao (Corresponding author)\n",
    "Documentation: https://easyewas-tutorial.github.io/\n",
    "Citation:https://doi.org/10.1093/bioadv/vbaf026\n"
  )

  packageStartupMessage(msg)
}

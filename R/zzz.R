.onAttach <- function(libname, pkgname) {
  
  suppressMessages({
    options('download.file.method.GEOquery' = 'auto')
    options('GEOquery.inmemory.gpl' = FALSE)
  })

  packageStartupMessage("Welcome! Thank you for using the R package easyEWAS! Detailed guidance can be found at https://easyewas.github.io/.\n
                        Please cite us:\n
                        Wang Y, Jiang M, Niu S, Gao X. easyEWAS: a flexible and user-friendly R package for epigenome-wide association study[J]. Bioinformatics Advances, 2025, 5(1): vbaf026.")
}

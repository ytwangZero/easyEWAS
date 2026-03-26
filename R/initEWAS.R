#' @title  Initialize the EWAS module
#' @description This function generates an R6 class for storing EWAS analysis
#' data and results. By default, results are kept in memory and no files are
#' written. If file export is requested, users must supply an explicit output
#' directory.
#'
#' @param outpath Optional path to an existing directory. When `export = TRUE`,
#' a subdirectory named `"EWASresult"` is created under `outpath` for exported
#' result files.
#' @param export Logical. If `TRUE`, create an output folder and export result
#' files. If `FALSE` (default), do not write files to disk; keep results in the
#' returned object.
#' @return input, an R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom tictoc tic toc
#' @importFrom lubridate now
#' @examples
#' res <- initEWAS(export = FALSE)
#' res
initEWAS <- function(outpath = NULL, export = FALSE){
  tictoc::tic()
  if (!is.logical(export) || length(export) != 1 || is.na(export)) {
    stop("'export' must be TRUE or FALSE.")
  }
  if (!is.null(outpath) && (!is.character(outpath) || length(outpath) != 1 || !nzchar(outpath))) {
    stop("'outpath' must be NULL or a non-empty character string.")
  }

  if (isTRUE(export)) {
    if (is.null(outpath)) {
      stop("'outpath' must be provided when export = TRUE.")
    }
    if (!dir.exists(outpath)) {
      stop("Invalid path: '", outpath, "'. Please provide an existing directory.")
    }
    target_outpath <- file.path(outpath, "EWASresult")
  } else {
    if (!is.null(outpath)) {
      message("'outpath' is ignored because export = FALSE.")
    }
    target_outpath <- NULL
  }

  # create R6 class input ----------------------------------------------------
  input <- R6::R6Class(
    "input",
    public = list(
      outpath = NULL,
      export_output = TRUE,
      Data = list(Expo = NULL,
                  Methy = NULL),

      initialize = function(outpath, export_output) {
        self$outpath = outpath
        self$export_output = export_output
      }
    ),

    lock_class = FALSE,
    lock_objects = FALSE
  )

  input <- input$new(outpath = target_outpath, export_output = export)

  if (isTRUE(export)) {
    dir.create(input$outpath, recursive = TRUE, showWarnings = FALSE)
    lubridate::now() -> NowTime
    message("EWAS module has been successfully initialized.\n",
            "Results will be saved to: ", input$outpath, "\n",
            "Initialization time: ", NowTime, "\n")
  } else {
    lubridate::now() -> NowTime
    message("EWAS module has been successfully initialized.\n",
            "File export is disabled. Results will be stored in the returned object only.\n",
            "Initialization time: ", NowTime, "\n")
  }

  tictoc::toc()

  return(input)
}

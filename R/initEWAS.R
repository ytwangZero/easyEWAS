#' @title  Initialize the EWAS module
#' @description This function is designed to generate an R6 class for storing all the data
#' and results of EWAS analysis, and also to create a folder on the local computer for
#' storing the analysis results.
#' @usage initEWAS(outpath = "default")
#'
#' @param outpath The user-specified path is used to store a generated folder named "EWASresult"
#' that contains all the analysis results. If "default" is specified, the folder will be generated
#' in the current working directory.
#' @return input, an R6 class object integrating all information.
#' @export
#' @import dplyr
#' @import stringr
#' @import tictoc
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' }
initEWAS <- function(outpath = "default"){
  tictoc::tic()

  # create R6 class input ----------------------------------------------------
  input <- R6::R6Class(
    "input",
    public = list(
      outpath = NULL,
      Data = list(Expo = NULL,
                  Methy = NULL),

      initialize = function(outpath) {
        self$outpath = outpath}
    ),

    lock_class = FALSE,
    lock_objects = FALSE
  )

  input <- input$new(outpath = outpath)

  # set output file path---------------
  if(outpath == "default"){
    stringr::str_c(getwd(),"/EWASresult") -> input$outpath
    dir.create(input$outpath)
  }else if(outpath != "default" & file.exists(outpath)){
    stringr::str_c(outpath,"/EWASresult") -> input$outpath
    dir.create(input$outpath)
  }else if(outpath != "default" & !file.exists(outpath)){
    print("No such file or directory! Please enter the correct path.")
  }


  lubridate::now() -> NowTime
  message("Complete initializing the EWAS module! Your EWAS results will be stored in ",input$outpath, "\n",
          NowTime, "\n")


  tictoc::toc()

  return(input)
}

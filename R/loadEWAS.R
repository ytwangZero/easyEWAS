#' @title  Load all data files for EWAS module
#' @description Upload sample data and methylation data for EWAS analysis.
#' @usage loadEWAS(input, ExpoPath = NULL, MethyPath = NULL, ExpoData = "default",
#' MethyData = "default")
#'
#' @param input An R6 class integrated with all the information obtained from the initEWAS function.
#' @param ExpoPath The path to store the user's sample data. Each row represents a sample, and each
#' column represents a variable (exposure variable or covariate). Both .csv and .xlsx file types are
#' supported. The first column must be the sample ID, which must be consistent with the IDs in the
#' methylation data.
#' @param MethyPath The path to store the user's methylation data. Each row represents a CpG site, and
#' each column represents a sample. Both .csv and .xlsx file types are supported. The first column must
#' be the CpG probes. The sample IDs must be consistent with the IDs in the sample data.
#' @param ExpoData The data.frame of the user-supplied sample data that has been loaded into
#' the R environment. If default, the example data inside the package is used. The first column
#' must be the sample name.
#' @param MethyData The data.frame of the user-supplied methylation data that has been loaded into
#' the R environment. If default, an example of methylation data inside the package is loaded.
#' The first column must be the CpG site name.
#'
#' @return input, an R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom tictoc tic toc
#' @importFrom ddpcr quiet
#' @importFrom lubridate now
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' }
loadEWAS <- function(input,
                     ExpoPath = NULL,
                     MethyPath = NULL,
                     ExpoData = "default",
                     MethyData = "default"

){


  tictoc::tic()


  condition1 = is.null(ExpoPath) & is.null(MethyPath)
  condition2 = is.character(ExpoData) && is.character(MethyData)
  if(!condition1){
    #### external environment data----
    ExpoPath -> input$ExpoPath
    MethyPath -> input$MethyPath

    flag01 = file.exists(input$ExpoPath)
    flag02 = file.exists(input$MethyPath)
    if(all(flag01,flag02)){
      #read data ---------------------------------------------------------------
      ddpcr::quiet(
        if(substr(ExpoPath,nchar(ExpoPath)-3,nchar(ExpoPath)) == "xlsx"){
          readxl::read_xlsx(ExpoPath) %>% as.data.frame() -> input$Data$Expo
        }else{
          vroom::vroom(ExpoPath, delim = ",",show_col_types = F) %>%
            as.data.frame() -> input$Data$Expo
        })


      ddpcr::quiet(
        if(substr(ExpoPath,nchar(ExpoPath)-3,nchar(ExpoPath)) == "xlsx"){
          readxl::read_xlsx(MethyPath) %>% as.data.frame() -> input$Data$Methy
        }else{
          vroom::vroom(MethyPath, delim = ",",show_col_types = F) %>%
            as.data.frame() -> input$Data$Methy
        }
      )
      sample_names <- input$Data$Expo[[1]]
      probe_names <- input$Data$Methy[[1]]

      input$Data$Methy <- input$Data$Methy[, sample_names, drop = FALSE]
      rownames(input$Data$Methy) <- probe_names

      lubridate::now() -> NowTime
      message("All data files have been successfully loaded.\n",
              "Timestamp: ", NowTime, "\n")


      tictoc::toc()

      return(input)

    }else{
      lubridate::now() -> NowTime
      stop("No such file or directory! Please enter the correct path. \n",
           "Timestamp: ", NowTime, "\n")

      tictoc::toc()

      return(input)
    }

  }else if(condition1){
    if(condition2){
      #### example data------
      if(ExpoData == "default" & MethyData == "default"){

        data("sampledata", package = "easyEWAS", envir = environment())
        as.data.frame(sampledata) -> input$Data$Expo

        data("methydata", package = "easyEWAS", envir = environment())
        as.data.frame(methydata) -> input$Data$Methy

        sample_names <- input$Data$Expo[[1]]
        probe_names <- input$Data$Methy[[1]]

        input$Data$Methy <- input$Data$Methy[, sample_names, drop = FALSE]
        rownames(input$Data$Methy) <- probe_names

        lubridate::now() -> NowTime
        message("Example sample data and methylation data have been successfully loaded.\n",
                "Timestamp: ", NowTime, "\n")
        tictoc::toc()

        return(input)

      }else{
        lubridate::now() -> NowTime
        stop("Invalid value detected in 'ExpoData' or 'MethyData'. Please check your input. \n",
             "Timestamp: ", NowTime, "\n")
        tictoc::toc()

        return(input)
      }


    }else{

      #### R environment data------
      as.data.frame(ExpoData) -> input$Data$Expo
      as.data.frame(MethyData) -> input$Data$Methy

      sample_names <- input$Data$Expo[[1]]
      probe_names <- input$Data$Methy[[1]]

      input$Data$Methy <- input$Data$Methy[, sample_names, drop = FALSE]
      rownames(input$Data$Methy) <- probe_names

      lubridate::now() -> NowTime
      message("All data files have been successfully loaded.\n",
              "Timestamp: ", NowTime, "\n")

      tictoc::toc()

      return(input)

    }


  }



}




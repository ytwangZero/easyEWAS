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
#' @import magrittr
#' @importFrom tictoc tic toc
#' @importFrom ddpcr quiet
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
        if(substr(ExpoPath,nchar(ExpoPath)-3,nchar(ExpoPath))){
          readxl::read_xlsx(MethyPath) %>% as.data.frame() -> input$Data$Methy
        }else{
          vroom::vroom(MethyPath, delim = ",",show_col_types = F) %>%
            as.data.frame() -> input$Data$Methy
        }
      )

      lubridate::now() -> NowTime
      message("Complete loading All files! \n", NowTime, "\n")

      tictoc::toc()

      return(input)

    }else{
      lubridate::now() -> NowTime
      message("Error: No such file or directory! Please enter the correct path. \n", NowTime, "\n")

      tictoc::toc()

      return(input)
    }

  }else if(condition1){
    if(condition2){
      #### example data------
      if(ExpoData == "default" & MethyData == "default"){

        data("sampledata", package = "easyEWAS", envir = environment())
        sampledata -> input$Data$Expo
        # rm("sampledata")

        data("methydata", package = "easyEWAS", envir = environment())
        methydata -> input$Data$Methy
        # rm("methydata")

        lubridate::now() -> NowTime
        message("Sample data and methylation data for example have been loaded! \n", NowTime, "\n")
        tictoc::toc()

        return(input)

      }else{
        lubridate::now() -> NowTime
        message("Error: Unrecognized characters were entered at ExpoData or MethyData! \n", NowTime, "\n")
        tictoc::toc()

        return(input)
      }


    }else{

      #### R environment data------
      ExpoData -> input$Data$Expo
      MethyData -> input$Data$Methy

      lubridate::now() -> NowTime
      message("Complete loading All files! \n", NowTime, "\n")
      tictoc::toc()

      return(input)

    }


  }



}




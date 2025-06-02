#' @title Convert variable type of sample data
#' @description Transform the variable types of sample data to the types specified by users.
#' @usage transEWAS(input, Vars = "default", TypeTo = "factor")
#' @param input An R6 class integrated with all the information obtained from the loadEWAS function.
#' @param Vars Variable names that the user wants to convert types for, with each variable name
#' separated by a comma. Ensure there are no spaces. e.g. "var1,var2,var3".
#' @param TypeTo The type of variable that the function allows to be converted, including numeric and factor.
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom tictoc tic toc
#' @importFrom lubridate now
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' }
transEWAS <- function(input,
                      Vars = "default",
                      TypeTo = "factor" #numeric, factor
){
  tictoc::tic()

  lubridate::now() -> NowTime

  if (!TypeTo %in% c("numeric", "factor")) {
    stop("'TypeTo' must be either 'numeric' or 'factor'.")
  }


  VarName = unlist(strsplit(Vars,","))

  missing_vars <- setdiff(VarName, colnames(input$Data$Expo))
  if (length(missing_vars) > 0) {
    stop("Error: The following variables are not found in the exposure data: ",
         paste(missing_vars, collapse = ", "))
  }


  if(all(VarName == "default")){
    VarName = "cov1"
  }
  #check the loading data
  if(length(input$Data$Expo) == 0 |
     length(input$Data$Methy) == 0){

    tictoc::toc()

    stop("Exposure data or methylation data is missing. Please run loadEWAS() first.")

    return(input)

  }else{

    # Define data -----------------------------------------------------------------------------
    input$Data$Expo  -> df.all

    # func_data.type -----------------------------------------------------------------------------
    switch(TypeTo,
           "numeric" = {
             df.all %>%
               dplyr::mutate(dplyr::across(all_of(VarName), as.numeric)) -> input$Data$Expo
           },
           "factor" = {
             df.all %>%
               dplyr::mutate(dplyr::across(all_of(VarName), as.factor)) -> input$Data$Expo
           }
    )

    rm(df.all)
    message("Variable type conversion completed successfully.\n",
            "Converted variables: ", paste(VarName, collapse = ", "), "\n",
            "Target type: ", TypeTo, "\n",
            "Timestamp: ", NowTime, "\n")


    tictoc::toc()

    return(input)
  }
}

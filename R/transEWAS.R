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
#' @import stringr
#' @import tictoc
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


  VarName = unlist(strsplit(Vars,","))
  if(all(VarName == "default")){
    VarName = "cov1"
  }
  #check the loading data
  if(length(input$Data$Expo) == 0 |
     length(input$Data$Methy) == 0){

    tictoc::toc()

    print("Error: No exposure data or methylation data were found!")

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
    message("Variable types have been transformed successfully!")

    tictoc::toc()

    return(input)
  }
}

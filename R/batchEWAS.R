#' @title Perform batch effect correction
#' @description Perform batch effect correction based on the function \code{\link[sva]{ComBat}}
#' form R package \pkg{sva}. It requires that the "batches" in the data set are known. It uses
#' either parametric or non-parametric empirical Bayes frameworks for adjusting data for batch effects.
#' @usage batchEWAS(input,adjustVar = NULL,batch = NULL, plot = TRUE, par.prior = TRUE,
#' mean.only = FALSE,ref.batch = NULL)
#'
#' @param input An R6 class integrated with all the information.
#' @param batch Name of the batch variable.
#' @param plot (Optional) TRUE give prior plots with black as a kernel estimate of the empirical
#' batch effect density and red as the parametric.
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates
#' non-parametric adjustments will be used.
#' @param mean.only (Optional) Default to FALSE. If TRUE ComBat only corrects the mean of the
#' batch effect (no scale adjustment).
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch
#' adjustment.
#' @param adjustVar (Optional) Names of the variate of interest and other covariates besides batch,
#' with each name separated by a comma. Ensure that when correcting for batch effects, the effects
#' of other factors are appropriately considered and adjusted for.Ensure there are no space. e.g.
#' "cov1,cov2".
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom ddpcr quiet
#' @import stringr
#' @import tictoc
#' @import sva
#' @import BiocParallel
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- batchEWAS(input = res, batch = "batch", par.prior=TRUE, ref.batch = NULL)
#' }
batchEWAS = function(input,
                   adjustVar = NULL,
                   batch = NULL,
                   plot = TRUE,
                   par.prior = TRUE,
                   mean.only = FALSE,
                   ref.batch = NULL
                   ){

  tictoc::tic()

  dfcpg = as.matrix(input$Data$Methy[,-1])
  rownames(dfcpg) = input$Data$Methy[[1]]
  batch = input$Data$Expo[[batch]]

  if(is.null(adjustVar)){
    mod = NULL
  }else{
    var = strsplit(adjustVar, ",")
    ff = as.formula(paste0("~",paste(c(var[[1]]), collapse = " + ")))
    mod <- model.matrix(ff, data = input$Data$Expo)
  }

  message("Start to adjust batch effects, please be patient ...")
  if(plot){
    pdf(file = paste0(input$outpath,"/combat_plots.pdf"), width = 8, height = 8)
    ddpcr::quiet(combat_data <- ComBat(dat = dfcpg,
                                       batch = batch,
                                       mod = mod,
                                       par.prior=par.prior,
                                       prior.plots=TRUE,
                                       mean.only = mean.only,
                                       ref.batch = ref.batch,
                                       BPPARAM = bpparam("SerialParam")))
    dev.off()
  }else{
    ddpcr::quiet(combat_data <- ComBat(dat = dfcpg,
                                       batch = batch,
                                       mod = mod,
                                       par.prior=par.prior,
                                       prior.plots=FALSE,
                                       mean.only = mean.only,
                                       ref.batch = ref.batch,
                                       BPPARAM = bpparam("SerialParam")))
  }


  combat_data = combat_data %>%
    as.data.frame() %>%
    mutate(probe = input$Data$Methy[[1]]) %>%
    dplyr::select(probe, everything())
  rownames(combat_data) = NULL

  input$Data$Methy = combat_data

  lubridate::now()  -> NowTime
  message(paste0("Batch effect adjustment has been completed !\nYou can find the adjusted  methylation data in input$Data$Methy.\n", NowTime))

  tictoc::toc()

  return(input)


}





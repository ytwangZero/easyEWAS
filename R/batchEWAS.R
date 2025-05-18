#' @title Perform batch effect correction
#' @description Perform batch effect correction based on the function \code{\link[sva]{ComBat}}
#' form R package \pkg{sva}. It requires that the "batches" in the data set are known. It uses
#' either parametric or non-parametric empirical Bayes frameworks for adjusting data for batch effects.
#' @usage batchEWAS(input,adjustVar = NULL,batch = NULL, plot = TRUE, par.prior = TRUE,
#' mean.only = FALSE,ref.batch = NULL)
#'
#' @param input An R6 class integrated with all the information.
#' @param batch Name of the batch variable.
#' @param plot Logical. TRUE give prior plots with black as a kernel estimate of the empirical
#' batch effect density and red as the parametric.
#' @param par.prior Logical. TRUE indicates parametric adjustments will be used, FALSE indicates
#' non-parametric adjustments will be used.
#' @param mean.only Logical. Default to FALSE. If TRUE, ComBat only corrects the mean of the
#' batch effect (no scale adjustment).
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch
#' adjustment.
#' @param adjustVar (Optional) Names of the variate of interest and other covariates besides batch,
#' with each name separated by a comma. Ensure that when correcting for batch effects, the effects
#' of other factors are appropriately considered and adjusted for.Ensure there are no space. e.g.
#' "cov1,cov2".
#' @param parallel Logical. Whether to enable parallel computing during batch effect correction. 
#' Default is \code{FALSE}.
#' @param cores Integer. Number of CPU cores to use if \code{parallel = TRUE}. Default is \code{NULL}.
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom ddpcr quiet
#' @importFrom magrittr %>%
#' @importFrom tictoc tic toc
#' @importFrom sva ComBat
#' @import BiocParallel bpparam SnowParam MulticoreParam
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
                   ref.batch = NULL,
                   parallel = FALSE,
                   cores = NULL
                   ){

  tictoc::tic()

  dfcpg = as.matrix(input$Data$Methy[,-1])
  rownames(dfcpg) = input$Data$Methy[[1]]
  
  if (is.null(batch) || !(batch %in% colnames(input$Data$Expo))) {
    stop("Error: 'batch' must be a valid column name in the sample data.")
  }
  batch = input$Data$Expo[[batch]]

  if(is.null(adjustVar)){
    mod = NULL
  }else{
    var = unlist(strsplit(adjustVar, ","))
    missing_vars <- setdiff(var, colnames(input$Data$Expo))
    if (length(missing_vars) > 0) {
      stop("Error: These variables are not found in sample data: ",
           paste(missing_vars, collapse = ", "))
    }
    ff = as.formula(paste0("~",paste(var, collapse = " + ")))
    mod <- model.matrix(ff, data = input$Data$Expo)
  }

  message("Start to adjust batch effects, please be patient ...")
  if(plot){
    pdf(file = file.path(input$outpath, "combat_plots.pdf"), width = 8, height = 8)
    
    if(parallel){
      
      if (.Platform$OS.type == "windows") {
        BPPARAM <- BiocParallel::SnowParam(workers = cores)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
      }
      message("Running ComBat in parallel using ", cores, " core(s).")
      
      ddpcr::quiet(combat_data <- sva::ComBat(dat = dfcpg,
                                         batch = batch,
                                         mod = mod,
                                         par.prior=par.prior,
                                         prior.plots=TRUE,
                                         mean.only = mean.only,
                                         ref.batch = ref.batch,
                                         BPPARAM = BPPARAM))
      
    }else{
      
      ddpcr::quiet(combat_data <- sva::ComBat(dat = dfcpg,
                                              batch = batch,
                                              mod = mod,
                                              par.prior=par.prior,
                                              prior.plots=TRUE,
                                              mean.only = mean.only,
                                              ref.batch = ref.batch))
      
    }
    
    dev.off()
    message("A diagnostic plot was saved to: ", file.path(input$outpath, "combat_plots.pdf"))
  }else{
    
    if(parallel){
      
      if (.Platform$OS.type == "windows") {
        BPPARAM <- BiocParallel::SnowParam(workers = cores)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(workers = cores)
      }
      message("Running ComBat in parallel using ", cores, " core(s).")
      
      ddpcr::quiet(combat_data <- sva::ComBat(dat = dfcpg,
                                              batch = batch,
                                              mod = mod,
                                              par.prior=par.prior,
                                              prior.plots=FALSE,
                                              mean.only = mean.only,
                                              ref.batch = ref.batch,
                                              BPPARAM = BPPARAM))
      
    }else{
      
      ddpcr::quiet(combat_data <- sva::ComBat(dat = dfcpg,
                                              batch = batch,
                                              mod = mod,
                                              par.prior=par.prior,
                                              prior.plots=FALSE,
                                              mean.only = mean.only,
                                              ref.batch = ref.batch))
      
    }
  }


  combat_data = combat_data %>%
    as.data.frame() %>%
    mutate(probe = input$Data$Methy[[1]]) %>%
    dplyr::select(probe, everything())
  rownames(combat_data) = NULL

  input$Data$Methy = combat_data

  lubridate::now()  -> NowTime
  message("Batch effect adjustment completed. Adjusted data is saved in: input$Data$Methy.\n", NowTime)
  tictoc::toc()

  return(input)


}





#' @title  Perform EWAS Analysis
#' @description Perform EWAS analysis to obtain the coefficient value, standard
#' deviation and significance p value (or adjust p value) of each site.
#' @usage startEWAS(input, filename ="default", model = "lm", expo = "default",
#' cov = NULL,random = NULL, time = NULL, status = NULL, chipType = "EPICV2",
#' adjustP = TRUE, core = "default")
#'
#' @param input An R6 class integrated with all the information obtained from the loadEWAS or
#' transEWAS function.
#' @param filename filename Name of the output CSV file to store EWAS results. If set to "default", 
#' the file will be named "ewasresult.csv" and saved in the specified output directory.
#' @param model Statistical model to use for EWAS analysis. Options include:
#'   - "lm": Linear regression (default)
#'   - "lmer": Linear mixed-effects model
#'   - "cox": Cox proportional hazards model
#' @param expo Name of the exposure variable used in the EWAS analysis.
#' @param cov Comma-separated list of covariate variable names to include in the model (e.g., 
#' "age,sex,bmi"). Do not include spaces between names. Optional.
#' @param random Name of the grouping variable for the random intercept, required only when using 
#' the "lmer" model.
#' @param time Name of the time-to-event variable, required only when using the "cox" model.
#' @param status Name of the event/censoring indicator variable, required only when using the "cox" model.
#' @param adjustP Logical. If TRUE (default), adjusts p-values using both FDR (Benjamini-Hochberg) and 
#' Bonferroni correction methods.
#' @param chipType Illumina array platform used for DNA methylation measurement. Available options:
#'   - "27K"
#'   - "450K"
#'   - "EPICV1"
#'   - "EPICV2" (default)
#'   - "MSA"
#' @param core Number of CPU cores to use for parallel processing. If set to "default", uses the number 
#' of available physical cores minus one.
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom tictoc tic toc
#' @importFrom lubridate now
#' @import parallel detectCores makeCluster stopCluster
#' @import foreach foreach %dopar%
#' @import doParallel registerDoParallel stopImplicitCluster
#' @importFrom vroom vroom_write
#' @importFrom survival coxph
#' @importFrom lmerTest lmer
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, filename = "default", chipType = "EPICV2", model = "lm",
#' expo = "var", cov = "cov1,cov2",adjustP = TRUE, core = "default")
#' }
startEWAS = function(input,
                     filename ="default",
                     model = "lm",
                     expo = NULL,
                     cov = NULL,
                     random = NULL,
                     time = NULL,
                     status = NULL,
                     adjustP = TRUE,
                     chipType = "EPICV2",
                     core = "default"
){
  tictoc::tic()
  
  message("It will take some time, please be patient...")
  
  # select number of cores
  if(core == "default"){
    no_cores = detectCores(logical=F) - 1
  }else{
    no_cores = core
  }
  
  expo <- if (is.null(expo) || expo == "default") "var" else expo
  
  # peform EWAS analysis--------------------------------------------------------
  ## calculate factor number
  expo_vec <- unlist(input$Data$Expo[, expo])
  facnum <- if (is.factor(expo_vec)) length(levels(expo_vec)) else 2
  if (!is.factor(expo_vec) && model %in% c("lm", "lmer") && length(unique(expo_vec)) <= 4) {
    warning(paste0("Exposure variable '", expo, "' is not a factor but has only ",
                   length(unique(expo_vec)), " unique values. ",
                   "If this is a categorical variable, please use transEWAS() to convert it."))
  }
  
  
  ## define ewas function-----
  ewasfun <- function(cg, ff, cov) {
    cov$cpg <- as.vector(t(cg))
    res <- tryCatch({
      if (model == "lm") {
        out <- summary(lm(ff, data = cov))
        unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 4)]))
      } else if (model == "lmer") {
        out <- summary(lmer(ff, data = cov))
        unlist(lapply(2:facnum, function(i) out$coefficients[i, c(1, 2, 5)]))
      } else if (model == "cox") {
        out <- summary(coxph(ff, data = cov))
        c(as.vector(out$conf.int[1, c(1, 3, 4)]), out$coefficients[1, 5])
      } else {
        stop("Unsupported model: ", model)
      }
    }, error = function(e) {
      if (model %in% c("lm", "lmer")) return(rep(NA_real_, 3 * (facnum - 1)))
      if (model == "cox") return(rep(NA_real_, 4))
    })
    return(res)
  }
  model -> input$model
  
  ## set ewas parameter
  VarCov <- if (!is.null(cov)) unlist(strsplit(cov, ",")) else character(0)
  required_cols <- switch(model,
                          "lm"   = c(expo, VarCov),
                          "lmer" = c(expo, VarCov, random),
                          "cox"  = c(time, status, VarCov)
  )
  covdata <- input$Data$Expo[, required_cols, drop = FALSE]
  
  
  if (model == "lmer") {
    colnames(covdata)[match(random, colnames(covdata))] <- "random"
    input$random <- random
  } else if (model == "cox") {
    colnames(covdata)[match(time, colnames(covdata))] <- "time"
    colnames(covdata)[match(status, colnames(covdata))] <- "status"
    input$time <- time
    input$status <- status
  }
  input$covdata <- covdata
  
  
  formula <- switch (model,
                     "lm" = as.formula(paste0("cpg ~ ",paste(colnames(covdata), collapse = " + "))),
                     "lmer" = as.formula(paste0("cpg ~ ",paste(colnames(covdata)[-random_index], collapse = " + "), " + (1 | random)")),
                     "cox" = {
                       if(length(VarCov) > 0){
                         as.formula(paste0("Surv(time, status) ~ cpg + ", paste(VarCov, collapse = "+")))
                       }else{
                         as.formula("Surv(time, status) ~ cpg")
                       }
                     }
  )
  formula -> input$formula
  
  sample_names <- input$Data$Expo[[1]]
  df_beta <- input$Data$Methy[, sample_names, drop = FALSE]
  rownames(df_beta) <- input$Data$Methy[[1]]
  
  ## set parallel parameters----
  message("Start the EWAS analysis...")
  len = nrow(df_beta)
  chunk.size <- ceiling(len/no_cores)
  result_cols <- switch(model,
                        "lm" = 3 * (facnum - 1),
                        "lmer" = 3 * (facnum - 1),
                        "cox" = 4
  ) # identify number of columns that each model returns
  
  
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("ewasfun", "formula", "covdata", "df_beta", "facnum"), envir = environment())
  
  if (model == "lmer") clusterEvalQ(cl, library(lmerTest))
  if (model == "cox") clusterEvalQ(cl, library(survival))
  
  ## parallel computing-------
  if(model %in% c("lm","lmer")){
    
    system.time(
      modelres <- foreach(i=1:no_cores, .combine='rbind') %dopar%
        {  # local data for results
          restemp <- matrix(0, nrow=min(chunk.size, len-(i-1)*chunk.size), ncol=3*(facnum-1))
          for(x in ((i-1)*chunk.size+1):min(i*chunk.size, len)) {
            restemp[x - (i-1)*chunk.size,] <- as.numeric(base::t(ewasfun(df_beta[x,],formula,covdata)))
          }
          # return local results
          restemp
        }
    )
    
    
  }else{
    system.time(
      
      modelres <- foreach(i=1:no_cores, .combine='rbind') %dopar%
        {  # local data for results
          restemp <- matrix(0, nrow=min(chunk.size, len-(i-1)*chunk.size), ncol=4)
          for(x in ((i-1)*chunk.size+1):min(i*chunk.size, len)) {
            restemp[x - (i-1)*chunk.size,] <- as.numeric(base::t(ewasfun(df_beta[x,],formula,covdata)))
          }
          # return local results
          restemp
        }
    )
    
  }
  
  
  stopImplicitCluster()
  stopCluster(cl)
  modelres = modelres[1:len,]
  
  if((model %in% c("lmer","lm")) & class(unlist(input$Data$Expo[,expo])) == "factor"){
    
    ### categorical variable-----
    modelres %>%
      as.data.frame() %>%
      purrr::set_names(paste(rep(c("BETA","SE","PVAL"),facnum-1),rep(1:(facnum-1),each = 3), sep = "_")) %>%
      mutate(probe = rownames(df_beta)) %>%
      dplyr::select(probe, everything()) -> modelres
    
    #### FDR---
    if(adjustP){
      FDRname = paste(rep(c("FDR","Bonfferoni"),each = (facnum-1)),rep(1:(facnum-1),2), sep = "_")
      pindex = grep("PVAL",colnames(modelres))
      FDR = matrix(0,nrow = nrow(modelres),ncol = length(FDRname))
      for(i in pindex){
        FDR[,(i-1)/3] = p.adjust(modelres[[i]],method = "BH")
        FDR[,((i-1)/3)+(facnum-1)] = p.adjust(modelres[[i]],method = "bonferroni")
      }
      FDR = as.data.frame(FDR)
      FDR <- FDR %>%
        purrr::set_names(paste(rep(c("FDR","Bonfferoni"),each = (facnum-1)),rep(1:(facnum-1),2), sep = "_"))
      modelres = cbind(modelres,FDR)
      
    }
    
    
  }else if((model %in% c("lmer","lm")) & class(unlist(input$Data$Expo[,expo])) != "factor"){
    
    ### continuous variable-----
    modelres %>%
      as.data.frame() %>%
      purrr::set_names("BETA","SE","PVAL") %>%
      mutate(probe = rownames(df_beta))-> modelres
    
    #### per SD & IQR---
    modelres$BETA_perSD = (modelres$BETA)*(sd(covdata[[expo]],na.rm = TRUE))
    modelres$BETA_perIQR = (modelres$BETA)*(IQR(covdata[[expo]],na.rm = TRUE))
    modelres$SE_perSD = (modelres$SE)*(sd(covdata[[expo]],na.rm = TRUE))
    modelres$SE_perIQR = (modelres$SE)*(IQR(covdata[[expo]],na.rm = TRUE))
    modelres %>%
      dplyr::select(probe,BETA,BETA_perSD,BETA_perIQR,SE,SE_perSD,SE_perIQR,PVAL) -> modelres
    
    #### FDR---
    if(adjustP){
      message("Start multiple testing correction ...\n")
      modelres$FDR = p.adjust(modelres$PVAL, method = "BH")
      modelres$Bonfferoni = p.adjust(modelres$PVAL,method = "bonferroni")
      
    }
    
    
  }else if(model == "cox"){
    
    modelres %>%
      as.data.frame() %>%
      purrr::set_names("HR","LOWER_95%","UPPER_95%","PVAL") %>%
      mutate(probe = rownames(df_beta)) %>%
      dplyr::select(probe, everything())-> modelres
    
    #### FDR---
    if(adjustP){
      modelres$FDR = p.adjust(modelres$PVAL, method = "BH")
      modelres$Bonfferoni = p.adjust(modelres$PVAL,method = "bonferroni")
      
    }
  }
  
  message("Start CpG sites annotation ...\n")
  #### set annotation file----
  if(!is.null(chipType) & chipType == "EPICV2"){
    colnames(annotationV2) = c("probe","chr","pos","relation_to_island","gene","location")
    modelres %>%
      left_join(annotationV2, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "EPICV1"){
    colnames(annotationV1) = c("probe","chr","pos","relation_to_island","gene","location")
    modelres %>%
      left_join(annotationV1, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "450K"){
    colnames(annotation450K) = c("probe","chr","pos","relation_to_island","gene","location")
    modelres %>%
      left_join(annotation450K, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "27K"){
    colnames(annotation27K) = c("probe","chr","pos","gene")
    modelres %>%
      left_join(annotation27K, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "MSA"){
    colnames(annotationMSA) = c("probe","chr","pos","relation_to_island","gene")
    modelres %>%
      left_join(annotationMSA, by = "probe") -> modelres
  }
  
  
  # save model result-----
  modelres -> input$result
  if(filename == "default"){
    vroom::vroom_write(modelres, paste0(input$outpath, "/ewasresult.csv"), ",")
  }else{
    vroom::vroom_write(modelres, paste0(input$outpath, "/",filename, ".csv"), ",")
  }
  
  lubridate::now() -> NowTime
  message(paste0("EWAS analysis has been completed! \nYou can find results in ",input$outpath, ".\n", NowTime))
  
  tictoc::toc()
  
  return(input)
  
}

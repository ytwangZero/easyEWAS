#' @title  Perform EWAS Analysis
#' @description Perform EWAS analysis to obtain the coefficient value, standard
#' deviation and significance p value (or adjust p value) of each site.
#' @usage startEWAS(input, filename ="default", model = "lm", expo = "default",
#' cov = NULL,random = NULL, adjustP = TRUE, core = "default")
#'
#' @param input An R6 class integrated with all the information obtained from the loadEWAS or
#' transEWAS function.
#' @param filename User-customized .csv file name for storing EWAS results. If
#' "default" is chosen, it will be named as "ewasresult".
#' @param model The statistical models used for EWAS analysis include "lm" (general linear
#' regression), "lmer" (linear mixed-effects model), and "cox" (Cox proportional hazards model).
#' The default model is "lm".
#' @param expo Name of the exposure variable used in the EWAS analysis. The default
#' is a variable in the example data.
#' @param cov Name(s) of covariate(s) used in the EWAS analysis, with each  name separated by
#' a comma. Ensure there are no space. e.g. "cov1,cov2,cov3".
#' @param random Random intercept item name, used only when selecting the "lmer" model.
#' @param time When the user selects the Cox proportional risk model, the name
#' of the time variable needs to be specified.
#' @param status When the user selects the Cox proportional risk model, the name
#' of the status variable needs to be specified.
#' @param adjustP Whether to calculate adjusted p-values(FDR and Bonferroni correction). The default
#' is set to TRUE.
#' @param chipType The Illumina chip versions for user measurement of methylation data,
#' including "27K", 450K ","EPICV1", "EPICV2", and "MSA". The default is "EPICV2".
#' @param core The number of cores used during parallel computation. If set to default, it calculates
#' the maximum number of available physical cores minus 1 and treats this as an operational kernel.
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import vroom
#' @import stringr
#' @importFrom survival coxph
#' @import tictoc
#' @importFrom lmerTest lmer
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, filename = "default", chipType = "EPICV2", model = "lm",
#' expo = "default", adjustP = TRUE, core = "default")
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

  if(exists("expo")){
    if(!is.null(expo)){
      if(expo == "default"){
        expo = "var"
      }
    }
  }

  # peform EWAS analysis
  ## calculate factor number
  if(class(unlist(input$Data$Expo[,expo])) == "factor"){
    facnum = length(levels(unlist(input$Data$Expo[,expo])))
  }else{
    facnum = 2
  }

  ## define ewas function-----
  ewasfun <- switch (model,
                     "lm" = function(cg,ff,cov){
                       cov$cpg = as.vector(t(cg))
                       out <- base::summary(lm(ff, data = cov))
                       temp = c()
                       for (i in 2:facnum) {
                         temp = append(temp,as.vector(out$coefficients[i,c(1,2,4)]))
                       }
                       return(temp)
                     },
                     "lmer" = function(cg,ff,cov){
                       cov$cpg = as.vector(t(cg))
                       out <- base::summary(lmer(ff, data = cov))
                       temp = c()
                       for (i in 2:facnum) {
                         temp = append(temp,as.vector(out$coefficients[i,c(1,2,5)]))
                       }
                       return(temp)
                     },
                     "cox" = function(cg,ff,cov){
                       cov$cpg = as.vector(t(cg))
                       out <- base::summary(coxph(ff, data = cov))
                       return(c(as.vector(out$conf.int[1,c(1,3,4)]),
                                out$coefficients[1,5]))
                     }
  )
  model -> input$model

  ## set ewas parameter
  if(!is.null(cov)){
    VarCov = unlist(strsplit(cov,","))
    covdata = switch(model,
                     "lm" = input$Data$Expo[,c(expo,VarCov)],
                     "lmer" = input$Data$Expo[,c(expo,VarCov,random)],
                     "cox" = input$Data$Expo[,c(time,status,VarCov)],

    )
  }else{
    covdata = switch(model,
                     "lm" = input$Data$Expo[expo],
                     "lmer" = input$Data$Expo[c(expo,random)],
                     "cox" = input$Data$Expo[c(time,status)],

    )
  }


  if(model == "lmer"){
    random_index = which(colnames(covdata) == random)
    colnames(covdata)[random_index] = "random"
    input$random = random

  }else if(model == "cox"){
    time_index = which(colnames(covdata) == time)
    colnames(covdata)[time_index] = "time"
    status_index = which(colnames(covdata) == status)
    colnames(covdata)[status_index] = "status"

  }
  covdata -> input$covdata

  formula <- switch (model,
                     "lm" = as.formula(paste0("cpg ~ ",paste(colnames(covdata), collapse = " + "))),
                     "lmer" = as.formula(paste0("cpg ~ ",paste(colnames(covdata)[-random_index],
                                                                collapse = " + "), " + (1 | random)")),
                     "cox" = as.formula(paste0("Surv(time, status) ~ cpg + ",
                                                paste(VarCov, collapse = " + ")))
  )
  formula -> input$formula

  input$Data$Methy %>%
    as.data.frame() %>%
    dplyr::select(input$Data$Expo[[1]]) -> df_beta
  rownames(df_beta) = input$Data$Methy[[1]]
  ## set cores number----
  cl <- makeCluster(no_cores)
  registerDoParallel(cl, cores=no_cores)
  if(model == "lmer"){
    clusterEvalQ(cl, library(lmerTest))
  }
  if(model == "cox"){
    clusterEvalQ(cl, library(survival))
  }


  ## parallel computing-------
  message("Start the EWAS analysis...")
  len = nrow(df_beta)
  chunk.size <- ceiling(len/no_cores)
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
    colnames(annotationV2) = c("probe","chr","pos","relation_to_island","gene")
    modelres %>%
      left_join(annotationV2, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "EPICV1"){
    colnames(annotationV1) = c("probe","chr","pos","relation_to_island","gene")
    modelres %>%
      left_join(annotationV1, by = "probe") -> modelres
  }
  if(!is.null(chipType) & chipType == "450K"){
    colnames(annotation450K) = c("probe","chr","pos","relation_to_island","gene")
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

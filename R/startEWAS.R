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
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom vroom vroom_write
#' @importFrom survival coxph
#' @importFrom lmerTest lmer
#' @importFrom stats lm
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

  message("Starting EWAS data preprocessing ...")
  preprocess_start_time <- Sys.time()

  # -------------------------------------------
  # Set number of cores for parallel processing
  # -------------------------------------------
  if(core == "default"){
    no_cores = detectCores(logical=F) - 1
  }else{
    no_cores = core
  }

  # -----------------------------
  # Set exposure variable name
  # -----------------------------
  expo <- if (is.null(expo) || expo == "default") "var" else expo

  # ----------------------------------------------------------
  # Check if exposure is a factor and compute number of levels
  # ----------------------------------------------------------
  expo_vec <- unlist(input$Data$Expo[, expo])
  facnum <- if (is.factor(expo_vec)) length(levels(expo_vec)) else 2
  if (!is.factor(expo_vec) && model %in% c("lm", "lmer") && length(unique(expo_vec)) <= 4) {
    warning(paste0("Exposure variable '", expo, "' is not a factor but has only ",
                   length(unique(expo_vec)), " unique values. ",
                   "If this is a categorical variable, please use transEWAS() to convert it."))
  }


  # ----------------------------------------------------------
  # Define EWAS model fitting function based on selected model
  # ----------------------------------------------------------

  model -> input$model

  # -----------------------------
  # Prepare covariate data
  # -----------------------------
  VarCov <- if (!is.null(cov)) unlist(strsplit(cov, ",")) else character(0)
  required_cols <- switch(model,
                          "lm"   = c(expo, VarCov),
                          "lmer" = c(expo, VarCov, random),
                          "cox"  = c(time, status, VarCov)
  )
  covdata <- input$Data$Expo[, required_cols, drop = FALSE]

  # Rename variables based on model needs
  if (model == "lmer") {
    colnames(covdata)[match(random, colnames(covdata))] <- "random"
    random_index <- which(colnames(covdata) == "random")
    input$random <- random
  } else if (model == "cox") {
    colnames(covdata)[match(time, colnames(covdata))] <- "time"
    colnames(covdata)[match(status, colnames(covdata))] <- "status"
    input$time <- time
    input$status <- status
  }
  input$covdata <- covdata

  # -----------------------------
  # Build model formula
  # -----------------------------
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

  # --------------------------------
  # Extract methylation beta matrix
  # --------------------------------
  sample_names <- input$Data$Expo[[1]]
  df_beta <- input$Data$Methy[, sample_names, drop = FALSE]
  rownames(df_beta) <- input$Data$Methy[[1]]

  preprocess_end_time <- Sys.time()
  message("EWAS data preprocessing completed in ", round(preprocess_end_time - preprocess_start_time, 2), " seconds.\n")

  # -----------------------------
  # Set up parallel computation
  # -----------------------------
  message("Starting parallel computation setup ...")
  len = nrow(df_beta)
  chunk.size <- ceiling(len/no_cores)
  result_cols <- switch(model,
                        "lm" = 3 * (facnum - 1),
                        "lmer" = 3 * (facnum - 1),
                        "cox" = 4)

  setup_time <- system.time({

    cl <- makeCluster(no_cores)
    registerDoParallel(cl)

    if (model == "lm") {
      clusterExport(cl, varlist = "ewasfun_lm", envir = asNamespace("easyEWAS"))
    } else if (model == "lmer") {
      clusterEvalQ(cl, library(lmerTest))
      clusterExport(cl, varlist = "ewasfun_lmer", envir = asNamespace("easyEWAS"))
    } else if (model == "cox") {
      clusterEvalQ(cl, library(survival))
      clusterExport(cl, varlist = "ewasfun_cox", envir = asNamespace("easyEWAS"))
    }


  })["elapsed"]

  message("Parallel setup completed in ", round(setup_time, 2), " seconds.\n")

  # --------------------------------
  # Run parallel EWAS model fitting
  # --------------------------------
  message("Running parallel EWAS model fitting ...")
  ewas_start_time <- Sys.time()

  modelres <- foreach(i=1:no_cores, .combine='rbind', .errorhandling = "pass") %dopar%
    {
      restemp <- matrix(0, nrow=min(chunk.size, len-(i-1)*chunk.size), ncol=result_cols)
      for(x in ((i-1)*chunk.size+1):min(i*chunk.size, len)) {
        restemp[x - (i-1)*chunk.size,] <- as.numeric(

          t(if (model == "lm") {
            ewasfun_lm(df_beta[x, ], formula, covdata, facnum)
          } else if (model == "lmer") {
            ewasfun_lmer(df_beta[x, ], formula, covdata, facnum)
          } else if (model == "cox") {
            ewasfun_cox(df_beta[x, ], formula, covdata)
          })
        )
      }
      restemp
    }

  stopImplicitCluster()
  stopCluster(cl)
  modelres = as.data.frame(modelres[1:len,])

  ewas_end_time <- Sys.time()
  message(sprintf("Parallel EWAS model fitting completed in %.2f seconds.\n", as.numeric(difftime(ewas_end_time, ewas_start_time, units = "secs"))))

  # -----------------------------
  # Post-processing results
  # -----------------------------
  if((model %in% c("lmer","lm")) & class(unlist(input$Data$Expo[,expo])) == "factor"){

    ## categorical variable-----
    colnames(modelres)[1:(3 * (facnum - 1))] <- paste0(
      rep(c("BETA", "SE", "PVAL"), facnum - 1), "_", rep(1:(facnum - 1), each = 3)
    )
    modelres <- cbind(probe = rownames(df_beta), modelres)

    ## FDR pr Bonferroni adjustment---
    if(adjustP){
      FDRname = paste(rep(c("FDR","Bonfferoni"),each = (facnum-1)),rep(1:(facnum-1),2), sep = "_")
      pindex = grep("PVAL",colnames(modelres))
      FDR = matrix(0,nrow = nrow(modelres),ncol = length(FDRname))
      for(i in pindex){
        FDR[,(i-1)/3] = p.adjust(modelres[[i]],method = "BH")
        FDR[,((i-1)/3)+(facnum-1)] = p.adjust(modelres[[i]],method = "bonferroni")
      }
      FDR <- as.data.frame(FDR)
      colnames(FDR) <- paste(rep(c("FDR", "Bonfferoni"), each = (facnum - 1)),
                             rep(1:(facnum - 1), 2), sep = "_")

      modelres = cbind(modelres,FDR)
      message("Multiple testing correction completed!\n")

    }


  }else if((model %in% c("lmer","lm")) & class(unlist(input$Data$Expo[,expo])) != "factor"){

    ## continuous variable-----
    names(modelres)[1:3] <- c("BETA", "SE", "PVAL")
    modelres <- cbind(probe = rownames(df_beta), modelres)

    ## per SD & IQR---
    modelres$BETA_perSD = (modelres$BETA)*(sd(covdata[[expo]],na.rm = TRUE))
    modelres$BETA_perIQR = (modelres$BETA)*(IQR(covdata[[expo]],na.rm = TRUE))
    modelres$SE_perSD = (modelres$SE)*(sd(covdata[[expo]],na.rm = TRUE))
    modelres$SE_perIQR = (modelres$SE)*(IQR(covdata[[expo]],na.rm = TRUE))
    modelres <- modelres[, c("probe", "BETA", "BETA_perSD", "BETA_perIQR",
                             "SE", "SE_perSD", "SE_perIQR", "PVAL")]


    ## FDR pr Bonferroni adjustment---
    if(adjustP){
      modelres$FDR = p.adjust(modelres$PVAL, method = "BH")
      modelres$Bonfferoni = p.adjust(modelres$PVAL,method = "bonferroni")
      message("Multiple testing correction completed!\n")

    }


  }else if(model == "cox"){

    modelres <- as.data.frame(modelres)
    colnames(modelres) <- c("HR", "LOWER_95%", "UPPER_95%", "PVAL")
    modelres$probe <- rownames(df_beta)
    modelres <- modelres[, c("probe", "HR", "LOWER_95%", "UPPER_95%", "PVAL")]

    ##  FDR pr Bonferroni adjustment---
    if(adjustP){
      modelres$FDR = p.adjust(modelres$PVAL, method = "BH")
      modelres$Bonfferoni = p.adjust(modelres$PVAL,method = "bonferroni")
      message("Multiple testing correction completed!\n")
    }
  }


  # ---------------------------------
  # Annotate CpG sites with chip info
  # ---------------------------------
  message("Start CpG sites annotation ...")

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

  if (!is.null(chipType)) {
    chipGenome <- switch(chipType,
                         "EPICV2" = "hg38 (GRCh38)",
                         "EPICV1" = "hg19 (GRCh37)",
                         "450K"   = "hg19 (GRCh37)",
                         "27K"    = "hg19 (GRCh37)",
                         "MSA"    = "hg19 (GRCh37)",
                         "Unknown genome"
    )
    message(sprintf("Using annotation for chip type: %s (Genome: %s)\n", chipType, chipGenome))
  }


  # -----------------------------
  # Save results to CSV
  # -----------------------------
  modelres -> input$result
  output_path <- if (filename == "default") {
    file.path(input$outpath, "ewasresult.csv")
  } else {
    file.path(input$outpath, paste0(filename, ".csv"))
  }
  vroom::vroom_write(modelres, output_path, delim = ",")

  lubridate::now() -> NowTime
  message(paste0("EWAS analysis has been completed! You can find results in ",input$outpath, ".\n", NowTime))

  tictoc::toc()

  return(input)

}

#' @title  Perform Bootstrap-based Internal Validation
#' @description Users can perform internal validation of the identified differently methylated
#' sites based on the bootstrap method.
#'
#' @usage bootEWAS(input, filterP = "PVAL", cutoff = 0.05, CpGs = NULL, times = 500,
#' bootCI = "perc",filename = "default")
#'
#' @param input An R6 class integrated with all the information obtained from the startEWAS or
#' plotEWAS function.
#' @param filterP The name of the p value columns such as "PVAL", "FDR", and "Bonfferoni." Users use
#' this P-value to screen for significance sites and further conduct internal validation.
#' @param cutoff The cutoff value of the P-value used to filter for further internal validation.
#' The default is 0.05.
#' @param CpGs The name of the methylation site specified by the user for bootstrap analysis,
#' separated by commas. Be careful not to have spaces, such as "cpg1,cpg2".
#' @param times Number of bootstrap times specified by the user. The default value is 100 times.
#' @param filename User-customized .csv file name for storing bootstrap results. If "default", it will
#' be named as "bootresult".
#' @param bootCI A vector of character strings representing the type of interval to base the test on.
#' The value should be one of "norm", "basic", "stud", "perc" (the default), and "bca".
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom tictoc tic toc
#' @importFrom ddpcr quiet
#' @importFrom boot boot boot.ci
#' @importFrom survival coxph
#' @importFrom lmerTest lmer
#'
#'
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, chipType = "EPICV2", model = "lm", expo = "default", adjustP = TRUE)
#' res <- bootEWAS(input = res, filterP = "PVAL", cutoff = 0.05, times = 100)
#' }
bootEWAS = function(input,
                    filterP = "PVAL", # FDR, Bonfferoni
                    cutoff = 0.05,
                    CpGs = NULL,
                    times = 500,
                    bootCI = "perc",
                    filename = "default"){

  tictoc::tic()
  if (is.null(input$result)) {
    stop("No EWAS result found in 'input$result'. Please run startEWAS() first.")
  }
  if (!filterP %in% colnames(input$result)) {
    stop("The specified filterP column '", filterP, "' was not found in the EWAS results.")
  }
  if (!bootCI %in% c("norm", "perc", "basic", "stud", "bca")) {
    stop("Invalid 'bootCI' value. Must be one of 'norm', 'perc', 'basic', 'stud', 'bca'.")
  }



  message("Starting bootstrap-based internal validation...")
  if (is.null(CpGs)) {
    sig_probes <- input$result$probe[input$result[[filterP]] < cutoff]
    if (length(sig_probes) == 0) {
      stop("No CpG sites meet the filtering criteria (", filterP, " < ", cutoff, ").")
    }
    message(length(sig_probes), " CpG sites selected by ", filterP, " < ", cutoff)
  } else {
    sig_probes <- unlist(strsplit(CpGs, ","))
    missing_probes <- setdiff(sig_probes, input$result$probe)
    if (length(missing_probes) > 0) {
      stop("Error: The following CpG names are not found in EWAS results:\n",
           paste(missing_probes, collapse = ", "), "\nTime: ", lubridate::now())
    }
    message(length(sig_probes), " CpG sites specified by user.")
  }


  df_beta <- input$Data$Methy[rownames(input$Data$Methy) %in% sig_probes, input$Data$Expo[[1]]]
  df <- as.data.frame(cbind(input$covdata, t(df_beta)))

  #### peform bootstrap analysis---------------------
  input$result %>%
    filter(probe %in% sig_probes) %>%
    dplyr::select(2) -> original


  set.seed(123)
  if (length(sig_probes) == 0) {
    stop("No CpG sites meet the filtering criteria.\nTime: ", lubridate::now())
  }else{

    message(length(sig_probes), " CpG sites selected for bootstrap analysis.")
    message("Bootstrap confidence interval method: ", bootCI)
    if(input$model == "lm"){
      coef_function <- function(data, formula, indices) {
        d <- data[indices,] #allows boot to select sample
        fit <- lm(formula, data = d)
        return(summary(fit)$coefficients[2,1])
      }
    }else if(input$model == "lmer"){
      coef_function <- function(data, formula, indices) {
        d <- data[indices,] #allows boot to select sample
        fit <- lmer(formula, data = d)
        return(summary(fit)$coefficients[2,1])
      }
    }else if(input$model == "cox"){
      coef_function <- function(data, formula, indices) {
        d <- data[indices,] #allows boot to select sample
        fit <- coxph(formula, data = d)
        return(summary(fit)$conf.int[1,1])
      }
    }

    if(input$model %in% c("lmer")){
      random_index = which(colnames(df) == input$random)
      colnames(df)[random_index] = "random"
    }
    if(input$model %in% c("lm","lmer","cox")){
      lower_CI <- c()
      upper_CI <- c()
      pval <- c()
      for(cpg in sig_probes){
        message("Bootstrap analysis for ", cpg, "...")
        if(input$model %in% c("lm","lmer")){
          formula=as.formula(paste0(cpg, "~", as.character(input$formula)[3]))
        }else if(input$model %in% c("cox")){
          formula = as.formula(paste0("Surv(time, status) ~ ", cpg,
                                      substr(as.character(input$formula)[3], 4,
                                             nchar(as.character(input$formula)[3]))))
        }

        ddpcr::quiet(boot.res <- boot(data=df, statistic=coef_function, R=times, formula=formula))
        ddpcr::quiet(boot.ci(boot.res, type=bootCI) -> temp_CI)
        # ddpcr::quiet(boot.pval(boot.res, type = "perc") -> temp_P)


        ci_vals <- switch(bootCI,
                          "norm"  = temp_CI$normal[2:3],
                          "perc"  = temp_CI$percent[4:5],
                          "basic" = temp_CI$basic[4:5],
                          "stud"  = temp_CI$student[4:5],
                          "bca"   = temp_CI$bca[4:5]
        )

        lower_temp <- ci_vals[1]
        upper_temp <- ci_vals[2]

        lower_CI = c(lower_CI, lower_temp)
        upper_CI = c(upper_CI, upper_temp)
      }



      input$bootres <- tibble(
        probe = sig_probes,
        original = original[[1]],
        lower_CI = lower_CI,
        upper_CI = upper_CI
      )


      if(filename == "default"){
        vroom::vroom_write(input$bootres, paste0(input$outpath, "/bootresult.csv"), ",")
      }else{
        vroom::vroom_write(input$bootres, file.path(input$outpath, paste0(filename, ".csv")), ",")
      }

      lubridate::now()  -> NowTime
      message(paste0("Bootstrap for internal validation has been completed!\nYou can find results in ",input$outpath, ".\n", NowTime))

      tictoc::toc()

      return(input)

    }

  }


}



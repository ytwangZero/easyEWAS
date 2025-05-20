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
#' @import tictoc
#' @importFrom ddpcr quiet
#' @import boot
#' @import boot.pval
#' @import survival
#' @importFrom lmerTest lmer
#' @importFrom lubridate now
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
  # message("It will take some time, please be patient...")
  #### method 1: p value filter--------

  if(is.null(CpGs)){
    base::subset(input$result, input$result[filterP] < cutoff, select = probe) -> cpgname

    # colnames(input$Data$Methy)[1]
    input$Data$Methy %>%
      filter(input$Data$Methy[[1]] %in% cpgname$probe) %>%
      as.data.frame() %>%
      dplyr::select(colnames(input$Data$Methy)[1], input$Data$Expo[[1]]) -> df_beta
    rownames(df_beta) = df_beta[[1]]
    df_beta[1] = NULL
    tdf_beta = t(df_beta)

    df = cbind(input$covdata, tdf_beta)
    input$cpgnames <- cpgname$probe

  }

  #### method 2: cpg name filter-------
  if(!is.null(CpGs)){

    cpgname = unlist(strsplit(CpGs,","))
    if(!all(cpgname %in% input$result$probe)){
      lubridate::now() -> NowTime
      message("Error: Not all CpG names are in the EWAS results! Please enter the correct names. \n",
              NowTime, "\n")

      tictoc::toc()

      return(input)
      stop("")
    }else{
      colnames(input$Data$Methy)[1] = "probe"
      input$Data$Methy %>%
        filter(probe %in% cpgname) %>%
        as.data.frame() %>%
        dplyr::select(colnames(input$Data$Methy)[1], input$Data$Expo[[1]]) -> df_beta
      rownames(df_beta) = df_beta[[1]]
      df_beta[1] = NULL
      tdf_beta = t(df_beta)

      df = cbind(input$covdata, tdf_beta)
      input$cpgnames = cpgname
    }

  }

  #### peform bootstrap analysis---------------------
  input$result %>%
    filter(probe %in% input$cpgnames) %>%
    dplyr::select(2) -> original


  set.seed(123)
  if(length(input$cpgnames) == 0){
    lubridate::now()  -> NowTime
    message("Error: No CpG sites meeting the filtering criteria were found! \n",
            NowTime, "\n")

    tictoc::toc()

    return(input)
    stop("")
  }else{
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
      for(cpg in input$cpgnames){
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


        if(bootCI == "norm"){
          lower_temp = temp_CI$normal[2]
          upper_temp = temp_CI$normal[3]

        }else if(bootCI == "perc"){
          lower_temp = temp_CI$percent[4]
          upper_temp = temp_CI$percent[5]

        }else if(bootCI == "basic"){
          lower_temp = temp_CI$basic[4]
          upper_temp = temp_CI$basic[5]

        }else if(bootCI == "stud"){
          lower_temp = temp_CI$student[4]
          upper_temp = temp_CI$student[5]

        }else if(bootCI == "bca"){
          lower_temp = temp_CI$bca[4]
          upper_temp = temp_CI$bca[5]

        }

        lower_CI = c(lower_CI, lower_temp)
        upper_CI = c(upper_CI, upper_temp)
        # pval = c(pval, temp_P)
      }



      input$bootres <- tibble(
        probe = input$cpgnames,
        original = original[[1]],
        lower_CI = lower_CI,
        upper_CI = upper_CI,

      )


      if(filename == "default"){
        vroom::vroom_write(input$bootres, paste0(input$outpath, "/bootresult.csv"), ",")
      }else{
        vroom::vroom_write(input$bootres, paste0(input$outpath, "/",filename, ".csv"), ",")
      }

      lubridate::now()  -> NowTime
      message(paste0("Bootstrap for internal validation has been completed !\nYou can find results in ",input$outpath, ".\n", NowTime))

      tictoc::toc()

      return(input)

    }

  }


}



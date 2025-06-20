#' @title Perform Differentially Methylated Region analysis
#' @description Perform differential methylation analysis based on the R package \pkg{DMRcate}.
#' Computes a kernel estimate against a null comparison to identify significantly DMRs.
#' @usage dmrEWAS(input, chipType = "EPICV2", what = "Beta", expo = NULL, cov = NULL, genome = "hg38",
#' lambda=1000, C = 2, filename = "default",fdrCPG = 0.05, pcutoff = "fdr", min.cpgs = 2,
#' epicv2Filter = "mean")
#'
#' @param input An R6 class integrated with all the information.
#' @param chipType The Illumina chip versions for user measurement of methylation data,
#' including "450K","EPICV1", and "EPICV2". The default is "EPICV2".
#' @param what Types of methylation values, including "Beta" and "M". Default to "Beta".
#' @param genome Reference genome for annotating DMRs.
#' Must be consistent with the array platform used:
#'
#' - Use `"hg38"` for EPICV2 arrays.
#' - Use `"hg19"` for 450K and EPICV1 arrays.
#'
#' Note: the `genome` argument does not currently affect the internal behavior of `extractRanges()`
#' (i.e., no liftover is performed).
#' @param lambda  If the distance between two significant CpG sites is greater than or equal
#' to lambda, they will be considered as belonging to different DMRs. The default value is 1000
#' nucleotides, meaning that if the distance between two significant CpG sites exceeds 1000
#' nucleotides, they will be separated into different DMRs.
#' @param filename User-customized .csv file name for storing DMR results. If "default", it will
#' be named as "DMRresult".
#' @param expo Name of the exposure variable used in the DMR analysis.
#' @param cov Name(s) of covariate(s) used in the DMR analysis, with each  name separated by
#' a comma. Ensure there are no space. e.g. "cov1,cov2,cov3".
#' @param C Scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma.
#' Empirical testing shows for both Illumina and bisulfite sequencing data that, when lambda=1000,
#' near-optimal prediction of sequencing-derived DMRs is obtained when C is approximately 2,
#' i.e. 1 standard deviation of Gaussian kernel = 500 base pairs. Cannot be < 0.2.
#' @param fdrCPG Used to individually assess the significance of each CpG site. If the FDR-adjusted
#' p-value of a CpG site is below the specified fdrCPG threshold, the site will be marked as significant.
#' The default value is 0.05.
#' @param pcutoff Used to determine the threshold for DMRs. It is strongly recommended to use the
#' default (fdr), unless you are confident about the risk of Type I errors (false positives).
#' @param min.cpgs Minimum number of consecutive CpGs constituting a DMR. Default to 2.
#' @param epicv2Filter Strategy for filtering probe replicates that map to the same CpG site.
#' "mean" takes the mean of the available probes; "sensitivity" takes the available probe most
#' sensitive to methylation change; "precision" either selects the available probe with the
#' lowest variation from the consensus value (most precise), or takes the mean if that confers
#' the lowest variation instead, "random" takes a single probe at random from each replicate group.
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom ddpcr quiet
#' @importFrom vroom vroom_write
#' @importFrom tictoc tic toc
#' @importFrom DMRcate cpg.annotate dmrcate extractRanges
#' @importFrom lubridate now
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- dmrEWAS(input = res, filename = "default", chipType = "EPICV2", what = "Beta", expo = "var",
#' cov = "cov1,cov2", genome = "hg38",)
#' }
#'

dmrEWAS = function(input,
                   chipType = "EPICV2",
                   what = "Beta",
                   epicv2Filter = "mean",
                   expo = NULL,
                   cov = NULL,
                   genome = "hg38",
                   fdrCPG = 0.05,
                   pcutoff = "fdr",
                   lambda=1000,
                   C = 2,
                   min.cpgs = 2,
                   filename = "default"
                   ){

  if (!chipType %in% c("EPICV2", "EPICV1", "450K")) {
    stop("Invalid 'chipType'. Must be one of: 'EPICV2', 'EPICV1', or '450K'.")
  }
  if (!genome %in% c("hg19", "hg38")) {
    stop("Invalid 'genome'. Must be 'hg19' or 'hg38'.")
  }
  if (!epicv2Filter %in% c("mean", "sensitivity", "precision", "random")) {
    stop("Invalid 'epicv2Filter'. Choose from: mean, sensitivity, precision, random.")
  }
  if (is.null(expo) || !expo %in% colnames(input$Data$Expo)) {
    stop("Exposure variable is missing or not found in input$Data$Expo.")
  }


  cov_list <- if (!is.null(cov)) strsplit(cov, ",")[[1]] else character(0)
  if (!is.null(cov)) {
    missing_covs <- setdiff(cov_list, colnames(input$Data$Expo))
    if (length(missing_covs) > 0) {
      stop("Covariates not found in input$Data$Expo: ", paste(missing_covs, collapse = ", "))
    }
  }



  tictoc::tic()
  array_map <- c("EPICV2" = "EPICv2", "EPICV1" = "EPICv1", "450K" = "450K")
  arraytype <- array_map[[chipType]]

  # Filter out position replicates from an EPICv2 beta- or M-matrix ---
  input$Data$Methy %>%
    as.data.frame() %>%
    dplyr::select(input$Data$Expo[[1]]) %>%
    as.matrix() -> dfcpg
  # rownames(dfcpg) = input$Data$Methy[[1]]

  formula <- reformulate(termlabels = c(expo, cov_list))
  design <- model.matrix(formula, data = input$Data$Expo)

  # DMR analysis ---
  message("Starting the differentially methylated region analysis. Please be patient...")
  if (arraytype == "EPICv2") {
    if (is.null(epicv2Filter) || epicv2Filter == "") {
      stop("Error: 'epicv2Filter' must be specified when chipType is 'EPICV2'.")
    }
    message("EPICv2 detected. Filtering position replicates using method: '", epicv2Filter, "'...")
  }

  ddpcr::quiet(myannotation <- cpg.annotate("array", dfcpg, arraytype = arraytype, what = what,
                               analysis.type="differential",
                               design=design, coef=2, fdr = fdrCPG,
                               epicv2Remap = TRUE, epicv2Filter = epicv2Filter))

  ddpcr::quiet(dmrcoutput <- dmrcate(myannotation, lambda=lambda, C=C, pcutoff=pcutoff))
  ddpcr::quiet(results.ranges <- extractRanges(dmrcoutput, genome = genome))
  input$dmrres = as.data.frame(results.ranges)

  message(length(myannotation@ranges), " CpGs used for DMR analysis.")
  message(length(results.ranges), " DMRs detected.")


  # save result ---
  outfile <- if (filename == "default") "DMRresult.csv" else paste0(filename, ".csv")
  vroom::vroom_write(input$dmrres, file.path(input$outpath, outfile), delim = ",")


  lubridate::now()  -> NowTime
  message("DMR result has been saved to: ", file.path(input$outpath, outfile))
  message("DMR analysis completed successfully.")
  message("Timestamp: ", NowTime)
  tictoc::toc()

  return(input)

}

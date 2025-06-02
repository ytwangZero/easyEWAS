#' @title  Enrichment analyses
#' @description Perform GO or KEGG enrichment analysis based on the \pkg{clusterProfiler} package.
#' @usage enrichEWAS(input, filename = "default", method = "GO", ont = "MF", pool = FALSE,
#' filterP = "PVAL", cutoff = 0.05, plot = TRUE, plotType = "dot", plotcolor = "pvalue",
#' showCategory= NULL, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
#'
#' @param input An R6 class integrated with all the information obtained from the startEWAS or plotEWAS
#' or bootEWAS function.
#' @param method Methods of enrichment analysis, including "GO" and "KEGG".
#' @param filterP The name of the p value columns such as "PVAL", "FDR", and "Bonfferoni." Users use
#' this P-value to screen for significance sites and further conduct enrichment analysis.
#' @param cutoff The cutoff value of the P-value used to filter for further enrichment analysis.
#' The default is 0.05.
#' @param filename User-customized .xlsx file name for storing EWAS results. If
#' "default" is chosen, it will be named as "enrichresult".
#' @param plot Whether the results of enrichment analysis need to be visualized, the default is TRUE
#' @param plotType Whether to draw a bar plot ("bar") or a dot plot ("dot"), the default is "dot".
#' @param plotcolor It is the vertical axis of the picture of the enrichment analysis results. Users can choose
#' "pvalue" or "p.adjust" or "qvalue". The default is "p.adjust".
#' @param showCategory The number of categories which will be displayed in the plots. Default to 10.
#' @param pvalueCutoff The p-value threshold used to filter enrichment results. Only results that pass
#' the p-value test (i.e., those smaller than this value) will be reported. This value refers to the
#' p-value before adjustment. The p-value represents the probability of observing the current level of enrichment
#' under the assumption of no enrichment. The smaller the p-value, the more significant the enrichment result.
#' @param pAdjustMethod The p-value adjustment method used for multiple hypothesis testing, aimed at reducing false
#' positives caused by multiple comparisons. One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. The q-value is the result
#' of controlling the false discovery rate (FDR) and represents the proportion of false positives that
#' may occur when conducting multiple tests.Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff
#' on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported. The default is 0.2.
#' @param ont When choosing GO enrichment analysis, select the GO sub-ontology for which the enrichment analysis
#' will be performed. One of "BP", "MF", and "CC" sub-ontologies, or "ALL" for all three. Default to "BP".
#' @param pool If ont='ALL', whether pool three GO sub-ontologies.
#' @param x Character string specifying the variable to be used on the x-axis of the plot.
#' Common options are "GeneRatio" or "Count".
#' - "GeneRatio": ratio of input genes annotated to a given term.
#' - "Count": the number of input genes annotated to the term.
#' @param width Width of the PDF output in inches. Default is 11.
#' @param height Height of the PDF output in inches. Default is 7.
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @importFrom tictoc tic toc
#' @importFrom clusterProfiler bitr enrichGO enrichKEGG
#' @import org.Hs.eg.db
#' @importFrom vroom vroom_write
#' @importFrom R.utils setOption
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, chipType = "EPICV2", model = "lm", expo = "default", adjustP = TRUE)
#' res <- plotEWAS(input = res, pval = "PVAL")
#' res <- bootEWAS(input = res, filterP = "PVAL", cutoff = 0.05, times = 100)
#' res <- enrichEWAS(input = res, method = "GO", filterP = "PVAL", cutoff = 0.05, pAdjustMethod = "BH")
#' }

enrichEWAS <- function(input,
                       filename = "default",
                       method = "GO",
                       filterP = "PVAL",
                       cutoff = 0.05,
                       ont = "BP",
                       pool = FALSE,
                       plot = TRUE,
                       plotType = "dot", # bar
                       plotcolor = "p.adjust",
                       x = "GeneRatio",
                       showCategory=10,
                       width = 11,
                       height = 7,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2
){
  tictoc::tic()
  lubridate::now() -> NowTime

  R.utils::setOption("clusterProfiler.download.method","auto")

  if (is.null(input$result)) {
    stop("No EWAS result found in 'input$result'.\nTime: ", NowTime)
  }
  if (!method %in% c("GO", "KEGG")) {
    stop("Invalid method. Must be one of: 'GO', 'KEGG'")
  }


  subset(input$result, input$result[filterP] < cutoff, select = gene) %>%
    as.data.frame() %>%
    mutate(genename = sub(";.*$", "", gene)) %>%
    dplyr::select(genename) %>%
    filter(genename != "") -> enrichdata
  if (nrow(enrichdata) == 0) {
    stop("No genes passed the filtering threshold (", filterP, " < ", cutoff, ").")
  }


  message("Converting gene symbols to Entrez IDs...")
  ddpcr::quiet(
    gene.df <- bitr(enrichdata$genename,
                    fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
  )
  message("Gene symbol conversion completed, and ", nrow(gene.df), " genes mapped.")

  gene<-gene.df$ENTREZID
  message("Start enrichment analysis ...")
  if(method == "GO"){

    enres <- enrichGO(gene = gene,
                      OrgDb = org.Hs.eg.db,
                      ont=ont,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = pAdjustMethod,
                      qvalueCutoff = qvalueCutoff,
                      readable = TRUE)
    enres@result -> input$enrichres
    input$enrichres$GeneRatio = as.character(input$enrichres$GeneRatio)
    input$enrichres$BgRatio = as.character(input$enrichres$BgRatio)

  }else if(method == "KEGG"){
    enres <- enrichKEGG(gene = gene,
                        organism = "hsa",
                        pvalueCutoff = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        qvalueCutoff = qvalueCutoff)
    enres@result -> input$enrichres
  }

  outfile <- if (filename == "default") "enrichresult.csv" else paste0(filename, ".csv")
  vroom::vroom_write(input$enrichres, file.path(input$outpath, outfile), delim = ",")


  if (plot) {
    message("Start result visualization ...")

    # Check if there are results to plot
    if (nrow(enres@result) == 0) {
      message("Warning: No enrichment results available for plotting.")
    } else {

      # Determine output file name suffix based on plot type
      suffix <- switch(plotType,
                       "dot" = "enrichdot.pdf",
                       "bar" = "enrichbar.pdf",
                       {
                         warning("Invalid plotType. Choose 'dot' or 'bar'. Skipping plot.")
                         return(input)
                       })

      # Compose the full file name
      file_name <- file.path(input$outpath,
                             if (filename == "default") suffix else paste0(filename, ".pdf"))

      # Generate the appropriate plot
      pdf(file = file_name, width = width, height = height)
      p <- switch(plotType,
                  "dot" = dotplot(enres,
                                  x = x,
                                  color = plotcolor,
                                  showCategory = showCategory),
                  "bar" = barplot(enres,
                                  x = x,
                                  color = plotcolor,
                                  showCategory = showCategory)
      )
      print(p)
      dev.off()

      message("Plot saved to: ", file_name)
    }
  }


  lubridate::now()  -> NowTime
  message(paste0("Enrichment analysis has been completed !\nYou can find results in ",input$outpath, ".\n", NowTime))

  tictoc::toc()

  return(input)

}

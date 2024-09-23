#' @title  Enrichment analyses
#' @description Perform GO or KEGG enrichment analysis based on the \pkg{clusterProfiler} package.
#' @usage enrichEWAS(input, filename = "default", method = "GO", ont = "MF", pool = FALSE, filterP = "PVAL",
#' cutoff = 0.05, plot = TRUE, plotType = "dot", plotcolor = "pvalue", showCategory= NULL, pAdjustMethod = "BH",
#' qvalueCutoff = 0.2)
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
#' "pvalue" or "p.adjust" or "qvalue". The default is "pvalue".
#' @param showCategory The number of categories which will be displayed in the plots.
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
#' will be performed. One of "BP", "MF", and "CC" sub-ontologies, or "ALL" for all three. Default to "MF".
#' @param pool If ont='ALL', whether pool 3 GO sub-ontologies.
#'
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @import stringr
#' @import tictoc
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import writexl
#' @importFrom ggplot2 ggsave
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
                       ont = "MF",
                       pool = FALSE,
                       plot = TRUE,
                       plotType = "dot", # bar
                       plotcolor = "pvalue",
                       showCategory=NULL,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2
){
  tictoc::tic()
  lubridate::now() -> NowTime

  R.utils::setOption("clusterProfiler.download.method","auto")

  if(is.null(input$result)){
    message("Error: No EWAS result file found.\n", NowTime)
  }else{
    subset(input$result, input$result[filterP] < cutoff, select = gene) %>%
      as.data.frame() %>%
      mutate(genename = sub(";.*$", "", gene)) %>%
      dplyr::select(genename) %>%
      filter(genename != "") -> enrichdata

    message("It will take some time, please be patient...")
    ddpcr::quiet(
      gene.df <- bitr(enrichdata$genename,
                      fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)
    )
    gene<-gene.df$ENTREZID
    message("Start enrichment analysis ...")
    if(method == "GO"){

      enres <- enrichGO(gene = gene,
                        OrgDb = org.Hs.eg.db,
                        ont="BP",
                        pvalueCutoff = pvalueCutoff,
                        qvalueCutoff = qvalueCutoff,
                        readable = TRUE)
      enres@result -> input$enrichres
      input$enrichres$GeneRatio = as.character(input$enrichres$GeneRatio)
      input$enrichres$BgRatio = as.character(input$enrichres$BgRatio)

    }else if(method == "KEGG"){
      enres <- enrichKEGG(gene = gene,
                          organism = "hsa",
                          pvalueCutoff = pvalueCutoff,
                          qvalueCutoff = qvalueCutoff)
      enres@result -> input$enrichres
    }

    if(filename == "default"){
      writexl::write_xlsx(input$enrichres, paste0(input$outpath, "/enrichresult.xlsx"))
    }else{
      writexl::write_xlsx(input$enrichres,paste0(input$outpath, "/",filename, ".xlsx"))
    }

    if (plot) {
      message("Start result visualization ...")

      if (plotType == "dot") {
        file_name <- if (filename == "default") {
          paste0(input$outpath, "/enrichdot.pdf")
        } else {
          paste0(input$outpath, "/", filename, ".pdf")
        }

        p <- dotplot(enres, x = "GeneRatio",
                     color = plotcolor,
                     decreasing = TRUE,
                     showCategory = showCategory)

        ggsave(filename = file_name, plot = p)
      }

      if (plotType == "bar") {
        file_name <- if (filename == "default") {
          paste0(input$outpath, "/enrichbar.pdf")
        } else {
          paste0(input$outpath, "/", filename, ".pdf")
        }

        p <- barplot(enres, x = "Count",
                     color = plotcolor,
                     showCategory = showCategory)

        ggsave(filename = file_name, plot = p)
      }
    }

    lubridate::now()  -> NowTime
    message(paste0("Enrichment analysis has been completed !\nYou can find results in ",input$outpath, ".\n", NowTime))

    tictoc::toc()

    return(input)

  }


}

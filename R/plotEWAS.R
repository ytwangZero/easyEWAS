#' @title  Visualize the results of EWAS analysis
#' @description Visualize EWAS results based on the \code{\link[CMplot]{CMplot}} package, including Manhattan
#' plots, QQ plots, etc. Please note that this function only supports plotting a single-layer circular
#' Manhattan plot. Additionally, the meaning of each parameter in this function is exactly the same as
#' in \code{\link[CMplot]{CMplot}} For more detailed information or to create multi-layer circular Manhattan
#' plots, please refer to \code{\link[CMplot]{CMplot}} (https://cran.r-project.org/web/packages/CMplot/index.html).
#' @usage plotEWAS(input, p = "PVAL", threshold=NULL, file=c("jpg","pdf","tiff","png"),
#' col=c("#4197d8","#f8c120","#413496","#495226","#d60b6f","#e66519","#d581b7","#83d3ad",
#' "#7c162c","#26755d"),LOG10=TRUE,pch=19,type="p",band=1,axis.cex=1,axis.lwd=1.5,
#' lab.cex=1.5,lab.font=2,plot.type=c("m","c","q","d"),r=0.3,cex=c(0.5,1,1),ylab="",
#' ylab.pos=3,xticks.pos=1,threshold.col="red", threshold.lwd=1,threshold.lty=2,
#' amplify=FALSE,signal.cex=1.5,signal.pch=19,signal.col=NULL,signal.line=2, highlight=NULL,
#' highlight.cex=1,highlight.pch=19,highlight.type="p",highlight.col="red",highlight.text=NULL,
#' highlight.text.col="black",highlight.text.cex=1,highlight.text.font=3,chr.labels=NULL,
#' chr.border=FALSE,chr.labels.angle=0,cir.axis=TRUE,cir.axis.col="black",cir.axis.grid=TRUE,
#' conf.int=TRUE,conf.int.col=NULL, file.name="",dpi=300,height=NULL,width=NULL,main="",
#' main.cex=1.5,main.font=2,box=FALSE,verbose=FALSE)
#'
#' @param input An R6 class integrated with all the information obtained from the startEWAS function.
#' @param p The user needs to specify the name of the p value selected for the result
#' visualization.
#' @param threshold The significant threshold.If threshold = 0 or NULL, then the threshold line will
#' not be added.
#' @param file The format of the output image file, including "jpg","pdf","tiff", and "png".
#' @param col A vector specifies the colors for the chromosomes. If the length of col is shorter than
#' the number of chromosomes, the colors will be applied cyclically.
#' @param LOG10 logical, whether to change the p-value into log10(p-value) scale.
#' @param pch a integer, the shape for the points, is the same with "pch" in \code{\link[plot]{plot}}.
#' @param type a character, could be "p" (point), "l" (cross line), "h" (vertical lines) and so on,
#' is the same with "type" in \code{\link[plot]{plot}}.
#' @param band a number, the size of space between chromosomes, the default is 1.
#' @param axis.cex a number, controls the size of ticks labels of X/Y-axis and the ticks labels of axis
#' for circle plot.
#' @param axis.lwd a number, controls the thickness of X/Y-axis lines and the thickness of axis for circle plot.
#' @param lab.cex a number, controls the size of labels of X/Y-axis and the labels of chromosomes for circle plot.
#' @param lab.font a number, controls the font of labels of all axis.
#' @param plot.type a character or vector, only "d", "c", "m", "q" can be used. if plot.type="d",
#' CpG density will be plotted; if plot.type="c", only circle-Manhattan plot will be plotted; if
#' plot.type="m",only Manhattan plot will be plotted; if plot.type="q",only Q-Q plot will be plotted;
#' if plot.type=c("m","q"), Both Manhattan and Q-Q plots will be plotted.
#' @param r a number, the radius for the circle (the inside radius), the default is 1.
#' @param cex a number or a vector, the size for the points, is the same with "size" in \code{\link[plot]{plot}}, and if
#' it is a vector, the first number controls the size of points in circle plot(the default is 0.5), the
#' second number controls the size of points in Manhattan plot (the default is 1), the third number controls
#' the size of points in Q-Q plot (the default is 1)
#' @param ylab a character, the labels for y axis.
#' @param ylab.pos the distance between ylab and yaxis.
#' @param xticks.pos 	the distance between labels of x ticks and x axis.
#' @param threshold.col a character or vector, the color for the line of threshold levels, it can also
#' control the color of the diagonal line of QQplot.
#' @param threshold.lwd a number or vector, the width for the line of threshold levels, it can also
#' control the thickness of the diagonal line of QQplot.
#' @param threshold.lty a number or vector, the type for the line of threshold levels, it can also
#' control the type of the diagonal line of QQplot
#' @param amplify logical, CMplot can amplify the significant points, if TRUE, then the points bigger
#' than the minimal significant level will be amplified, the default: amplify=TRUE.
#' @param signal.cex a number, if amplify=TRUE, users can set the size of significant points.
#' @param signal.pch a number, if amplify=TRUE, users can set the shape of significant points.
#' @param signal.col a character, if amplify=TRUE, users can set the colour of significant points,
#' if signal.col=NULL, then the colors of significant points will not be changed.
#' @param signal.line a number, the thickness of the lines of significant CpGs cross the circle.
#' @param highlight a vector, names of CpGs which need to be highlighted.
#' @param highlight.cex a vector, the size of points for CpGs which need to be highlighted.
#' @param highlight.pch a vector, the pch of points for CpGs which need to be highlighted.
#' @param highlight.type a vector, the type of points for CpGs which need to be highlighted.
#' @param highlight.col a vector, the col of points for CpGs which need to be highlighted.
#' @param highlight.text a vector, the text which would be added around the highlighted CpGs.
#' @param highlight.text.col a vector, the color for added text.
#' @param highlight.text.cex a value, the size for added text.
#' @param highlight.text.font text font for the highlighted CpGs
#' @param chr.labels a vector, the labels for the chromosomes of density plot and Manhattan plot.
#' @param chr.border a logical, whether to plot the dot line between chromosomes.
#' @param chr.labels.angle a value, rotate tick labels of x-axis for Manhattan plot (-90 < chr.labels.angle < 90).
#' @param cir.axis a logical, whether to add the axis of circle Manhattan plot.
#' @param cir.axis.col a character, the color of the axis for circle.
#' @param cir.axis.grid logical, whether to add axis grid line in circles.
#' @param conf.int logical, whether to plot confidence interval on QQ-plot.
#' @param conf.int.col character or vector, the color of confidence interval of QQplot.
#' @param file.name a character or vector, the names of output files.
#' @param dpi a number, the picture resolution for '.jpg', '.npg', and '.tiff' files. The default is 300.
#' @param height the height of output files.
#' @param width the width of output files.
#' @param main character of vector, the title of the plot for manhattan plot and qqplot.
#' @param main.cex size of title.
#' @param main.font font of title.
#' @param box logical, this function draws a box around the current plot.
#' @param verbose whether to print the log information.
#' @param bin.size a integer, the size of bin in bp for marker density plot.
#' @param bin.breaks a vector, set the breaks for the legend of density plot, e.g., seq(min, max, step),
#' the windows in which the number of markers is out of the this range will be plotted in the same colors
#' with the min or max value.
#' @param ylim vector (c(min, max)), CMplot will only plot the points among this interval.
#' @param outward logical, if TRUE, all points will be plotted from inside to outside for circular Manhattan plot.
#' @param mar the size of white gaps around the plot, 4 values should be provided, indicating the direction
#' of bottom, left, up, and right.
#' @param chr.den.col a character or vector or NULL, the colour for the CpG density. If the length of parameter
#' 'chr.den.col' is bigger than 1, CpG density that counts the number of CpG within given size ('bin.size')
#' will be plotted around the circle. If chr.den.col=NULL, the density bar will not be attached on the bottom
#' of manhattan plot.
#' @param chr.pos.max logical, whether the physical positions of each chromosome contain the maximum length
#' of the chromosome.
#' @param cir.chr logical, a boundary that represents chromosomes will be plotted on the periphery of a circle,
#' the default is TRUE.
#' @param cir.chr.h a number, the width for the boundary, if cir.chr=FALSE, then this parameter will be useless.
#' @param file.output a logical, users can choose whether to output the plot results.
#'
#' @return The updated input object, including CMplot-ready data stored in input$CMplot.
#' @export
#' @import dplyr
#' @importFrom CMplot CMplot
#' @importFrom tictoc tic toc
#' @importFrom lubridate now
#' @importFrom withr with_dir
#'
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, chipType = "EPICV2", model = "lm", expo = "var", adjustP = TRUE)
#' res <- plotEWAS(input = res, p = "PVAL")
#' }
plotEWAS <- function(input,
                     p = "PVAL",
                     threshold=NULL,
                     file=c("jpg","pdf","tiff","png"),
                     col=c("#4197d8","#f8c120","#413496","#495226","#d60b6f",
                           "#e66519","#d581b7","#83d3ad","#7c162c","#26755d"),
                     bin.size=1e6,
                     bin.breaks=NULL,
                     LOG10=TRUE,
                     pch=19,
                     type="p",
                     band=1,
                     H=1.5,
                     ylim=NULL,
                     axis.cex=1,
                     axis.lwd=1.5,
                     lab.cex=1.5,
                     lab.font=2,
                     plot.type=c("m","c","q","d"),
                     multracks=FALSE,
                     multracks.xaxis=FALSE,
                     multraits=FALSE,
                     points.alpha=100L,
                     r=0.3,
                     cex=c(0.5,1,1),
                     outward=FALSE,
                     ylab=expression(-log[10](italic(p))),
                     ylab.pos=3,
                     xticks.pos=1,
                     mar=c(3,6,3,3),
                     mar.between=0,
                     threshold.col="red",
                     threshold.lwd=1,
                     threshold.lty=2,
                     amplify=FALSE,
                     signal.cex=1.5,
                     signal.pch=19,
                     signal.col=NULL,
                     signal.line=2,
                     highlight=NULL,
                     highlight.cex=1,
                     highlight.pch=19,
                     highlight.type="p",
                     highlight.col="red",
                     highlight.text=NULL,
                     highlight.text.col="black",
                     highlight.text.cex=1,
                     highlight.text.font=3,
                     chr.labels=NULL,
                     chr.border=FALSE,
                     chr.labels.angle=0,
                     chr.den.col="black",
                     chr.pos.max=FALSE,
                     cir.band=1,
                     cir.chr=TRUE,
                     cir.chr.h=1.5,
                     cir.axis=TRUE,
                     cir.axis.col="black",
                     cir.axis.grid=TRUE,
                     conf.int=TRUE,
                     conf.int.col=NULL,
                     file.output=TRUE,
                     file.name="",
                     dpi=300,
                     height=NULL,
                     width=NULL,
                     main="",
                     main.cex=1.5,
                     main.font=2,
                     legend.ncol=NULL,
                     legend.cex=1,
                     legend.pos=c("left","middle","right"),
                     box=FALSE,
                     verbose=FALSE){

  tictoc::tic()
  if (is.null(input$result)) {
    stop("No EWAS result found in 'input'. Please run startEWAS() first.")
  }
  if (!p %in% colnames(input$result)) {
    stop(sprintf("The p-value column '%s' does not exist in input$result. Please check the input or specify a correct column name.", p))
  }
  if (file.name == "") {
    message("No file name provided. Output images will use default names from CMplot.")
  }
  if (is.null(threshold)) {
    threshold <- 0.05
    message("No significance threshold provided. Default threshold (0.05) will be used.")
  }




  input$result %>%
    dplyr::select("probe","chr","pos", all_of(p)) %>%
    arrange(chr) -> df
  sx = factor(df$chr,levels = c("1","2","3","4","5","6","7","8","9","10",
                                "11","12","13","14","15","16","17","18",
                                "19","20","21","22","X","Y"), ordered = T)
  df[order(sx),] -> input$CMplot

  #plot the ewas result----------
  withr::with_dir(input$outpath, {
    CMplot::CMplot(input$CMplot,
                   col=col,
                   bin.size=bin.size, bin.breaks=bin.breaks, LOG10=LOG10, pch=pch, type=type, band=band,
                   H=H, ylim=ylim, axis.cex=axis.cex, axis.lwd=axis.lwd, lab.cex=lab.cex, lab.font=lab.font,
                   plot.type=plot.type, multracks=multracks,
                   multracks.xaxis=multracks.xaxis, multraits=multraits, points.alpha=points.alpha, r=r,
                   cex=cex, outward=outward, ylab=ylab,
                   ylab.pos=ylab.pos, xticks.pos=xticks.pos, mar=mar, mar.between=mar.between, threshold=threshold,
                   threshold.col=threshold.col, threshold.lwd=threshold.lwd, threshold.lty=threshold.lty, amplify=amplify,
                   signal.cex=signal.cex, signal.pch=signal.pch, signal.col=signal.col, signal.line=signal.line,
                   highlight=highlight, highlight.cex=highlight.cex, highlight.pch=highlight.pch, highlight.type=highlight.type,
                   highlight.col=highlight.col, highlight.text=highlight.text, highlight.text.col=highlight.text.col,
                   highlight.text.cex=highlight.text.cex, highlight.text.font=highlight.text.font, chr.labels=chr.labels,
                   chr.border=chr.border, chr.labels.angle=chr.labels.angle, chr.den.col=chr.den.col,
                   chr.pos.max=chr.pos.max, cir.band=cir.band, cir.chr=cir.chr, cir.chr.h=cir.chr.h,
                   cir.axis=cir.axis, cir.axis.col=cir.axis.col, cir.axis.grid=cir.axis.grid, conf.int=conf.int,
                   conf.int.col=conf.int.col, file.output=file.output, file.name=file.name, file=file, dpi=dpi,
                   height=height, width=width, main=main, main.cex=main.cex,
                   main.font=main.font, legend.ncol=legend.ncol, legend.cex=legend.cex,
                   legend.pos=legend.pos, box=box, verbose=verbose)
  })


  lubridate::now() -> NowTime
  message("EWAS Visualization has completed! \nYou can find results in ",input$outpath, ".\n", NowTime, "\n")

  tictoc::toc()

  return(input)
}

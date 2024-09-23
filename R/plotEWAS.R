#' @title  Visualize the results of EWAS analysis
#' @description Visualize EWAS results based on the CMplot package, including Manhattan
#' plots, QQ plots, etc. For other detailed drawing parameters, please refer to
#' \code{\link[CMplot]{CMplot}}.
#' @usage plotEWAS(input, p = "PVAL", threshold=NULL, threshold=NULL, file="jpg")
#' @param input An R6 class integrated with all the information obtained from the startEWAS function.
#' @param p The user needs to specify the name of the p value selected for the result
#' visualization.
#' @param threshold The significant threshold.If threshold = 0 or NULL, then the threshold line will
#' not be added.
#' @param file The format of the output image file, including "jpg","pdf","tiff", and "png".
#' @return input, An R6 class object integrating all information.
#' @export
#' @import dplyr
#' @import CMplot
#' @import tictoc
#'
#' @examples \dontrun{
#' res <- initEWAS(outpath = "default")
#' res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
#' res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
#' res <- startEWAS(input = res, chipType = "EPICV2", model = "lm", expo = "default", adjustP = TRUE)
#' res <- plotEWAS(input = res, p = "PVAL")
#' }
plotEWAS <- function(input,
                     p = "PVAL",
                     threshold=NULL,
                     file=c("jpg","pdf","tiff","png"),
                     col=c("#4197d8","#f8c120","#413496","#495226",
                           "#d60b6f","#e66519","#d581b7","#83d3ad","#7c162c","#26755d"),
                     bin.size=1e6,bin.breaks=NULL,LOG10=TRUE,pch=19,type="p",band=1,
                     H=1.5,ylim=NULL,axis.cex=1,axis.lwd=1.5,lab.cex=1.5,lab.font=2,
                     plot.type=c("m","c","q","d"),multracks=FALSE,
                     multracks.xaxis=FALSE,multraits=FALSE,points.alpha=100L,r=0.3,
                     cex=c(0.5,1,1),outward=FALSE,ylab=expression(-log[10](italic(p))),
                     ylab.pos=3,xticks.pos=1,mar=c(3,6,3,3),mar.between=0,
                     threshold.col="red",threshold.lwd=1,threshold.lty=2,amplify=FALSE,
                     signal.cex=1.5,signal.pch=19,signal.col=NULL,signal.line=2,
                     highlight=NULL,highlight.cex=1,highlight.pch=19,highlight.type="p",
                     highlight.col="red",highlight.text=NULL,highlight.text.col="black",
                     highlight.text.cex=1,highlight.text.font=3,chr.labels=NULL,
                     chr.border=FALSE,chr.labels.angle=0,chr.den.col="black",
                     chr.pos.max=FALSE,cir.band=1,cir.chr=TRUE,cir.chr.h=1.5,
                     cir.axis=TRUE,cir.axis.col="black",cir.axis.grid=TRUE,conf.int=TRUE,
                     conf.int.col=NULL,file.output=TRUE,file.name="",
                     dpi=300,height=NULL,width=NULL,main="",main.cex=1.5,
                     main.font=2,legend.ncol=NULL,legend.cex=1,
                     legend.pos=c("left","middle","right"),box=FALSE,verbose=FALSE){

  tictoc::tic()
  input$result %>%
    dplyr::select("probe","chr","pos", all_of(p)) %>%
    arrange(chr) -> df
  sx = factor(df$chr,levels = c("1","2","3","4","5","6","7","8","9","10",
                                "11","12","13","14","15","16","17","18",
                                "19","20","21","22","X","Y"), ordered = T)
  df[order(sx),] -> input$CMplot

  #plot the ewas result----------
  path = getwd()
  setwd(input$outpath)
  CMplot::CMplot(input$CMplot,col=col,
                 bin.size=bin.size,bin.breaks=bin.breaks,LOG10=LOG10,pch=pch,type=type,band=band,
                 H=H,ylim=ylim,axis.cex=axis.cex,axis.lwd=axis.lwd,lab.cex=lab.cex,lab.font=lab.font,
                 plot.type=plot.type,multracks=multracks,
                 multracks.xaxis=multracks.xaxis,multraits=multraits,points.alpha=points.alpha,r=r,
                 cex=cex,outward=outward,ylab=ylab,
                 ylab.pos=ylab.pos,xticks.pos=xticks.pos,mar=mar,mar.between=mar.between,threshold=threshold,
                 threshold.col=threshold.col,threshold.lwd=threshold.lwd,threshold.lty=threshold.lty,amplify=amplify,
                 signal.cex=signal.cex,signal.pch=signal.pch,signal.col=signal.col,signal.line=signal.line,
                 highlight=highlight,highlight.cex=highlight.cex,highlight.pch=highlight.pch,highlight.type=highlight.type,
                 highlight.col=highlight.col,highlight.text=highlight.text,highlight.text.col=highlight.text.col,
                 highlight.text.cex=highlight.text.cex,highlight.text.font=highlight.text.font,chr.labels=chr.labels,
                 chr.border=chr.border,chr.labels.angle=chr.labels.angle,chr.den.col=chr.den.col,
                 chr.pos.max=chr.pos.max,cir.band=cir.band,cir.chr=cir.chr,cir.chr.h=cir.chr.h,
                 cir.axis=cir.axis,cir.axis.col=cir.axis.col,cir.axis.grid=cir.axis.grid,conf.int=conf.int,
                 conf.int.col=conf.int.col,file.output=file.output,file.name=file.name,file=file,dpi=dpi,
                 height=height,width=width,main=main,main.cex=main.cex,
                 main.font=main.font,legend.ncol=legend.ncol,legend.cex=legend.cex,
                 legend.pos=legend.pos,box=box,verbose=verbose)
  setwd(path)


  lubridate::now() -> NowTime
  message("Ewas Visualization has completed! \nYou can find results in ",input$outpath, ".\n", NowTime, "\n")




  tictoc::toc()

  return(input)
}

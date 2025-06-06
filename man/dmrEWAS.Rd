% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmrEWAS.R
\name{dmrEWAS}
\alias{dmrEWAS}
\title{Perform Differentially Methylated Region analysis}
\usage{
dmrEWAS(input, chipType = "EPICV2", what = "Beta", expo = NULL, cov = NULL, genome = "hg38",
lambda=1000, C = 2, filename = "default",fdrCPG = 0.05, pcutoff = "fdr", min.cpgs = 2,
epicv2Filter = "mean")
}
\arguments{
\item{input}{An R6 class integrated with all the information.}

\item{chipType}{The Illumina chip versions for user measurement of methylation data,
including "450K","EPICV1", and "EPICV2". The default is "EPICV2".}

\item{what}{Types of methylation values, including "Beta" and "M". Default to "Beta".}

\item{epicv2Filter}{Strategy for filtering probe replicates that map to the same CpG site.
"mean" takes the mean of the available probes; "sensitivity" takes the available probe most
sensitive to methylation change; "precision" either selects the available probe with the
lowest variation from the consensus value (most precise), or takes the mean if that confers
the lowest variation instead, "random" takes a single probe at random from each replicate group.}

\item{expo}{Name of the exposure variable used in the DMR analysis.}

\item{cov}{Name(s) of covariate(s) used in the DMR analysis, with each  name separated by
a comma. Ensure there are no space. e.g. "cov1,cov2,cov3".}

\item{genome}{Reference genome for annotating DMRs.
Must be consistent with the array platform used:
\itemize{
\item Use \code{"hg38"} for EPICV2 arrays.
\item Use \code{"hg19"} for 450K and EPICV1 arrays.
}

Note: the \code{genome} argument does not currently affect the internal behavior of \code{extractRanges()}
(i.e., no liftover is performed).}

\item{fdrCPG}{Used to individually assess the significance of each CpG site. If the FDR-adjusted
p-value of a CpG site is below the specified fdrCPG threshold, the site will be marked as significant.
The default value is 0.05.}

\item{pcutoff}{Used to determine the threshold for DMRs. It is strongly recommended to use the
default (fdr), unless you are confident about the risk of Type I errors (false positives).}

\item{lambda}{If the distance between two significant CpG sites is greater than or equal
to lambda, they will be considered as belonging to different DMRs. The default value is 1000
nucleotides, meaning that if the distance between two significant CpG sites exceeds 1000
nucleotides, they will be separated into different DMRs.}

\item{C}{Scaling factor for bandwidth. Gaussian kernel is calculated where lambda/C = sigma.
Empirical testing shows for both Illumina and bisulfite sequencing data that, when lambda=1000,
near-optimal prediction of sequencing-derived DMRs is obtained when C is approximately 2,
i.e. 1 standard deviation of Gaussian kernel = 500 base pairs. Cannot be < 0.2.}

\item{min.cpgs}{Minimum number of consecutive CpGs constituting a DMR. Default to 2.}

\item{filename}{User-customized .csv file name for storing DMR results. If "default", it will
be named as "DMRresult".}
}
\value{
input, An R6 class object integrating all information.
}
\description{
Perform differential methylation analysis based on the R package \pkg{DMRcate}.
Computes a kernel estimate against a null comparison to identify significantly DMRs.
}
\examples{
\dontrun{
res <- initEWAS(outpath = "default")
res <- loadEWAS(input = res, ExpoData = "default", MethyData = "default")
res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")
res <- dmrEWAS(input = res, filename = "default", chipType = "EPICV2", what = "Beta", expo = "var",
cov = "cov1,cov2", genome = "hg38",)
}

}

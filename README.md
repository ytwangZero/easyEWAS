
# easyEWAS

<!-- badges: start -->
  <!-- badges: end -->
  
  easyEWAS is a flexible and user-friendly R package that systematically performs EWAS analyses under various study designs, along with downstream analyses and result visualization. It can be easily integrated into various DNA methylation microarray detected by Illumina HumanMethylation Bead Chip (27K, 450K, EPICV1, EPICV2, and MSA), significantly enhancing the accessibility of EWAS analysis.

## Installation

You can install the development version of easyEWAS like so:
  
``` r
devtools::install_github("ytwangZero/easyEWAS")
```

💡 **Tip for new users in bioinformatics:**  
If you're not familiar with installing multiple packages manually, we recommend using the provided installation script.  

Since `easyEWAS` depends on many Bioconductor packages, the easiest and safest way to install everything correctly is to download the script from the following link:  

🔗 [Install_easyEWAS.R](https://github.com/ytwangZero/easyEWAS_materials)  

Save this file to your current R working directory, then run:

```r
source("Install_easyEWAS.R")
``` 


## Example

This is an example of performing an EWAS analysis using internal sample data and methylation data with easyEWAS.
  
``` r
library(easyEWAS)
getwd()

# prepare the data file ------
res <- initEWAS(outpath = "default")
res <- loadEWAS(input = res,
                ExpoData = sampledata,
                MethyData = methydata)
res <- transEWAS(input = res, Vars = "cov1", TypeTo = "factor")

# remove batch effect -----
res <- batchEWAS(res,
                 adjustVar = "cov1,cov2",
                 batch = "batch",
                 plot = TRUE,
                 par.prior = TRUE,
                 mean.only = FALSE,
                 ref.batch = NULL)
                 
# perform the EWAS analysis ------
res <- startEWAS(input = res,
                model = "lm",
                expo = "var",
                cov = "cov1,cov2",
                core = "default")

# visualize the EWAS result ------
res <- plotEWAS(input = res,
                file = "jpg",
                p = "PVAL",
                threshold = 0.05)

# internal validation based on the bootstrap method ------
res <- bootEWAS(input = res,
                filterP = "PVAL",
                cutoff = 0.001,
                bootCI = "perc",
                times = 100)

# conduct enrichment analysis ------
res <- enrichEWAS(input = res,
                  method = "GO",
                  filterP = "PVAL",
                  cutoff = 0.05,
                  plot = TRUE,
                  plotType = "dot",
                  plotcolor = "pvalue",
                  showCategory = 10)
                  
# DMR analysis -----
res <- dmrEWAS(input = res,
               chipType = "EPICV2",
               what = "Beta",
               expo = "var",
               cov = "cov1,cov2",
               genome = "hg38",
               lambda=1000,
               C = 2,
               filename = "default",
               pcutoff = 0.05,
               epicv2Filter = "mean")

```


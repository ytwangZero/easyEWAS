# easyEWAS

<div class="home-hero">
  <p class="hero-kicker">Epigenome-Wide Association Study in R</p>
  <p class="hero-lead">easyEWAS is an R package for conducting Epigenome-Wide Association Study (EWAS) in a unified and reproducible way. It supports Illumina methylation array platforms including 27K, 450K, EPIC v1, EPIC v2, and MSA, and provides an end-to-end analysis workflow covering association modeling, batch-effect handling, result visualization, bootstrap-based internal validation, enrichment analysis, and optional DMR discovery.</p>
</div>

## Website
For tutorial, you can refer to https://easyewas-tutorial.github.io/.

## Installation

Before installing `easyEWAS`, we recommend pre-installing optional dependencies used by advanced modules.

### Core dependencies by function

| Function | Required packages | Notes |
|---|---|---|
| `batchEWAS()` | `sva` | Required for ComBat batch correction |
| `batchEWAS(..., parallel = TRUE)` | `sva`, `BiocParallel` | `BiocParallel` is needed only for parallel mode |
| `enrichEWAS()` | `clusterProfiler`, `org.Hs.eg.db` | Required for ID conversion and GO/KEGG enrichment |
| `enrichEWAS(plot = TRUE, plotType = "dot")` | `enrichplot` | Required for dotplot visualization |
| `dmrEWAS()` | `DMRcate` | Required for DMR analysis |

### Recommended pre-install

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "sva",
  "BiocParallel",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "DMRcate"
), ask = FALSE, update = TRUE)
```

Install from GitHub:

```r
remotes::install_github("ytwangZero/easyEWAS")
```

Load package:

```r
library(easyEWAS)
```

If you prefer installing only when needed:

- `batchEWAS()`: `BiocManager::install(c("sva", "BiocParallel"))`
- `enrichEWAS()`: `BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))`
- `dmrEWAS()`: `BiocManager::install("DMRcate")`

## Update: Annotation Data Loading

Large chip annotation tables are now downloaded on demand and cached locally, instead of being bundled inside the main package tarball. This makes installation lighter, faster, and more stable for GitHub users.

Download annotation (one-time per chip type):

```r
downloadAnnotEWAS(chipType = "EPICV2")
```


## Citation

Wang Y, Jiang M, Niu S, Gao X. easyEWAS: a flexible and user-friendly R package for epigenome-wide association study[J]. Bioinformatics Advances, 2025, 5(1): vbaf026.

# easyEWAS

<div class="home-hero">
  <p class="hero-kicker">Epigenome-Wide Association Study in R</p>
  <p class="hero-lead">easyEWAS is an R package for conducting Epigenome-Wide Association Study (EWAS) in a unified and reproducible way. It supports Illumina methylation array platforms including 27K, 450K, EPIC v1, EPIC v2, and MSA, and provides an end-to-end analysis workflow covering association modeling, batch-effect handling, result visualization, bootstrap-based internal validation, enrichment analysis, and optional DMR discovery.</p>
</div>

## Installation

Install from GitHub:

```r
remotes::install_github("ytwangZero/easyEWAS")
```

Load package:

```r
library(easyEWAS)
```

### Optional dependency for DMR analysis

`dmrEWAS()` depends on `DMRcate`, which is in `Suggests` and is not installed automatically in minimal setups.

```r
BiocManager::install("DMRcate")
```

## Update: Annotation Data Loading

Large chip annotation tables are now downloaded on demand and cached locally, instead of being bundled inside the main package tarball. This makes installation lighter, faster, and more stable for GitHub users.

Download annotation (one-time per chip type):

```r
downloadAnnotEWAS(chipType = "EPICV2")
```

If needed, set a custom annotation host:

```r
options(easyEWAS.annotation_base_url = "https://github.com/ytwangZero/easyEWAS_materials/raw/main/annotation")
```

## Citation

Wang Y, Jiang M, Niu S, Gao X. easyEWAS: a flexible and user-friendly R package for epigenome-wide association study[J]. Bioinformatics Advances, 2025, 5(1): vbaf026.

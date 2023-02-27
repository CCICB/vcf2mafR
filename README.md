
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vcf2mafR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/vcf2mafR)](https://CRAN.R-project.org/package=vcf2mafR)
<!-- badges: end -->

**vcf2mafR** converts vcf files and vcf-like dataframes to MAF files

This is an **unofficial** R implementation of VCF to MAF mapping, and is
inspired by the official perl implementation:
[vcf2maf](https://github.com/mskcc/vcf2maf).

**vcf2mafR** implements only a subset of the functionality available in
[vcf2maf](https://github.com/mskcc/vcf2maf), so we recommend any users
who don’t require a native R solution to use the original!

## Installation

You can install the development version of vcf2mafR like so:

``` r
# install.packages('remotes')
remotes::install_github('CCICB/vcf2mafR')
```

## Conversions

-   `df2maf` : Converts VCF-like data.frames to MAF format
-   `vcf2maf`: Converts VCF files to MAF format

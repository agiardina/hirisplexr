# hirisplexr

Build HIrisPlex / HIrisPlex-S CSV files directly from PLINK 1.9 BED/BIM/FAM.

## Installation (from source)

```r
# in R
install.packages("BEDMatrix")
install.packages("data.table")

# build & install this package (using R CMD build/check/INSTALL or devtools)
```

## Usage

```r
library(hirisplexr)

# Example prefix chosen for documentation purposes
# (provide your own PLINK files at this path)
prefix <- "inst/extdata/testprefix"

write_hirisplex_csv(prefix, panel = "hirisplexs", out = "hirisplexs.csv")
```

## How it works

- Genotypes are read on-demand from `.bed` using [BEDMatrix], which returns
  dosages of **A1** (first allele in `.bim`).
- The HIrisPlex(-S) webtool requires a count (0/1/2) of a specific **input allele**
  for each SNP. We map `A1`/`A2` to that input allele; for strand issues we can
  also consider complements.

## Panels

Panel definitions (rsID and required input allele, in the precise order expected
by the webtool) are packaged in `inst/extdata/hirisplex_panels.csv`. They are
sourced from the official webtool and manual (see citations below).

## Citations

- HIrisPlex-S webtool: https://hirisplex.erasmusmc.nl/ (columns and upload format)
- HIrisPlex-S user manual (2018), Erasmus MC (instructions and caveats)
  - Notes include: count of 0/1/2 per SNP; use `NA` if missing; `rs6119471` may
    require strand conversion from NCBI forward.


# SNPopt (SNPoptimizer)

A tool to choose a minimal subset of SNP to discriminate samples. SNPopt is an R/Shiny application that uses a Genetic Algorithm (GA) to select a minimal subset of SNPs able to maximally discriminate samples in large genotypic matrices (e.g., for cultivar identification, panel design, and downstream assay development).

## Key features
- GA-based SNP panel optimization to maximize sample discrimination (e.g., uniqueness / diversity metrics).
- Interactive Shiny interface to explore results, duplicates, and discriminatory power.
- Optional second-round optimization on unresolved duplicates to add extra SNPs and obtain a fully discriminant final matrix.
- Export utilities (selected SNP list, final genotype matrix, duplicated sample list, summary reports).

## Installation
### Option A — Run locally (recommended)
```r
# install dependencies
install.packages(c("shiny", "data.table", "ggplot2"))
# + other dependencies used in your project

# run
shiny::runApp("app")


## Demo and benchmark datasets

Due to size constraints, full HapMap (HMP) benchmark datasets
used in the manuscript are not stored directly in this repository.

They are publicly available on Zenodo at:
**10.5281/zenodo.17975237**


# SNPOptimizer (SNPOpt) 

A tool to choose a minimal subset of SNP to discriminate samples. SNPopt is an R/Shiny application that uses a Genetic Algorithm (GA) to select a minimal subset of SNPs able to maximally discriminate samples in large genotypic matrices (e.g., for cultivar identification, panel design, and downstream assay development).

## Key features
- GA-based SNP panel optimization to maximize sample discrimination (e.g., uniqueness / diversity metrics).
- Interactive Shiny interface to explore results, duplicates, and discriminatory power.
- Optional second-round optimization on unresolved duplicates to add extra SNPs and obtain a fully discriminant final matrix.
- Export utilities (selected SNP list, final genotype matrix, duplicated sample list, summary reports).

## Demo and benchmark datasets

Due to size constraints, full HapMap (HMP) benchmark datasets
used in the manuscript are not stored directly in this repository. They are publicly available on Zenodo at:
**10.5281/zenodo.17975237**

## Run locally

### Requirements
- R (>= 4.2)  
- (Recommended) RStudio  
- On Windows: Rtools may be required for some packages.

### 1) Download the repository
Clone or download this repository and open the project folder.

### 2) Install dependencies and run
```r
source("install.R")

source("run.R")

If you keep the app inside a subfolder (e.g. app/), run:

shiny::runApp("app")

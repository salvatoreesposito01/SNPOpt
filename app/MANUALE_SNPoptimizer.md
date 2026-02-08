% SNPoptimizer - User Manual
% Updated on 2026-02-08

# 1. Overview
SNPoptimizer is a Shiny app to:
- select a minimal SNP set that maximizes sample discrimination;
- evaluate an already defined SNP set;
- generate flanking sequences and Primer3 input files.

The app supports **HapMap** and **VCF** inputs (biallelic variants with `GT`).

# 2. Start the app
From the project folder:

```r
shiny::runApp()
```

Alternative:

```r
source("run.R")
```

Upload size limit is set to **1 GB**.

# 3. Discovery section
## 3.1 Input
- Upload HapMap (`.txt/.tsv/.hmp`) or VCF (`.vcf/.vcf.gz`).
- Set SNP heterozygosity threshold.
- Set `k` (number of SNPs to select).

## 3.2 Selection method
Two methods are available:
- **Hybrid (greedy + local search)**: usually faster.
- **GA**: classical genetic algorithm with `popSize` and `maxiter`.

## 3.3 Reproducibility
You can enable a fixed seed (`useSeed`) for reproducible results.

## 3.4 Additional run (optional)
If duplicates remain, you can run a second search on duplicate samples only to add extra SNPs.

# 4. Discovery outputs
The app shows:
- KPIs (filtered SNPs, selected SNPs, samples, fitness, duplicates);
- method used (Hybrid/GA);
- SNP discrimination curve;
- selected SNP table;
- remaining duplicate sample groups.

Available downloads:
- selected SNP metadata;
- selected genotype matrix;
- discovery summary;
- curve plot;
- duplicated sample list.

# 5. Evaluation (given SNPs) section
This section evaluates a SNP list (CHR, POS) on:
- a dedicated evaluation dataset, or
- the dataset already loaded in Discovery.

Supported SNP-list formats:
- file with `CHR` and `POS` columns;
- free text, e.g. `chr1:12345` or `1 12345`.

Output:
- requested/found SNP counts;
- list of missing SNPs;
- downloadable subset matrix.

# 6. Post-discovery: Flanking and Primer3
## 6.1 Supported FASTA
- `.fa`, `.fasta`, `.fa.gz`

## 6.2 Performance
If `samtools` is available and FASTA is not `.gz`, the app uses `faidx` for fast extraction.

## 6.3 Chromosome matching
CHR matching is robust for common naming variants (e.g. `1`, `chr1`, `SL4.0ch01`).
If matching fails, the app shows a warning with FASTA header preview.

## 6.4 New skipped-SNP diagnostics
During flanking generation, some SNPs may be excluded (e.g. unresolved CHR, out-of-range POS).
The app now:
- shows a notification with skipped SNP count;
- provides `Download skipped SNPs` (TSV) with skip reasons.

This prevents silent SNP loss in outputs.

# 7. Common issues and fixes
## 7.1 "Maximum upload size exceeded"
- Use the updated app version (1 GB limit);
- verify you are launching the app from the correct folder.

## 7.2 "No requested chromosomes were found in the FASTA"
- check assembly and chromosome naming consistency;
- verify required chromosomes exist in the FASTA;
- if available, use indexed FASTA (`samtools faidx`).

## 7.3 Flanking output has fewer SNPs than expected
- download `Download skipped SNPs` to see which SNPs were excluded and why.

# 8. Requirements
Main R packages:
- shiny
- GA
- data.table
- ggplot2

Optional external tool:
- samtools (recommended for large FASTA files)

# 9. Method note (paper)
- **GA** mode remains available for GA-based workflows and reporting.
- **Hybrid** mode is a practical faster alternative for many datasets.
- In the paper, report which optimizer was used for each analysis.

# 10. Main project files
- `app.R`
- `ui.R`
- `server.R`
- `run.R`
- `install.R`


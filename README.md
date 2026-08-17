# SNPoptimizer

SNPoptimizer is a Shiny app to discover and evaluate minimal SNP panels for genotype discrimination.

## Key Features

- **Discovery mode** for selecting a minimal informative SNP set.
- **Two optimizers**:
  - `Hybrid (greedy + local search)` (usually faster)
  - `GA` (genetic algorithm)
- **Optional additional run** on unresolved duplicate sample groups.
- **Evaluation mode** for checking a user-provided SNP list (`CHR`, `POS`).
- **Post-discovery flanking export** from reference FASTA (`.fa`, `.fasta`, `.fa.gz`).
- **Primer3 export** generation.
- **Robust chromosome-name matching** across common naming styles (e.g. `1`, `chr1`, `SL4.0ch01`).
- **Skipped-SNP diagnostics** during flanking extraction with downloadable report.

## Input Support

- **HapMap** (`.txt`, `.tsv`, `.hmp`)
- **VCF** (`.vcf`, `.vcf.gz`) with biallelic SNPs and `GT`

## Requirements

R packages:

- `shiny`
- `GA`
- `data.table`
- `ggplot2`

Optional (recommended for large FASTA):

- `samtools` (for `faidx`-based fast extraction)

## Run the App

From the repository root:

```r
shiny::runApp("app")
```

or:

```r
source("app/run.R")
```

Upload limit is configured to **1 GB**.

## Discovery Outputs

- Discovery KPI summary
- SNP discrimination curve
- Selected SNP table
- Duplicate sample groups
- Downloads:
  - selected SNP metadata
  - selected genotype matrix
  - discovery summary
  - curve plot
  - duplicate sample list

## Flanking / Primer3 Outputs

After Discovery:

- `Download flanking FASTA`
- `Download Primer3 file`
- `Download skipped SNPs` (TSV with reasons, e.g. unresolved chromosome, out-of-range position)

## Notes for Manuscripts

- `GA` mode is available for GA-based reporting.
- `Hybrid` mode is a practical faster alternative on many datasets.
- Report the optimizer used for each analysis.


## Citation

If you use SNPoptimizer, please cite the associated manuscript and this software repository:

Esposito, S. (2026). Molecular Breeding

## Project Files

- `app.R`
- `ui.R`
- `server.R`
- `run.R`
- `install.R`
- `MANUALE_SNPoptimizer.md`
- `MANUALE_SNPoptimizer.pdf`

## License

See `LICENSE` in the repository root.

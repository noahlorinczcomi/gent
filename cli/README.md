# gent CLI

Command-line interface for **gent** — Gene-based Association Testing using GWAS summary statistics.

## TL;DR

```bash
git clone https://github.com/noahlorinczcomi/gent.git
cd cli
pixi install
pixi run fetch-ld EUR ./ld_matrices
pixi run gent \
  --gwas ~/data/my_gwas.tsv \
  --ld-dir ~/data/ld_matrices \
  --ld-pop EUR \
  --build grch37 \
  --snp MarkerName \
  --chr Chromosome \
  --pos Position \
  --ea Effect_allele \
  --z z \
  --out results/gent_output.csv
```

## Requirements

- [Git](https://git-scm.com/)
- [pixi](https://pixi.sh/) — cross-platform package manager (installs all R dependencies automatically)

## Installation

Clone the repository and navigate to the CLI directory:

```bash
git clone https://github.com/noahlorinczcomi/gent.git
cd gent/cli
pixi install
```

`pixi install` creates an isolated environment containing R and all required packages (`data.table`, `dplyr`, `tidyr`, `mvnfast`, `susieR`, `ggplot2`, `ggrepel`, `optparse`). No system-wide R installation is needed.

## LD Reference Data

`gent_genomewide` requires pre-computed gene-specific LD matrices (see `fetch-ld` below). After extraction the directory will have this structure:

```
ld_matrices/
  EUR/
    chr1/
      GENE1.Rds
      GENE2.Rds
      ...
    chr2/
    ...
  AFR/
  EAS/
  SAS/
  AMR/
```

## Tasks

| Task | Description |
|------|-------------|
| `pixi run fetch-ld <pop> [dir]` | Download and extract LD reference matrices for a population |
| `pixi run gent [options]` | Run genome-wide gene-based association testing |

## Usage

All commands are run from within the `gent/cli/` directory.

### `fetch-ld`

```bash
pixi run fetch-ld <population> [output_dir]
```

Downloads and extracts LD reference matrices from Zenodo. `output_dir` defaults to the current directory.

```bash
pixi run fetch-ld EUR ./ld_matrices     # single population
pixi run fetch-ld AFR ./ld_matrices     # AFR (two-part download, merged automatically)
pixi run fetch-ld ALL ./ld_matrices     # all populations
```

Supported populations: `EUR`, `EAS`, `SAS`, `AMR`, `AFR`, `ALL` (case-insensitive).

### `gent`

```bash
pixi run gent --gwas <file> --ld-dir <dir> --ld-pop <pop> [options]
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--test` | `gent` | Test type: `gent` (single-trait) or `mugent` (multi-trait/ancestry) |
| `--gwas` | — | Path to GWAS file. For `mugent`, comma-separated list (e.g. `eur.tsv,afr.tsv`) |
| `--ld-dir` | — | Path to LD reference directory |
| `--ld-pop` | `EUR` | LD population: `EUR`, `AFR`, `EAS`, `SAS`, `AMR`. For `mugent`, comma-separated |
| `--build` | `grch37` | Genomic build: `grch37` (hg19) or `grch38` (hg38) |
| `--kb-window` | `50` | Kilobase window around each gene boundary |
| `--snp` | `rsid` | Column name for SNP rsID |
| `--chr` | `chr` | Column name for chromosome |
| `--pos` | `position` | Column name for base pair position |
| `--ea` | `effect_allele` | Column name for effect allele |
| `--z` | `z` | Column name for Z-statistic |
| `--out` | `gent_results.csv` | Output CSV file path |
| `--verbose` | `TRUE` | Print progress |

GWAS files are read with `data.table::fread()`, which auto-detects tab/comma separators.

### Single-trait genome-wide GenT

```bash
pixi run gent \
  --gwas ~/data/my_gwas.tsv \
  --ld-dir ~/data/ld_matrices \
  --ld-pop EUR \
  --build grch37 \
  --snp MarkerName \
  --chr Chromosome \
  --pos Position \
  --ea Effect_allele \
  --z z \
  --out results/gent_output.csv
```

Output columns: `pval`, `shape`, `rate`, `mu_h0`, `sigma2_h0`, `gene`, `m`, `chr`, `gene_start`, `window_start`, `gene_end`, `window_end`.

### Multi-ancestry genome-wide MuGenT

```bash
pixi run gent \
  --test mugent \
  --gwas ~/data/eur_gwas.tsv,~/data/afr_gwas.tsv \
  --ld-dir ~/data/ld_matrices \
  --ld-pop EUR,AFR \
  --build grch37 \
  --out results/mugent_output.csv
```

Three output files are written:
- `results/mugent_output_MuGenT.csv` — multi-ancestry association test
- `results/mugent_output_MuGenT-PH.csv` — population heterogeneity test
- `results/mugent_output_MuGenT-Pleiotropy.csv` — pleiotropy test

## Help

```bash
pixi run gent --help
```

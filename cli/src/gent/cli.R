library(optparse)

option_list <- list(
  make_option(
    c("--test"), type = "character", default = "gent",
    help = "Test to run: 'gent' (single-trait genome-wide) or 'mugent' (multi-trait/ancestry genome-wide) [default: %default]"
  ),
  make_option(
    c("--gwas"), type = "character", default = NULL,
    help = "Path(s) to GWAS summary statistics file(s). For mugent, supply comma-separated paths (e.g. eur.tsv,afr.tsv). Files are read by data.table::fread() and auto-detect separators."
  ),
  make_option(
    c("--ld-dir"), type = "character", default = NULL,
    help = "Path to the LD reference directory. Must contain population sub-directories (e.g. EUR/, AFR/) each with per-chromosome sub-directories (e.g. chr1/) containing gene-specific .Rds LD matrices."
  ),
  make_option(
    c("--ld-pop"), type = "character", default = "EUR",
    help = "LD reference population(s): EUR, AFR, EAS, SAS, or AMR. For mugent, supply comma-separated values matching the order of --gwas [default: %default]"
  ),
  make_option(
    c("--build"), type = "character", default = "grch37",
    help = "Genomic build of the GWAS data: grch37 (hg19) or grch38 (hg38) [default: %default]"
  ),
  make_option(
    c("--kb-window"), type = "integer", default = 50,
    help = "Kilobase window added around each gene boundary when assigning SNPs [default: %default]"
  ),
  make_option(
    c("--snp"), type = "character", default = "rsid",
    help = "Column name for SNP rsID in GWAS file(s) [default: %default]"
  ),
  make_option(
    c("--chr"), type = "character", default = "chr",
    help = "Column name for chromosome in GWAS file(s) [default: %default]"
  ),
  make_option(
    c("--pos"), type = "character", default = "position",
    help = "Column name for base pair position in GWAS file(s) [default: %default]"
  ),
  make_option(
    c("--ea"), type = "character", default = "effect_allele",
    help = "Column name for effect allele in GWAS file(s) [default: %default]"
  ),
  make_option(
    c("--z"), type = "character", default = "z",
    help = "Column name for Z-statistic in GWAS file(s) [default: %default]"
  ),
  make_option(
    c("--out"), type = "character", default = "gent_results.csv",
    help = "Output file path for results (CSV). For mugent, three files are written with the suffixes _MuGenT.csv, _MuGenT-PH.csv, and _MuGenT-Pleiotropy.csv [default: %default]"
  ),
  make_option(
    c("--verbose"), action = "store_true", default = TRUE,
    help = "Print progress to the console [default: TRUE]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = paste(
    "gent CLI: Gene-based Association Testing",
    "",
    "Runs genome-wide gene-based association tests using GWAS summary statistics",
    "and pre-computed LD reference matrices. See the README for details on the",
    "required LD directory structure and how to download reference data.",
    sep = "\n"
  ),
  usage = "pixi run run-gent --gwas <file> --ld-dir <dir> --ld-pop <pop> [options]"
)

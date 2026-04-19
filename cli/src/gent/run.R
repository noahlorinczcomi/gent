#!/usr/bin/env Rscript

# Resolve the directory containing this script so paths stay correct regardless
# of the working directory from which `pixi run run-gent` is invoked.
args0 <- commandArgs(FALSE)
file_arg <- args0[grep("--file=", args0)]
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", file_arg[1])))
} else {
  script_dir <- getwd()
}

gent_root     <- normalizePath(file.path(script_dir, "../../../"))
gent_r_dir    <- file.path(gent_root, "R")
gent_data_dir <- file.path(gent_root, "data")

source(file.path(script_dir, "cli.R"))
source(file.path(script_dir, "utils.R"))

# Install a data() override in the global env so that calls like
# data(EnsemblHg19GenePos) inside functions.r resolve from local .rda files.
data <- make_data_loader(gent_data_dir)

load_gent_functions(gent_r_dir)

# ── Parse arguments ────────────────────────────────────────────────────────────
opt <- parse_args(opt_parser)

if (is.null(opt$gwas)) {
  stop("--gwas is required. Run with --help for usage.")
}
if (is.null(opt[["ld-dir"]])) {
  stop("--ld-dir is required. Run with --help for usage.")
}

test_type <- tolower(opt$test)
if (!test_type %in% c("gent", "mugent")) {
  stop("--test must be 'gent' or 'mugent'.")
}

# ── Run GenT (single-trait genome-wide) ───────────────────────────────────────
if (test_type == "gent") {
  gwas_paths <- split_csv_arg(opt$gwas)
  if (length(gwas_paths) != 1) {
    stop("--test gent expects exactly one --gwas file path.")
  }
  ld_pops <- split_csv_arg(opt[["ld-pop"]])
  if (length(ld_pops) != 1) {
    stop("--test gent expects exactly one --ld-pop value.")
  }

  message("Reading GWAS data from: ", gwas_paths)
  gwas <- read_gwas(gwas_paths)

  results <- gent_genomewide(
    gwas           = gwas,
    KbWindow       = opt[["kb-window"]],
    ld_population  = ld_pops,
    ld_directory   = opt[["ld-dir"]],
    build          = opt$build,
    snp            = opt$snp,
    chromosome     = opt$chr,
    position       = opt$pos,
    effect_allele  = opt$ea,
    z_statistic    = opt$z,
    verbose        = opt$verbose
  )

  write_results(results, opt$out, "gent")
}

# ── Run MuGenT (multi-trait/ancestry genome-wide) ─────────────────────────────
if (test_type == "mugent") {
  gwas_paths <- split_csv_arg(opt$gwas)
  ld_pops    <- split_csv_arg(opt[["ld-pop"]])

  if (length(gwas_paths) < 2) {
    stop("--test mugent requires at least two comma-separated --gwas file paths.")
  }
  if (length(ld_pops) != length(gwas_paths)) {
    stop("--ld-pop must supply one population per --gwas file (comma-separated, same order).")
  }

  message("Reading GWAS data from: ", paste(gwas_paths, collapse = ", "))
  gwas_list <- lapply(gwas_paths, read_gwas)
  names(gwas_list) <- ld_pops

  col_names <- list(
    snp            = opt$snp,
    chromosome     = opt$chr,
    position       = opt$pos,
    effect_allele  = opt$ea,
    z_statistic    = opt$z
  )
  make_list <- function(val) setNames(lapply(seq_along(ld_pops), function(i) val), ld_pops)

  results <- mugent_genomewide(
    gwas_list            = gwas_list,
    ld_population_list   = as.list(setNames(ld_pops, ld_pops)),
    ld_directory         = opt[["ld-dir"]],
    build                = opt$build,
    KbWindow             = opt[["kb-window"]],
    snp_list             = make_list(opt$snp),
    chromosome_list      = make_list(opt$chr),
    position_list        = make_list(opt$pos),
    effect_allele_list   = make_list(opt$ea),
    z_statistic_list     = make_list(opt$z),
    verbose              = opt$verbose
  )

  write_results(results, opt$out, "mugent")
}

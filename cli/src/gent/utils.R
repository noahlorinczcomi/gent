# gent_root and gent_data_dir are set by run.R before this file is sourced.

load_gent_functions <- function(r_dir) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(tibble)
  })
  source(file.path(r_dir, "helpers.r"), local = FALSE)
  source(file.path(r_dir, "functions.r"), local = FALSE)
}

# Override base::data() so that calls inside functions.r (e.g. data(EnsemblHg19GenePos))
# load from local .rda files rather than requiring the gent package to be installed.
make_data_loader <- function(data_dir) {
  function(..., list = character(), envir = .GlobalEnv) {
    names <- as.character(match.call(expand.dots = FALSE)$`...`)
    names <- c(names, list)
    for (nm in names) {
      rda_path <- file.path(data_dir, paste0(nm, ".rda"))
      if (file.exists(rda_path)) {
        load(rda_path, envir = envir)
      } else {
        base::data(list = nm, envir = envir)
      }
    }
    invisible(NULL)
  }
}

read_gwas <- function(path) {
  as.data.frame(data.table::fread(path, data.table = FALSE))
}

split_csv_arg <- function(x) {
  trimws(strsplit(x, ",")[[1]])
}

write_results <- function(results, out_path, test_type) {
  if (test_type == "gent") {
    write.csv(results$result, out_path, row.names = FALSE, quote = FALSE)
    message("Results written to: ", out_path)
  } else {
    base_path <- sub("\\.csv$", "", out_path, ignore.case = TRUE)
    for (nm in names(results)) {
      fp <- paste0(base_path, "_", nm, ".csv")
      write.csv(results[[nm]], fp, row.names = FALSE, quote = FALSE)
      message("Results written to: ", fp)
    }
  }
}

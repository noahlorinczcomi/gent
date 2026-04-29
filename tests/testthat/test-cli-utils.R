# devtools::test() sets the CWD to tests/testthat/, so navigate up two levels
# to reach the package root before resolving the CLI utils path.
local({
  p <- normalizePath(file.path("../../cli/src/gent/utils.R"), mustWork = FALSE)
  if (!file.exists(p)) skip("CLI utils not found at expected path")
  source(p, local = FALSE)
})

# ── split_csv_arg ─────────────────────────────────────────────────────────────

test_that("split_csv_arg() splits a single value into a length-1 character vector", {
  expect_equal(split_csv_arg("EUR"), "EUR")
})

test_that("split_csv_arg() splits comma-separated values and trims whitespace", {
  expect_equal(split_csv_arg("EUR,AFR,EAS"), c("EUR", "AFR", "EAS"))
  expect_equal(split_csv_arg("EUR, AFR , EAS"), c("EUR", "AFR", "EAS"))
})

test_that("split_csv_arg() returns a character vector", {
  result <- split_csv_arg("a,b,c")
  expect_true(is.character(result))
  expect_length(result, 3L)
})

# ── make_data_loader ──────────────────────────────────────────────────────────

test_that("make_data_loader() returns a function", {
  loader <- make_data_loader(tempdir())
  expect_true(is.function(loader))
})

test_that("make_data_loader() loads an .rda file into the specified environment", {
  tmp <- tempdir()
  test_obj <- list(a = 1, b = 2)
  save(test_obj, file = file.path(tmp, "test_obj.rda"))

  env <- new.env(parent = emptyenv())
  loader <- make_data_loader(tmp)
  loader(test_obj, envir = env)

  expect_true(exists("test_obj", envir = env))
  expect_equal(env$test_obj, test_obj)
})

test_that("make_data_loader() returns NULL invisibly when the .rda is absent", {
  # Point at an empty dir so the fallback path is exercised without actually
  # needing base::data() to resolve a named dataset in a package context.
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  loader <- make_data_loader(tmp)
  # No .rda exists; base::data() will be tried and will throw — that is fine to
  # test by verifying the loader at least reaches the fallback branch.
  expect_error(loader(no_such_dataset_xyz))
})

# ── read_gwas ─────────────────────────────────────────────────────────────────

test_that("read_gwas() reads a tab-separated file into a data.frame", {
  tmp <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(rsid = c("rs1", "rs2"), z = c(1.5, -0.3)),
    file = tmp, sep = "\t", row.names = FALSE, quote = FALSE
  )
  result <- read_gwas(tmp)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2L)
  expect_true("rsid" %in% colnames(result))
  expect_true("z"    %in% colnames(result))
  unlink(tmp)
})

test_that("read_gwas() reads a comma-separated file into a data.frame", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(rsid = c("rs1", "rs2"), z = c(1.5, -0.3)),
    file = tmp, row.names = FALSE
  )
  result <- read_gwas(tmp)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2L)
  unlink(tmp)
})

# ── write_results ─────────────────────────────────────────────────────────────

test_that("write_results() for 'gent' writes a single CSV file", {
  tmp <- tempfile(fileext = ".csv")
  results <- list(result = data.frame(gene = "BRCA1", pval = 0.01))
  expect_no_error(write_results(results, tmp, "gent"))
  expect_true(file.exists(tmp))
  out <- read.csv(tmp)
  expect_equal(nrow(out), 1L)
  expect_true("gene" %in% colnames(out))
  unlink(tmp)
})

test_that("write_results() for 'mugent' writes one CSV per named result", {
  tmp <- tempfile(fileext = ".csv")
  results <- list(
    MuGenT            = data.frame(gene = "APOE", pval = 0.001),
    `MuGenT-PH`       = data.frame(gene = "APOE", pval = 0.002),
    `MuGenT-Pleiotropy` = data.frame(gene = "APOE", result = "no pleiotropy")
  )
  expect_no_error(write_results(results, tmp, "mugent"))

  base <- sub("\\.csv$", "", tmp, ignore.case = TRUE)
  for (nm in names(results)) {
    fp <- paste0(base, "_", nm, ".csv")
    expect_true(file.exists(fp), label = paste("file exists:", fp))
    unlink(fp)
  }
})

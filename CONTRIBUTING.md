# Contributing to gent

Thank you for your interest in improving gent. This document covers everything you need to go from a fresh clone to a passing CI run.

## Table of contents

1. [Dev environment setup](#dev-environment-setup)
2. [Running the tests](#running-the-tests)
3. [Code style](#code-style)
4. [Submitting a pull request](#submitting-a-pull-request)
5. [Reporting bugs and requesting features](#reporting-bugs-and-requesting-features)

---

## Dev environment setup

### Prerequisites

| Tool | Minimum version | Notes |
|------|----------------|-------|
| R | 4.3 | Install from [r-project.org](https://www.r-project.org/) |
| renv | 1.0+ | `install.packages("renv")` |
| devtools | any | `install.packages("devtools")` |

### First-time setup

```r
# 1 – Clone the repo then open R inside it
# 2 – Restore the exact package versions recorded in renv.lock
renv::restore()

# 3 – Load the package in development mode
devtools::load_all()
```

`renv::restore()` installs all dependencies listed in `renv.lock` into a project-local library so your global R installation is not affected.

### Updating dependencies

If you add or remove a package dependency, record the change before committing:

```r
# Add the package to DESCRIPTION Imports/Suggests first, then:
renv::snapshot()
```

Commit both `DESCRIPTION` and the updated `renv.lock` together.

---

## Running the tests

The test suite uses [testthat](https://testthat.r-lib.org/) edition 3.

```r
# Run all tests
devtools::test()

# Run a single test file
devtools::test(filter = "gent")     # matches tests/testthat/test-gent.R
devtools::test(filter = "helpers")  # matches tests/testthat/test-helpers.R
```

### R CMD check

Before opening a PR, run the full package check locally:

```r
devtools::check(args = c("--no-manual", "--no-vignettes"))
```

The CI workflow fails on any WARNING or ERROR, so confirm there are none locally first.

---

## Code style

### General

- Follow the existing conventions in each file you touch. Consistency within a file beats stylistic perfection across files.
- Keep lines under 100 characters where practical.
- Use `<-` for assignment, not `=` (except inside function argument lists).

### Functions

- Exported functions live in `R/functions.r`. Internal helpers belong in `R/helpers.r`.
- Every exported function must have a roxygen2 doc block with `@param`, `@return`, `@export`, and at least one `@examples` entry.
- Avoid `library()` or `require()` inside package functions. Use `::` for calls to other packages (e.g. `dplyr::filter()`), and list the package in `DESCRIPTION` under `Imports`.

### Tests

- One test file per source file: `test-gent.R` for `functions.r`, `test-helpers.R` for `helpers.r`, etc.
- Test names should read as plain English sentences: `"gent() p-value is in [0, 1]"`.
- Use `set.seed()` when a test depends on specific random values.
- Internal functions are accessed via `gent:::fn()` in tests.

### CLI (`cli/`)

- CLI source files follow the layout described in `cli/CLAUDE.md`.
- New CLI options go in `cli/src/gent/cli.R`; new utilities go in `cli/src/gent/utils.R`.
- Tests for CLI utilities live in `tests/testthat/test-cli-utils.R`.

---

## Submitting a pull request

1. **Fork** the repo and create a branch off `main` (`git checkout -b my-feature`).
2. Make your changes and add tests covering every new or changed behaviour.
3. Verify `devtools::test()` and `devtools::check(args = "--as-cran")` both pass with no warnings or errors.
4. Update `renv.lock` if you changed any dependencies (`renv::snapshot()`).
5. Open a PR against `main`. The PR template will prompt you for the key information.

The CI workflow runs `R CMD check --as-cran` on Ubuntu and macOS on every push and PR. A red CI badge means the PR is not mergeable.

---

## Reporting bugs and requesting features

Use the GitHub issue tracker:

- **Bug** — choose the *Bug report* template. Include a minimal reproducible example.
- **Feature request** — choose the *Feature request* template. Explain the motivation before the solution.

For security-sensitive issues, email the maintainer directly (see `DESCRIPTION`) rather than opening a public issue.

# Goal

`gent` currently exists as an R package installable using the `devtools` or `remotes` R packages, but I want users to be able to run the `gent` functionality as a CLI tool.

## Implementation

Generate a `pixi` environment in `gent/cli` that contains `r-base` and all dependent R packages in `gent`. Use the `optparse` package to create a `run.R` script that simply calls the pre-existing `gent/R/*.R` scripts to perform the analysis.

## Technical details

Use the `optparse` package and lay out the `cli/` directory like this:

- `README.md`: instructions for how to use the `gent` CLI
    - i.e., requiring cloning the repo and installing the `pixi` environment
- `src/gent/run.R`: running the CLI
- `src/gent/cli.R`: arguments passed at the command line
- `src/gent/utils.R`: any utils scripts you need to load functions from `gent/R/*R`

 Use a `pixi` task named `run-gent` to run the CLI.
 

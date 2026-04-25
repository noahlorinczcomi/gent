---
name: Bug report
about: Something isn't working as expected
title: '[Bug] '
labels: bug
assignees: ''
---

## Describe the bug

A clear description of what the bug is.

## Minimal reproducible example

```r
# Paste the smallest R snippet that triggers the bug
library(gent)

LD <- diag(5)
zs <- c(1, 2, 3, 4, 5)
gent(zs, LD)
```

## Expected behaviour

What did you expect to happen?

## Actual behaviour

What happened instead? Include the full error message or unexpected output.

## Environment

- **gent version** (run `packageVersion("gent")`):
- **R version** (run `R.version.string`):
- **OS**:
- **Installed via** (devtools/remotes/CRAN):

## Additional context

Any other context, screenshots, or data that helps explain the problem.

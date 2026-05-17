[![DOI](https://zenodo.org/badge/839937647.svg)](https://doi.org/10.5281/zenodo.17449854)

> [!NOTE]
> `gent` can be run from the command line using the [`cli/`](cli/) repo. Follow the [link](cli/) for more details.

## TL;DR

```R
library(gent)
LD_matrix = 0.5^toeplitz(0:9)
z = MASS::mvrnorm(n=1, mu=rep(0,10), Sigma=LD_matrix)
result = gent(
    zs=z,         # vector of variant Z-statistics
    LD=LD_matrix  # LD correlated matrix allele-harmonized to ``z_statistics_vector`` 
)
```

## Summary
The methods contained in this R package test the association between a set of SNPs and a phenotype. These methods only require GWAS summary statistics and an LD reference panel. If you are performing xGenT (xQTL-weighted gene-based testing), you also need xQTL effect sizes. When SNP sets are gene-specific (e.g., containing SNPs near a gene), we refer to our set of methods as 'gene-based association' tests. 

[Pre-computed results for 50+ phenotypes (Shiny app)](https://nlorinczcomi.shinyapps.io/gent/)

## R installation
```
remotes::install_github('noahlorinczcomi/gent')
```

## How to use `gent`
Jump to examples for testing a **single gene**:
* [GenT](https://github.com/noahlorinczcomi/gent/wiki/Single%E2%80%90gene-GenT) `gent::gent()`
* [xGenT](https://github.com/noahlorinczcomi/gent/wiki/Single%E2%80%90gene-xGenT) `gent::gent()`
* [MuGenT](https://github.com/noahlorinczcomi/gent/wiki/Single%E2%80%90gene-MuGenT) `gent::mugent()`
* [MuGenT Population heterogeneity](https://github.com/noahlorinczcomi/gent/wiki/Single%E2%80%90gene-MuGenT%E2%80%90PH) `gent::mugent_ph()`
* [MuGenT Pleiotropy](https://github.com/noahlorinczcomi/gent/wiki/Single%E2%80%90gene-MuGenT%E2%80%90Pleiotropy) `gent::mugent_pleio()`

Jump to examples for testing **genome-wide**:
* [GenT](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-GenT) `gent::gent_genomewide()`
* [xGenT](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-xGenT) `gent::xgent_genomewide()`
* [MuGenT + Post-hoc Heterogeneity and Pleiotropy tests](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-MuGenT-&--post%E2%80%90hoc-tests) `gent::mugent_genomewide()`
* [Fine-mapping gene-based test statistics](https://github.com/noahlorinczcomi/gent/wiki/Fine%E2%80%90mapping-gene-statistics) `gent::gent_finemap()`

# References

- Lorincz-Comi, N., Song, W., Chen, X., Paz, I. R., Hou, Y., Zhou, Y., ... & Cheng, F. (2026). Combining xQTL and genome-wide association studies from diverse populations improves druggable gene discovery. Nature Communications.

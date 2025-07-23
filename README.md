## Summary
The methods contained in this R package test the association between a set of SNPs and a phenotype. These methods only require GWAS summary statistics and an LD reference panel. When SNP sets are gene-specific (e.g., containing SNPs near a gene), we refer to our set of methods as 'gene-based association' tests. 

[Our preprint](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5080346)

[Pre-computed results for 50+ phenotypes](https://nlorinczcomi.shinyapps.io/gent/)

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
* [MuGenT](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-MuGenT) `gent::mugent_genomewide()`
* [MuGenT Population heterogeneity](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-MuGenT%E2%80%90PH) `gent::mugent_ph_genomewide()`
* [MuGenT Pleiotropy](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-MuGenT%E2%80%90Pleiotropy) `gent::mugent_pleio_genomewide()`
* [Fine-mapping gene-based test statistics](https://github.com/noahlorinczcomi/gent/wiki/Genome%E2%80%90wide-GenT%E2%80%90Finemap) `gent::gent_pips_genomewide()`


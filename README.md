# Installation
```
devtools::install_github('noahlorinczcomi/gent')
# or
remotes::install_github('noahlorinczcomi/gent')
```

We also have a Shiny application for GenT applied to 50+ complex diseases here: [https://nlorinczcomi.shinyapps.io/gent/](https://nlorinczcomi.shinyapps.io/gent/).

# Summary
The methods contained in this R package perform joint tests of association on sets of SNPs which may be in linkage disequilibrium (LD) with each other. When SNP sets are gene-specific (e.g., containing SNPs near a gene), we refer to our set of methods as 'gene-based association' tests. When LD between SNPs in each set is known, these methods are exact. When LD is estimated without bias, these methods have controlled false positive and negative rates in simulations (see our [preprint](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5080346)) though no theoretical guarantee. The methods available in this R package are
* `gent()`: Standard gene-based association test, or an xQTL-integrated gene-based association test (weights are xQTL effect sizes).
* `mugent()`: Multi-ancestry gene-based association test.
* `mugent_ph()`: Gene-based association test of heterogeneity in effect size between multiple populations.
* `mugent_pleio()`: Gene-based test that the gene is associated with the phenotype in all populations.

All methods use only GWAS summary statistics and an LD reference panel. The dosage allele in GWAS must match the dosage allele in the LD reference panel.

# (GenT) Gene-based association test
This tests the null hypothesis that no SNPs in the the gene-specific set are associated with the disease trait. The image below shows an overview of the statistical and computational frameworks of the GenT method.

![](https://github.com/noahlorinczcomi/gent/blob/main/gent_overview.png)

We now show how to load our example data for the *SYK* gene and perform a gene-based association test with GenT. Z-statistics are from the Alzheimer's disease (AD) GWAS by Bellenguez et al. (2022) and the LD matrix is estimated using 1000 Genomes Phase 3 European samples.
```
library(gent)
data=readRDS('example_GenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data/example_GenT_data.Rds
z=data$z # (vector) AD Z-statistics for 805 SNPs corresponding to SYK gene
LD=data$LD # (matrix) LD matrix for the 805 SNPs (allele-harmonized)
results=gent(z,LD)

print(results)
$pval
[1] 1.96426e-06

$shape
[1] 12.70002

$rate
[1] 0.01577642

$mu_h0
[1] 805

$sigma2_h0
[1] 51025.51

$mu_h1
[1] 3106.22

$sigma2_h1
[1] 316919.3
```

# (xGenT) Gene-based association test integrating xQTLs
This tests the null hypothesis that no SNPs in the gene-specific set are associated with the disease trait using an xQTL-weighted statistic. The image below shows an overview of the statistical and computational frameworks of the xGenT method.

![](https://github.com/noahlorinczcomi/gent/blob/main/xgent_overview.png)

We now show how to load our example data for the *RIPK2* gene and perform a gene-based association test integrating brain eQTLs with xGenT. Disease Z-statistics are from the Alzheimer's disease GWAS by Bellenguez et al. (2022); eQTL Z-statistics are from GTEx v8 (https://gtexportal.org/home/) in cerebellum, spinal cord, frontal cortex, cortex, and hippocampal tissues; the LD matrix is estimated using 1000 Genomes Phase 3 European samples.
```
library(gent)
data=readRDS('example_xGenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
ad_z=data$ad_z # (vector) Z-statistics for 58 SNPs from the AD GWAS for RIPK2 gene
eqtl_z=data$eqtl_z # (matrix) Z-statistics for same 58 SNPs (allele-harmonized) for gene expression association
LD=data$LD # (matrix) LD matrix for 58 SNPs (allele-harmonized)
results=gent(ad_z,LD,xqtl_Z=eqtl_z) # adding `xqtl_Z` weighting matrix makes it xGenT

print(results)
$pval
[1] 5.126952e-05

$shape
[1] 1.906503

$rate
[1] 0.03287074

$mu_h0
[1] 58

$sigma2_h0
[1] 1764.487

$mu_h1
[1] 430.9083

$sigma2_h1
[1] 11732.78
```

# (MuGenT) Multi-ancestry gene-based association test
This tests the null hypothesis that no SNPs in the gene-specific set for any population are associated with the disease trait. The image below shows an overview of the statistical and computational frameworks of the MuGenT method and also shows how to perform the MuGenT-PH (MuGenT population heterogeneity) test (see the **MuGenT-PH** subsection below).

![](https://github.com/noahlorinczcomi/gent/blob/main/mugent_overview.png)

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of association with type 2 diabetes (T2D) integrating GWAS data from African American (AFA), East and South Asian (EAS, SAS), European (EUR), and Hispanic (HIS) populations. T2D Z-statistics are from GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_MuGenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS (allele-harmonized) for PPP3CA gene
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs (allele-harmonized)
results=mugent(t2d_z,LD_list)

print(results)
$pval
[1] 4.606173e-47

$shape
[1] 8.437202

$rate
[1] 1.68744

$mu_h0
[1] 5

$sigma2_h0
[1] 2.963067
```

# (MuGenT-PH) Multi-ancestry gene-based test of association heterogeneity
This tests the null hypothesis that the gene is not associated with the disease trait at varying magnitudes across populations.

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of association heterogeneity with type 2 diabetes (T2D) across multiple populations. T2D  Z-statistics are from the GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_mugenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS (allele-harmonized) for PPP3CA gene
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs (allele-harmonized)
results=mugent_ph(t2d_z,LD_list)

print(results)
$pval
[1] 0.007877747

$shape
[1] 20.19291

$rate
[1] 0.02692388

$mu_h0
[1] 750

$sigma2_h0
[1] 27856.31
```

# (MuGenT-Pleio) Multi-ancestry gene-based test of pleiotropic association with all populations
This tests the null hypothesis that the gene is not associated with the disease trait in all populations.

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of pleiotropy across multiple  populations for type 2 diabetes (T2D). T2D  Z-statistics are from the GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_mugenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS (allele-harmonized) for PPP3CA gene
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs (allele-harmonized)
results=mugent_pleio(t2d_z,LD_list)

print(results)
$result
[1] "pleiotropy"

$adjusted_significance_quantile
[1] 1.478305
```
This test simply compares a test statistic to a critical value and so returns only whether the null hypothesis was rejected (i.e., `$result == 'pleiotropy'`) or not rejected (i.e., `$result == 'no pleiotropy'`).

## References
1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526(7571), 68.

Bellenguez, C., Küçükali, F., Jansen, I. E., Kleineidam, L., Moreno-Grau, S., Amin, N., ... & Goldhardt, O. (2022). New insights into the genetic etiology of Alzheimer’s disease and related dementias. Nature genetics, 54(4), 412-436.

Suzuki, K., Hatzikotoulas, K., Southam, L., Taylor, H. J., Yin, X., Lorenz, K. M., ... & Kamanu, F. K. (2024). Genetic drivers of heterogeneity in type 2 diabetes pathophysiology. Nature, 627(8003), 347-357.





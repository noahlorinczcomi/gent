# Installation
```
devtools::install_github('noahlorinczcomi/gent')
# or
remotes::install_github('noahlorinczcomi/gent')
```
# Preliminaries
All gene-based association test methods [```gent()```, ```mugent()```, ```mugent_ph()```, ```mugent_pleio()```,```mugent_sel()```] require estimated LD matrices for a set of SNPs. The effect allele in your GWAS must be the effective allele allele in the reference file (e.g., the ```a1``` allele in PLINK-formatted .bim files). If it is not for a particular SNP, simply multiple the GWAS beta or Z-statistic by -1 (negative one).

# (GenT) Gene-based association test
$\color{blue}{\textsf{This tests the alternative hypothesis that the gene is associated with the disease trait.}}$

We show how to load our example data for the *SYK* gene and perform a gene-based association test with GenT. Z-statistics are from the Alzheimer's disease (AD) GWAS by Bellenguez et al. (2022) and the LD matrix is estimated using 1000 Genomes Phase 3 European samples.
```
library(gent)
data=readRDS('example_GenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
z=data$z # (vector) AD Z-statistics for 805 SNPs corresponding to SYK gene
LD=data$LD # (matrix) LD matrix for the 805 SNPs
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
$\color{blue}{\textsf{This tests the alternative hypothesis that the gene is associated with the disease trait and genetically}}$
$\color{blue}{\textsf{correlated with local xQTLs, which is equivalent to Mendelian Randomization if its assumptions hold.}}$

We show how to load our example data for the *RIPK2* gene and perform a gene-based association test integrating brain eQTLs with xGenT. Disease Z-statistics are from the Alzheimer's disease GWAS by Bellenguez et al. (2022); eQTL Z-statistics are from GTEx v8 (https://gtexportal.org/home/) in cerebellum, spinal cord, frontal cortex, cortex, and hippocampal tissues; the LD matrix is estimated using 1000 Genomes Phase 3 European samples.
```
library(gent)
data=readRDS('example_xGenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
ad_z=data$ad_z # (vector) Z-statistics for 58 SNPs from the AD GWAS
eqtl_z=data$eqtl_z # (matrix) Z-statistics for same 58 SNPs for gene expression association
LD=data$LD # (matrix) LD matrix for 58 SNPs
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
$\color{blue}{\textsf{This tests the alternative hypothesis that the gene is associated with the disease trait in any population.}}$

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of association with type 2 diabetes (T2D) integrating GWAS data from African American (AFA), East and South Asian (EAS, SAS), European (EUR), and Hispanic (HIS) populations. T2D Z-statistics are from GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_mugenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs
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
$\color{blue}{\textsf{This tests the alternative hypothesis that the gene is associated with the disease trait at varying magnitudes}}$
$\color{blue}{\textsf{across populations.}}$

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of association heterogeneity with type 2 diabetes (T2D) across multiple populations. T2D  Z-statistics are from the GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_mugenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs
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
$\color{blue}{\textsf{This tests the alternative hypothesis that the gene is associated with the disease trait in all populations.}}$

We show how to load our example data for the *PPP3CA* gene and perform a gene-based test of pleiotropy across multiple  populations for type 2 diabetes (T2D). T2D  Z-statistics are from the GWAS by Suzuki et al. (2024) and the LD matrices are population-specific and estimated using 1000 Genomes Phase 3 samples.
```
library(gent)
data=readRDS('example_mugenT_data.Rds') # https://github.com/noahlorinczcomi/gent/tree/main/example_data 
t2d_z=data$t2d_z # (matrix) Z-statistics for 150 SNPs (rows) from the EUR, AFA, SAS, EAS, and HIS (columns) GWAS
LD_list=data$LD_list # list of population-specific LD matrices for 150 SNPs
results=mugent_pleio(t2d_z,LD_list)

print(results)
$result
[1] "pleiotropy"

$adjusted_significance_quantile
[1] 1.478305
```

## References
1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526(7571), 68.

Bellenguez, C., Küçükali, F., Jansen, I. E., Kleineidam, L., Moreno-Grau, S., Amin, N., ... & Goldhardt, O. (2022). New insights into the genetic etiology of Alzheimer’s disease and related dementias. Nature genetics, 54(4), 412-436.

Suzuki, K., Hatzikotoulas, K., Southam, L., Taylor, H. J., Yin, X., Lorenz, K. M., ... & Kamanu, F. K. (2024). Genetic drivers of heterogeneity in type 2 diabetes pathophysiology. Nature, 627(8003), 347-357.





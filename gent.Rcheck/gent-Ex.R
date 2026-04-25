pkgname <- "gent"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('gent')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AFRGenTStatLD")
### * AFRGenTStatLD

flush(stderr()); flush(stdout())

### Name: AFRGenTStatLD
### Title: Gene-based test statistic correlations (African 1000 Genomes
###   Phase 3)
### Aliases: AFRGenTStatLD
### Keywords: datasets

### ** Examples

data(AFRGenTStatLD)
str(AFRGenTStatLD)



cleanEx()
nameEx("AMRGenTStatLD")
### * AMRGenTStatLD

flush(stderr()); flush(stdout())

### Name: AMRGenTStatLD
### Title: Gene-based test statistic correlations (Admixed American 1000
###   Genomes Phase 3)
### Aliases: AMRGenTStatLD
### Keywords: datasets

### ** Examples

data(AMRGenTStatLD)
str(AMRGenTStatLD)



cleanEx()
nameEx("EASGenTStatLD")
### * EASGenTStatLD

flush(stderr()); flush(stdout())

### Name: EASGenTStatLD
### Title: Gene-based test statistic correlations (East Asian 1000 Genomes
###   Phase 3)
### Aliases: EASGenTStatLD
### Keywords: datasets

### ** Examples

data(EASGenTStatLD)
str(EASGenTStatLD)



cleanEx()
nameEx("EURGenTStatLD")
### * EURGenTStatLD

flush(stderr()); flush(stdout())

### Name: EURGenTStatLD
### Title: Gene-based test statistic correlations (European 1000 Genomes
###   Phase 3)
### Aliases: EURGenTStatLD
### Keywords: datasets

### ** Examples

data(EURGenTStatLD)
str(EURGenTStatLD)



cleanEx()
nameEx("EnsemblHg19GenePos")
### * EnsemblHg19GenePos

flush(stderr()); flush(stdout())

### Name: EnsemblHg19GenePos
### Title: Ensembl hg19 coordinates of gene start and end base pair
###   positions
### Aliases: EnsemblHg19GenePos
### Keywords: datasets

### ** Examples

data(EnsemblHg19GenePos)
str(EnsemblHg19GenePos)



cleanEx()
nameEx("EnsemblHg38GenePos")
### * EnsemblHg38GenePos

flush(stderr()); flush(stdout())

### Name: EnsemblHg38GenePos
### Title: Ensembl hg38 coordinates of gene start and end base pair
###   positions
### Aliases: EnsemblHg38GenePos
### Keywords: datasets

### ** Examples

data(EnsemblHg38GenePos)
str(EnsemblHg38GenePos)



cleanEx()
nameEx("SASGenTStatLD")
### * SASGenTStatLD

flush(stderr()); flush(stdout())

### Name: SASGenTStatLD
### Title: Gene-based test statistic correlations (South Asian 1000 Genomes
###   Phase 3)
### Aliases: SASGenTStatLD
### Keywords: datasets

### ** Examples

data(SASGenTStatLD)
str(SASGenTStatLD)



cleanEx()
nameEx("gene_clump")
### * gene_clump

flush(stderr()); flush(stdout())

### Name: gene_clump
### Title: Clumping of gene-based test statistics
### Aliases: gene_clump

### ** Examples

gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
gene_clump(gent_results,'EUR')



cleanEx()
nameEx("gent")
### * gent

flush(stderr()); flush(stdout())

### Name: gent
### Title: Gene-based association test (GenT)
### Aliases: gent

### ** Examples

# Example for 5 SNPs
LD=cov2cor(rWishart(1,100,diag(5))[,,1])
z=c(mvnfast::rmvn(1,rep(0,5),diag(5)))
gent(z,LD)



cleanEx()
nameEx("gent_finemap")
### * gent_finemap

flush(stderr()); flush(stdout())

### Name: gent_finemap
### Title: Fine-mapping gene-based association test statistics
### Aliases: gent_finemap

### ** Examples

# Example for Alzheimer's disease GWAS data
gent_finemap(
  gent_results=gent_results,
  ld_population='EUR',
  gwas_n=50000,
  index_genes=c('CR1','MYCL'),
  window_kb_width=2000,
  verbose=TRUE)



cleanEx()
nameEx("gent_genomewide")
### * gent_genomewide

flush(stderr()); flush(stdout())

### Name: gent_genomewide
### Title: Genome-wide gene-based association test (GenT)
### Aliases: gent_genomewide

### ** Examples

# Example for Alzheimer's disease GWAS data
gent_genomewide(
  gwas=ad_gwas,
  ld_population='EUR',
  ld_directory='ld_matrices',
  build='grch37',
  KbWindow=50,
  snp='MarkerName',
  chromosome='Chromosome',
  position='Position',
  effect_allele='Effect_allele',
  z_statistic='z',
  verbose=TRUE,
  return_snp_gene_pairs=FALSE)



cleanEx()
nameEx("gent_manhattan")
### * gent_manhattan

flush(stderr()); flush(stdout())

### Name: gent_manhattan
### Title: Gene-based test Manhattan plot
### Aliases: gent_manhattan

### ** Examples

gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
gent_manhattan(gent_results)



cleanEx()
nameEx("gent_qq")
### * gent_qq

flush(stderr()); flush(stdout())

### Name: gent_qq
### Title: Gene-based test QQ plot
### Aliases: gent_qq

### ** Examples

gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
gent_qq(gent_results)



cleanEx()
nameEx("ld")
### * ld

flush(stderr()); flush(stdout())

### Name: ld
### Title: Return LD matrix for gene
### Aliases: ld

### ** Examples

# Example for 10 SNPs
rs=paste0('rs',1:10)
ea=sample(c('A','C','T','G'),10,replace=T)
ld_ref='~/1kg.v3/EUR'
wd=getwd()
ld(rs,ea,ld_ref,wd,plink_exec='plink',verbose=TRUE)



cleanEx()
nameEx("mugent")
### * mugent

flush(stderr()); flush(stdout())

### Name: mugent
### Title: Multi-ancestry gene-based association test (MuGenT)
### Aliases: mugent

### ** Examples

# Example for 5 SNPs and 2 populations
ldlist=rWishart(2,100,diag(5))
ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
mugent(Z,ldlist)



cleanEx()
nameEx("mugent_genomewide")
### * mugent_genomewide

flush(stderr()); flush(stdout())

### Name: mugent_genomewide
### Title: Genome-wide multi-trait/ancestry gene-based association test
###   (MuGenT)
### Aliases: mugent_genomewide

### ** Examples

# Example for COVID-19 GWAS data in EUR and AFR populations
mugent_genomewide(
  gwas_list = list(EUR=eur_gwas, AFR=afr_gwas),
  ld_population_list = list(EUR='EUR', AFR='AFR'),
  ld_directory = 'ld_matrices',
  build 'grch37',
  KbWindow = 50,
  snp_list = list(EUR='rsid', AFR='rsid'),
  chromosome_list = list(EUR='#CHR', AFR='#CHR'),
  position_list = list(EUR='POS', AFR='POS'),
  effect_allele_list = list(EUR='ALT', AFR='ALT'),
  z_statistic_list = list(EUR='z', AFR='z'),
  verbose = TRUE.
  return_snp_gene_pairs = FALSE)



cleanEx()
nameEx("mugent_ph")
### * mugent_ph

flush(stderr()); flush(stdout())

### Name: mugent_ph
### Title: Multi-ancestry gene-based association heterogeneity test
###   (MuGenT-PH)
### Aliases: mugent_ph

### ** Examples

# Example for 5 SNPs and 2 populations
ldlist=rWishart(2,100,diag(5))
ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
mugent_ph(Z,ldlist)



cleanEx()
nameEx("mugent_pleio")
### * mugent_pleio

flush(stderr()); flush(stdout())

### Name: mugent_pleio
### Title: Multi-ancestry joint gene-based associated test (MuGenT-Pleio)
### Aliases: mugent_pleio

### ** Examples

# Example for 5 SNPs and 2 populations
ldlist=rWishart(2,100,diag(5))
ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
mugent_pleio(Z,ldlist)



cleanEx()
nameEx("mugent_sel")
### * mugent_sel

flush(stderr()); flush(stdout())

### Name: mugent_sel
### Title: Pairwise test of only one population with a gene association
###   (MuGenT-Sel)
### Aliases: mugent_sel

### ** Examples

# Example for 5 SNPs and 2 populations
ldlist=rWishart(2,100,diag(5))
ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
mugent_sel(Z,ldlist)



cleanEx()
nameEx("multipop_anova")
### * multipop_anova

flush(stderr()); flush(stdout())

### Name: multipop_anova
### Title: SNP-based ANOVA test for population heterogeneity
### Aliases: multipop_anova

### ** Examples

# Example for 5 SNPs and 2 populations
z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
b=z/sqrt(30000)
multipop_anova(b,b/z)



cleanEx()
nameEx("wgent_genomewide")
### * wgent_genomewide

flush(stderr()); flush(stdout())

### Name: wgent_genomewide
### Title: Genome-wide inverse MAF-weighted gene-based association test
### Aliases: wgent_genomewide

### ** Examples

# Example for Alzheimer's disease GWAS data
gent_genomewide(
  gwas=ad_gwas,
  snp_weights='sqrt_inverse_mafs',
  ld_population='EUR',
  ld_directory='ld_matrices',
  build='grch37',
  KbWindow=50,
  snp='MarkerName',
  chromosome='Chromosome',
  position='Position',
  effect_allele='Effect_allele',
  effect_size='beta',
  standard_error='se',
  verbose=TRUE,
  return_snp_gene_pairs=FALSE)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

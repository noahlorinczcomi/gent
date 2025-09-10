#' Gene-based association test (GenT)
#'
#' This function performs a single gene-based association test.
#' @param zs vector of SNP-specific Z-statistics from GWAS whose indices correspond to those of \code{LD}
#' @param LD LD matrix whose row/column indices correspond to the indices of \code{zs}
#' @param A A weight matrix which forms the quadratic form \code{zs}'(A^2)\code{zs}. Equivalent to setting zs<-A%*%zs. Weights SNPs in gene-based testing.
#' @param xqtl_Z An m x p matrix of xQTL effect sizes (e.g., Z-statistics). Rows correspond to SNPs; columns correspond to tissues or xQTL types.
#' @param chisquares SNP chi-square statistics. If you supply a value for this and leave \code{A} and \code{xqtl_Z} empty, it is equivalent to supplying \code{zs} and \code{LD}.
#' @return A list with these components:
#' \itemize{
#'  \item \code{pval}: P-value for testing H0: gene is not associated with trait.
#'  \item \code{shape}: shape parameter of null (Gamma) distribution.
#'  \item \code{rate}: rate parameter of null (Gamma) distribution.
#'  \item \code{mu_h0}: expectation of null (Gamma) distribution.
#'  \item \code{sigma2_h0}: variance of null (Gamma) distribution.
#' }
#'
#' @export
#' @examples
#' # Example for 5 SNPs
#' LD=cov2cor(rWishart(1,100,diag(5))[,,1])
#' z=c(mvnfast::rmvn(1,rep(0,5),diag(5)))
#' gent(z,LD)
gent=function (zs = NULL, LD, A = NULL, xqtl_Z = NULL, chisquares = NULL) {
  if (is.null(A) & is.null(xqtl_Z)) {
    mu = tr(LD)
    trASAS = tr(LD %*% LD)
    y = sum(zs^2)
  }
  else if (!is.null(A) & is.null(xqtl_Z)) {
    mu = tr(A %*% LD)
    trASAS = tr(A %*% LD %*% A %*% LD)
    y = c(t(zs) %*% A %*% zs)
  }
  else if (is.null(A) & !is.null(xqtl_Z)) {
    xqtl_Z = as.matrix(xqtl_Z)
    m = nrow(xqtl_Z)
    p = ncol(xqtl_Z)
    L = diag(c(rowSums(xqtl_Z^2)))/sqrt(m * p)
    mu = tr(L %*% LD)
    trASAS = tr(L %*% LD %*% L %*% LD)
    y = c(t(zs) %*% L %*% zs)
  }
  if (!is.null(chisquares)) y = sum(chisquares)
  sigma2 = 2 * trASAS
  beta = mu/sigma2
  alpha = beta * mu
  pval = pgamma(y, shape = alpha, rate = beta, lower.tail = FALSE)
  out=list(pval = pval, shape = alpha, rate = beta, mu_h0 = mu, sigma2_h0 = sigma2)
  out=lapply(out,c)
  return(out)
}

#' Multi-ancestry gene-based association test (MuGenT)
#'
#' This function performs a gene-based association test using multiple populations.
#' @param Z matrix of Z-statistics from GWAS. Rows are SNPs and columns are populations.
#' @param ldlist list of population-specific LD matrices whose ordering corresponds to the column ordering of \code{Z}
#' @return
#' \itemize {
#' \item \code{pval}: P-value for testing H0: gene is not associated with trait.
#' \item \code{shape}: shape parameter of null (Gamma) distribution.
#' \item \code{rate}: rate parameter of null (Gamma) distribution.
#' \item \code{mu_h0}: expectation of null (Gamma) distribution.
#' \item \code{sigma2_h0}: variance of null (Gamma) distribution.
#' }
#' @export
#' @examples
#' # Example for 5 SNPs and 2 populations
#' ldlist=rWishart(2,100,diag(5))
#' ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
#' Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
#' mugent(Z,ldlist)
mugent=function(Z,ldlist) {
  if(is.matrix(ldlist)) ldlist=lapply(1:ncol(Z), function(h) ldlist)
  Z=as.matrix(Z)
  p=ncol(Z);m=nrow(Z);j=rep(1,p)
  stat=c(t(j)%*%(t(Z)%*%Z/m)%*%j)
  EZ=t(j)%*%diag(ncol(Z))%*%j # under H0
  Kj=kronecker(t(j),t(j))
  VZ=tr(Kj%*%varmat(p,ldlist)%*%t(Kj))/m^2
  beta=EZ/VZ
  alpha=EZ*beta
  p=pgamma(stat,shape=alpha,rate=beta,lower.tail=FALSE)
  out=list(pval=p,shape=alpha,rate=beta,mu_h0=EZ,sigma2_h0=VZ)
  lapply(out,c)
}

#' SNP-based ANOVA test for population heterogeneity
#'
#' This function performs a single ANOVA test for SNP across populations/phenotypes.
#' @param effect_size_matrix Matrix of GWAS SNP effect sizes. Rows are SNPs. Columns are populations/phenotypes.
#' @param standard_error_matrix Matrix of GWAS standard errors corresponding to \code{effect_size_matrix}.
#' @return matrix whose columns respectively are the ANOVA P-value and test statistics
#' @export
#' @examples
#' # Example for 5 SNPs and 2 populations
#' z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
#' b=z/sqrt(30000)
#' multipop_anova(b,b/z)
multipop_anova=function(effect_size_matrix,standard_error_matrix) {
  # assuming population independence (under H0) and effect_size_matrix and standard_error_matrix are gene-specific matrices
  BS=as.matrix(effect_size_matrix);VS=as.matrix(standard_error_matrix)^2
  if(any(dim(BS)!=dim(VS))) stop('dimensions of `effect_size_matrix` and `standard_error_matrix` must be the same')
  p=stat=c()
  for(i in 1:nrow(BS)) {
    bs=BS[i,]
    vs=VS[i,]
    w=1/vs # inverse-variance weights to estimate grand mean
    mu=sum(w*bs)/sum(w) # estimate of grand mean
    varmu=(1/sum(w))^2*sum(w^2*vs) # estimate of variance of grand mean
    taus=bs-mu # esimates of taus
    vartaus=varmu+vs-2/sum(w) # estimate of variance of taus
    # perform joint test
    stat[i]=t(taus)%*%diag(1/vartaus)%*%taus
    p[i]=c(pchisq(stat[i],length(taus),lower.tail=FALSE)) # p-value to test H0: all true betas are equal
  }
  out=cbind(p,stat)
  colnames(out)=c('ANOVA_pval','ANOVA_chisq')
  return(out)
}

#' Multi-ancestry gene-based association heterogeneity test (MuGenT-PH)
#'
#' This function performs a gene-based test of association heterogeneity across multiple populations/phenotypes.
#' @param Z matrix of Z-statistics from GWAS. Rows are SNPs and columns are populations.
#' @param ldlist list of population-specific LD matrices whose ordering corresponds to the column ordering of \code{Z}
#' @return
#' \itemize {
#' \item \code{pval}: P-value for testing H0: gene is not associated with trait.
#' \item \code{shape}: shape parameter of null (Gamma) distribution.
#' \item \code{rate}: rate parameter of null (Gamma) distribution.
#' \item \code{mu_h0}: expectation of null (Gamma) distribution.
#' \item \code{sigma2_h0}: variance of null (Gamma) distribution.
#' }
#' @export
#' @examples
#' # Example for 5 SNPs and 2 populations
#' ldlist=rWishart(2,100,diag(5))
#' ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
#' Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
#' mugent_ph(Z,ldlist)
mugent_ph=function(Z,ldlist) {
  # calculate SNP-level ANOVA P-values
  Z=as.matrix(Z);p=ncol(Z)
  chisquares=multipop_anova(Z,matrix(1,nrow(Z),p))[,2]
  # does NOT all LD matrices are the same across populations
  chisquares=c(chisquares)
  m=length(chisquares)
  R1=matrix(0,nr=m,nc=m)
  for(ll in 1:p) R1=R1+ldlist[[ll]]^2; R1=R1/p
  D=diag(sqrt(2*p),m)
  j=rep(1,m)
  mu=m*p
  variance=c(t(j)%*%D%*%R1%*%D%*%j) # actually same as sum(D%*%R1%*%D)
  beta=mu/variance
  alpha=beta*mu
  p=pgamma(sum(chisquares),shape=alpha,rate=beta,lower.tail=FALSE)
  out=list(pval=p,shape=alpha,rate=beta,mu_h0=mu,sigma2_h0=variance)
  lapply(out,c)
}

#' Multi-ancestry joint gene-based associated test (MuGenT-Pleio)
#'
#' This function performs a test of H1 that the gene is associated with in all populations/for all phenotypes.
#' @param Z matrix of Z-statistics from GWAS. Rows are SNPs and columns are populations.
#' @param ldlist list of population-specific LD matrices whose ordering corresponds to the column ordering of \code{Z}
#' @param alpha the nominal (uncorrected) significance threshold for testing a single gene. If testing all genes genome-wide, this should be less than the threshold for testing a single gene.
#' @return
#' \itemize {
#' \item \code{result}: An indication if H0 (the gene is not association in all populations/with all traits) can be rejected. 'Pleiotropy' if it is rejected and 'no pleiotropy' otherwise.
#' \item \code{adjusted_significance_quantile}: The chi-square quantile of the adjusted nominal significance threshold for inferring H1 given the number of SNPs used and populations/traits tested.
#' }
#' @export
#' @examples
#' # Example for 5 SNPs and 2 populations
#' ldlist=rWishart(2,100,diag(5))
#' ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
#' Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
#' mugent_pleio(Z,ldlist)
mugent_pleio=function(Z,ldlist,alpha=0.05) {
  p=ncol(Z);meff=f_meff(ldlist)
  q0=qchisq(1-(1-(1-alpha)^(1/(2*meff)))^(1/p),1)
  boo=apply(Z^2,1,function(h) all(h>q0))
  out=list(result=ifelse(any(boo),'pleiotropy','no pleiotropy'),adjusted_significance_quantile=c(q0))
  out
}

#' Pairwise test of only one population with a gene association (MuGenT-Sel)
#'
#' This function returns the posterior probability that a gene is only associated with one population/trait in a pair, pairwise for all available populations/traits.
#' @param Z matrix of Z-statistics from GWAS. Rows are SNPs and columns are populations.
#' @param ldlist list of population-specific LD matrices whose ordering corresponds to the column ordering of \code{Z}
#' @return This function returns a matrix of pairwise comparisons between populations/traits. Non-missing values in this matrix represent posterior probabilities that the gene is only associated with one of the two traits/populations.
#' @export
#' @examples
#' # Example for 5 SNPs and 2 populations
#' ldlist=rWishart(2,100,diag(5))
#' ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
#' Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
#' mugent_sel(Z,ldlist)
mugent_sel=function(Z,ldlist,verbose=T) {
  p=ncol(Z)
  if(!is.list(ldlist)) lapply(1:p,function(h) ldlist)
  P=matrix(nr=p,nc=p)
  colnames(P)=rownames(P)=colnames(Z)
  for(i in 1:p) {
    for(j in 1:p) {
      if(i<=j) next
      ldi=list(ldlist[[i]],ldlist[[j]])
      P[i,j]=mixp(Z[,c(i,j)],ldi)
    }
  }
  colnames(P)=rownames(P)=paste0('Pop./trait ',1:p)
  if(verbose) cat('P(only one associated population/trait) \n(pairwise)\n\n')
  P
}

#' Return LD matrix for gene
#'
#' This function can externally use PLINK to estimate the matrix of LD correlations between a set of SNPs.
#' @param rsids vector of rsIDs for the SNPs you want to estimate LD between.
#' @param effect_alleles effect alleles of the rsIDs from GWAS.
#' @param ld_ref_file filename (without extension) of LD reference panel in PLINK (.bim/.bed/.fam) format.
#' @param writeable_directory a directory in which you have read-write access. Files will be created in this directory then deleted from it.
#' @param plink_exec command to execute PLINK from the command line, typically just 'plink'.
#' @param verbose should a statement about your effect alleles and those in the LD reference be printed to the console?
#' @return
#' \itemize {
#' \item \code{ld}: LD matrix. Rows and columns are in the original order provided to this function, filtered to only those present in the LD reference.
#' \item \code{signmat}: This is a vector of +1s and -1s which indicate if your effect alleles did match (+1) or did not match (-1) the effect alleles in the LD reference panel. Effect alleles must match for later inference. If they don't simply multiply your Z-statistics, for SNPs present in the LD reference panel (see row/column names of \code{ld}), by \code{signment} element-wise.
#' }
#' @export
#' @examples
#' # Example for 10 SNPs
#' rs=paste0('rs',1:10)
#' ea=sample(c('A','C','T','G'),10,replace=T)
#' ld_ref='~/1kg.v3/EUR'
#' wd=getwd()
#' ld(rs,ea,ld_ref,wd,plink_exec='plink',verbose=TRUE)
ld=function(rsids,effect_alleles,ld_ref_file=NULL,writeable_directory=getwd(),plink_exec='plink',chromosome=NULL,need_bim=TRUE,verbose=T,temp_file_prefix='') {
  library(data.table);library(dplyr)
  if(is.null(ld_ref_file)) stop('you must specify an LD reference in PLINK format without file extension')
  if(need_bim) {
    bf=paste0(ld_ref_file,'.bim')
    if(!file.exists(bf)) stop(paste0(bf,' does not exist but it must'))
    bim=fread(bf) %>% as.data.frame();
    colnames(bim)=c('chr','rsid','x','position','a1','a2')
    if(!is.null(chromosome)) bim=bim %>% filter(chr %in% chromosome)
    ri=rsids[rsids %in% bim$rsid]
    ea=effect_alleles[rsids %in% bim$rsid]
    bim=bim %>% left_join(data.frame(rsid=ri,effect_allele=ea)) %>% arrange(rsid)
    ix=sapply(ri,function(h) which(bim$rsid==h))
    bim=bim %>% arrange(position)
    # directly calculate LD
    setwd(writeable_directory)
    rd=temp_file_prefix
    write.table(bim$rsid,paste0(rd,'.txt'),row.names=F,quote=F,col.names=F)
    cmd=paste0(plink_exec,' --bfile ',ld_ref_file,' --r square --extract ',rd,'.txt --out ',rd)
    oo=system(cmd,intern=T)
    ld=fread=paste0(rd,'.ld') %>% as.matrix()
    ld=ld[ix,ix]
    bim=bim[ix,]
    rownames(ld)=colnames(ld)=bim$rsid
    signmat=ifelse(bim$ea==bim$a1,1,-1)
    oo=system(paste0('rm ',rd,'*'),intern=T)
    out=list(ld=ld,signmat=signmat)
    if(verbose) cat('`signmat` is a vector of constants to multiple your Z-statistics by since their effect alleles did not match the A1 alleles in the .bim file')
    out
  } else {
    # directly calculate LD
    setwd(writeable_directory)
    rd=temp_file_prefix
    write.table(rsids,paste0(rd,'.txt'),row.names=F,quote=F,col.names=F)
    cmd=paste0(plink_exec,' --bfile ',ld_ref_file,' --r square --extract ',rd,'.txt --out ',rd)
    oo=system(cmd,intern=T)
    ld=fread(paste0(rd,'.ld')) %>% as.matrix()
    oo=system(paste0('rm ',rd,'*'),intern=T)
    out=list(ld=ld,signmat=rep(1,nrow(ld)))
    colnames(out$ld)=rownames(out$ld)=rsids
  }
  return(out)
}

#' Genome-wide gene-based association test (GenT)
#'
#' This function performs gene-based association testing genome-wide.
#' @param gwas GWAS summary statistics
#' @param KbWindow Kilobase window used to assign SNPs to genes
#' @param ld_population LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR'
#' @param ld_directory Relative or absolute filepath to LD directory (see Github Wiki)
#' @param snp Column name of SNP rsID in GWAS data
#' @param chromosome Column name of SNP chromosome in GWAS data
#' @param position Column name of SNP hg19 base pair position in GWAS data
#' @param effect_allele Column name of SNP effect allele in GWAS data
#' @param z_statistic Column name of SNP Z-statistic in GWAS data
#' @param index Gene index file. see \code{data(EnsemblHg19GenePos)} for the expected format.
#' @param verbose TRUE if progress should be printed to the console, FALSE otherwise
#' @param return_snp_gene_pairs If TRUE, will return a chromosome-specific list of tested gene-specific SNP sets annotated by gene
#' @return A dataframe with these components:
#' \itemize{
#'  \item `pval`: P-value testing the gene-based null hypothesis
#'  \item `shape`: Shape parameter of the Gamma null distribution
#'  \item `rate`: Rate parameter of the Gamma null distribution
#'  \item `mu_h0`: Expected value of the gene-based test statistic under the null hypothesis
#'  \item `sigma2_h0`: Variance of the gene-based test statistic under the null hypothesis
#'  \item `gene`: Gene symbol
#'  \item `m`: Number of SNPs tested in gene-specific set
#'  \item `chr`: Chromosome of gene
#'  \item `gene_start`: Start position of gene (Ensembl hg19)
#'  \item `window_start`: Start position of window in which SNPs were captured for this gene (hg19)
#'  \item `gene_end`: End position of gene (Ensembl hg19)
#'  \item `window_end`: End position of window in which SNPs were captured for this gene (hg19)
#' }
#'
#' @export
#' @examples
#' # Example for Alzheimer's disease GWAS data
#' gent_genomewide(
#'   gwas=ad_gwas,
#'   ld_population='EUR',
#'   ld_directory='ld_matrices',
#'   KbWindow=50,
#'   snp='MarkerName',
#'   chromosome='Chromosome',
#'   position='Position',
#'   effect_allele='Effect_allele',
#'   z_statistic='z',
#'   verbose=TRUE,
#'   return_snp_gene_pairs=FALSE)
gent_genomewide=function(gwas,
                         KbWindow=50,
                         ld_population='EUR',
                         ld_directory='ld_directory',
                         snp='rsid',
                         chromosome='chr',
                         position='position',
                         effect_allele='effect_allele',
                         z_statistic='z',
                         index=NULL,
                         verbose=TRUE,
                         return_snp_gene_pairs=FALSE) {
  if(is.null(index)) {data(EnsemblHg19GenePos);index=EnsemblHg19GenePos}
  setwd(ld_directory)
  gwas=gwas %>% na.omit()
  gwas=gwas %>% rename(rsid=!!sym(snp), chr=!!sym(chromosome), position=!!sym(position), effect_allele=!!sym(effect_allele), z=!!sym(z_statistic))
  chrs=gwas %>% select(chr) %>% pull() %>% unique() %>% as.numeric() %>% na.omit() %>% sort()
  chrs=intersect(1:22,chrs)
  sgp=list()
  rdf=data.frame()
  for(cc in 1:length(chrs)) {
    setwd(ld_directory)
    setwd(ld_population)
    setwd(paste0('chr',cc))
    if(verbose) cat('Starting chromosome', chrs[cc], '\n')
    gwas_chr=gwas %>% filter(chr==chrs[cc])
    index_chr=index %>% filter(chr==cc)
    genes=unique(index_chr$symbol)
    #pb=txtProgressBar(min=0,max=length(genes),style=3)
    k=0
    chrlist=list()
    for(i in 1:length(genes)) {
      #if(verbose) setTxtProgressBar(pb, i)
      tryCatch(
        {
          # check up front that the data contains at least one SNP for this gene
          starti=index_chr %>% filter(symbol==genes[i]) %>% select(start) %>% head(.,1) %>% pull()
          endi=index_chr %>% filter(symbol==genes[i]) %>% select(end) %>% head(.,1) %>% pull()
          starti=starti-KbWindow*1e3
          endi=endi+KbWindow*1e3
          gwas_chri=gwas_chr %>% filter(position>starti,position<endi)
          if(nrow(gwas_chri)==0) stop()
          # then proceed to load LD
          genefp=paste0(genes[i],'.Rds')
          if(!file.exists(genefp)) stop()
          ld=readRDS(genefp)
          ldrn=rownames(ld)
          spp=sapply(ldrn,\(.) unlist(strsplit(.,'_')))
          ld_df=t(spp) %>% as.data.frame() %>% rename(rsid=V1,a1=V2)
          ldrn=ld_df$rsid
          usesnps=intersect(ldrn,gwas_chri$rsid)
          gwas_chri=gwas_chri %>% filter(rsid %in% usesnps) %>% distinct(rsid,.keep_all=TRUE) %>% arrange(position)
          z=gwas_chri %>% left_join(ld_df,by='rsid') %>% mutate(z=ifelse(effect_allele==a1,z,-z)) %>% pull(z)
          # make sure data contains SNPs for this gene
          if(length(z)==0) stop()
          rownames(ld)=colnames(ld)=ldrn
          ld=ld[gwas_chri$rsid,gwas_chri$rsid]
          result=gent(z,ld)
          toadd=as.data.frame(result) %>% mutate(gene=genes[i],m=nrow(ld),chr=chrs[cc],gene_start=starti+KbWindow*1e3,window_start=starti,gene_end=endi-KbWindow*1e3,window_end=endi)
          rdf=rbind(rdf,toadd)
          # record SNP-gene pairs
          k=k+1
          if(return_snp_gene_pairs) {
            chrlist[[k]]=paste(gwas_chri$rsid,collapse=',')
            names(chrlist)[k]=genes[i]
          }
        },
        error=function(x) NA
      )
    }
    #close(pb)
    if(return_snp_gene_pairs) {
      sgp[[cc]]=chrlist
      names(sgp)[[cc]]=paste0('chr',chrs[cc])
    }
  }
  if(return_snp_gene_pairs) {
    return(list(result=as_tibble(rdf),snp_sets=sgp))
  } else {
    return(as_tibble(rdf))
  }
}

#' Genome-wide multi-trait/ancestry gene-based association test (MuGenT)
#'
#' This function performs multi-trait/ancestry gene-based association testing genome-wide.
#' @param gwas List of population/trait-specific GWAS summary statistics
#' @param KbWindow Kilobase window used to assign SNPs to genes
#' @param ld_population_list List of LD populations (see Github Wiki)
#' @param ld_directory Relative or absolute filepath to LD directory (see Github Wiki)
#' @param snp_list List of column names of SNP rsIDs in GWAS datasets
#' @param chromosome_list List of column name of SNP chromosomes in GWAS datasets
#' @param position_list List of column name of SNP hg19 base pair positions in GWAS datasetes
#' @param effect_allele_list List of column name of SNP effect alleles in GWAS datasets
#' @param z_statistic_list List of column name of SNP Z-statistics in GWAS datasets
#' @param mugentpleio_alpha MuGenT-Pleio nominal Type I error rate (uncorrected significance threshold)
#' @param index Gene index file. see \code{data(EnsemblHg19GenePos)} for the expected format.
#' @param verbose TRUE if progress should be printed to the console, FALSE otherwise
#' @param return_snp_gene_pairs If TRUE, will return a chromosome-specific list of tested gene-specific SNP sets annotated by gene
#' @return A list of dataframes with these components:
#' \itemize{
#' \itemize{
#'  \code{$MuGenT}
#'  \item `test`: Type of test (MuGenT)
#'  \item `gene`: Gene symbol
#'  \item `chr`: Chromosome of gene
#'  \item `pval`: P-value testing the gene-based null hypothesis
#'  \item `shape`: Shape parameter of the Gamma null distribution
#'  \item `rate`: Rate parameter of the Gamma null distribution
#'  \item `mu_h0`: Expected value of the gene-based test statistic under the null hypothesis
#'  \item `sigma2_h0`: Variance of the gene-based test statistic under the null hypothesis
#' }
#' \itemize{
#' \code{$MuGenT-PH}
#'  \item `test`: Type of test (MuGenT-PH)
#'  \item `gene`: Gene symbol
#'  \item `chr`: Chromosome of gene
#'  \item `pval`: P-value testing the gene-based null hypothesis
#'  \item `shape`: Shape parameter of the Gamma null distribution
#'  \item `rate`: Rate parameter of the Gamma null distribution
#'  \item `mu_h0`: Expected value of the gene-based test statistic under the null hypothesis
#'  \item `sigma2_h0`: Variance of the gene-based test statistic under the null hypothesis
#' }
#' \itemize{
#' \code{$MuGenT-Pleiotropy}
#'  \item `test`: Type of test (MuGenT-PH)
#'  \item `gene`: Gene symbol
#'  \item `chr`: Chromosome of gene
#'  \item `result`: Test result: 'pleitoropy' for evidence of pleiotropy, 'no pleiotropy' otherwise
#'  \item `adjusted_significance_quantile`: Quantile of the null Gamma distribution corrected for multiple testing
#' }
#' }
#'
#' @export
#' @examples
#' # Example for COVID-19 GWAS data in EUR and AFR populations
#' mugent_genomewide(
#'   gwas_list = list(EUR=eur_gwas, AFR=afr_gwas),
#'   ld_population_list = list(EUR='EUR', AFR='AFR'),
#'   ld_directory = 'ld_matrices',
#'   KbWindow = 50,
#'   snp_list = list(EUR='rsid', AFR='rsid'),
#'   chromosome_list = list(EUR='#CHR', AFR='#CHR'),
#'   position_list = list(EUR='POS', AFR='POS'),
#'   effect_allele_list = list(EUR='ALT', AFR='ALT'),
#'   z_statistic_list = list(EUR='z', AFR='z'),
#'   verbose = TRUE.
#'   return_snp_gene_pairs = FALSE)
mugent_genomewide=function(
    gwas_list,
    ld_population_list,
    ld_directory='ld_directory',
    KbWindow=50,
    snp_list=lapply(1:length(gwas_list),'rsid'),
    chromosome_list=lapply(1:length(gwas_list),'chr'),
    position_list=lapply(1:length(gwas_list),'position'),
    effect_allele_list=lapply(1:length(gwas_list),'effect_allele'),
    z_statistic_list=lapply(1:length(gwas_list),'z'),
    mugentpleio_alpha=0.05/12727,
    index=NULL,
    verbose=TRUE,
    return_snp_gene_pairs=FALSE) {
  if(is.null(index)) {data(EnsemblHg19GenePos);index=EnsemblHg19GenePos}
  # clean each GWAS
  k=length(gwas_list)
  used_snps=c()
  for(i in 1:k) {
    gwas=gwas_list[[i]]
    gwas=gwas %>% tidyr::drop_na()
    gwas=gwas %>% select(
      rsid=!!sym(snp_list[[i]]),
      chr=!!sym(chromosome_list[[i]]),
      position=!!sym(position_list[[i]]),
      effect_allele=!!sym(effect_allele_list[[i]]),
      z=!!sym(z_statistic_list[[i]]))
    if(i==1) used_snps=gwas$rsid else used_snps=intersect(used_snps,gwas$rsid)
    gwas_list[[i]]=gwas %>% filter(rsid %in% used_snps)
  }
  gwas_list=lapply(gwas_list, function(h) h %>% filter(rsid %in% used_snps))
  chrs=lapply(gwas_list,function(h) h %>% pull(chr) %>% unique())
  tt=table(unlist(chrs))
  chrs=as.numeric(names(tt[tt==k]))
  chrs=intersect(1:22,chrs)
  sgp=list()
  rdf1=rdf2=rdf3=data.frame()
  for(cc in 1:length(chrs)) {
    k=0
    chrlist=list()
    if(verbose) cat('Chromosome',chrs[cc],'\n')
    # subset GWAS to this chromosome
    gwas_chr=lapply(gwas_list,function(h) h %>% filter(chr==chrs[cc]))
    ## load LD matrices (about 0.5 gb per chromsome per population)
    # setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices/chr_specific_files')
    # ld_list=lapply(1:k,\(.) readRDS(paste0(ld_population_list[[.]],'/chr',cc,'.Rds')))
    # ld_list is a length-k list. each entry is a list whose entries are gene-specific LD matrices
    index_chr=index %>% filter(chr==chrs[cc])
    genes=unique(index_chr$symbol)
    #pb=txtProgressBar(min=0,max=length(genes),style=3)
    for(i in 1:length(genes)) {
      #if(verbose) setTxtProgressBar(pb, i)
      beep=tryCatch(
        {
          # load gene-specific LD matrices (more memory efficient than loading all genes at once)
          ldi=lapply(names(ld_population_list),function(h) {
            setwd(ld_directory)
            setwd(h)
            setwd(paste0('chr',chrs[cc]))
            fps=dir()
            fps=sapply(fps,function(o) unlist(strsplit(o,'[.]'))[1])
            ix=which(fps==genes[i])
            if(length(ix)==0) stop()
            as.matrix(readRDS(paste0(fps[ix],'.Rds')))
          })
          # dataframes of SNP-effect allele pairs
          ld_df=lapply(ldi,function(h) {
            spp=sapply(rownames(h), function(.) strsplit(.,'_'))
            do.call(rbind,spp) %>% as.data.frame() %>% rename(rsid=V1,a1=V2)
          })
          ## overlapping SNPs in LD matrices
          #   ldi=lapply(ld_list,function(h) h[[which(names(h)==genes[i])]])
          ldsnps=list(); for(o in 1:k) ldsnps[[o]]=unname(sapply(rownames(ldi[[o]]),function(h) unlist(strsplit(h,'_'))[1]))
          tt=table(unlist(ldsnps))
          ldsnps=names(tt[tt==k])
          # use SNPs in both LD and GWAS (all GWAS already share the same SNPs because of code above)
          usesnps=intersect(ldsnps,gwas_chr[[1]]$rsid)
          ldi=lapply(ldi,function(h) {
            h=as.matrix(h)
            rownames(h)=colnames(h)=sapply(rownames(h),function(hh) unlist(strsplit(hh,'_'))[1])
            h[usesnps,usesnps]
          })
          # match GWAS SNP to LD SNPs' order
          gwas_chri=lapply(gwas_chr,function(h) {
            x=h %>% filter(rsid %in% usesnps) %>% as.data.frame()
            rownames(x)=x$rsid
            x[rownames(ldi[[1]]),]
          })
          # harmonize SNPs to their respective effect alleles
          for(o in 1:k) {
            eao=ld_df[[o]]
            gwas_chri[[o]]=gwas_chri[[o]] %>%
              left_join(eao,by='rsid') %>%
              mutate(z=ifelse(effect_allele==a1,z,-z)) %>%
              select(-effect_allele) %>%
              rename(effect_allele=a1)
          }
          # Z-matrix
          Z=lapply(gwas_chri,function(h) h %>% pull(z)) %>% do.call(cbind,.)
          mugent_result=mugent(Z,ldi)
          mugentph_result=mugent_ph(Z,ldi)
          mugentpleio_result=mugent_pleio(Z,ldi,alpha=mugentpleio_alpha)
          # store results
          l1=data.frame(test='MuGenT',gene=genes[i],chr=chrs[cc]) %>% bind_cols(as.data.frame(mugent_result))
          l2=data.frame(test='MuGenT-PH',gene=genes[i],chr=chrs[cc]) %>% bind_cols(as.data.frame(mugentph_result))
          l3=data.frame(test='MuGenT-Pleiotropy',gene=genes[i],chr=chrs[cc]) %>% bind_cols(as.data.frame(mugentpleio_result))
          rdf1=rbind(rdf1,l1)
          rdf2=rbind(rdf2,l2)
          rdf3=rbind(rdf3,l3)
          # record SNP-gene pairs
          k=k+1
          if(return_snp_gene_pairs) {
            chrlist[[k]]=paste(gwas_chri$rsid,collapse=',')
            names(chrlist)[k]=genes[i]
          }
        },
        error=function(x) NA
      )
      # if(is.logical(beep)) cat(i,'\n')
    }
    #close(pb)
    if(return_snp_gene_pairs) {
      sgp[[cc]]=chrlist
      names(sgp)[[cc]]=paste0('chr',chrs[cc])
    }
  }
  outlist=list('MuGenT'=rdf1, 'MuGenT-PH'=rdf2, 'MuGenT-Pleiotropy'=rdf3)
  if(return_snp_gene_pairs) outlist$snp_sets=sgp

}

#' Fine-mapping gene-based association test statistics
#'
#' This function performs fine-mapping of gene-based association tests statistics.
#' @param gent_results Direct output of \code{gent_genomwide()} or \code{mugent_genomewide()$MuGenT}
#' @param ld_population LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR'
#' @param gwas_n GWAS sample size
#' @param index_genes Vector of index genes to perform fine-mapping around. Will be found automatically using clumping if NULL
#' @param chromosome Column name of gene chromosome in \code{gent_results}
#' @param gene_start Column name of gene start position (hg19) in \code{gent_results}
#' @param gene_symbol Column name of gene symbol in \code{gent_results}
#' @param pval Column name of gene-based test P-value in \code{gent_results}
#' @param null_mean Column name of gene-based test null mean (number of tested SNPs) in \code{gent_results}
#' @param null_variance Column name of gene-based test null variance in \code{gent_results}
#' @param window_kb_width Kb of window in which fine-mapping should be performed (total window will be this size)
#' @param R_ridge_penalty Ridge penalty to add to gene-gene correlation matrix
#' @param clump_p P-value threshold to use when finding loci in which to perform fine-mapping (only invoked if \code{index_genes} is NULL)
#' @param clump_r2 Gene-gene squared correlation upper threshold to use when finding loci in which to perform fine-mapping (only envoked if \code{index_genes} is NULL). Genes correlated below the square root of this threshold will be considered as tagging separate loci.
#' @param verbose TRUE if progress should be printed to the console, FALSE otherwise
#' @param index Gene index file. see \code{data(EnsemblHg19GenePos)} for the expected format.
#' @return A dataframe with these components:
#' \itemize{
#'  \item `gene`: Gene symbol
#'  \item `locus_index_gene`: This gene is at the center of the fine-mapped locus and was either the most significant gene or was set by you
#'  \item `chr`: Chromosome location
#'  \item `gene_start`: Start position of gene using Ensembl hg19
#'  \item `finemapping_pip`: SuSiE posterior inclusion probability (PIP)
#'  \item `finemapping_CS`: SuSiE credible set (NA if not in any credible set; set numbering is arbitrary)
#' }
#'
#' @export
#' @examples
#' # Example for Alzheimer's disease GWAS data
#' gent_finemap(
#'   gent_results=gent_results,
#'   ld_population='EUR',
#'   gwas_n=50000,
#'   index_genes=c('CR1','MYCL'),
#'   window_kb_width=2000,
#'   verbose=TRUE)
gent_finemap=function(
    gent_results,
    ld_population,
    gwas_n,
    index_genes=NULL,
    chromosome='chr',
    gene_start='gene_start',
    gene_symbol='gene',
    pval='pval',
    null_mean='mu_h0',
    null_variance='sigma2_h0',
    window_kb_width=2000,
    R_ridge_penalty=0,
    clump_p=0.05/12727,
    clump_r2=0.01,
    verbose=TRUE,
    index=NULL, ...) { # ...'s are for susieR::susie_rss()
  ## steps:
  # 1) Find index genes
  # 2) perform fine-mapping across index genes
  # load gent statistic correlations for this population
  if(toupper(ld_population)=='EUR') {data(EURGenTStatLD);gent_ld=EURGenTStatLD}
  if(toupper(ld_population)=='AFR') {data(AFRGenTStatLD);gent_ld=AFRGenTStatLD}
  if(toupper(ld_population)=='EAS') {data(EASGenTStatLD);gent_ld=EASGenTStatLD}
  if(toupper(ld_population)=='SAS') {data(SASGenTStatLD);gent_ld=SASGenTStatLD}
  if(toupper(ld_population)=='AMR') {data(AMRGenTStatLD);gent_ld=AMRGenTStatLD}
  # halve the total window size because I will go left then right later
  window_kb_width=round(window_kb_width/2)
  # load gene positional index
  if(is.null(index)) {data(EnsemblHg19GenePos);index=EnsemblHg19GenePos}
  # clean gent results and perform asymptotic transformation
  gent_results=gent_results %>%
    rename(gene=!!sym(gene_symbol),
           chr=!!sym(chromosome),
           gene_start=!!sym(gene_start),
           pval=!!sym(pval),
           null_mean=!!sym(null_mean),
           null_variance=!!sym(null_variance)) %>%
    mutate(rate=null_mean/null_variance, shape=null_mean*rate) %>%
    mutate(Q=qgamma(pval,shape=shape,rate=rate,,lower.tail=FALSE)) %>%
    mutate(a=(Q/null_mean-1)*null_mean/sqrt(null_variance))
  notused=c()
  if(!is.null(index_genes)) {
    ug=paste(index_genes[1:3],collapse=', ')
    if(verbose) message(paste0('\nchoosing loci around ',ug,' ...'))
    # use genes that the user gave
    # subset gene positional index to just genes to be analyzed
    index=index %>% filter(symbol %in% index_genes)
    use_genes=index$symbol
  } else {
    # find index genes using clump procedure
    if(verbose) message('\nclumping genes to identify loci')
    clumps=gene_clump(
      gent_results,
      ld_population,
      chromosome=chromosome,
      gene_start=gene_start,
      gene_symbol=gene_symbol,
      pval=pval,
      clump_p=clump_p,
      clump_kb=window_kb_width,
      clump_r2=clump_r2,
      verbose=FALSE)
    use_genes=names(clumps)
  }
  # 2) now perform fine-mapping using selected genes as input loci
  rdf=data.frame()
  for(i in 1:length(use_genes)) {
    # subset everything to this gene's chromosome
    chri=index %>% filter(symbol==use_genes[i]) %>% head(.,1) %>% pull(chr)
    gent_resultsi=gent_results %>% filter(chr==chri)
    gent_ld_chr=as.matrix(gent_ld[[paste0('chr',chri)]])
    if(!(use_genes[i] %in% gent_resultsi$gene)) {notused=c(notused,use_genes[i]);next}
    indexbp=gent_results %>% filter(gene==use_genes[i]) %>% head(.,1) %>% pull(gene_start)
    indexbp=c(indexbp,indexbp)+c(-window_kb_width,window_kb_width)*1e3
    gent_resultsi=gent_results %>% filter(data.table::between(gene_start,indexbp[1],indexbp[2]))
    # correlation matrix for these genes
    ldgenes=intersect(gent_resultsi$gene,rownames(gent_ld_chr))
    gent_resultsi=gent_resultsi %>% filter(gene %in% ldgenes)
    gent_ldi=gent_ld_chr[gent_resultsi$gene,gent_resultsi$gene]
    m=nrow(gent_ldi)
    gent_ldi=(1-R_ridge_penalty)*gent_ldi+R_ridge_penalty*diag(m)
    # apply fine-mapping here
    fit=suppressWarnings(susieR::susie_rss(z=gent_resultsi$a,R=gent_ldi,n=gwas_n,verbose=FALSE,...))
    # fit=susieR::susie_rss(z=gent_resultsi$a,R=gent_ldi,n=gwas_n,verbose=FALSE,...)
    # assign credible sets to each SNP
    cs=rep(NA,m)
    if(!is.null(fit$sets$cs)) for(o in 1:length(fit$sets$cs)) cs[fit$sets$cs[[o]]]=o
    gent_resultsi=gent_resultsi %>%
      mutate(finemapping_pip=unname(fit$pip), finemapping_CS=cs) %>%
      mutate(locus_index_gene=use_genes[i]) %>%
      select(gene,locus_index_gene,chr,gene_start,finemapping_pip,finemapping_CS)
    rdf=rbind(rdf,gent_resultsi)
  }
  # print messages?
  if(verbose & length(notused)>0) {
    notused=paste(notused,collapse=', ')
    message(paste0('(',notused,') had insufficient SNPs to estimate gene-level correlations'))
  }
  return(rdf)
}

#' Genome-wide inverse MAF-weighted gene-based association test
#'
#' This function performs gene-based association testing genome-wide using weights.
#' @param gwas GWAS summary statistics
#' @param snp_weights Column name in \code{gwas} indicating SNP weights (e.g., square root of inverse of MAF)
#' @param KbWindow Kilobase window used to assign SNPs to genes
#' @param ld_population LD population, one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR'
#' @param ld_directory Relative or absolute filepath to LD directory (see Github Wiki)
#' @param snp Column name of SNP rsID in GWAS data
#' @param chromosome Column name of SNP chromosome in GWAS data
#' @param position Column name of SNP hg19 base pair position in GWAS data
#' @param effect_allele Column name of SNP effect allele in GWAS data
#' @param effect_size Column name of SNP estimated effect size in GWAS data
#' @param standard_error Column name of SNP standard error in GWAS data
#' @param index Gene index file. see \code{data(EnsemblHg19GenePos)} for the expected format.
#' @param verbose TRUE if progress should be printed to the console, FALSE otherwise
#' @param return_snp_gene_pairs If TRUE, will return a chromosome-specific list of tested gene-specific SNP sets annotated by gene
#' @return A dataframe with these components:
#' \itemize{
#'  \item `pval`: P-value testing the gene-based null hypothesis
#'  \item `shape`: Shape parameter of the Gamma null distribution
#'  \item `rate`: Rate parameter of the Gamma null distribution
#'  \item `mu_h0`: Expected value of the gene-based test statistic under the null hypothesis
#'  \item `sigma2_h0`: Variance of the gene-based test statistic under the null hypothesis
#'  \item `gene`: Gene symbol
#'  \item `m`: Number of SNPs tested in gene-specific set
#'  \item `chr`: Chromosome of gene
#'  \item `gene_start`: Start position of gene (Ensembl hg19)
#'  \item `window_start`: Start position of window in which SNPs were captured for this gene (hg19)
#'  \item `gene_end`: End position of gene (Ensembl hg19)
#'  \item `window_end`: End position of window in which SNPs were captured for this gene (hg19)
#' }
#'
#' @export
#' @examples
#' # Example for Alzheimer's disease GWAS data
#' gent_genomewide(
#'   gwas=ad_gwas,
#'   snp_weights='sqrt_inverse_mafs',
#'   ld_population='EUR',
#'   ld_directory='ld_matrices',
#'   KbWindow=50,
#'   snp='MarkerName',
#'   chromosome='Chromosome',
#'   position='Position',
#'   effect_allele='Effect_allele',
#'   effect_size='beta',
#'   standard_error='se',
#'   verbose=TRUE,
#'   return_snp_gene_pairs=FALSE)
wgent_genomewide=function(gwas,
                          snp_weights='weights',
                          KbWindow=50,
                          ld_population='EUR',
                          ld_directory='ld_directory',
                          snp='rsid',
                          chromosome='chr',
                          position='position',
                          effect_allele='effect_allele',
                          effect_size='beta',
                          standard_error='se',
                          index=NULL,
                          verbose=TRUE,
                          return_snp_gene_pairs=FALSE) {
  if(is.null(index)) {data(EnsemblHg19GenePos);index=EnsemblHg19GenePos}
  setwd(ld_directory)
  gwas=gwas %>% na.omit()
  gwas=gwas %>%
    rename(rsid=!!sym(snp),
           chr=!!sym(chromosome),
           position=!!sym(position),
           effect_allele=!!sym(effect_allele),
           effect_size=!!sym(effect_size),
           standard_error=!!sym(standard_error),
           weight=!!sym(snp_weights))
  chrs=gwas %>% select(chr) %>% pull() %>% unique() %>% as.numeric() %>% na.omit() %>% sort()
  chrs=intersect(1:22,chrs)
  sgp=list()
  rdf=data.frame()
  for(cc in 1:length(chrs)) {
    setwd(ld_directory)
    setwd(ld_population)
    setwd(paste0('chr',cc))
    if(verbose) cat('Starting chromosome', chrs[cc], '\n')
    gwas_chr=gwas %>% filter(chr==chrs[cc])
    index_chr=index %>% filter(chr==cc)
    genes=unique(index_chr$symbol)
    #pb=txtProgressBar(min=0,max=length(genes),style=3)
    k=0
    chrlist=list()
    for(i in 1:length(genes)) {
      #if(verbose) setTxtProgressBar(pb, i)
      tryCatch(
        {
          # check up front that the data contains at least one SNP for this gene
          starti=index_chr %>% filter(symbol==genes[i]) %>% select(start) %>% head(.,1) %>% pull()
          endi=index_chr %>% filter(symbol==genes[i]) %>% select(end) %>% head(.,1) %>% pull()
          starti=starti-KbWindow*1e3
          endi=endi+KbWindow*1e3
          gwas_chri=gwas_chr %>% filter(position>starti,position<endi)
          if(nrow(gwas_chri)==0) stop()
          # then proceed to load LD
          genefp=paste0(genes[i],'.Rds')
          if(!file.exists(genefp)) stop()
          ld=readRDS(genefp)
          ldrn=rownames(ld)
          spp=sapply(ldrn,\(.) unlist(strsplit(.,'_')))
          ld_df=t(spp) %>% as.data.frame() %>% rename(rsid=V1,a1=V2)
          ldrn=ld_df$rsid
          usesnps=intersect(ldrn,gwas_chri$rsid)
          gwas_chri=gwas_chri %>%
            filter(rsid %in% usesnps) %>%
            distinct(rsid,.keep_all=TRUE) %>%
            arrange(position) %>%
            left_join(ld_df,by='rsid') %>%
            mutate(effect_size=ifelse(effect_allele==a1,effect_size,-effect_size))
          betahat=gwas_chri %>% pull(effect_size)
          se=gwas_chri %>% pull(standard_error)
          w=gwas_chri %>% pull(weight)
          # make sure data contains SNPs for this gene
          if(length(betahat)==0) stop()
          rownames(ld)=colnames(ld)=ldrn
          ld=ld[gwas_chri$rsid,gwas_chri$rsid]
          # transform LD to be covariance matrix of estimated effect sizes (the math works out)
          D=diag(c(se))
          W=diag(c(w))
          result=gent(zs=betahat,LD=D%*%ld%*%D,A=W^2) # squared bc weighted version is W%*%betahat
          result_unweighted=gent(zs=betahat/se,LD=ld)
          toadd=as.data.frame(result) %>%
            mutate(gene=genes[i],
                   m=nrow(ld),
                   chr=chrs[cc],
                   gene_start=starti+KbWindow*1e3,
                   window_start=starti,
                   gene_end=endi-KbWindow*1e3,
                   window_end=endi) %>%
            mutate(pval_unweighted=result_unweighted$pval)
          rdf=rbind(rdf,toadd)
          # record SNP-gene pairs
          k=k+1
          if(return_snp_gene_pairs) {
            chrlist[[k]]=paste(gwas_chri$rsid,collapse=',')
            names(chrlist)[k]=genes[i]
          }
        },
        error=function(x) NA
      )
    }
    #close(pb)
    if(return_snp_gene_pairs) {
      sgp[[cc]]=chrlist
      names(sgp)[[cc]]=paste0('chr',chrs[cc])
    }
  }
  if(return_snp_gene_pairs) {
    return(list(result=as_tibble(rdf),snp_sets=sgp))
  } else {
    return(as_tibble(rdf))
  }
}

#' Clumping of gene-based test statistics
#'
#' This function performs PLINK-style clumping of gene-based test statistics
#' @param gentres Gene-based association test statistic results for multiple genes (e.g., output of <gent/mugent/xgent>_genomewide())
#' @param ld_population one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR' from 1000 Genomes Phase 3. Used for internal loading of GenT correlation matrices.
#' @param chromosome Chromosome variable name in \code{genedf}
#' @param gene_start Gene start position variable name in \code{genedf}
#' @param symbol Gene symbol variable name in \code{genedf}
#' @param pval Gene-based test statistic P-value variable name in \code{genedf}
#' @param clump_p Only genes with a P-value less than this threshold may index a locus for clumping
#' @param clump_kb Kilobase size of the entire clumping window. Left and right windows from the index gene will be half the size of \code{clump_kb}
#' @param clump_r2 Only genes correlated with lead genes beyond this threshold may be clumped to other genes
#' @param verbose TRUE if progress should be printed to the console, FALSE otherwise
#' @return A vector of clumped genes
#' @export
#' @examples
#' gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
#' gene_clump(gent_results,'EUR')
gene_clump=function(gentres,
                    ld_population,
                    chromosome='chr',
                    gene_start='gene_start',
                    gene_symbol='symbol',
                    pval='pval',
                    clump_p=0.05/12727,
                    clump_kb=1000,
                    clump_r2=0.01,
                    verbose=TRUE) {
  # load correlations
  if(toupper(ld_population)=='EUR') {data(EURGenTStatLD);gent_ld=EURGenTStatLD}
  if(toupper(ld_population)=='AFR') {data(AFRGenTStatLD);gent_ld=AFRGenTStatLD}
  if(toupper(ld_population)=='EAS') {data(EASGenTStatLD);gent_ld=EASGenTStatLD}
  if(toupper(ld_population)=='SAS') {data(SASGenTStatLD);gent_ld=SASGenTStatLD}
  if(toupper(ld_population)=='AMR') {data(AMRGenTStatLD);gent_ld=AMRGenTStatLD}
  gentres=gentres %>%
    dplyr::rename(symbol=!!sym(gene_symbol),
                  gene_start=!!sym(gene_start),
                  pval=!!sym(pval),
                  chr=!!sym(chromosome))
  df=gentres %>%
    dplyr::select(symbol,gene_start,pval,chr) %>%
    dplyr::filter(pval<clump_p) %>%
    dplyr::arrange(pval)
  if(nrow(df)==0) {if(verbose) {cat('no significant clumps; returning NA\n')};return(NA)}
  clump_kb=round(clump_kb/2)
  clumps=list()
  skipped_genes=c()
  k=0
  # loop over each chromosome
  chrs=unique(df$chr)
  for(cc in 1:length(chrs)) {
    gent_ldchr=gent_ld[[paste0('chr',chrs[cc])]] %>% as.matrix()
    df_chr=df %>% dplyr::filter(chr==chrs[cc])
    for(i in 1:nrow(df_chr)) {
      allclumps=unlist(clumps);allclumps=c(names(clumps),allclumps)
      allclumps=unname(allclumps)
      if(df_chr$symbol[i] %in% allclumps) next
      k=k+1
      genesaround=df_chr %>%
        dplyr::filter(abs(gene_start-df_chr$gene_start[i])<(clump_kb*1e3)) %>%
        dplyr::filter(symbol!=df_chr$symbol[i])
      if(nrow(genesaround)==0) {
        clumps[[k]]=NA
        names(clumps)[k]=df_chr$symbol[i]
        next
      }
      # attach LD with this SNP
      # if gene not in LD reference, use closest gene as proxy
      if(!(df_chr$symbol[i] %in% rownames(gent_ldchr))) {
        posi=df_chr$gene_start[i]
        closest_gene=df_chr %>%
          dplyr::filter(symbol!=df_chr$symbol[i]) %>%
          dplyr::mutate(d=abs(gene_start-posi)) %>%
          dplyr::arrange(d) %>%
          head(.,1) %>%
          dplyr::pull(symbol)
        if(closest_gene %in% names(clumps) || closest_gene %in% unlist(clumps)) next
        df_chr$symbol[i]=closest_gene
        df_chr=df_chr %>% dplyr::distinct()
      }
      ix=which(rownames(gent_ldchr)==df_chr$symbol[i])
      lddf=data.frame(symbol=rownames(gent_ldchr),r2=gent_ldchr[ix,]^2) %>%
        filter(symbol %in% genesaround$symbol) %>%
        filter(r2>clump_r2)
      if(nrow(lddf)==0) {
        # need to have nearby SNPs in LD to be clumps
        clumps[[k]]=NA
        names(clumps)[k]=df_chr$symbol[i]
        next
      }
      clumps[[k]]=lddf %>% dplyr::filter(!(symbol %in% allclumps)) %>% pull(symbol)
      names(clumps)[k]=df_chr$symbol[i]
    }
  }
  clumps=lapply(clumps,function(h) if(length(h)==0) NA else h)
  if(length(skipped_genes)>0) cat(skipped_genes,'not in LD reference\n')
  clumps
}

#' Gene-based test Manhattan plot
#'
#' This function generates a Manhattan plot of gene-based test statistics
#' @param gentres Gene-based association test statistic results for multiple genes (e.g., output of <gent/mugent/xgent>_genomewide())
#' @param chromosome Column name of gene chromosome in \code{gentres} data
#' @param gene_position Column name of gene base pair position (e.g., TSS, start position) in \code{gentres} data
#' @param gene_pvalue Column name of gene-based test P-values in  \code{gentres} data
#' @param significance_threshold Genome-wide significance threshold. Used to draw a reference line and for annotating genes if you set \code{label_genes=TRUE}
#' @param plot_threshold Lower P-value threshold at which to truncate gene-based test P-values.
#' @param label_genes If TRUE, clumped genes (see \code{gene_clump()}) will be annotated at their position
#' @param gene_label Column name of gene identifier (e.g., gene symbol) in \code{gentres} data
#' @param ld_population LD population, one of 'EUR','AFR','EAS','SAS', or 'AMR'
#' @return A \code{ggplot2} object which is a Manhattan plot of gene-based association test P-values
#' @export
#' @examples
#' gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
#' gent_manhattan(gent_results)
gent_manhattan=function(gentres,
                        chromosome='chr',
                        gene_position='gene_start',
                        gene_pvalue='pval',
                        significance_threshold=0.05/12727,
                        plot_threshold=1e-100,
                        label_genes=FALSE,
                        gene_label='gene',
                        ld_population='EUR') {
  if(label_genes) require(ggrepel)
  chrdf=gentres %>%
    dplyr::mutate(!!sym(chromosome):=as.numeric(!!sym(chromosome)),
                  !!sym(gene_position):=as.numeric(!!sym(gene_position)),
                  !!sym(gene_pvalue):=as.numeric(!!sym(gene_pvalue))) %>%
    dplyr::mutate(!!sym(gene_pvalue):=ifelse(!!sym(gene_pvalue)<1e-300,1e-300,!!sym(gene_pvalue))) %>%
    dplyr::group_by(!!sym(chromosome)) %>%
    dplyr::mutate(newbp=!!sym(chromosome)+!!sym(gene_position)/max(!!sym(gene_position),na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(chrcol=!!sym(chromosome) %in% seq(1,22,2)) %>%
    dplyr::mutate(!!sym(gene_pvalue):=ifelse(!!sym(gene_pvalue)<plot_threshold,plot_threshold,!!sym(gene_pvalue)))
  axisdf=chrdf %>%
    dplyr::group_by(!!sym(chromosome)) %>%
    dplyr::summarise(newbp=median(newbp,na.rm=TRUE)) %>%
    dplyr::ungroup()
  manplot=chrdf %>%
    ggplot(aes(newbp,-log10(!!sym(gene_pvalue)),color=chrcol)) +
    geom_point(pch=19) +
    scale_color_manual(values=c('skyblue','royalblue')) +
    geom_hline(yintercept=-log10(significance_threshold)) +
    scale_x_continuous(breaks=axisdf$newbp,labels=axisdf$chr) +
    theme_bw() +
    theme(legend.position='none',
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    labs(x='chromosome',
         y=expression('-log'[10]*'(gene-based P-value)'))
  manplot
  if(plot_threshold>0) manplot=manplot + labs(subtitle=paste('P-values truncated to',as.character(plot_threshold)))
  if(label_genes) {
    clumps=gene_clump(chrdf,
                      ld_population=ld_population,
                      chromosome=chromosome,
                      gene_start=gene_position,
                      gene_symbol=gene_label,
                      pval=gene_pvalue,
                      clump_p=significance_threshold,
                      clump_kb=1000,
                      clump_r2=0.01,
                      verbose=TRUE)
    if(length(clumps)>0) {
      labdf=chrdf %>% dplyr::filter(!!sym(gene_label) %in% names(clumps))
      manplot=manplot +
        geom_text_repel(aes(newbp,-log10(!!sym(gene_pvalue)),label=!!sym(gene_label)),
                        max.overlaps=Inf,
                        ylim=c(-log10(significance_threshold),-log10(plot_threshold)),
                        force=10,
                        data=labdf,
                        color='black',
                        min.segment.length=0)
    }
  }
  return(manplot)
}

#' Gene-based test QQ plot
#'
#' This function generates a quantile-quantile plot of gene-based test statistics
#' @param gentres Gene-based association test statistic results for multiple genes (e.g., output of <gent/mugent/xgent>_genomewide())
#' @param gene_pvalue Column name of gene-based test P-values in  \code{gentres} data
#' @return A \code{ggplot2} object which is a QQ plot of gene-based association test P-values
#' @export
#' @examples
#' gent_results=gent_genomewide(gwas_data,50,'EUR','ld_directory')
#' gent_qq(gent_results)
gent_qq=function(gentres,gene_pvalue='pval') {
  pv=as.data.frame(gentres)[,gene_pvalue]
  pv=sort(na.omit(pv))
  pp=ppoints(length(pv))
  ppdf=data.frame(e=pp,o=pv)
  qqplot=ppdf %>%
    ggplot(aes(-log10(e),-log10(o))) +
    geom_point(pch=19,color='royalblue') +
    geom_abline(intercept=0,slope=1,lwd=2/3) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    labs(x=expression('expected -log'[10]*'(gene-based P-value)'),
         y=expression('observed -log'[10]*'(gene-based P-value)'))
  return(qqplot)
}

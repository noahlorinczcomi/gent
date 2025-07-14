#' Gene-based association test (GenT)
#'
#' This function performs a single gene-based association test.
#' @param zs vector of SNP-specific Z-statistics from GWAS whose indices correspond to those of \code{LD}
#' @param LD LD matrix whose row/column indices correspond to the indices of \code{zs}
#' @param mafs Minor allele frequencies of the SNPs whose Z-statistics are in \code{zs} (in corresponding order). If non-null, inverse 2MAF(1-MAF) weights will be aplied.
#' @param chisquares if you instead give chi-square statistics (squared Z-statistics), supply them here and leave \code{zs} empty
#' @return A list with these components:
#' \itemize{
#'  \item \code{pval}: P-value for testing H0: gene is not associated with trait.
#'  \item \code{shape}: shape parameter of null (Gamma) distribution.
#'  \item \code{rate}: rate parameter of null (Gamma) distribution.
#'  \item \code{mu_h0}: expectation of null (Gamma) distribution.
#'  \item \code{sigma2_h0}: variance of null (Gamma) distribution.
#'  \item \code{mu_h1}: observed statistic used to test H0.
#'  \item \code{sigma2_h1}: variance of observed statistic used to test H0.
#' }
#'
#' @export
#' @examples
#' # Example for 5 SNPs
#' LD=cov2cor(rWishart(1,100,diag(5))[,,1])
#' z=c(mvnfast::rmvn(1,rep(0,5),diag(5)))
#' gent(z,LD)
gent=function(zs=NULL,LD,mafs=NULL,xqtl_Z=NULL,chisquares=NULL) {
    # find null distribution
    if(is.null(mafs) & is.null(xqtl_Z)) {
      mu=nrow(LD)
      trASAS=tr(LD%*%LD)
    } else if(!is.null(mafs) & is.null(xqtl_Z)){
     # if the analysis is MAF-weighted, return early
      mafs=c(mafs)
      if(length(mafs)!=length(zs)) stop('length of MAF vector not equal to length of Z-statistic vector')
      mafs=2*mafs*(1-mafs)
      A=diag(1/mafs)
      mu=sum(diag(A%*%LD))
      trASAS=tr(A%*%LD%*%A%*%LD)
      mu_1=sum(diag(A%*%LD))+t(zs)%*%A%*%zs
      sigma2_h1=2*trASAS+4*t(zs)%*%A%*%LD%*%A%*%zs
      beta=mu/(2*trASAS)
      alpha=beta*mu
      y=c(t(zs)%*%A%*%zs)
      pval=pgamma(y,shape=alpha,rate=beta,lower.tail=FALSE)
      return(list(pval=pval,shape=alpha,rate=beta,mu_h0=mu,sigma2_h0=2*trASAS,mu_h1=mu_h1,sigma2_h1=sigma2_h1))
    } else if(is.null(mafs) & !is.null(xqtl_Z)) {
      xqtl_Z=as.matrix(xqtl_Z);m=nrow(xqtl_Z);p=ncol(xqtl_Z)
      L=matrix(0,m,m);for(o in 1:p) L=L+xqtl_Z[,o]%*%t(xqtl_Z[,o])
      L=L/sqrt(m*p)
      mu=sum(diag(L%*%LD))
      trASAS=tr(L%*%LD%*%L%*%LD)
    }
    sigma2=2*trASAS
    beta=(mu/trASAS)/2
    alpha=beta*mu
    # P-values for zs (should be in LD)
    if(!is.null(chisquares)) y=sum(chisquares) else y=sum(zs^2)
    mu_h1=mu+y
    if(!is.null(zs)) sigma2_h1=sigma2+4*t(zs)%*%LD%*%zs else sigma2_h1=NULL
    pval=pgamma(y,shape=alpha,rate=beta,lower.tail=FALSE)
    out=list(pval=pval,shape=alpha,rate=beta,mu_h0=mu,sigma2_h0=sigma2,mu_h1=mu_h1,sigma2_h1=sigma2_h1)
    lapply(out,c)
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






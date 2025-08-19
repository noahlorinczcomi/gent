# GenT
gent_genomewide=function(gwas,
                         KbWindow=50,
                         ld_population='EUR',
                         ld_directory='ld_directory',
                         snp='rsid',
                         chromosome='chr',
                         position='position',
                         effect_allele='effect_allele',
                         z='z',
                         index=NULL,
                         verbose=TRUE) {
  # expects:
  ## 1) gwas: GWAS summary statistics with at least rsid, chr, position, effect_allele, z
  ## 2) ld_population: population of LD. data will be list of LD matrices. each entry corresponds to one chromosome and is a list. each entry of each list entry is a gene-specific LD matrix of +/-1Mb 1kg v3 SNPs
  ## 3) KbWindow: Kb window size of SNPs to use
  ## 4) index: dataframe of (gene,chr,start,end) positions
  ## 5) savefp: directory to save results in (will also return them)
  if(is.null(index)) {data(EnsemblHg19GenePos);index=EnsemblHg19GenePos}
  setwd(ld_directory)
  gwas=gwas %>% na.omit()
  gwas=gwas %>% rename(rsid=!!sym(snp), chr=!!sym(chromosome), position=!!sym(position), effect_allele=!!sym(effect_allele), z=!!sym(z))
  chrs=gwas %>% select(chr) %>% pull() %>% unique() %>% as.numeric() %>% na.omit() %>% sort()
  rdf=data.frame()
  for(cc in 1:length(chrs)) {
    setwd(ld_directory)
    setwd(ld_population)
    setwd(paste0('chr',cc))
    if(verbose) cat('Starting chromosome', chrs[cc], '\n')
    gwas_chr=gwas %>% filter(chr==chrs[cc])
    index_chr=index %>% filter(chr==cc)
    genes=unique(index_chr$symbol)
    pb=txtProgressBar(min=0,max=length(genes),style=3)
    for(i in 1:length(genes)) {
      if(verbose) setTxtProgressBar(pb, i)
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
        },
        error=function(x) NA
          )
    }
    close(pb)
  }
  rdf %>% as_tibble()
}

# MuGenT
mugent_genomewide=function(
    gwas_list,
    ld_population_list,
    ld_directory='ld_directory',
    KbWindow=50,
    snp_list=lapply(1:length(gwas_list),'rsid'),
    chromosome_list=lapply(1:length(gwas_list),'chr'),
    position_list=lapply(1:length(gwas_list),'position'),
    effect_allele_list=lapply(1:length(gwas_list),'effect_allele'),
    z_list=lapply(1:length(gwas_list),'z'),
    # must use Z-stats
    mugentpleio_alpha=0.05/12727,index=NULL,verbose=TRUE) {
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
        z=!!sym(z_list[[i]]))
    if(i==1) used_snps=gwas$rsid else used_snps=intersect(used_snps,gwas$rsid)
    gwas_list[[i]]=gwas %>% filter(rsid %in% used_snps)
  }
  gwas_list=lapply(gwas_list, function(h) h %>% filter(rsid %in% used_snps))
  chrs=lapply(gwas_list,function(h) h %>% pull(chr) %>% unique())
  tt=table(unlist(chrs))
  chrs=as.numeric(names(tt[tt==k]))
  rdf1=rdf2=rdf3=data.frame()
  for(cc in 1:length(chrs)) {
    if(verbose) cat('Chromosome',chrs[cc],'\n')
    # subset GWAS to this chromosome
    gwas_chr=lapply(gwas_list,function(h) h %>% filter(chr==chrs[cc]))
    ## load LD matrices (about 0.5 gb per chromsome per population)
    # setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices/chr_specific_files')
    # ld_list=lapply(1:k,\(.) readRDS(paste0(ld_population_list[[.]],'/chr',cc,'.Rds')))
    # ld_list is a length-k list. each entry is a list whose entries are gene-specific LD matrices
    index_chr=index %>% filter(chr==chrs[cc])
    genes=unique(index_chr$symbol)
    pb=txtProgressBar(min=0,max=length(genes),style=3)
    for(i in 1:length(genes)) {
        if(verbose) setTxtProgressBar(pb, i)
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
            },
            error=function(x) NA
        )
        # if(is.logical(beep)) cat(i,'\n')
    }
    close(pb)
  }
  return(list('MuGenT'=rdf1, 'MuGenT-PH'=rdf2, 'MuGenT-Pleiotropy'=rdf3))
}


# gent_finemap()
gent_finemap=function(
    gent_results,
    ld_population,
    gwas_n,
    index_genes=NULL,
    chromosome='chr',
    gene_start='gene_start',
    symbol='gene',
    pval='pval',
    null_mean='mu_h0',
    null_variance='sigma2_h0',
    pval_threshold=0.05/12727,
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
    rename(gene=!!sym(symbol),
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
      symbol=symbol,
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

# function to perform PLINK-style clumping of gene-based test statistics
gene_clump=function(genedf,ld_population,chromosome='chr',gene_start='gene_start',symbol='symbol',pval='pval',
                    clump_p=0.05/12727,clump_kb=1000,clump_r2=0.01,verbose=TRUE) {
  # genedf: gene-based association test statistic results for multiple genes (e.g., output of <gent/mugent/xgent>_genomewide())
  # ld_population: one of 'EUR', 'AFR', 'EAS', 'SAS', or 'AMR' from 1000 Genomes Phase 3. Used for internal loading of GenT correlation matrices.
  # chromosome: Chromosome variable name in `genedf`
  # gene_start: Gene start position variable name in `genedf`
  # symbol: Gene symbol variable name in `genedf`
  # pval: Gene-based test statistic P-value variable name in `genedf`
  # clump_p: Only genes with a P-value less than this threshold may index a locus for clumping
  # clump_kb: Kilobase size of the entire clumping window. Left and right windows from the index gene will be half the size of `clump_kb`
  # clump_r2: Only genes correlated with lead genes beyond this threshold may be clumped to other genes
  # verbose: TRUE if progress should be printed to the console, FALSE otherwise
  
  # load correlations
  if(toupper(ld_population)=='EUR') {data(EURGenTStatLD);gent_ld=EURGenTStatLD}
  if(toupper(ld_population)=='AFR') {data(AFRGenTStatLD);gent_ld=AFRGenTStatLD}
  if(toupper(ld_population)=='EAS') {data(EASGenTStatLD);gent_ld=EASGenTStatLD}
  if(toupper(ld_population)=='SAS') {data(SASGenTStatLD);gent_ld=SASGenTStatLD}
  if(toupper(ld_population)=='AMR') {data(AMRGenTStatLD);gent_ld=AMRGenTStatLD}
  genedf=genedf %>% rename(symbol=!!sym(symbol),gene_start=!!sym(gene_start),pval=!!sym(pval),chr=!!sym(chromosome))
  df=genedf %>%
    select(symbol,gene_start,pval,chr) %>%
    filter(pval<clump_p) %>%
    arrange(pval)
  if(nrow(df)==0) {if(verbose) {cat('no significant clumps; returning NA\n')};return(NA)}
  clump_kb=round(clump_kb/2)
  clumps=list()
  k=0
  # loop over each chromosome
  chrs=unique(df$chr)
  for(cc in 1:length(chrs)) {
    gent_ldchr=gent_ld[[paste0('chr',chrs[cc])]] %>% as.matrix()
    df_chr=df %>% filter(chr==chrs[cc])
    for(i in 1:nrow(df)) {
      allclumps=unlist(clumps);allclumps=c(names(clumps),allclumps)
      allclumps=unname(allclumps)
      if(df$symbol[i] %in% allclumps) next
      k=k+1
      genesaround=df %>%
        filter(abs(gene_start-df$gene_start[i])<(clump_kb*1e3)) %>%
        filter(symbol!=df$symbol[i])
      if(nrow(genesaround)==0) {
        clumps[[k]]=NA
        names(clumps)[k]=df$symbol[i]
        next
      }
      # attach LD with this SNP
      ix=which(rownames(gent_ldchr)==df$symbol[i])
      lddf=data.frame(symbol=rownames(gent_ldchr),r2=gent_ldchr[ix,]^2) %>%
        filter(symbol %in% genesaround$symbol) %>%
        filter(r2>clump_r2)
      if(nrow(lddf)==0) {
        # need to have nearby SNPs in LD to be clumps
        clumps[[k]]=NA
        names(clumps)[k]=df$symbol[i]
        next
      }
      # toadd=lddf %>% filter(!(symbol %in% allclumps)) %>% pull(symbol)
      # clumps[[k]]=ifelse(length(toadd)==0,NA,toadd)
      clumps[[k]]=lddf %>% filter(!(symbol %in% allclumps)) %>% pull(symbol)
      names(clumps)[k]=df$symbol[i]
    }
  }
  clumps=lapply(clumps,function(h) if(length(h)==0) NA else h)
  clumps
}



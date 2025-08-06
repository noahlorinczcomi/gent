
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# genome-wide GenT
# this is much more memory efficient than loading all genes into memory at once
gent_genomewide=function(gwas,KbWindow=50,ld_population='EUR',snp='rsid',chromosome='chr',position='position',effect_allele='effect_allele',z='z',beta=NULL,se=NULL,index=NULL,savefp=NULL,verbose=TRUE) {
  # expects:
  ## 1) gwas: GWAS summary statistics with at least rsid, chr, position, effect_allele, z
  ## 2) ld_population: population of LD. data will be list of LD matrices. each entry corresponds to one chromosome and is a list. each entry of each list entry is a gene-specific LD matrix of +/-1Mb 1kg v3 SNPs
  ## 3) KbWindow: Kb window size of SNPs to use
  ## 4) index: dataframe of (gene,chr,start,end) positions
  ## 5) savefp: directory to save results in (will also return them)
  if(is.null(z) & !is.null(beta) & !is.null(se)) gwas=gwas %>% mutate(z=!!sym(beta)/!!sym(se))
  if(is.null(index)) index=data(hg19genepos)  else index=fread('~/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt') # NEED TO ADD HERE
  setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices')
  ld_population=tail(unlist(strsplit(ld_population,'/')),1)
  if(toupper(ld_population)=='HIS') ldpop='AMR' else ldpop=toupper(ld_population)
  gwas=gwas %>% na.omit()
  gwas=gwas %>% rename(rsid=!!sym(snp), chr=!!sym(chromosome), position=!!sym(position), effect_allele=!!sym(effect_allele), z=!!sym(z))
  chrs=gwas %>% select(chr) %>% pull() %>% unique() %>% as.numeric() %>% na.omit() %>% sort()
  for(cc in 1:length(chrs)) {
    setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices')
    setwd(ldpop)
    setwd(paste0('chr',cc))
    if(verbose) cat('Starting chromosome ',chrs[cc],'\n')
    gwas_chr=gwas %>% filter(chr==chrs[cc])
    index_chr=index %>% filter(chr==cc)
    genes=unique(index_chr$symbol)
    pb=txtProgressBar(min=0,max=length(genes),style=3)
    for(i in 1:length(genes)) {
        if(verbose) setTxtProgressBar(pb, i)
        tryCatch(
            {
              ld=readRDS(paste0(genes[i],'.Rds'))
              ldrn=rownames(ld)
              spp=sapply(ldrn,\(.) unlist(strsplit(.,'_')))
              ld_df=do.call(rbind,spp) %>% as.data.frame() %>% rename(rsid=V1,a1=V2)
              ldrn=lapply(spp,\(.) .[1]
              starti=index_chr %>% filter(symbol==genes[i]) %>% select(start) %>% head(.,1) %>% pull()
              endi=index_chr %>% filter(symbol==genes[i]) %>% select(end) %>% head(.,1) %>% pull()
              starti=starti-KbWindow*1e3
              endi=endi+KbWindow*1e3
              gwas_chri=gwas_chr %>% filter(position>starti,position<endi)
              usesnps=intersect(ldrn,gwas_chri$rsid)
              gwas_chri=gwas_chri %>% filter(rsid %in% usesnps) %>% distinct(rsid,.keep_all=TRUE) %>% arrange(position)
              z=gwas_chri %>% left_join(ld_df,by='rsid') %>% mutate(z=ifelse(effect_allele==a1,z,-z)) %>% pull(z)
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
  if(!is.null(savefp)) saveRDS(rdf,savefp)
  rdf %>% as_tibble()
}


gent_manhattan=function(result) {
  chrdf=result %>%
    mutate(chr=as.numeric(chr)) %>%
    na.omit() %>%
    group_by(chr) %>%
    mutate(newbp=chr+gene_start/max(gene_start)) %>%
    ungroup() %>%
    mutate(chrcol=chr %in% seq(1,22,2)) %>%
    mutate(pval=ifelse(pval<1e-300,1e-300,pval)) %>%
    mutate(tp=-log10(pval))
  axisdf=chrdf %>% group_by(chr) %>% summarise(newbp=median(newbp)) %>% ungroup()
  ggplot(chrdf,aes(newbp,tp,color=chrcol)) +
    geom_point(alpha=2/3,show.legend=FALSE) +
    scale_color_manual(values=c('skyblue','cornflowerblue')) +
    geom_hline(yintercept=-log10(0.05/12727),lwd=1/2,linetype='dashed') +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
    labs(x='chromosome',y=expression('-log'[10]*'(GenT P-value)')) +
    scale_x_continuous(breaks=axisdf$newbp,labels=axisdf$chr)
}

gent_qq=function(result) {
    result %>% 
        mutate(pval=ifelse(pval<1e-300,1e-300,pval)) %>%
        arrange(pval) %>%
        mutate(pp=ppoints(n)) %>%
        ggplot(aes(-log10(pval),-log10(pp))) +
        geom_point(alpha=2/3,color='gold') +
        geom_abline(intercept=0,slope=1,lwd=1/2) +
        theme_bw() +
        theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
        labs(x=expression('expected -log'[10]*'(GenT P-value)'),y=expression('observed -log'[10]*'(GenT P-value)'))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MuGenT, MuGenT-PH, MuGenT-Pleio
library(data.table)
library(gent)
library(dplyr)
eur_gwas=fread('HOSPITALIZED_vs_POPULATION_EUR.tsv.gz') # European COVID GWAS
afr_gwas=fread('HOSPITALIZED_vs_POPULATION_AFR.tsv.gz') # African COVID GWAS
# add Z-scores
eur_gwas=eur_gwas %>% mutate(z=-all_inv_var_meta_beta/all_inv_var_meta_sebeta)
afr_gwas=afr_gwas %>% mutate(z=all_inv_var_meta_beta/all_inv_var_meta_sebeta)


# Genome-wide multi-population gene-based association testing (MuGenT)
covid_result=mugent_genomewide(
    gwas_list = list(EUR=eur_gwas %>% filter(`#CHR`==22), AFR=afr_gwas %>% filter(`#CHR`==22)),                                  # full GWAS summary statistics in each population
    # gwas_list = list(EUR=eur_gwas, AFR=afr_gwas),                                  # full GWAS summary statistics in each population
    ld_population_list = list(EUR='EUR', AFR='AFR'),                               # populations of the LD references to use (matching GWAS order)
    KbWindow = 50,                                                            # Kb window size of SNP-gene set assignment
    snp_list = list(EUR='rsid', AFR='rsid'),                                       # column names of unique SNP identifiers
    chromosome_list = list(EUR='#CHR', AFR='#CHR'),                                # column names of chromosomes
    position_list = list(EUR='POS', AFR='POS'),                                    # column names of SNP base pair positions (hg19)
    effect_allele_list = list(EUR='ALT', AFR='ALT'),                               # effect allele column names
    z_list = list(EUR='z', AFR='z'),                                                                 # Z-statistic column is not present in either GWAS so set NULL
    verbose = TRUE)         

## recall mugent() looks like this:
# results=mugent(t2d_z,LD_list)

library(gent,lib='~/isilon/Cheng-Noah/Rpkgs')
mugent_genomewide=function(
    gwas_list,
    ld_population_list,
    KbWindow=50,
    snp_list=lapply(1:length(gwas_list),'rsid'),
    chromosome_list=lapply(1:length(gwas_list),'chr'),
    position_list=lapply(1:length(gwas_list),'position'),
    effect_allele_list=lapply(1:length(gwas_list),'effect_allele'),
    z_list=lapply(1:length(gwas_list),'z'),
    # must use Z-stats
    mugentpleio_alpha=0.05/12727,index=NULL,verbose=TRUE) {
  if(is.null(index)) index=data.table::fread('~/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt') # NEED TO ADD HERE
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
                setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices/')
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














## this is less memory efficient because it loads the entire chromosome's gene LD matrices into memory at once
gent_genomewide=function(gwas,KbWindow=50,ld_population='EUR',snp='rsid',chromosome='chr',position='position',effect_allele='effect_allele',z='z',beta=NULL,se=NULL,index=NULL,savefp=NULL,verbose=TRUE) {
  ## now uses entire chromosome-specific files which are lists of LD matrices
  # expects:
  ## 1) gwas: GWAS summary statistics with at least rsid, chr, position, effect_allele, z
  ## 2) ld_population: population of LD. data will be list of LD matrices. each entry corresponds to one chromosome and is a list. each entry of each list entry is a gene-specific LD matrix of +/-1Mb 1kg v3 SNPs
  ## 3) KbWindow: Kb window size of SNPs to use
  ## 4) index: dataframe of (gene,chr,start,end) positions
  ## 5) savefp: directory to save results in (will also return them)
  if(is.null(z) & !is.null(beta) & !is.null(se)) gwas=gwas %>% mutate(z=!!sym(beta)/!!sym(se))
  if(is.null(index)) index=fread('~/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt') # NEED TO ADD HERE
  setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices')
  ld_population=tail(unlist(strsplit(ld_population,'/')),1)
  if(toupper(ld_population)=='HIS') ldpop='AMR' else ldpop=toupper(ld_population)
  gwas=gwas %>% na.omit()
  gwas=gwas %>% rename(rsid=!!sym(snp), chr=!!sym(chromosome), position=!!sym(position), z=!!sym(z))
  chrs=gwas %>% select(chr) %>% pull() %>% unique() %>% as.numeric() %>% na.omit() %>% sort()
  rdf=data.frame()
  for(cc in 1:length(chrs)) {
    setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices')
    setwd(ldpop)
    setwd(paste0('chr',cc))
    if(verbose) cat('Starting chromosome ',chrs[cc],'\n')
    gwas_chr=gwas %>% filter(chr==chrs[cc])
    index_chr=index %>% filter(chr==cc)
    # load chromosome-specific LD matrices for EUR
    setwd('~/isilon/Cheng-Noah/reference_data/ld_matrices/chr_specific_files/EUR')
    ld_chr=readRDS(paste0('chr',cc,'.Rds'))
    # genes=unique(index_chr$symbol)
    genes=names(ld_chr)
    pb=txtProgressBar(min=0,max=length(genes),style=3)
    for(i in 1:length(genes)) {
        if(verbose) setTxtProgressBar(pb, i)
        tryCatch(
            {
              ld=as.matrix(ld_chr[[i]])
            #   ld=readRDS(paste0(genes[i],'.Rds'))
              ldrn=rownames(ld)
              ldrn=sapply(ldrn,\(.) unlist(strsplit(.,'_'))[1])
              starti=index_chr %>% filter(symbol==genes[i]) %>% select(start) %>% head(.,1) %>% pull()
              endi=index_chr %>% filter(symbol==genes[i]) %>% select(end) %>% head(.,1) %>% pull()
              starti=starti-KbWindow*1e3
              endi=endi+KbWindow*1e3
              gwas_chri=gwas_chr %>% filter(position>starti,position<endi)
              usesnps=intersect(ldrn,gwas_chri$rsid)
              gwas_chri=gwas_chri %>% filter(rsid %in% usesnps) %>% distinct(rsid,.keep_all=TRUE) %>% arrange(position)
              z=gwas_chri$z
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
  if(!is.null(savefp)) saveRDS(rdf,savefp)
  rdf %>% as_tibble()
}

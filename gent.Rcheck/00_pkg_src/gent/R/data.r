#' Gene-based test statistic correlations (European 1000 Genomes Phase 3)
#'
#' Each entry in this list is a chromosome-specific matrix of correlations between GenT test statistics under the gene-based null hypothesis (-+/50Kb SNP-gene assignment).
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs in gene-specific sets. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by within-chromosome base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 European reference panel
#' @examples
#' data(EURGenTStatLD)
#' str(EURGenTStatLD)
"EURGenTStatLD"


#' Gene-based test statistic correlations (African 1000 Genomes Phase 3)
#'
#' Each entry in this list is a chromosome-specific matrix of correlations between GenT test statistics under the gene-based null hypothesis (-+/50Kb SNP-gene assignment).
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs in gene-specific sets. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by within-chromosome base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 African reference panel
#' @examples
#' data(AFRGenTStatLD)
#' str(AFRGenTStatLD)
"AFRGenTStatLD"


#' Gene-based test statistic correlations (East Asian 1000 Genomes Phase 3)
#'
#' Each entry in this list is a chromosome-specific matrix of correlations between GenT test statistics under the gene-based null hypothesis (-+/50Kb SNP-gene assignment).
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs in gene-specific sets. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by within-chromosome base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 East Asian reference panel
#' @examples
#' data(EASGenTStatLD)
#' str(EASGenTStatLD)
"EASGenTStatLD"


#' Gene-based test statistic correlations (South Asian 1000 Genomes Phase 3)
#'
#' Each entry in this list is a chromosome-specific matrix of correlations between GenT test statistics under the gene-based null hypothesis (-+/50Kb SNP-gene assignment).
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs in gene-specific sets. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by within-chromosome base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 South Asian reference panel
#' @examples
#' data(SASGenTStatLD)
#' str(SASGenTStatLD)
"SASGenTStatLD"


#' Gene-based test statistic correlations (Admixed American 1000 Genomes Phase 3)
#'
#' Each entry in this list is a chromosome-specific matrix of correlations between GenT test statistics under the gene-based null hypothesis (-+/50Kb SNP-gene assignment).
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs in gene-specific sets. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by within-chromosome base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 Admixed American reference panel
#' @examples
#' data(AMRGenTStatLD)
#' str(AMRGenTStatLD)
"AMRGenTStatLD"


#' Ensembl hg19 coordinates of gene start and end base pair positions
#'
#' Each row in this data set is a unique pair of HGNC gene symbols and Ensembl IDs. Their chromosomes, start, end, and base pair midpoint positions are provided.
#'
#' @format Each row in this data set is a unique pair of HGNC gene symbols and Ensembl IDs. Their chromosomes, start, end, and base pair midpoint positions are provided
#' @source Ensembl hg19 coordinates downloaded from the Biomart web tool.
#' @examples
#' data(EnsemblHg19GenePos)
#' str(EnsemblHg19GenePos)
"EnsemblHg19GenePos"

#' Ensembl hg38 coordinates of gene start and end base pair positions
#'
#' Each row in this data set is a unique pair of HGNC gene symbols and Ensembl IDs. Their chromosomes, start, end, and base pair midpoint positions are provided.
#'
#' @format Each row in this data set is a unique pair of HGNC gene symbols and Ensembl IDs. Their chromosomes, start, end, and base pair midpoint positions are provided
#' @source Ensembl hg38 coordinates downloaded from the Biomart web tool.
#' @examples
#' data(EnsemblHg38GenePos)
#' str(EnsemblHg38GenePos)
"EnsemblHg38GenePos"

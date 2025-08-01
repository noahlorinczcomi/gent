% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.r
\name{mugent}
\alias{mugent}
\title{Multi-ancestry gene-based association test (MuGenT)}
\usage{
mugent(Z, ldlist)
}
\arguments{
\item{Z}{matrix of Z-statistics from GWAS. Rows are SNPs and columns are populations.}

\item{ldlist}{list of population-specific LD matrices whose ordering corresponds to the column ordering of \code{Z}}
}
\value{
\itemize {
\item \code{pval}: P-value for testing H0: gene is not associated with trait.
\item \code{shape}: shape parameter of null (Gamma) distribution.
\item \code{rate}: rate parameter of null (Gamma) distribution.
\item \code{mu_h0}: expectation of null (Gamma) distribution.
\item \code{sigma2_h0}: variance of null (Gamma) distribution.
}
}
\description{
This function performs a gene-based association test using multiple populations.
}
\examples{
# Example for 5 SNPs and 2 populations
ldlist=rWishart(2,100,diag(5))
ldlist=lapply(1:2,function(h) cov2cor(ldlist[,,h]))
Z=t(mvnfast::rmvn(2,rep(0,5),diag(5)))
mugent(Z,ldlist)
}

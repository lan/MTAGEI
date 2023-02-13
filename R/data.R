#' Example1: individual-level data
#'
#' An example of calculating summary statistics in MAGENTA given the individual-level data set.
#' @docType data
#' @usage data(rawdata)
#' @format A list with 2 sublists:
#' \describe{
#'   \item{geno.dat}{a matrix of genotype data, each row denotes subject \eqn{i} and column represents SNP \eqn{j} in the gene, the first column is the unique subject IDs}
#'   \item{dat}{a matrix of traits, covariates and environmental/treatment variable data, the first column is the unique subject IDs.
#' There are 3 traits, 86 variants (both rare and common). Specifically, there are 1000 subjects with three trait measurements. The environmental variable is treatment.}}
"rawdata"

#' Example2: the score summary statistics and their variance for continuous traits
#'
#' An example of calculating summary statistics in MAGENTA given the score summary statistics U and its variance V.
#' @docType data
#' @format A list with 3 sublists:
#' \describe{
#'   \item{sumstats}{a list of score summary statistics for all the SNPs in the gene for each trait in each environmental group}
#'   \item{common.sumstats}{a list of score summary statistics for a large number of independent common variants for each trait in each environmental group}
#'   \item{MAF}{a vector of minor allele frequency for each SNP in the gene.}
#'   \item{R.C}{a SNP correlation matrix for all the SNPs in the gene}
#'   \item{KA}{a \eqn{K \times K} matrix of genetic correlation matrix for \eqn{K} traits}}
"sumstats.dat"


#' Example2: the score summary statistics and their variance for binary traits
#'
#' An example of calculating summary statistics in MAGENTA given the score summary statistics U and its variance V.
#' @docType data
#' @format A list with 3 sublists:
#' \describe{
#'   \item{sumstats}{a list of score summary statistics for all the SNPs in the gene for each trait in each environmental group}
#'   \item{common.sumstats}{a list of score summary statistics for a large number of independent common variants for each trait in each environmental group}
#'   \item{MAF}{a vector of minor allele frequency for each SNP in the gene.}
#'   \item{R.C}{a SNP correlation matrix for all the SNPs in the gene}
#'   \item{KA}{a \eqn{K \times K} matrix of genetic correlation matrix for \eqn{K} traits}
#'   \item{sumstats.SPA}{a list of SPA-adjusted score summary statistics for all the SNPs in the gene for each trait in each environmental group}}
"sumstats.dat.bin"

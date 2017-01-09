#' Beta values of selected CpG sites
#'
#' A matrix of beta values corresponding to 151 probes and 78 samples.
#'
#' @format A matrix with 151 rows and 78 columns
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51954}
"sunexp_data"

#' Genotype of the simulated individuals
#'
#' A data frame whose rows are markers and columns are individuals. The marker 101 corresponds to the 'all' variable of the 'indiv' data frame, it is the parental allele linked with the phenotype. From this marker 101, genotype is propagated (upstream and downstream) with recombination factor of 5%. 
#'
#' @format A data frame with 78 rows and 4 variables:
#' \describe{
#'   \item{sex}{Character string describing sex of the individual (male or female).}
#'   \item{age}{Character string describing the age of thne individual.}
#'   \item{histo}{Character string describing the histological type of the sample (dermis or epidermis).}
#'   \item{expo}{Character string describing the sun exposure (protected or exposed).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51954}
"sunexp_design"

#' Description of selected CpG sites
#'
#' A data frame describing CpG positions.
#'
#' @format A data frame with 151 rows and 2 variables:
#' \describe{
#'   \item{CHR}{Character string corresponding to chromosome (hg19).}
#'   \item{MAPINFO}{Integer corresponding to position on the chromosome (hg19).}
#'   \item{UCSC_CpG_Islands_Name}{Character string describing CpG island position (hg19).}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51954}
"sunexp_platform"



#' Description of selected genes
#'
#' A data frame describing selected genes.
#'
#' @format A data frame with 9 rows and 9 variables:
#' \describe{
#'   \item{chrom}{Character string corresponding to chromosome (hg19).}
#'   \item{chromStart}{Integer corresponding to lowest bound of the gene (hg19).}
#'   \item{chromEnd}{Integer corresponding  to upperest bound of the gene (hg19).}
#'   \item{name}{Character string corresponding to the name of the gene.}
#'   \item{score}{Integer corresponding to nothing (mandatory in bed format).}
#'   \item{strand}{Character corresponding to the strand of the gene.}
#' }
"genes"
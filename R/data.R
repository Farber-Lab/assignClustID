#’ DO mice genotypes
#'
#' Dataset containing genotypes of selected DO mice
#'
#' @format A data.frame with 112913 rows and 9 variables
#' \describe{
#'   \item{marker}{GigaMUGA markers}
#'   \item{chr}{chromosome}
#'   \item{position}{position of the marker on chromosome: mm10 reference}
#'   \item{Xnum}{X followed by ID of a genotyped DO mouse}
#' }
#' @name DO_mm10
#' @source <sftp://rivanna.hpc.virginia.edu/sfs/qumulo/qproject/farber_lab/Users/Current_users/Kristyna_K/annotations/DO_mice>
"DO_mm10"


#’ Comparison table to match dinucleotide pairs
#'
#' Matrix containing all possible combinations of what could be dinucleotide
#' matches
#'
#' @format A matrix with 16 rows and 16 variables
#' \describe{
#'   \item{dinucleotide}
#' }
#' @name comparisonTable
"comparisonTable"



#’ Comparison table to match dinucleotide pairs
#'
#' Matrix containing dinucleotide matches excluding reverse complements
#'
#' @format A matrix with 16 rows and 16 variables
#' \describe{
#'   \item{dinucleotide}
#' }
#' @name comparisonTableNoRevComp
"comparisonTableNoRevComp"

#' Simulated dataset and counts
#'
#' I created the `RNAseqCovarImpute_data` dataset to show how to use the package functions. 
#' The exact code used to generate these data are found in the Example_Data_for_RNAseqCovarImpute
#' vignette. In short, `example_data` contains 500 rows with data for variables x, y, and z, which
#' are continuous normally distributed, and a and b, which are binary variables. Missigness
#' was simulated for all variables other than x such that a complete case analysis would drop
#' 24.2% of participants. `example_DGE` contains random count data from the Poisson distribution
#' for 500 made up genes, ENS1-ENS500
#'
#' \describe{
#' \item{example_data}{random covariate data for 5 variables with 24.2% missingness}
#' \item{example_DGE}{random count data for 500 made up genes in a DGElist}
#' ...
#' }
#' @return Tibble with 500 rows of data for variables x, y, and z
#' @docType data
#' @keywords datasets
#' @name RNAseqCovarImpute_data
#' @usage data(RNAseqCovarImpute_data)
#' @format Two objects: A data frame with 500 rows and 5 variables and A DGElist with 500 genes and 500 samples
NULL

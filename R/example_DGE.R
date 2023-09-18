#' Simulated counts in DGE list
#'
#' The exact code used to generate these data are found in the Example_Data_for_RNAseqCovarImpute
#' vignette. In short, `example_data` contains 500 rows with data for variables x, y, and z, which
#' are continuous normally distributed, and a and b, which are binary variables. Missigness
#' was simulated for all variables other than x such that a complete case analysis would drop
#' 24.2% of participants. `example_DGE` contains random count data from the Poisson distribution
#' for 500 made up genes, ENS1-ENS500
#' 
#' @format ## `example_DGE`
#' A DGElist with 500 genes and 500 samples
#' 
#' @return DGElist for 500 made up genes, ENS1-ENS500
#' @docType data
#' @keywords datasets
#' @name example_DGE
#' @usage data(example_DGE)
#' @examples
#' data(example_DGE)
"example_DGE"

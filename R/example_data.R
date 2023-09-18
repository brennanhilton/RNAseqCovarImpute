#' Simulated dataset
#'
#' The exact code used to generate these data are found in the Example_Data_for_RNAseqCovarImpute
#' vignette. In short, `example_data` contains 500 rows with data for variables x, y, and z, which
#' are continuous normally distributed, and a and b, which are binary variables. Missigness
#' was simulated for all variables other than x such that a complete case analysis would drop
#' 24.2% of participants. `example_DGE` contains random count data from the Poisson distribution
#' for 500 made up genes, ENS1-ENS500
#'
#' @format ## `example_data`
#' data frame with 500 rows and 5 variables
#' \describe{
#' \item{x}{continuous normally distributed}
#' \item{y}{continuous normally distributed}
#' \item{z}{continuous normally distributed}
#' \item{a}{binary}
#' \item{b}{binary}
#' ...
#' }
#' @return Tibble with 500 rows of data for variables x, y, and z
#' @docType data
#' @keywords datasets
#' @name example_data
#' @usage data(example_data)
#' @examples
#' data(example_data)
"example_data"

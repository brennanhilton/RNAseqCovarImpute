#' limmavoom_imputed_data_list
#'
#' Loops through the imputed data list (output from "impute_by_gene_bin" function)
#' and runs limma-voom RNA seq analysis.
#' @return A dataframe with coefficient, standard error, sigma, and residual degrees of freedom values from limma-voom gene expression analysis. One row per gene and one set of values per imputed dataset.
#' @param gene_intervals Output from get_gene_bin_intervals function. A dataframe where each row contains the start (first col) and end (second col) values for each gene bin interval.
#' @param DGE A DGEList object.
#' @param imputed_data_list Output from impute_by_gene_bin.
#' @param m Number of imputed data sets.
#' @param voom_formula Formula for design matrix.
#' @param predictor Independent variable of interest. Must be a variable in voom_formula
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select rename as_tibble all_of bind_cols as_tibble
#' @importFrom foreach %do% foreach
#' @importFrom edgeR cpm
#' @importFrom mice complete
#' @importFrom stats model.matrix as.formula
#' @importFrom limma lmFit normalizeBetweenArrays
#' @importFrom rlang .data
#' @importFrom dplyr left_join
#'
#' @examples
#' data(RNAseqCovarImpute_data)
#' intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
#' gene_bin_impute <- impute_by_gene_bin(example_data,
#'     intervals,
#'     example_DGE,
#'     m = 2,
#'     param = SerialParam()
#' )
#' coef_se <- limmavoom_imputed_data_list(
#'     gene_intervals = intervals,
#'     DGE = example_DGE,
#'     imputed_data_list = gene_bin_impute,
#'     m = 2,
#'     voom_formula = "~x + y + z + a + b",
#'     predictor = "x",
#'     param = SerialParam()
#' )
#'
#' final_res <- combine_rubins(
#'     DGE = example_DGE,
#'     model_results = coef_se,
#'     voom_formula = "~x + y + z + a + b"
#' )
#' @export


limmavoom_imputed_data_list <- function(gene_intervals, DGE, imputed_data_list, m, voom_formula, predictor, param = bpparam()) {
  # Validity tests
  if (!class(DGE) %in% "DGEList") {
    stop("Input 'DGE' is not a valid DGEList object.")
  }
  if (!(class(m) %in% c("numeric"))) {
    stop("Input 'm' must be numeric.")
  }
  if (!class(as.formula(voom_formula))%in% c("formula")) {
    stop()
  }
  if (!(class(predictor) %in% c("character"))) {
    stop("Input 'predictor' must be a character")
  }
  # Get mean-variance curve from all genes across all M imputations
  sx_sy <- lowess_all_gene_bins(gene_intervals, DGE, imputed_data_list, m, voom_formula, predictor)
  
  # Parallelize the loop using bplapply
  all_coefs_se <- bplapply(seq(length(imputed_data_list)),
                           gene_intervals = gene_intervals,
                           m = m,
                           imputed_data = imputed_data_list, 
                           DGE = DGE, 
                           voom_formula = voom_formula, 
                           predictor = predictor, 
                           sx_sy = sx_sy,
                           FUN = limmavoom_imputed_data_list_helper, 
                           BPPARAM = param)
  
  all_coefs_se <- do.call(rbind, all_coefs_se)
  return(all_coefs_se)
}


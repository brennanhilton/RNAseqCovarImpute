#' lowess_all_gene_bins
#'
#' Loops through all bins and all M imputations, prepares DGE and design to run
#' voom_sx_sy, which fits gene-wise linear models and extracts log count size (sx)
#' and sqrt resudual standard deviations (sy) to make the lowess curve
#' @return All sx and sy values for lowess function across all M imputation.
#' 
#' @examples
#' data(example_data)
#' data(example_DGE)
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
#' 
#' @export
#' @keywords internal
lowess_all_gene_bins <- function(gene_intervals, DGE, imputed_data_list, m, voom_formula, predictor) {
  lowess_all_gene_bins <- foreach(i = seq(m), .combine = "rbind") %do% {
    
    foreach(gene_bin = seq(length(imputed_data_list)), .combine = "rbind") %do% {
      # Get imputed data for this gene bin and this M imputed dataset
      imputed_data <- complete(imputed_data_list[[gene_bin]],i)
      # Get DGE list for this gene interval
      alldg_bin <- DGE[as.numeric(gene_intervals[gene_bin, 1]):as.numeric(gene_intervals[gene_bin, 2]), ]
      design1 <- model.matrix(as.formula(voom_formula), imputed_data)
      # voom_sx_sy is the first half of the voom function to get sx and sy for lowess curve
      out <- voom_sx_sy(alldg_bin, design1, lib.size.all = DGE$samples$lib.size*DGE$samples$norm.factors)
    }}}

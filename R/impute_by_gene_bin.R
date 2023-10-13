#' impute_by_gene_bin
#'
#' Loops through DGE list using the gene bin intervals from the
#' "get_gene_bin_intervals" function and makes imputed datasets. For instance,
#' if n = 100 and intervals contains 200 gene bin intervals, output will be a list of
#' 200 sets of 100 imputed datasets. Each of the 200 sets are imputed
#' using only the genes in one gene bin.
#' @return A list of sets of n imputed datasets, one per gene bin.
#' @param data Sample data with one row per sample. Sample row order should match the col order in the DGEList.
#' @param intervals Output from get_gene_bin_intervals function. A dataframe where each row contains the start (first col) and end (second col) values for each gene bin interval.
#' @param DGE A DGEList object.
#' @param m Number of imputed data sets.
#' @param maxit Used by mice function.
#' @param BPPARAM A BiocParallelParam object
#'
#' @include impute_by_gene_bin_helper.R
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate bind_cols as_tibble
#' @importFrom foreach %do% foreach
#' @importFrom edgeR cpm
#' @importFrom mice mice quickpred
#' @importFrom rlang .data
#'
#' @examples
#' data(example_data)
#' data(example_DGE)
#' intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
#' gene_bin_impute <- impute_by_gene_bin(example_data,
#'     intervals,
#'     example_DGE,
#'     m = 2
#' )
#' coef_se <- limmavoom_imputed_data_list(
#'     gene_intervals = intervals,
#'     DGE = example_DGE,
#'     imputed_data_list = gene_bin_impute,
#'     m = 2,
#'     voom_formula = "~x + y + z + a + b"
#' )
#'
#' final_res <- combine_rubins(
#'     DGE = example_DGE,
#'     model_results = coef_se,
#'     predictor = "x"
#' )
#' @export


# Define the function to perform parallel imputation using bplapply
impute_by_gene_bin <- function(data, intervals, DGE, m, maxit = 10, BPPARAM = bpparam()) {
    # Validity tests
    if (!class(DGE) %in% "DGEList") {
        stop("Input 'DGE' is not a valid DGEList object.")
    }
    if (!any((class(data) %in% c("tbl_df", "tbl", "data.frame")))) {
        stop("Input 'data' is not a valid data.frame, tbl, or tbl_df object.")
    }
    if (!(class(m) %in% c("numeric"))) {
        stop("Input 'm' must be numeric.")
    }
    if (!(class(maxit) %in% c("numeric"))) {
        stop("Input 'maxit' must be numeric.")
    }
    # get the counts per million for all genes in DGE
    cpm_all <- cpm(DGE, log = TRUE, prior.count = 5)

    gene_bin_impute <- bplapply(seq(nrow(intervals)),
        intervals = intervals,
        cpm_all = cpm_all,
        data = data,
        m = m,
        maxit = maxit,
        FUN = impute_gene_bin_helper,
        BPPARAM = BPPARAM
    )

    return(gene_bin_impute)
}

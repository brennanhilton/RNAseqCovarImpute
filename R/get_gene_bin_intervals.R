#' get_gene_bin_intervals
#'
#' Creates gene bins. Input DGE list, sample data,
#' and 'n' number of individuals per genes. By default, number of bins and
#' genes per bin are set so that each bin has approximately 1 gene per
#' 10 individuals in the data.
#' @return Data frame with one row per gene bin. Columns indicate the start and end positions and the number of genes of each bin.
#' @param DGE A DGEList object.
#' @param data Sample data with one row per sample. Sample row order should match the col order in the DGEList.
#' @param n Genes per bin are set so that each bin has approximately 1 gene per n individuals in the data.
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

get_gene_bin_intervals <- function(DGE, data, n = 10) {
    # Validity tests
    if (!class(DGE) %in% "DGEList") {
        stop("Input 'DGE' is not a valid DGEList object.")
    }
    if (!any((class(data) %in% c("tbl_df", "tbl", "data.frame")))) {
        stop("Input 'data' is not a valid data.frame, tbl, or tbl_df object.")
    }
    if (!(class(n) %in% c("numeric"))) {
        stop("Input 'n' must be numeric.")
    }

    gene_number <- as.numeric(nrow(mat))
    genes_per_bin <- floor(nrow(ampute_i)/n)
    num_bins <- floor(gene_number/genes_per_bin)
    remainder = gene_number - (genes_per_bin*num_bins)
    
    # Initialize the bins with the lowest possible genes_per_bin
    bins <- rep(genes_per_bin, num_bins)
    
    # Distribute the remainder evenly across the bins
    i <- 1
    while (remainder > 0) {
        bins[i] <- bins[i] + 1
        remainder <- remainder - 1
        i <- i %% num_bins + 1 #bump i to the next bin, can also go back to bin 1
    }
  
    # Calculate the start and end genes for each bin
    start_genes <- c(1, cumsum(bins[-num_bins]) + 1)
    end_genes <- cumsum(bins)
    output <- data.frame(start = start_genes, end = end_genes, number = bins)
}

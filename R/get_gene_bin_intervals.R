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
#' intervals <- get_gene_bin_intervals(DGE, data, n = 10)
#' gene_bin_impute <- impute_by_gene_bin(data,
#'     intervals,
#'     DGE,
#'     n = 3
#' )
#' coef_se <- limmavoom_imputed_data_list(
#'     gene_intervals = intervals,
#'     DGE = DGE,
#'     imputed_data_list = gene_bin_impute,
#'     n = 3,
#'     voom_formula = "~x + y + z + a + b",
#'     predictor = "x"
#' )
#'
#' final_res <- combine_rubins(
#'     DGE = DGE,
#'     model_results = coef_se,
#'     voom_formula = "~x + y + z + a + b"
#' )
#' @export

get_gene_bin_intervals <- function(DGE, data, n = 10) {
    gene_number <- as.numeric(nrow(DGE))
    genes_per_bin <- floor(nrow(data) / n) # Recommended 10 to 1 individuals to genes ratio so default n is 10
    bins <- floor(gene_number / genes_per_bin)
    gene_remainder <- gene_number %% bins
    # When gene_number is not divisible by genes_per_bin, intervals1 will have all the larger gene bins,
    # and intervals2 will have all the smaller gene bins. e.g. if you have 11 genes, 5 genes per bin,
    # and 2 bins, intevals1 will have the first bin of 6 and intervals2 will have the last bin of 5.
    # When number of genes per bin evenly divides into gene_number, intervals1 will be empty
    intervals1 <- tibble(
        start = seq(from = 1, by = genes_per_bin + 1, length.out = gene_remainder),
        end = seq(from = genes_per_bin + 1, by = genes_per_bin + 1, length.out = gene_remainder)
    )

    intervals2 <- tibble(
        start = seq(from = 1 + gene_remainder * (genes_per_bin + 1), by = genes_per_bin, length.out = bins - gene_remainder),
        end = seq(from = genes_per_bin + gene_remainder * (genes_per_bin + 1), by = genes_per_bin, length.out = bins - gene_remainder)
    )

    intervals <- rbind(intervals1, intervals2) %>% mutate(number = 1 + end - start)
}

#' impute_by_gene_bin_helper
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
#' @param param Arguments passed to BiocParallel::bpparam()
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate bind_cols as_tibble
#' @importFrom foreach %do% foreach
#' @importFrom edgeR cpm
#' @importFrom mice mice quickpred
#' @importFrom rlang .data
#'
#' @keywords internal
# Define a function for imputation that will be called by bplapply
impute_gene_bin_helper <- function(i, intervals, cpm_all, data, m, maxit) {
    # Extract cpm for the current bin of genes based on the interval
    cpm_bin <- cpm_all[as.numeric(intervals[i, 1]):as.numeric(intervals[i, 2]), ] %>%
        t() %>%
        as_tibble()

    # Add the bin of genes to the covariate data
    data_mice <- data %>% bind_cols(cpm_bin)

    # Impute the data using mice
    imputed_data <- mice(data_mice, m = m, maxit = maxit, predictorMatrix = quickpred(data_mice))

    return(imputed_data)
}

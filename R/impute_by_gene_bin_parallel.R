#' impute_by_gene_bin_parallel
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
#' @param n Number of imputed data sets.
#' @param cores Number of cores to run in parallel using 'PSOCK' back-end implemented with doParallel.
#'
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom edgeR cpm
#' @importFrom mice mice
#' @importFrom mice quickpred
#' @importFrom dplyr bind_cols
#' @importFrom dplyr as_tibble
#' @importFrom rlang .data
#'
#' @examples
#' intervals <- get_gene_bin_intervals(DGE, data, n = 10)
#' gene_bin_impute <- impute_by_gene_bin_parallel(data,
#'     intervals,
#'     DGE,
#'     n = 3,
#'     cores = 2
#' )
#' coef_se <- limmavoom_imputed_data_list_parallel(
#'     gene_intervals = intervals,
#'     DGE = DGE,
#'     imputed_data_list = gene_bin_impute,
#'     n = 3,
#'     voom_formula = "~x + y + z + a + b",
#'     predictor = "x",
#'     cores = 2
#' )
#'
#' final_res <- combine_rubins(
#'     DGE = DGE,
#'     model_results = coef_se,
#'     voom_formula = "~x + y + z + a + b"
#' )
#' @export

impute_by_gene_bin_parallel <- function(data, intervals, DGE, n, cores) {
    # get the counts per million for all genes in DGE
    cpm_all <- cpm(DGE, log = TRUE, prior.count = 5)
    myCluster <- makeCluster(cores, # number of cores to use
        type = "PSOCK"
    ) # type of cluster
    registerDoParallel(myCluster)
    gene_bin_impute <- foreach(i = seq(nrow(intervals))) %dopar% {
        # extract cpm for just current bin of genes based on the interval
        cpm_bin <- cpm_all[as.numeric(intervals[i, 1]):as.numeric(intervals[i, 2]), ] %>%
          t() %>%
          as_tibble()
        # add the bin of genes to the covaraite data
        data_mice <- data %>% bind_cols(cpm_bin)
        imputed_data <- mice(data_mice, m = n, maxit = 10, seed = 2022, predictorMatrix = quickpred(data_mice))
    }
}

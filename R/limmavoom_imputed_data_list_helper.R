#' limmavoom_imputed_data_list_helper
#'
#' Loops through the imputed data list (output from "impute_by_gene_bin" function)
#' and runs limma-voom RNA seq analysis.
#' @return A dataframe with coefficient, standard error, sigma, and residual degrees of freedom values from limma-voom gene expression analysis. One row per gene and one set of values per imputed dataset.
#' @param gene_intervals Output from get_gene_bin_intervals function. A dataframe where each row contains the start (first col) and end (second col) values for each gene bin interval.
#' @param DGE A DGEList object.
#' @param imputed_data_list Output from impute_by_gene_bin.
#' @param m Number of imputed data sets.
#' @param voom_formula Formula for design matrix.
#' @param predictor Independent variable of interest. Must be a variable in voom_formula.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate bind_cols as_tibble left_join
#' @importFrom foreach %do% foreach
#' @importFrom edgeR cpm
#' @importFrom mice mice quickpred
#' @importFrom stats model.matrix
#' @importFrom rlang .data
#' @importFrom limma normalizeBetweenArrays lmFit
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
#' @keywords internal
limmavoom_imputed_data_list_helper <- function(gene_bin, gene_intervals, DGE, imputed_data_list, m, voom_formula, predictor, sx_sy) {
  # get mean-variance curve from all genes across all M imputations
  
        # get imputed data
        imputed_data <- imputed_data_list[[gene_bin]]
        # get dge list for this gene interval
        alldg_bin <- DGE[as.numeric(gene_intervals[gene_bin, 1]):as.numeric(gene_intervals[gene_bin, 2]), ]

        all_coef_se_within_bin <- foreach(i = seq(m), .combine = "left_join") %do% {
            # we have imputed data for a particular gene bin interval.
            # now we get the ith imputed data set within
            data_i <- complete(imputed_data, i)
          
            # run limmavoom
            design1 <- model.matrix(as.formula(voom_formula), data_i)
          
            voom1 <- voom_master_lowess(alldg_bin, design1, lib.size.all = DGE$samples$lib.size*DGE$samples$norm.factors, sx = sx_sy$sx, sy = sx_sy$sy)
            fit1 <- lmFit(voom1)
          
            # get coefficients unscaled SE, df residual, and sigma from fit1
            coef <- fit1$coefficients %>%
                as_tibble() %>%
                dplyr::select(all_of(starts_with(predictor))) %>%
                dplyr::rename(coef = all_of(starts_with(predictor)))
          
            SE_unscaled <- fit1$stdev.unscaled * fit1$sigma
            SE_unscaled <- as_tibble(SE_unscaled) %>%
                dplyr::select(all_of(starts_with(predictor))) %>%
                dplyr::rename(SE_unscaled = all_of(starts_with(predictor)))
          
            degrees_freedom_residual <- fit1$df.residual
          
            sigma <- fit1$sigma
          
            output1 <- coef %>%
                cbind(SE_unscaled) %>%
                mutate(probe = rownames(fit1),
                       sigma = sigma,
                       df_residual = degrees_freedom_residual)
            # rename fit values to include info on which imputed data they come from
            colnames(output1)[colnames(output1)!= "probe"] <- paste0(colnames(output1)[colnames(output1)!= "probe"],".",i)
            output1
        }
        return(all_coef_se_within_bin)
}

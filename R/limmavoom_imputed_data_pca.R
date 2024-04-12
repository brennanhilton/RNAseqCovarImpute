#' limmavoom_imputed_data_pca
#'
#' Combines results from each imputed dataset using Rubin's rules.
#' @return Dataframe with one row per gene. Columns contain coefficients, standard errors, and p-values from the limma-voom pipeline.
#' @param imp Imputed data from mice (mids object)
#' @param DGE A DGEList object.
#' @param voom_formula Formula for design matrix.
#' @param BPPARAM A BiocParallelParam object
#' 
#' @include limmavoom_imputed_data_pca_helper.R
#' @include combine_rubins_pca.R
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select everything
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom BiocGenerics Reduce
#' 
#' @examples
#' data(example_data)
#' data(example_DGE)
#' pca_data = limma::voom(example_DGE)$E
#' p = PCAtools::pca(pca_data)
#' pcs = p$rotated[,1:78]
#' example_data = cbind(example_data, pcs)
#' imp = mice::mice(example_data, m=3)
#' mi_pca_res = limmavoom_imputed_data_pca(
#'    imp = imp,
#'    DGE = example_DGE,
#'    voom_formula = "~x + y + z + a + b"
#'    )
#' @export
limmavoom_imputed_data_pca <- function(imp, DGE, voom_formula, BPPARAM = bpparam()) {
  results <- BiocParallel::bplapply(seq(1, imp$m),
                      FUN = limmavoom_imputed_data_pca_helper,
                      BPPARAM = BPPARAM,
                      imp=imp,
                      DGE=DGE,
                      voom_formula= voom_formula)
  
  results_final <- BiocParallel::bplapply(seq_len(nrow(DGE)),
                            BPPARAM = BPPARAM,
                            FUN = function(gene) {
                              BiocParallel::bplapply(seq_along(results),
                                       BPPARAM = BiocParallel::bpparam("SerialParam"),
                                       FUN = function(imp_data) {
                                         res <- rbind(results[[imp_data]]$coef[gene, ], results[[imp_data]]$SE[gene, ])
                                         rownames(res) <- c("coef", "SE")
                                         res
                                         })
                              })
  
  results_final <- lapply(results_final, function(x) {
    do.call(rbind, x)
  })
  
  dfs <- lapply(results, function(inner_list) inner_list$df)
  avg_df = BiocGenerics::Reduce(`+`, dfs) / length(dfs)
  
  final_res_pca <- 
    do.call(rbind,
            BiocParallel::bplapply(
              seq_along(results_final),
              FUN = combine_rubins_pca,
              results_final = results_final,
              avg_df = avg_df,
              m = imp$m,
              BPPARAM = BPPARAM))
  final_res_pca <- as.data.frame(final_res_pca) %>% 
    mutate(probe = rownames(DGE)) %>% 
    dplyr::select(probe, dplyr::everything())
}

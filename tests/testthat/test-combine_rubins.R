test_that("consistent results", {
    data(example_data)
    data(example_DGE)
    library(BiocParallel)
    intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
    gene_bin_impute <- impute_by_gene_bin(example_data,
        intervals,
        example_DGE,
        m = 2,
        BPPARAM = SerialParam(RNGseed = 2023)
    )
    coef_se <- limmavoom_imputed_data_list(
        gene_intervals = intervals,
        DGE = example_DGE,
        imputed_data_list = gene_bin_impute,
        m = 2,
        voom_formula = "~x + y + z + a + b",
        predictor = "x",
        BPPARAM = SerialParam(RNGseed = 2023)
    )
    final_res <- combine_rubins(
        DGE = example_DGE,
        model_results = coef_se,
        voom_formula = "~x + y + z + a + b"
    )
    expect_equal(final_res$coef_combined[1], -0.02099867, tolerance = 0.001)
})

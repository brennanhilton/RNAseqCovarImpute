test_that("impute_by_gene_bin outputs mids", {
  data(example_data)
  data(example_DGE)
  intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
  gene_bin_impute <- impute_by_gene_bin(example_data,
     intervals,
     example_DGE,
     m = 2,
     param = SerialParam()
  )
  expect_s3_class(gene_bin_impute[[1]], "mids")
})
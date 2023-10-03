test_that("Intervals correct type", {
    data(example_data)
    data(example_DGE)
    intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
    expect_type(intervals, "list")
})

test_that("Intervals includes all rows in data", {
    data(example_data)
    data(example_DGE)
    intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
    expect_equal(nrow(example_data), sum(intervals$number))
})

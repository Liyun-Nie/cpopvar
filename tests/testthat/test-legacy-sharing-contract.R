test_that("legacy sharing analysis retains its output contract", {
  sharing <- data.frame(
    gene = paste0("gene", seq_len(8)),
    num_species = c(1, 1, 1, 2, 2, 3, 4, 4)
  )

  result <- perform_sharing_statistical_test(
    hotspot_analysis_results = list(
      all_gene_sharing = sharing,
      parameters = list(total_species = 4)
    ),
    config = list(),
    task_params = list(alpha_level = 0.05)
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("tests", "comparison_table", "sharing_summary", "report", "test_metadata")
  )
  expect_true(is.data.frame(result$comparison_table))
  expect_identical(result$test_metadata$test_method, "stepwise_two_speed")
  expect_named(
    result$tests,
    c(
      "private_genes",
      "core_genes",
      "two_speed_model",
      "hypergeometric",
      "bipolarity_index"
    )
  )
})

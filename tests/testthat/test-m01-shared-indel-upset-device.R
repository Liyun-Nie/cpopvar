m01_shared_indel_matrix <- function() {
  mat <- matrix(
    c(
      1, 1, 1, 0,
      1, 1, 0, 0,
      1, 0, 1, 0,
      0, 1, 1, 1,
      1, 0, 0, 1,
      0, 1, 0, 1,
      1, 1, 1, 1,
      0, 0, 1, 1
    ),
    nrow = 8,
    ncol = 4,
    byrow = TRUE
  )
  colnames(mat) <- c("sA", "sB", "sC", "sD")
  rownames(mat) <- paste0("indel", seq_len(nrow(mat)))
  mat
}

m01_upset_config <- function() {
  yaml::read_yaml(system.file(
    "extdata",
    "synthetic",
    "config.yml",
    package = "cpopvar"
  ))
}

expect_valid_m01_upset_pdf <- function(path) {
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)
  if (requireNamespace("pdftools", quietly = TRUE)) {
    info <- pdftools::pdf_info(path)
    expect_gt(info[["pages"]], 0)
    expect_equal(info[["pages"]], 1)
    if (isTRUE(capabilities("cairo"))) {
      text_n <- nchar(trimws(paste(pdftools::pdf_text(path), collapse = " ")))
      expect_gt(text_n, 0)
    }
  }
}

test_that("M01 shared INDEL UpSet writes one nonblank intended PDF without Rplots.pdf", {
  skip_if_not_installed("UpSetR")
  skip_if_not_installed("withr")

  workdir <- withr::local_tempdir()
  withr::local_dir(workdir)

  intended <- file.path(workdir, "M01_shared_indel_upset.pdf")
  result <- cpopvar:::generate_shared_indel_upset_plot(
    upset_matrix = m01_shared_indel_matrix(),
    output_file = intended,
    config = m01_upset_config()
  )

  expect_identical(result, intended)
  expect_false(file.exists(file.path(workdir, "Rplots.pdf")))
  expect_valid_m01_upset_pdf(intended)
  expect_equal(unname(grDevices::dev.cur()), 1)
})

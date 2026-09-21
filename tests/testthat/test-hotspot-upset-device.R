hotspot_upset_candidates <- function(feature_col = "gene") {
  data <- data.frame(
    species = c("sA", "sB", "sA", "sC", "sB", "sC", "sA", "sB"),
    stringsAsFactors = FALSE
  )
  data[[feature_col]] <- c("g1", "g1", "g2", "g2", "g3", "g3", "g4", "g4")
  data
}

hotspot_upset_config <- function() {
  yaml::read_yaml(system.file(
    "extdata",
    "synthetic",
    "config.yml",
    package = "cpopvar"
  ))
}

expect_valid_upset_pdf <- function(path) {
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

test_that("UpSet helper writes the intended PDF and does not create Rplots.pdf", {
  skip_if_not_installed("UpSetR")
  skip_if_not_installed("withr")

  workdir <- withr::local_tempdir()
  withr::local_dir(workdir)

  intended <- file.path(workdir, "M03_candidate_upset.pdf")
  result <- generate_hotspot_upset_plot(
    candidate_hotspots = hotspot_upset_candidates("gene"),
    config = hotspot_upset_config(),
    task_params = list(),
    feature_col = "gene",
    output_path = intended
  )

  expect_false(is.null(result))
  expect_false(is.null(result$plot))
  expect_true(is.data.frame(result$intersection_data))
  expect_gt(nrow(result$intersection_data), 0)
  expect_false(file.exists(file.path(workdir, "Rplots.pdf")))
  expect_valid_upset_pdf(intended)
  expect_equal(unname(grDevices::dev.cur()), 1)
})

test_that("M04 poiGS UpSet helper uses the M04 filename contract without Rplots.pdf", {
  skip_if_not_installed("UpSetR")
  skip_if_not_installed("withr")

  workdir <- withr::local_tempdir()
  withr::local_dir(workdir)

  intended <- file.path(workdir, "M04_candidate_upset.pdf")
  result <- generate_hotspot_upset_plot(
    candidate_hotspots = hotspot_upset_candidates("poiGS_ID"),
    config = hotspot_upset_config(),
    task_params = list(),
    feature_col = "poiGS_ID",
    output_path = intended
  )

  expect_false(is.null(result))
  expect_false(is.null(result$plot))
  expect_true("poiGS_count" %in% names(result$intersection_data) ||
                "poiGS" %in% names(result$intersection_data))
  expect_false(file.exists(file.path(workdir, "Rplots.pdf")))
  expect_valid_upset_pdf(intended)
  expect_equal(unname(grDevices::dev.cur()), 1)
})

test_that("UpSet helper preserves a caller-owned graphics device", {
  skip_if_not_installed("UpSetR")
  skip_if_not_installed("withr")

  workdir <- withr::local_tempdir()
  withr::local_dir(workdir)

  caller_pdf <- file.path(workdir, "caller_device.pdf")
  cpopvar:::open_cpopvar_pdf(caller_pdf, width = 4, height = 3)
  graphics::plot(1:3, main = "CALLER_DEVICE")
  caller_dev <- grDevices::dev.cur()
  expect_gt(caller_dev, 1)

  intended <- file.path(workdir, "M03_candidate_upset.pdf")
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && caller_dev %in% open_devices) {
      grDevices::dev.off(which = caller_dev)
    }
  }, add = TRUE)

  result <- generate_hotspot_upset_plot(
    candidate_hotspots = hotspot_upset_candidates("gene"),
    config = hotspot_upset_config(),
    task_params = list(),
    output_path = intended
  )

  expect_false(is.null(result$plot))
  expect_equal(grDevices::dev.cur(), caller_dev)
  expect_true(caller_dev %in% grDevices::dev.list())
  expect_false(file.exists(file.path(workdir, "Rplots.pdf")))
  expect_valid_upset_pdf(intended)

  grDevices::dev.off(which = caller_dev)
  expect_true(file.exists(caller_pdf))
  expect_gt(file.info(caller_pdf)$size, 0)
  if (requireNamespace("pdftools", quietly = TRUE) && isTRUE(capabilities("cairo"))) {
    caller_text <- paste(pdftools::pdf_text(caller_pdf), collapse = " ")
    expect_match(caller_text, "CALLER_DEVICE", fixed = TRUE)
  }
})

test_that("M03 upset-only caller writes intended files without Rplots.pdf", {
  skip_if_not_installed("UpSetR")
  skip_if_not_installed("withr")

  workdir <- withr::local_tempdir()
  withr::local_dir(workdir)
  output_dir <- file.path(workdir, "m03_out")
  dir.create(output_dir)

  normalized <- data.frame(
    species = rep(c("sA", "sB", "sC"), each = 6),
    gene = rep(letters[1:6], 3),
    frequency_per_kb = rep(c(0, 1, 2, 8, 9, 10), 3),
    region_type = "CDS",
    var_type = "snp",
    stringsAsFactors = FALSE
  )

  result <- generate_M03_hotspot_plots(
    normalized_data = normalized,
    config = hotspot_upset_config(),
    output_dir = output_dir,
    task_name = "M03_upset_device",
    task_params = list(
      plot_types = "upset",
      frequency_percentile = 0.75
    )
  )

  expect_true(isTRUE(result$success))
  expect_true(file.exists(file.path(output_dir, "M03_candidate_upset.pdf")))
  expect_true(file.exists(file.path(output_dir, "M03_upset_intersections.csv")))
  expect_false(file.exists(file.path(workdir, "Rplots.pdf")))
  expect_false(file.exists(file.path(output_dir, "Rplots.pdf")))
  expect_valid_upset_pdf(file.path(output_dir, "M03_candidate_upset.pdf"))
  expect_equal(unname(grDevices::dev.cur()), 1)
})

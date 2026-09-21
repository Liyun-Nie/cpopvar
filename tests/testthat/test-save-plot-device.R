test_that("save_plot preserves PDF text with the cross-platform device", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pdftools")
  skip_if_not(capabilities("cairo"), "Cairo PDF support is unavailable")

  probe_text <- "CPOPVAR_TEXT_PROBE_A1B2"
  plot <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 1:3),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = probe_text,
      x = "Probe x axis",
      y = "Probe y axis"
    )

  base_path <- tempfile("cpopvar_pdf_probe_")
  config <- list(
    visualization_settings = list(
      output = list(default_width = 6, default_height = 4)
    )
  )

  saved_path <- save_plot(plot, base_path, config)

  expect_identical(saved_path, paste0(base_path, ".pdf"))
  expect_true(file.exists(saved_path))

  pdf_text <- paste(pdftools::pdf_text(saved_path), collapse = " ")
  expect_match(pdf_text, probe_text, fixed = TRUE)
  expect_match(pdf_text, "Probe x axis", fixed = TRUE)
  expect_match(pdf_text, "Probe y axis", fixed = TRUE)
})

test_that("base and grid callers can open the same PDF device", {
  skip_if_not(capabilities("cairo"), "Cairo PDF support is unavailable")

  output_file <- tempfile(fileext = ".pdf")
  cpopvar:::open_cpopvar_pdf(output_file, width = 4, height = 3)
  graphics::plot(1:3, main = "Base PDF probe")
  grDevices::dev.off()

  expect_true(file.exists(output_file))
  expect_gt(file.info(output_file)$size, 0)
})

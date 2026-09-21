test_that("conversion preflight names every missing raster dependency", {
  check <- cpopvar:::check_plot_conversion_dependencies(
    export_formats = c("png", "jpeg"),
    package_available = function(package) FALSE
  )

  expect_false(check$ok)
  expect_setequal(check$required, c("magick", "pdftools"))
  expect_setequal(check$missing, c("magick", "pdftools"))
  expect_match(check$message, "magick", fixed = TRUE)
  expect_match(check$message, "pdftools", fixed = TRUE)
  expect_match(check$message, "install.packages", fixed = TRUE)
  expect_match(check$message, "no raster conversion was attempted", fixed = TRUE)
})

test_that("conversion preflight does not require raster packages when unused", {
  package_check_called <- FALSE
  check <- cpopvar:::check_plot_conversion_dependencies(
    export_formats = character(0),
    package_available = function(package) {
      package_check_called <<- TRUE
      FALSE
    }
  )

  expect_true(check$ok)
  expect_length(check$required, 0)
  expect_false(package_check_called)
})

test_that("PDF to PNG conversion succeeds when optional dependencies exist", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("magick")
  skip_if_not_installed("pdftools")
  skip_if_not(capabilities("cairo"), "Cairo PDF support is unavailable")

  plots_dir <- tempfile("cpopvar_conversion_")
  dir.create(plots_dir)

  plot <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 1:3),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Conversion probe")

  config <- list(
    visualization_settings = list(
      output = list(
        export_formats = "png",
        png_dpi = 72,
        png_bg = "white"
      )
    )
  )

  pdf_path <- save_plot(plot, file.path(plots_dir, "probe"), config)
  result <- convert_plots_for_reporting(plots_dir, config)

  expect_true(file.exists(pdf_path))
  expect_true(result$success)
  expect_identical(result$total_failed, 0)
  expect_identical(result$total_converted, 1)
  expect_true(file.exists(file.path(plots_dir, "probe.png")))
})

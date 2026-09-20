test_that("smart grid keeps all ggpubr comparison plots under ggplot2 4", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("patchwork")

  set.seed(1)
  make_plot <- function(n_levels) {
    d <- data.frame(
      grp = factor(rep(letters[seq_len(n_levels)], each = 20)),
      y = abs(stats::rnorm(20 * n_levels))
    )
    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = .data$grp, y = .data$y, fill = .data$grp, color = .data$grp)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::annotate(
        "text",
        x = (n_levels + 1) / 2,
        y = max(d$y) * 1.12,
        label = "p = 0.01",
        size = 3.5
      )
    if (n_levels > 2) {
      comps <- utils::combn(levels(d$grp), 2, simplify = FALSE)
      p <- p + ggpubr::stat_compare_means(
        comparisons = comps,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = FALSE
      )
    }
    p
  }

  plots <- list(
    region_type = make_plot(3),
    genome_region = make_plot(3),
    phylogeny = make_plot(2),
    life_form = make_plot(2),
    status = make_plot(2)
  )

  arranged <- arrange_plots_smart_grid(
    plot_list = plots,
    max_cols = 3,
    auto_arrange = TRUE,
    common_title = "Single-Factor Analysis",
    title_size = 14
  )

  expect_s3_class(arranged, "patchwork")
  expect_gte(length(arranged), 5L)

  outfile <- tempfile(fileext = ".png")
  ggplot2::ggsave(outfile, arranged, width = 10, height = 8, dpi = 72)
  expect_true(file.exists(outfile))
  expect_gt(file.info(outfile)$size, 20000)
})

test_that("single-factor overall p-value uses annotate so ggplot2 4 PDFs are not blank", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")

  set.seed(1)
  d <- data.frame(
    grp = factor(rep(c("CDS", "IGS", "intron"), each = 40)),
    y = abs(stats::rnorm(120))
  )
  comps <- utils::combn(levels(d$grp), 2, simplify = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$grp, y = .data$y, fill = .data$grp)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::annotate("text", x = 2, y = max(d$y) * 1.12, label = "p < 2e-16", size = 3.5) +
    ggpubr::stat_compare_means(
      comparisons = comps,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    )

  outfile <- tempfile(fileext = ".png")
  ggplot2::ggsave(outfile, p, width = 6, height = 4, dpi = 100)
  expect_gt(file.info(outfile)$size, 8000)
})

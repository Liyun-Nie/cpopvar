#' Smart Grid Layout for Multiple Plots
#'
#' Automatically arranges multiple plots in an intelligent grid layout
#' with configurable maximum columns. Supports n x max_cols arrangement.
#'
#' @param plot_list List of ggplot objects to arrange
#' @param max_cols Maximum number of columns (default: 3)
#' @param auto_arrange Whether to automatically arrange plots (default: TRUE)
#' @param common_title Optional common title for the combined plot
#' @param title_size Size of the common title (default: 16)
#' @return Combined plot object or list of plots if auto_arrange = FALSE
#' arrange_plots_smart_grid
#' @export
arrange_plots_smart_grid <- function(plot_list,
                                   max_cols = 3,
                                   auto_arrange = TRUE,
                                   common_title = NULL,
                                   title_size = 16) {

  # Validate inputs
  if (length(plot_list) == 0) {
    log_message("Warning: Empty plot list provided to arrange_plots_smart_grid", level = "warning")
    return(ggplot2::ggplot() + ggplot2::ggtitle("No Plots Available"))
  }

  # Filter out NULL plots
  valid_plots <- plot_list[!sapply(plot_list, is.null)]
  n_plots <- length(valid_plots)

  if (n_plots == 0) {
    log_message("Warning: No valid plots found after filtering NULLs", level = "warning")
    return(ggplot2::ggplot() + ggplot2::ggtitle("No Valid Plots Available"))
  }

  log_message(sprintf("Arranging %d plots in smart grid (max_cols = %d)", n_plots, max_cols))

  # If auto_arrange is FALSE, just return the list
  if (!auto_arrange) {
    return(valid_plots)
  }

  # Calculate optimal grid dimensions
  n_cols <- min(n_plots, max_cols)
  n_rows <- ceiling(n_plots / n_cols)

  log_message(sprintf("Grid layout: %d rows x %d columns", n_rows, n_cols))

  # ggplot2 4.x cannot convert ggpubr::stat_compare_means(comparisons = ...)
  # layers through ggplotGrob()/gridExtra/cowplot. Try backends independently
  # so a grob conversion failure does not collapse a multi-plot grid to one panel.
  try_backend <- function(name, builder) {
    if (!requireNamespace(name, quietly = TRUE)) {
      return(NULL)
    }
    result <- tryCatch(
      builder(),
      error = function(e) {
        log_message(
          sprintf("%s arrangement failed: %s", name, e$message),
          level = "warning"
        )
        NULL
      }
    )
    if (!is.null(result)) {
      log_message(sprintf("[SUCCESS] Successfully arranged plots using %s", name))
    }
    result
  }

  combined_plot <- try_backend("patchwork", function() {
    arranged <- patchwork::wrap_plots(valid_plots, nrow = n_rows, ncol = n_cols)
    if (!is.null(common_title)) {
      arranged <- arranged +
        patchwork::plot_annotation(
          title = common_title,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
              size = title_size,
              hjust = 0.5,
              face = "bold"
            )
          )
        )
    }
    arranged
  })
  if (!is.null(combined_plot)) {
    return(combined_plot)
  }

  combined_plot <- try_backend("cowplot", function() {
    arranged <- cowplot::plot_grid(
      plotlist = valid_plots,
      nrow = n_rows,
      ncol = n_cols,
      align = "hv"
    )
    if (!is.null(common_title)) {
      title_plot <- cowplot::ggdraw() +
        cowplot::draw_label(common_title, fontface = "bold", size = title_size)
      arranged <- cowplot::plot_grid(
        title_plot,
        arranged,
        ncol = 1,
        rel_heights = c(0.1, 1)
      )
    }
    arranged
  })
  if (!is.null(combined_plot)) {
    return(combined_plot)
  }

  combined_plot <- try_backend("gridExtra", function() {
    if (!is.null(common_title)) {
      title_grob <- grid::textGrob(
        common_title,
        gp = grid::gpar(fontsize = title_size, fontface = "bold")
      )
      return(gridExtra::grid.arrange(
        grobs = valid_plots,
        nrow = n_rows,
        ncol = n_cols,
        top = title_grob
      ))
    }
    gridExtra::grid.arrange(
      grobs = valid_plots,
      nrow = n_rows,
      ncol = n_cols
    )
  })
  if (!is.null(combined_plot)) {
    return(combined_plot)
  }

  log_message(
    "No plot arrangement backend succeeded; returning the first plot only",
    level = "warning"
  )
  fallback_plot <- valid_plots[[1]]
  if (!is.null(common_title)) {
    fallback_plot <- fallback_plot +
      ggplot2::ggtitle(paste(common_title, "(Showing 1 of", n_plots, "plots)"))
  }
  fallback_plot
}

#' Calculate Optimal Grid Dimensions
#'
#' Helper function to calculate the optimal number of rows and columns
#' for arranging a given number of plots with a maximum column constraint.
#'
#' @param n_plots Number of plots to arrange
#' @param max_cols Maximum number of columns allowed
#' @return List with nrow and ncol values
calculate_grid_dimensions <- function(n_plots, max_cols = 3) {
  if (n_plots <= 0) {
    return(list(nrow = 1, ncol = 1))
  }

  n_cols <- min(n_plots, max_cols)
  n_rows <- ceiling(n_plots / n_cols)

  return(list(nrow = n_rows, ncol = n_cols))
}

#' Create Plot Grid Summary
#'
#' Creates a summary information plot showing the arrangement details
#'
#' @param plot_names Names of the plots being arranged
#' @param grid_dims Grid dimensions (from calculate_grid_dimensions)
#' @return Summary ggplot object
create_arrangement_summary <- function(plot_names, grid_dims) {
  summary_text <- paste(
    sprintf("Plot Arrangement Summary"),
    sprintf("Total plots: %d", length(plot_names)),
    sprintf("Grid: %d rows x %d columns", grid_dims$nrow, grid_dims$ncol),
    sprintf("Plots: %s", paste(head(plot_names, 6), collapse = ", ")),
    if (length(plot_names) > 6) "..." else "",
    sep = "\n"
  )

  summary_plot <- ggplot2::ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = summary_text,
             hjust = 0.5, vjust = 0.5, size = 4) +
    theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = "gray90"),
      plot.margin = margin(20, 20, 20, 20)
    )

  return(summary_plot)
}

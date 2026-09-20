#' Get Internal Visual Standards for All Plotting Systems
#'
#' Provides a centralized, internal list of aesthetic defaults,
#' structured by plotting system (e.g., ggplot2, pheatmap).
#' This avoids scattering hardcoded values across the codebase
#' and keeps the user-facing config file simple.
#'
#' @return A nested list of default aesthetic values.
#' get_visual_standards
#' @export
get_visual_standards <- function() {
  list(
    # Defaults for the ggplot2 System
    ggplot2 = list(
      point_size = 1.5, point_alpha = 0.8, point_outlier_shape = NA,
      jitter_width = 0.2,
      violin_alpha = 0.7, violin_scale = "width",
      boxplot_width = 0.1, boxplot_alpha = 1.0,
      bar_fill = "skyblue", bar_color = "white", bar_linewidth = 0.2,
      line_size = 1.0, line_ref_linetype = "dashed", line_ref_alpha = 0.5,
      density_alpha = 0.8, density_linewidth = 1.2,
      text_size = 2, text_color = "black"
    ),

    # Defaults for the pheatmap Package
    pheatmap = list(
      border_color = "grey60",
      color = colorRampPalette(c("#2166AC", "#FFFFBF", "#D73027"))(100),
      fontsize = 10, fontsize_row = 8,
      cellwidth = NA, # Let pheatmap decide by default
      cellheight = NA # Let pheatmap decide by default
    )
  )
}
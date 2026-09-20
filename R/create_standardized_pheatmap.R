#' Create standardized pheatmap visualization
#' 
#' @title Create standardized pheatmap visualization
#' @description Creates a standardized heatmap using pheatmap with consistent styling
#' @param data_matrix Numeric matrix to visualize as heatmap
#' @param title Plot title (optional)
#' @param color_palette Color palette for the heatmap
#' @param cluster_rows Whether to cluster rows (default: TRUE)
#' @param cluster_cols Whether to cluster columns (default: TRUE)
#' @param display_numbers Whether to display values in cells (default: TRUE)
#' @param breaks Numeric vector specifying color breaks (optional)
#' @param ... Additional arguments passed to pheatmap
#' @return pheatmap object
#' @export
create_standardized_pheatmap <- function(data_matrix, 
                                       title = NULL,
                                       color_palette,  # No default - must be provided via dependency injection
                                       cluster_rows = TRUE,
                                       cluster_cols = TRUE,
                                       display_numbers = TRUE,
                                       breaks = NULL,
                                       ...) {
  
  # Load required package for heatmap generation
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for heatmap generation but not available")
  }
  
  # Set default breaks if not provided
  if (is.null(breaks)) {
    data_range <- range(data_matrix, na.rm = TRUE)
    if (all(data_range >= -1 & data_range <= 1)) {
      # Correlation-like data - use symmetric breaks around 0
      breaks <- seq(-1, 1, length.out = 101)
    } else {
      # General data
      breaks <- seq(data_range[1], data_range[2], length.out = 101)
    }
  }
  
  # Create color palette with improved mapping for correlation data
  if (all(range(data_matrix, na.rm = TRUE) >= -1 & range(data_matrix, na.rm = TRUE) <= 1)) {
    # For correlation matrices, ensure proper color mapping
    colors <- grDevices::colorRampPalette(color_palette)(length(breaks) - 1)
  } else {
    colors <- grDevices::colorRampPalette(color_palette)(100)
  }
  
  # Create pheatmap with standardized parameters using namespace prefix
  pheatmap::pheatmap(
    data_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method = "complete",
    color = colors,
    breaks = breaks,
    border_color = "white",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 9,
    main = title,
    angle_col = 45,
    display_numbers = display_numbers,
    number_color = "black",
    fontsize_number = if(is.character(display_numbers)) 6 else 7,
    silent = TRUE,
    ...
  )
}

############################################################
#### Visualization Utilities - Common Helper Functions ####
############################################################
#
# Provides shared visualization helper functions for all modules
# Including color schemes, themes, and common data processing
#
############################################################

# Load required libraries

# get_color_palette function removed - using version from get_color_palette.R

#' Get unified ggplot theme
#' 
#' @param base_size Base font size
#' @param title_size Title font size
#' @return ggplot theme object
#' @importFrom ggplot2 theme_minimal theme element_text element_rect element_line
#' @export
get_main_theme <- function(base_size = 12, title_size = 14) {
  
  ggplot2::theme_minimal(base_size = base_size) +
  ggplot2::theme(
      # Title styling
      plot.title =ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle =ggplot2::element_text(size = base_size - 1, hjust = 0.5),
      
      # Axis titles and text
      axis.title =ggplot2::element_text(size = base_size, face = "bold"),
      axis.text =ggplot2::element_text(size = base_size - 1),
      
      # Legend
      legend.title =ggplot2::element_text(size = base_size, face = "bold"),
      legend.text =ggplot2::element_text(size = base_size - 1),
      legend.position = "bottom",
      
      # Panel
      panel.background =ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = "grey90", size = 0.5),
      panel.grid.minor = ggplot2::element_line(color = "grey95", size = 0.25),
      panel.border =ggplot2::element_rect(color = "grey70", fill = NA, size = 0.5),
      
      # Facet labels
      strip.text =ggplot2::element_text(size = base_size, face = "bold"),
      strip.background =ggplot2::element_rect(fill = "grey95", color = "grey70")
    )
}

# detect_outliers function removed - using version from detect_outliers.R

#' Generate standardized filename
#' @param module_name Module name (e.g. "M01", "M02", "M03")
#' @param plot_type Plot type (e.g. "boxplot", "heatmap", "density")
#' @param factors Vector of factor names
#' @param extension File extension
#' @return Standardized filename
#' @export
generate_filename <- function(module_name, plot_type, factors = NULL, extension = "png") {
  
  # Base filename
  filename_parts <- c(module_name)
  
  # Add factor information
  if (!is.null(factors) && length(factors) > 0) {
    filename_parts <- c(filename_parts, factors)
  }
  
  # Add plot type
  filename_parts <- c(filename_parts, plot_type)
  
  # Combine and clean filename
  filename <- paste(filename_parts, collapse = "_")
  filename <- gsub("[^A-Za-z0-9_.-]", "_", filename)  # Clean special characters
  filename <- paste0(filename, ".", extension)
  
  return(filename)
}

# safe_save_plot function removed - using version from viz_utils.R

#' Validate required columns in data frame
#' 
#' @param data Data frame
#' @param required_cols Vector of required column names
#' @return Logical value, TRUE if all required columns exist
validate_required_columns <- function(data, required_cols) {
  
  if (is.null(data) || !is.data.frame(data)) {
    return(FALSE)
  }
  
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    warning("Missing required columns: ", paste(missing_cols, collapse = ", "))
    return(FALSE)
  }
  
  return(TRUE)
}

#' Create data summary information
#' 
#' @param data Data frame
#' @return String containing data summary
create_data_summary <- function(data) {
  
  if (is.null(data) || !is.data.frame(data)) {
    return("No data available")
  }
  
  summary_info <- sprintf(
    "Data Summary: %d rows, %d columns\nColumns: %s",
    nrow(data),
    ncol(data),
    paste(names(data), collapse = ", ")
  )
  
  return(summary_info)
}

# Export all functions for use by other modules
# message("✓ Visualization utilities loaded successfully")
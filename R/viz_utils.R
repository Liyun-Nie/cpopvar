############################################################
#### visualization_utils.R - Unified Plot Output System ####
############################################################
#
# Unified high-quality plot output system for cpopvar.
# Provides a single, configuration-driven save_plot() function
# that handles all plot types (ggplot, pheatmap, etc.) with
# consistent DPI and sizing standards.
#
# Key Features:
# - Configuration-driven: Reads standards from config.yml
# - Type-aware: Auto-detects ggplot vs pheatmap objects
# - High-quality output: Enforces 300 DPI standard
# - Error isolation: Individual plot failures don't crash pipeline
# - Consistent interface: Single function for all plot saving
#
# Plot-dimension helpers:
# - get_plot_dimensions() reads width/height/dpi from config
# - hardcoded dimension values are avoided in favour of the config path
#
############################################################

# Load required libraries

# %||% operator defined in config_yaml_utils.R

#' Get plot dimensions from configuration
#' 
#' Helper to extract plot dimensions from config
#' 
#' @param config List. Global configuration object
#' @param width Numeric. Override width (optional)
#' @param height Numeric. Override height (optional)
#' @param dpi Numeric. Override dpi (optional)
#' @return List with width, height, and dpi
#' @export
get_plot_dimensions <- function(config, width = NULL, height = NULL, dpi = NULL) {
  # Read dimensions from the output configuration path
  viz_config <- config$visualization_settings$output %||% 
                config$visualization_settings$theme_and_sizing %||% 
                list()
  
  list(
    width = width %||% viz_config$default_width %||% 12,
    height = height %||% viz_config$default_height %||% 8,
    dpi = dpi %||% viz_config$default_dpi %||% 300
  )
}

#' Unified vectorized plot saving function
#' 
#' This is the single entry point for all plot saving in the system.
#' It automatically detects plot type and saves as high-quality vector PDF.
#' 
#' @param plot_object Plot object (ggplot, pheatmap, or other)
#' @param base_path Character. File path without extension (e.g., "/path/to/plot")
#' @param config List. Global configuration object containing visualization settings
#' @param width Numeric. Override width in inches (optional)
#' @param height Numeric. Override height in inches (optional)
#' @return Character. Path to saved file, or NULL if failed
#' 
#' @examples
#' if (require("ggplot2")) {
#'   # Create a dummy plot and config for the example
#'   my_ggplot <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
#'                ggplot2::geom_point()
#'   config <- list(visualization_settings = list(theme_and_sizing = list(
#'     default_width = 8, default_height = 6
#'   )))
#'
#'   # Define a temporary path for the example output
#'   temp_plot_path <- tempfile(fileext = ".pdf")
#'
#'   # Run the save_plot function
#'   save_plot(my_ggplot, sub("\\.pdf$", "", temp_plot_path), config)
#'
#'   # Clean up the created file
#'   if (file.exists(temp_plot_path)) {
#'     print(paste("Example plot saved to:", temp_plot_path))
#'     unlink(temp_plot_path)
#'   }
#' }
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices pdf dev.off dev.cur
#' @importFrom grid grid.draw
#' @importFrom tools file_path_sans_ext file_ext
#' @export
save_plot <- function(plot_object, base_path, config, 
                      width = NULL, height = NULL) {
  
  # Extract configuration standards
  viz_config <- config$visualization_settings$output %||% 
                config$visualization_settings$theme_and_sizing %||% 
                list()
  
  # Get standard parameters from config
  default_width <- width %||% viz_config$default_width %||% 12
  default_height <- height %||% viz_config$default_height %||% 8
  default_dpi <- viz_config$default_dpi %||% 300
  
  # Construct vectorized PDF file path
  file_path <- paste0(base_path, ".pdf")
  
  # Ensure output directory exists
  output_dir <- dirname(file_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  tryCatch({
    
    # Detect plot type and handle appropriately
    plot_class <- class(plot_object)
    
    if ("ggplot" %in% plot_class) {
      # Handle ggplot objects - vectorized PDF output
      ggplot2::ggsave(
        filename = file_path,
        plot = plot_object,
        width = default_width,
        height = default_height,
        device = "pdf",
        units = "in"
      )
      
    } else if ("pheatmap" %in% plot_class) {
      # Handle pheatmap objects - vectorized PDF output
      grDevices::pdf(
        file = file_path,
        width = default_width,
        height = default_height
      )
      grid::grid.draw(plot_object)
      grDevices::dev.off()
      
    } else if (any(c("recordedplot", "grob", "gtable") %in% plot_class)) {
      # Handle other grid/base graphics objects - vectorized PDF output
      grDevices::pdf(
        file = file_path,
        width = default_width,
        height = default_height
      )
      grid::grid.draw(plot_object)
      grDevices::dev.off()
      
    } else {
      # Fallback: try to draw the object using grid - vectorized PDF output
      grDevices::pdf(
        file = file_path,
        width = default_width,
        height = default_height
      )
      grid::grid.draw(plot_object)
      grDevices::dev.off()
    }
    
    # Verify file was created
    if (file.exists(file_path)) {
      # Log successful save
      file_size <- round(file.size(file_path) / 1024 / 1024, 2)
      cat(sprintf("[SUCCESS] Vectorized plot saved: %s (%.2f MB, %dx%d inches)\n", 
                  basename(file_path), file_size, default_width, default_height))
      return(file_path)
    } else {
      warning("Plot file was not created: ", file_path)
      return(NULL)
    }
    
  }, error = function(e) {
    # Ensure graphics device is closed on error
    if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    
    warning(sprintf("Failed to save plot %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}

#' Legacy wrapper for backward compatibility
#' 
#' This function maintains compatibility with existing safe_save_plot calls
#' while redirecting to the new vectorized save_plot() function.
#' 
#' @param plot ggplot object
#' @param filename Character. Full filename with extension
#' @param output_dir Character. Output directory path
#' @param width Numeric. Width in inches
#' @param height Numeric. Height in inches
#' @param dpi Numeric. DPI resolution (ignored for PDF output)
#' @param config List. Configuration object (added for new system)
#' @return Character. Path to saved file or NULL if failed
safe_save_plot <- function(plot, filename, output_dir, width = 10, height = 8, dpi = 300, config = NULL) {
  
  # Construct base path (remove extension from filename)
  base_filename <- tools::file_path_sans_ext(filename)
  base_path <- file.path(output_dir, base_filename)
  
  # Extract format from filename
  format <- tools::file_ext(filename)
  if (format == "") format <- "png"
  
  # Use new vectorized function if config is available
  if (!is.null(config)) {
    return(save_plot(plot, base_path, config, width = width, height = height))
  } else {
    # Fallback to old logic for backward compatibility
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    file_path <- file.path(output_dir, filename)
    
    tryCatch({
      ggplot2::ggsave(filename = file_path, plot = plot, 
                      width = width, height = height, dpi = dpi)
      
      cat("Plot saved:", file_path, "\n")
      return(file_path)
      
    }, error = function(e) {
      warning("Failed to save plot:", filename, " - Error:", e$message)
      return(NULL)
    })
  }
}

# Export notification
# cat("[SUCCESS] Vectorized visualization output system loaded successfully\n")
# cat("   Main function: save_plot(plot_object, base_path, config)\n") 
# cat("   Features: Auto-detection, config-driven, vector PDF output\n")
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

#' Cross-platform PDF graphics device
#'
#' Prefer Cairo when it can actually open. Some hosts (notably GitHub Actions
#' macOS) report `capabilities("cairo") == TRUE` but `cairo_pdf()` only emits
#' a warning such as `failed to load cairo DLL` and does not open a usable
#' device. Catch those warnings (and hard errors), then fall back to the
#' native PDF device. The native Windows PDF device can silently omit text
#' grobs under some R/ggplot2 combinations, so Cairo remains the preferred
#' path when it opens successfully.
#'
#' @param filename Output PDF path.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param ... Additional graphics-device arguments.
#' @return Invisibly returns the active graphics device number.
#' @keywords internal
cpopvar_pdf_device <- function(filename, width, height, ...) {
  cairo_opened <- FALSE
  if (isTRUE(capabilities("cairo"))) {
    cairo_opened <- isTRUE(tryCatch(
      {
        withCallingHandlers(
          {
            grDevices::cairo_pdf(
              filename = filename,
              width = width,
              height = height,
              ...
            )
          },
          warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("cairo|DLL", msg, ignore.case = TRUE)) {
              tryInvokeRestart("muffleWarning")
              stop(msg, call. = FALSE)
            }
          }
        )
        # Refuse a silent no-op: device must actually advance past null device.
        if (grDevices::dev.cur() <= 1L) {
          stop("cairo_pdf() returned without opening a graphics device", call. = FALSE)
        }
        TRUE
      },
      error = function(e) {
        warning(
          paste0(
            "Cairo PDF device failed to open (",
            conditionMessage(e),
            "); falling back to the native PDF device."
          ),
          call. = FALSE
        )
        # Close a half-opened Cairo device if one exists before native fallback.
        if (grDevices::dev.cur() > 1L) {
          try(grDevices::dev.off(), silent = TRUE)
        }
        FALSE
      }
    ))
  }

  if (!cairo_opened) {
    if (!isTRUE(capabilities("cairo"))) {
      warning(
        paste0(
          "Cairo PDF support is unavailable; falling back to the native PDF device. ",
          "Text rendering may be incomplete on Windows."
        ),
        call. = FALSE
      )
    }
    grDevices::pdf(
      file = filename,
      width = width,
      height = height,
      ...
    )
  }

  invisible(grDevices::dev.cur())
}

#' Open the cross-platform PDF device for base/grid drawing
#'
#' @param file Output PDF path.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param ... Additional graphics-device arguments.
#' @return Invisibly returns the active graphics device number.
#' @keywords internal
open_cpopvar_pdf <- function(file, width, height, ...) {
  cpopvar_pdf_device(
    filename = file,
    width = width,
    height = height,
    ...
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
  
  device_before_save <- grDevices::dev.cur()
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
        device = cpopvar_pdf_device,
        units = "in"
      )
      
    } else if ("pheatmap" %in% plot_class) {
      # Handle pheatmap objects - vectorized PDF output
      open_cpopvar_pdf(
        file = file_path,
        width = default_width,
        height = default_height
      )
      grid::grid.draw(plot_object)
      grDevices::dev.off()
      
    } else if (any(c("recordedplot", "grob", "gtable") %in% plot_class)) {
      # Handle other grid/base graphics objects - vectorized PDF output
      open_cpopvar_pdf(
        file = file_path,
        width = default_width,
        height = default_height
      )
      grid::grid.draw(plot_object)
      grDevices::dev.off()
      
    } else {
      # Fallback: try to draw the object using grid - vectorized PDF output
      open_cpopvar_pdf(
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
    # Close only a device opened by this save attempt. Do not close the
    # caller's RStudio or other interactive device after ggsave restores it.
    current_device <- grDevices::dev.cur()
    if (current_device > 1 && current_device != device_before_save) {
      grDevices::dev.off(which = current_device)
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
      if (identical(tolower(tools::file_ext(file_path)), "pdf")) {
        ggplot2::ggsave(
          filename = file_path,
          plot = plot,
          width = width,
          height = height,
          dpi = dpi,
          device = cpopvar_pdf_device
        )
      } else {
        ggplot2::ggsave(
          filename = file_path,
          plot = plot,
          width = width,
          height = height,
          dpi = dpi
        )
      }
      
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
############################################################
#### format_converter.R - Plot Format Conversion Module ####
############################################################
#
# Post-processing format converter
# Converts PDF plots to additional formats (PNG, SVG, etc.) after analysis completes
# Called before reporting module to ensure Markdown-compatible images are available
#
# Architecture: PDF as the primary figure format, with optional conversion
# - Analysis modules: Always output PDF (vectorized, high-quality)
# - Post-processing: Batch convert PDF to PNG/other formats
# - Reporting: Intelligently select best format for embedding
#
############################################################

#' Check optional dependencies needed by requested plot conversions
#'
#' @param export_formats Character vector of requested export formats.
#' @param package_available Function used to test package availability.
#' @return List containing `ok`, `required`, `missing`, and `message`.
#' @keywords internal
check_plot_conversion_dependencies <- function(
    export_formats,
    package_available = function(package) {
      requireNamespace(package, quietly = TRUE)
    }) {
  formats <- unique(tolower(as.character(export_formats)))
  raster_formats <- intersect(formats, c("png", "jpg", "jpeg"))
  required <- if (length(raster_formats) > 0) {
    c("magick", "pdftools")
  } else {
    character(0)
  }
  missing <- required[!vapply(required, package_available, logical(1))]

  if (length(missing) == 0) {
    return(list(
      ok = TRUE,
      required = required,
      missing = character(0),
      message = NULL
    ))
  }

  install_hint <- sprintf(
    "install.packages(c(%s))",
    paste(sprintf('"%s"', missing), collapse = ", ")
  )
  list(
    ok = FALSE,
    required = required,
    missing = missing,
    message = paste0(
      "Plot conversion cannot start because optional package(s) are missing: ",
      paste(missing, collapse = ", "),
      ". Install them in the R library used to run cpopvar with ",
      install_hint,
      ". PDF plots were preserved and no raster conversion was attempted."
    )
  )
}

#' Convert PDF plots to additional formats
#' 
#' Batch conversion of PDF plots to PNG or other formats
#' for reporting and presentation purposes.
#' 
#' @param plots_dir Path to plots directory
#' @param config Configuration object
#' @param session_id Session identifier for logging
#' @return List of conversion results
#' @export
convert_plots_for_reporting <- function(plots_dir, config, session_id = NULL) {
  
  log_message("=== Starting Plot Format Conversion ===")
  log_message(sprintf("Plots directory: %s", plots_dir))
  
  # Extract configuration
  viz_config <- config$visualization_settings$output %||% list()
  export_formats <- unique(tolower(
    as.character(viz_config$export_formats %||% character(0))
  ))
  
  # Check if conversion is needed
  if (length(export_formats) == 0) {
    log_message("No additional export formats specified, skipping conversion")
    return(list(
      success = TRUE,
      total_converted = 0,
      total_failed = 0,
      message = "No conversion needed"
    ))
  }

  dependency_check <- check_plot_conversion_dependencies(export_formats)
  if (!dependency_check$ok) {
    log_message(dependency_check$message, level = "error")
    return(list(
      success = FALSE,
      total_converted = 0,
      total_failed = 0,
      missing_packages = dependency_check$missing,
      error = dependency_check$message
    ))
  }

  log_message(sprintf("Export formats: %s", paste(export_formats, collapse = ", ")))
  
  # Find all PDF files
  if (!dir.exists(plots_dir)) {
    log_message(sprintf("Plots directory not found: %s", plots_dir), level = "warning")
    return(list(
      success = FALSE,
      total_converted = 0,
      total_failed = 0,
      error = "Plots directory not found"
    ))
  }
  
  pdf_files <- list.files(plots_dir, pattern = "\\.pdf$", 
                          recursive = TRUE, full.names = TRUE)
  
  if (length(pdf_files) == 0) {
    log_message("No PDF files found for conversion")
    return(list(
      success = TRUE,
      total_converted = 0,
      total_failed = 0,
      message = "No PDF files to convert"
    ))
  }
  
  log_message(sprintf("Found %d PDF files for format conversion", length(pdf_files)))
  
  # Initialize results tracking
  conversion_results <- list()
  total_converted <- 0
  total_failed <- 0
  
  # Process each format
  for (format in export_formats) {
    log_message(sprintf("Converting to %s format...", toupper(format)))
    
    format_results <- convert_to_format(
      pdf_files = pdf_files,
      target_format = format,
      config = viz_config
    )
    
    conversion_results[[format]] <- format_results
    total_converted <- total_converted + format_results$success_count
    total_failed <- total_failed + format_results$failure_count
  }
  
  # Summary
  log_message(sprintf("Conversion complete: %d files converted, %d failed", 
                     total_converted, total_failed))
  
  conversion_success <- total_failed == 0
  return(list(
    success = conversion_success,
    total_converted = total_converted,
    total_failed = total_failed,
    error = if (conversion_success) NULL else sprintf(
      "%d plot conversion(s) failed. See the session log for details.",
      total_failed
    ),
    by_format = conversion_results
  ))
}

#' Convert PDF files to a specific format
#' 
#' @param pdf_files Vector of PDF file paths
#' @param target_format Target format ("png", "svg", etc.)
#' @param config Visualization configuration
#' @return List with conversion statistics
convert_to_format <- function(pdf_files, target_format, config) {
  
  success_count <- 0
  failure_count <- 0
  converted_files <- character(0)
  
  for (pdf_path in pdf_files) {
    output_path <- sub("\\.pdf$", paste0(".", target_format), pdf_path)
    
    # Skip if already exists
    if (file.exists(output_path)) {
      log_message(sprintf("Skipping existing file: %s", basename(output_path)), 
                 level = "debug")
      success_count <- success_count + 1
      converted_files <- c(converted_files, output_path)
      next
    }
    
    result <- tryCatch({
      if (target_format == "png") {
        convert_pdf_to_png(pdf_path, output_path, config)
      } else if (target_format == "svg") {
        convert_pdf_to_svg(pdf_path, output_path, config)
      } else if (target_format == "jpg" || target_format == "jpeg") {
        convert_pdf_to_jpg(pdf_path, output_path, config)
      } else {
        stop(sprintf("Unsupported format: %s", target_format))
      }
      
      success_count <- success_count + 1
      converted_files <- c(converted_files, output_path)
      TRUE
      
    }, error = function(e) {
      log_message(sprintf("Failed to convert %s: %s", 
                         basename(pdf_path), e$message), 
                 level = "warning")
      failure_count <- failure_count + 1
      FALSE
    })
  }
  
  return(list(
    format = target_format,
    success_count = success_count,
    failure_count = failure_count,
    converted_files = converted_files
  ))
}

#' Convert PDF to PNG
#' 
#' @param pdf_path Input PDF file path
#' @param png_path Output PNG file path
#' @param config Visualization configuration
#' @return TRUE if successful
convert_pdf_to_png <- function(pdf_path, png_path, config) {
  
  # Get PNG-specific settings
  png_dpi <- config$png_dpi %||% 300
  png_bg <- config$png_bg %||% "white"
  
  # Read PDF as image
  img <- magick::image_read_pdf(pdf_path, density = png_dpi)
  
  # Set background (important for transparency)
  img <- magick::image_background(img, png_bg)
  
  # Write PNG
  magick::image_write(img, png_path, format = "png", quality = 100)
  
  # Verify file was created
  if (!file.exists(png_path)) {
    stop("PNG file was not created")
  }
  
  return(TRUE)
}

#' Convert PDF to SVG
#' 
#' @param pdf_path Input PDF file path
#' @param svg_path Output SVG file path
#' @param config Visualization configuration
#' @return TRUE if successful
convert_pdf_to_svg <- function(pdf_path, svg_path, config) {
  
  # Note: SVG conversion requires additional tools (pdf2svg or similar)
  # This is a placeholder for future implementation
  
  log_message("SVG conversion not yet implemented", level = "warning")
  return(FALSE)
}

#' Convert PDF to JPG
#' 
#' @param pdf_path Input PDF file path
#' @param jpg_path Output JPG file path
#' @param config Visualization configuration
#' @return TRUE if successful
convert_pdf_to_jpg <- function(pdf_path, jpg_path, config) {
  
  # Get JPG-specific settings
  jpg_dpi <- config$png_dpi %||% 300  # Reuse PNG DPI setting
  jpg_quality <- 95  # High quality for scientific figures
  
  # Read PDF as image
  img <- magick::image_read_pdf(pdf_path, density = jpg_dpi)
  
  # Set white background (JPG doesn't support transparency)
  img <- magick::image_background(img, "white")
  
  # Write JPG
  magick::image_write(img, jpg_path, format = "jpeg", quality = jpg_quality)
  
  # Verify file was created
  if (!file.exists(jpg_path)) {
    stop("JPG file was not created")
  }
  
  return(TRUE)
}

#' Find best plot format for reporting
#' 
#' Select the best available plot format for reporting
#' Priority: PNG > PDF (for better Markdown compatibility)
#' 
#' @param base_path Base path without extension
#' @return Full path to best available format, or NULL if not found
#' @export
find_best_plot_format <- function(base_path) {
  
  # Priority order for Markdown embedding
  format_priority <- c(".png", ".jpg", ".jpeg", ".pdf")
  
  for (ext in format_priority) {
    full_path <- paste0(base_path, ext)
    if (file.exists(full_path)) {
      return(full_path)
    }
  }
  
  return(NULL)
}


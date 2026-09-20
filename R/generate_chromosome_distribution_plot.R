############################################################
#### generate_chromosome_distribution_plot.R - Chromosome Ideogram Visualization ####
############################################################
#
# Generates publication-quality chromosome ideogram plots with variant density
# overlays using the RIdeogram package. Supports configurable color schemes
# and output formats.
#
# Part of M05 Chromosome Distribution Analysis Module
#
############################################################

#' Lighten a hex color by a specified amount
#' @param hex_color A hex color string (e.g., "#3B4992FF")
#' @param amount Lightening factor (0-1, where 1 is white)
#' @return Lightened hex color
#' @keywords internal
lighten_color <- function(hex_color, amount = 0.5) {
  # Convert hex to RGB
  rgb_vals <- col2rgb(hex_color)[,1]
  
  # Lighten by interpolating towards white (255, 255, 255)
  lightened_rgb <- rgb_vals + (255 - rgb_vals) * amount
  
  # Convert back to hex
  rgb(lightened_rgb[1], lightened_rgb[2], lightened_rgb[3], maxColorValue = 255)
}

#' Generate pan-species chromosome distribution ideogram plot with dual tracks
#'
#' Creates a chromosome ideogram visualization with dual variant density overlays
#' using the RIdeogram package. The function generates publication-quality
#' plots showing the spatial distribution of two variant types across multiple species.
#'
#' @param karyotype_data Data frame containing pan-species chromosome structure with columns:
#'   - Chr: Species chromosome identifier (character, from species_map)  
#'   - Start: Start position (numeric, typically 0)
#'   - End: End position/chromosome length (numeric)
#' @param track1_data Data frame containing first track density information with columns:
#'   - Chr: Species chromosome identifier (character)
#'   - Start: Start position of window (numeric)
#'   - End: End position of window (numeric)
#'   - Value: Density value for visualization (numeric)
#'   Additional columns are ignored
#' @param track2_data Data frame containing second track information. Format depends on label_type:
#'   For label_type="heatmap": columns Chr, Start, End, Value (same as track1_data)
#'   For label_type="marker": columns Type, Shape, Chr, Start, End, Label, color
#' @param label_type Character string specifying track2 display mode:
#'   "heatmap" for density heatmap overlay (default)
#'   "marker" for genomic region markers (IRa, IRb, LSC, SSC annotations)
#' @param config Configuration object containing visualization settings
#' @param task_params Task-specific parameters for configuration-driven plotting
#' @param output_path_base Base path for output files (without extension)
#'
#' @return Character string containing the path to the final generated plot file
#'
#' @details
#' This function creates chromosome ideogram visualizations by:
#' 1. Preparing data in RIdeogram-compatible format
#' 2. Calling RIdeogram::ideogram() with density overlay
#' 3. Converting SVG output to final format (PDF/PNG) using RIdeogram::convertSVG()
#' 4. Applying configuration-driven visual parameters
#'
#' The function supports configurable color schemes through the config object.
#' Default colors are provided if configuration is not available.
#'
#' Output files are generated in the same directory as output_path_base with
#' appropriate extensions (.svg, .pdf, .png).
#'
#' @examples
#' \dontrun{
#' # Example density data
#' density_data <- data.frame(
#'   Chr = c("1", "1", "2"),
#'   Window_Start = c(0, 2000, 0),
#'   Window_End = c(1999, 3999, 1999),
#'   Variant_Density = c(0.005, 0.008, 0.003)
#' )
#' 
#' # Example karyotype data
#' karyotype_data <- data.frame(
#'   Chr = c("1", "2"),
#'   Start = c(1, 1),
#'   End = c(10000, 8000)
#' )
#' 
#' # Generate plot
#' plot_path <- generate_chromosome_distribution_plot(
#'   density_data, karyotype_data, config, "/path/to/output"
#' )
#' }
#'
#' @importFrom RIdeogram ideogram convertSVG
#' @export
generate_chromosome_distribution_plot <- function(karyotype_data,
                                                 track1_data, 
                                                 track2_data,
                                                 label_type = "heatmap",
                                                 config,
                                                 task_params,
                                                 output_path_base) {
  
  # Input validation
  if (is.null(karyotype_data) || !is.data.frame(karyotype_data)) {
    stop("karyotype_data must be a non-null data frame")
  }
  
  if (is.null(track1_data) || !is.data.frame(track1_data)) {
    stop("track1_data must be a non-null data frame")
  }
  
  if (is.null(track2_data) || !is.data.frame(track2_data)) {
    stop("track2_data must be a non-null data frame")
  }
  
  if (is.null(output_path_base) || !is.character(output_path_base)) {
    stop("output_path_base must be a non-null character string")
  }
  
  # Check required columns in karyotype_data
  required_karyotype_cols <- c("Chr", "Start", "End")
  missing_karyotype_cols <- setdiff(required_karyotype_cols, names(karyotype_data))
  if (length(missing_karyotype_cols) > 0) {
    stop(sprintf("karyotype_data is missing required columns: %s", 
                 paste(missing_karyotype_cols, collapse = ", ")))
  }
  
  # Check required columns in track1_data (RIdeogram format)
  required_density_cols <- c("Chr", "Start", "End", "Value")
  missing_track1_cols <- setdiff(required_density_cols, names(track1_data))
  if (length(missing_track1_cols) > 0) {
    stop(sprintf("track1_data is missing required columns: %s", 
                 paste(missing_track1_cols, collapse = ", ")))
  }
  
  # Check required columns in track2_data based on label_type
  if (label_type == "marker") {
    # For marker mode, require marker-specific columns
    required_marker_cols <- c("Type", "Shape", "Chr", "Start", "End", "color")
    missing_track2_cols <- setdiff(required_marker_cols, names(track2_data))
    if (length(missing_track2_cols) > 0) {
      stop(sprintf("track2_data for marker mode is missing required columns: %s", 
                   paste(missing_track2_cols, collapse = ", ")))
    }
  } else {
    # For heatmap mode, require density-specific columns
    missing_track2_cols <- setdiff(required_density_cols, names(track2_data))
    if (length(missing_track2_cols) > 0) {
      stop(sprintf("track2_data for heatmap mode is missing required columns: %s", 
                   paste(missing_track2_cols, collapse = ", ")))
    }
  }
  
  # Ensure chromosome identifiers are character type
  karyotype_data$Chr <- as.character(karyotype_data$Chr)
  track1_data$Chr <- as.character(track1_data$Chr)
  track2_data$Chr <- as.character(track2_data$Chr)
  
  # Prepare RIdeogram-compatible data format
  # RIdeogram expects specific column names and data structure
  
  # Prepare karyotype data for RIdeogram (Chr, Start, End format)
  ideogram_karyotype <- karyotype_data[, c("Chr", "Start", "End")]
  
  # Prepare track 1 density data for RIdeogram overlay (overlaid parameter)
  # RIdeogram expects (Chr, Start, End, Value) format for overlaid data
  ideogram_density_1 <- track1_data[, c("Chr", "Start", "End", "Value")]
  
  # Prepare track 2 data for RIdeogram label (label parameter) based on label_type
  if (label_type == "marker") {
    # For marker mode, use the complete marker data structure
    ideogram_track2 <- track2_data
  } else {
    # For heatmap mode, use density data format  
    # RIdeogram expects (Chr, Start, End, Value) format for label data
    ideogram_track2 <- track2_data[, c("Chr", "Start", "End", "Value")]
  }
  
  # Data structure validation before RIdeogram call
  if (exists("log_message")) {
    log_message("Final data validation before RIdeogram visualization:")
    log_message(sprintf("Karyotype: %d species, Chr range: %s", 
                       nrow(ideogram_karyotype), 
                       paste(range(ideogram_karyotype$Chr), collapse = "-")))
    log_message(sprintf("Track 1 density: %d windows, Chr range: %s", 
                       nrow(ideogram_density_1), 
                       paste(range(ideogram_density_1$Chr), collapse = "-")))
    if (label_type == "marker") {
      log_message(sprintf("Track 2 markers: %d markers, Chr range: %s", 
                         nrow(ideogram_track2), 
                         paste(range(ideogram_track2$Chr), collapse = "-")))
    } else {
      log_message(sprintf("Track 2 density: %d windows, Chr range: %s", 
                         nrow(ideogram_track2), 
                         paste(range(ideogram_track2$Chr), collapse = "-")))
    }
    
    # Validate Chr consistency
    karyotype_chrs <- sort(unique(ideogram_karyotype$Chr))
    track1_chrs <- sort(unique(ideogram_density_1$Chr))
    track2_chrs <- sort(unique(ideogram_track2$Chr))
    
    log_message(sprintf("Chr consistency check:"))
    log_message(sprintf("  - Karyotype Chrs: %s", paste(karyotype_chrs, collapse = ", ")))
    log_message(sprintf("  - Track1 Chrs: %s", paste(track1_chrs, collapse = ", ")))
    log_message(sprintf("  - Track2 Chrs: %s", paste(track2_chrs, collapse = ", ")))
  }
  
  # Extract visualization parameters using configuration-driven approach
  # Dynamic color construction based on variant type from track1_data
  
  # Determine variant type from track1_data metadata (if available in task_params)
  variant_type_track1 <- NULL
  if (!is.null(task_params$current_plot_mode)) {
    # Extract variant type from plot mode (e.g., "snp_vs_indel" -> "snp" for track1)
    if (grepl("^snp_", task_params$current_plot_mode)) {
      variant_type_track1 <- "snp"
    } else if (grepl("^indel_", task_params$current_plot_mode)) {
      variant_type_track1 <- "INDEL"
    }
  }
  
  # Build gradient colors dynamically based on variant type using var_type colors
  var_type_colors <- get_color_palette(config, "var_type")
  
  if (!is.null(variant_type_track1) && !is.null(var_type_colors) && variant_type_track1 %in% names(var_type_colors)) {
    # Build 3-level gradient: white -> light color -> configured color
    base_color <- var_type_colors[[variant_type_track1]]
    colorset1 <- c("#f7f7f7", lighten_color(base_color, 0.5), base_color)
    log_message(sprintf("M05: Using variant-specific color for %s: %s", variant_type_track1, base_color))
    log_message(sprintf("M05: Track 1 gradient: %s", paste(colorset1, collapse = " -> ")))
  } else {
    # Fallback to default blue gradient
    colorset1 <- c("#f7f7f7", "#80B1D3", "#0000FF")
    log_message("M05: Using default track 1 colors (variant color not found in config)")
  }
  
  # colorset2 is only used for heatmap mode (track2 when label_type == "heatmap")
  if (label_type == "heatmap") {
    # Determine track2 variant type
    variant_type_track2 <- NULL
    if (!is.null(task_params$current_plot_mode) && grepl("snp_vs_indel", task_params$current_plot_mode)) {
      variant_type_track2 <- "INDEL"  # In snp_vs_indel, track2 is INDEL
    }
    
    # Build gradient colors using var_type colors for track2
    if (!is.null(variant_type_track2) && !is.null(var_type_colors) && variant_type_track2 %in% names(var_type_colors)) {
      base_color <- var_type_colors[[variant_type_track2]]
      colorset2 <- c("#f7f7f7", lighten_color(base_color, 0.5), base_color)
      log_message(sprintf("M05: Using variant-specific color for %s: %s", variant_type_track2, base_color))
      log_message(sprintf("M05: Track 2 gradient: %s", paste(colorset2, collapse = " -> ")))
    } else {
      # Fallback to default orange gradient
      colorset2 <- c("#f7f7f7", "#FB8072", "#FF4500")
      log_message("M05: Using default track 2 colors (variant color not found in config)")
    }
  } else {
    # For marker mode, colorset2 is not used
    colorset2 <- NULL
  }
  
  # Get plotting parameters using standard configuration approach
  # label_type is now passed as parameter and validated above
  width_param <- get_task_parameter(task_params, config, "width", 5.85)
  height_param <- get_task_parameter(task_params, config, "height", 8.26)
  dpi_param <- get_task_parameter(task_params, config, "dpi", 600)
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_path_base)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set working directory to output directory for RIdeogram
  # RIdeogram creates files in the current working directory
  original_wd <- getwd()
  setwd(output_dir)
  
  tryCatch({
    
    # Generate ideogram with track 1 (overlaid) and track 2 (label) visualization
    # Track 2 display mode is determined by label_type parameter
    RIdeogram::ideogram(
      karyotype = ideogram_karyotype,
      overlaid = ideogram_density_1,    # Track 1: First variant type density (e.g., SNP)
      label = ideogram_track2,          # Track 2: Second track data (density or markers)
      label_type = label_type,          # Display mode: "heatmap" or "marker"
      colorset1 = colorset1,            # Color scheme for track 1
      colorset2 = colorset2,            # Color scheme for track 2 (used for heatmap mode)
      Lx = 185,                         # Legend X position
      Ly = 20                           # Legend Y position
    )
    
    # Define output file paths
    svg_file <- "chromosome.svg"  # RIdeogram default output name
    final_output <- paste0(basename(output_path_base), ".pdf")
    
    # Convert SVG to final format (default to PDF)
    if (grepl("\\.pdf$", final_output) || !grepl("\\.(png|jpg|jpeg)$", final_output)) {
      # Convert to PDF format
      RIdeogram::convertSVG(
        svg_file, 
        device = "pdf",
        width = width_param,
        height = height_param,
        dpi = dpi_param
      )
      final_file <- paste0(tools::file_path_sans_ext(svg_file), ".pdf")
    } else {
      # Convert to PNG format
      RIdeogram::convertSVG(
        svg_file, 
        device = "png",
        width = width_param,
        height = height_param,
        dpi = dpi_param
      )
      final_file <- paste0(tools::file_path_sans_ext(svg_file), ".png")
    }
    
    # Rename to desired output name if different
    if (final_file != final_output) {
      if (file.exists(final_file)) {
        file.rename(final_file, final_output)
        final_file <- final_output
      }
    }
    
    # Construct full path for return
    final_path <- file.path(output_dir, final_file)
    
    # Log success
    if (exists("log_message")) {
      log_message(sprintf("Successfully generated pan-species dual-track chromosome ideogram: %s", final_path))
      log_message(sprintf("Track 1: %d density windows across %d species", 
                         nrow(ideogram_density_1), 
                         length(unique(ideogram_density_1$Chr))))
      if (label_type == "marker") {
        log_message(sprintf("Track 2: %d markers across %d species", 
                           nrow(ideogram_track2), 
                           length(unique(ideogram_track2$Chr))))
      } else {
        log_message(sprintf("Track 2: %d density windows across %d species", 
                           nrow(ideogram_track2), 
                           length(unique(ideogram_track2$Chr))))
      }
      log_message(sprintf("Total visualization contains %d species represented as chromosomes", 
                         length(unique(ideogram_karyotype$Chr))))
    }
    
    return(final_path)
    
  }, error = function(e) {
    
    # Log error details
    if (exists("log_message")) {
      log_message(sprintf("Error generating chromosome ideogram: %s", e$message), level = "error")
    }
    
    # Re-throw error for upstream handling
    stop(sprintf("Failed to generate chromosome ideogram: %s", e$message))
    
  }, finally = {
    
    # Always restore original working directory
    setwd(original_wd)
    
  })
}
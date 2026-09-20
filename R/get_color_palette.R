#' Get a color palette based on V4 architecture unified color system.
#'
#' This function implements the three-layer override model:
#' 1. Custom Palettes (highest priority): by_variable and by_concept
#' 2. Default Palette (global fallback): ggsci palette selection  
#' 3. ggplot2 Engine (basic elements): handled elsewhere
#'
#' @param config The full configuration object for the current run.
#' @param concept_or_variable Either a plot_type concept name (e.g., "global_heatmap") or variable name (e.g., "region_type").
#' @param variable_values A vector of the unique actual values from the data for the given variable (for variable mapping).
#' @return A named vector/list of colors suitable for ggplot2 scales.
#' @export
get_color_palette <- function(config, concept_or_variable, variable_values = NULL) {

  # V4 Architecture: Three-layer override model
  
  # Layer 1 (Highest Priority): Custom Palettes
  if (!is.null(config$visualization_settings$color_system$custom_palettes)) {
    custom_palettes <- config$visualization_settings$color_system$custom_palettes
    
    # V4.1 PRIORITY FIX: Check by_variable FIRST (highest priority for data-driven colors)
    if (!is.null(custom_palettes$by_variable) && 
        concept_or_variable %in% names(custom_palettes$by_variable)) {
      palette_def <- custom_palettes$by_variable[[concept_or_variable]]
      # Custom variable palette found
      
      if (palette_def$type == "manual") {
        if (!is.null(variable_values)) {
          defined_colors <- palette_def$values
          # Ensure all required values are present
          missing_keys <- setdiff(unique(variable_values), names(defined_colors))
          if (length(missing_keys) > 0) {
            warning(sprintf("Custom palette '%s' missing keys: %s", 
                          concept_or_variable, paste(missing_keys, collapse=", ")))
            fallback_colors <- setNames(rep("#808080", length(missing_keys)), missing_keys)
            final_colors <- c(defined_colors, fallback_colors)
            # Ensure we return named vector, not list
            if (is.list(final_colors)) {
              final_colors <- unlist(final_colors)
            }
            return(final_colors)
          }
          # Ensure we return named vector, not list
          if (is.list(defined_colors)) {
            defined_colors <- unlist(defined_colors)
          }
          return(defined_colors)
        }
        # Ensure we return named vector, not list
        final_values <- palette_def$values
        if (is.list(final_values)) {
          final_values <- unlist(final_values)
        }
        return(final_values)
      }
    }
    
    # Check by_concept second (for plot_types like "global_heatmap")
    if (!is.null(custom_palettes$by_concept) && 
        concept_or_variable %in% names(custom_palettes$by_concept)) {
      palette_def <- custom_palettes$by_concept[[concept_or_variable]]
      # Custom concept palette found
      # Return the palette definition as-is to preserve structure
      # For concept palettes, maintain list structure (e.g., list(gradient = c(...)))
      return(palette_def)
    }
  }

  # Layer 2 (Global Default): Default Palette
  if (!is.null(config$visualization_settings$color_system$default_palette)) {
    default_palette <- config$visualization_settings$color_system$default_palette$selected
    if (!is.null(default_palette)) {
      # Using default palette
      
      # For data variables, generate named colors
      if (!is.null(variable_values)) {
        n_colors <- length(unique(variable_values))
        colors <- get_ggsci_colors(default_palette, n_colors)
        
        # Safe attribute setting - check if colors is not NULL
        if (!is.null(colors) && length(colors) > 0) {
          names(colors) <- unique(variable_values)
          return(colors)
        } else {
          # Fallback if get_ggsci_colors returns NULL
          fallback_colors <- rainbow(n_colors)
          names(fallback_colors) <- unique(variable_values)
          return(fallback_colors)
        }
      }
      
      # For concepts, return a reasonable number of colors
      colors <- get_ggsci_colors(default_palette, 10)
      if (!is.null(colors) && length(colors) > 0) {
        return(colors)
      } else {
        # Fallback if get_ggsci_colors returns NULL
        return(rainbow(10))
      }
    }
  }

  # Fallback error
  stop(sprintf("V4 Color System Error: No color definition found for '%s' in custom_palettes or default_palette", concept_or_variable))
}

#' Generate enhanced pie chart matrix for variant type distribution
#' 
#' @title Generate enhanced pie chart matrix
#' @description Creates a matrix of pie charts showing variant type distribution across species and genomic regions
#' @param data Normalized frequency data with variant counts
#' @param config Configuration list containing visualization parameters
#' @param task_params Task-specific parameters for plot customization
#' @return List containing the pie matrix plot and metadata
#' @importFrom magrittr %>%
#' @export
generate_enhanced_pie_matrix <- function(data, config, task_params) {
  
  log_message("Generating enhanced pie chart matrix")
  
  tryCatch({
    # Get visual standards for consistent styling
    visual_standards <- get_visual_standards()
    geom_defaults <- visual_standards$ggplot2
    
    # Validate input data
    if (is.null(data) || nrow(data) == 0) {
      stop("No data provided for pie matrix")
    }
    
    # Check required columns
    required_cols <- c("species", "var_type")
    region_col <- "region_type"
    if (!region_col %in% names(data)) {
      if ("var_location" %in% names(data)) {
        region_col <- "var_location"
      } else {
        stop("Required region column not found")
      }
    }
    required_cols <- c(required_cols, region_col)
    
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Filter out "all_types" to focus on individual variant types
    if ("var_type" %in% names(data)) {
      pie_data <- data[data$var_type != "all_types", ]
    } else {
      pie_data <- data
    }
    
    # Trust the incoming data has already been correctly scoped upstream.
    # All data filtering and scoping is handled by the module orchestrator before this function is called.
    
    if (nrow(pie_data) == 0) {
      warning("No data available for enhanced pie matrix after filtering")
      return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("No Data Available"), 
                 data = data, metadata = list()))
    }
    
    # Aggregate data by species, var_type, and region for 4x17 layout
    # Each pie will show region composition within species-vartype combination
    if (!"variant_count" %in% names(pie_data)) {
      # Fallback if variant_count column is missing, count rows instead
      log_message("variant_count column not found, counting rows as a fallback.", level = "warning")
      pie_data_aggregated <- pie_data %>%
        dplyr::group_by(species, var_type, !!rlang::sym(region_col)) %>%
        dplyr::summarise(
          count = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::rename(region = !!rlang::sym(region_col))
    } else {
      # Standard path: sum the variant_count for accurate representation
      pie_data_aggregated <- pie_data %>%
        dplyr::group_by(species, var_type, !!rlang::sym(region_col)) %>%
        dplyr::summarise(
          count = sum(variant_count, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::rename(region = !!rlang::sym(region_col))
    }
    
    # V10 FIX: Complete factor combination matrix for dual-path architecture
    # This function receives pie_data which preserves all factor levels (including RNA)
    pie_summary <- pie_data_aggregated %>%
      dplyr::ungroup() %>%
      tidyr::complete(species, var_type, region, fill = list(count = 0))
    
    # V10 ENHANCEMENT: Preserve all factor levels for complete visualization
    # IMPORTANT: Unlike freq_data path, pie_data path must maintain ALL factor levels
    # including RNA to show complete species-variant-region combinations
    log_message("Converting character columns to factors for pie matrix plot.", level = "info")

    # Trust orchestrator's species ordering - only ensure it's a factor
    if (!is.factor(pie_summary$species)) {
      pie_summary$species <- factor(pie_summary$species, levels = unique(pie_summary$species))
    }
    # Preserve all region levels - critical for RNA display
    pie_summary$region_type <- factor(pie_summary$region, levels = unique(pie_summary$region))
    pie_summary$var_type <- factor(pie_summary$var_type, levels = unique(pie_summary$var_type))
    
    # Keep region column as factor with same levels for compatibility
    pie_summary$region <- factor(pie_summary$region, levels = unique(pie_summary$region))
    
    # Calculate total counts per species-vartype combination for percentages
    pie_summary <- pie_summary %>%
      dplyr::group_by(species, var_type) %>%
      dplyr::mutate(
        total_count_in_pie = sum(count),
        percentage = ifelse(total_count_in_pie == 0, 0, (count / total_count_in_pie) * 100)
      ) %>%
      dplyr::ungroup()
    
    # For zero-total combinations, create a single "No Data" slice for visualization
    pie_summary <- pie_summary %>%
      dplyr::group_by(species, var_type) %>%
      dplyr::mutate(
        # If all percentages are 0 for this pie, make the first region 100% for visualization
        percentage = ifelse(total_count_in_pie == 0 & dplyr::row_number() == 1, 100, 
                                 ifelse(total_count_in_pie == 0 & dplyr::row_number() > 1, 0, percentage))
      ) %>%
      dplyr::ungroup()
    
    # Calculate global count range for color gradient
    global_max_count <- max(pie_summary$total_count_in_pie, na.rm = TRUE)
    global_min_count <- min(pie_summary$total_count_in_pie, na.rm = TRUE)
    
    # Add normalized count for color mapping (0-1 scale)
    pie_summary <- pie_summary %>%
      dplyr::mutate(
        normalized_count = ifelse(global_max_count == global_min_count, 0,
                                       (total_count_in_pie - global_min_count) / (global_max_count - global_min_count))
      )
    
    # Get consistent colors for regions from the central manager
    region_colors <- get_color_palette(config, "region_type", levels(pie_summary$region))
    
    # Create the label text to show percentages instead of counts
    pie_summary <- pie_summary %>%
      dplyr::group_by(species, var_type) %>%
      dplyr::mutate(
        label_text = dplyr::case_when(
          total_count_in_pie == 0 & dplyr::row_number() == 1 ~ "No Data",  # Show "No Data" for empty pies
          total_count_in_pie == 0 ~ "",  # Hide labels for other slices of empty pies
          percentage >= 5 ~ paste0(round(percentage), "%"),  # Show percentage for significant slices (>5%)
          TRUE ~ ""  # Hide labels for small slices
        )
      ) %>%
      dplyr::ungroup()
    
    # Create the pie matrix using ggplot with polar coordinates
    # 4x17 layout: rows = variant types, columns = species, fill = region
    # Use percentage for y to ensure each pie is a complete circle
    pie_matrix_plot <- ggplot2::ggplot(pie_summary, ggplot2::aes(x = "", y = percentage, fill = region)) +
      ggplot2::geom_bar(stat = "identity", width = 1, color = geom_defaults$bar_color, linewidth = geom_defaults$bar_linewidth) +
      # Add text labels showing count values
      ggplot2::geom_text(ggplot2::aes(label = label_text), 
                        color = geom_defaults$text_color, size = geom_defaults$text_size, fontface = "bold",
                        position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::facet_grid(var_type ~ species, 
                          scales = "free",
                          labeller = ggplot2::labeller(
                            var_type = c("snp" = "SNP", "INDEL" = "INDEL", "complex" = "Complex", "mnp" = "MNP"),
                            species = function(x) gsub("_", " ", x)
                          )) +
      ggplot2::scale_fill_manual(values = region_colors, name = "Region") +
      ggplot2::labs(
        title = "Region Distribution Matrix: 4 Variant Types x 17 Species",
        subtitle = "Pie charts show region composition (CDS, IGS, intron) within each variant type-species combination",
        caption = "Colors represent genomic regions; percentages show relative distribution per region"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        # Remove axis elements for clean pie charts
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        
        # Style facet labels
        strip.text.x = ggplot2::element_blank(), # Remove species names from the top/bottom of pie charts
        strip.text.y =ggplot2::element_text(size = 10, face = "bold"),
        strip.background = ggplot2::element_blank(), # Remove facet label background
        
        # Style plot titles
        plot.title =ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle =ggplot2::element_text(hjust = 0.5, size = 12, color = "gray40"),
        plot.caption =ggplot2::element_text(hjust = 0.5, size = 10, color = "gray60"),
        
        # Clean panel appearance
        panel.grid = ggplot2::element_blank(),
        panel.background =ggplot2::element_rect(fill = "white", color = NA),
        plot.background =ggplot2::element_rect(fill = "white", color = NA),
        
        # Legend styling
        legend.position = "right",
        legend.title =ggplot2::element_text(size = 11, face = "bold"),
        legend.text =ggplot2::element_text(size = 10),
        legend.key.size = grid::unit(0.8, "cm"),
        plot.margin = ggplot2::margin(b = -10) # Reduce bottom margin to decrease gap to bar chart
      )
    
    # Calculate comprehensive statistics with crash fix
    matrix_stats <- pie_summary %>%
      dplyr::group_by(species, var_type) %>%
      dplyr::summarise(
        total_variants = sum(count),
        regions_present = dplyr::n_distinct(region[count > 0]),  # Only count regions with data
        dominant_region = if (all(is.na(count))) NA_character_ else region[which.max(count)],
        dominant_percentage = ifelse(sum(count, na.rm = TRUE) == 0, 0, max(percentage, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # Calculate global statistics
    global_stats <- list(
      total_species = dplyr::n_distinct(pie_summary$species),
      total_regions = dplyr::n_distinct(pie_summary$region),
      total_variant_types = dplyr::n_distinct(pie_summary$var_type),
      total_combinations = nrow(matrix_stats),
      overall_variants = sum(pie_summary$count),
      species_list = levels(pie_summary$species),
      regions_list = levels(pie_summary$region),
      variant_types_list = levels(pie_summary$var_type)
    )
    
    log_message(sprintf("Enhanced pie matrix generated: %dx%d = %d pie charts", 
                       global_stats$total_variant_types, global_stats$total_species, 
                       global_stats$total_combinations))
    
    # --- START: New Bar Chart for Magnitude Comparison ---
    
    # Create summary data for the bar chart view
    # Aggregate by species and region to show total variants per region (sum across all variant types)
    bar_chart_data <- pie_summary %>%
      dplyr::group_by(species, region) %>%
      dplyr::summarise(total_count = sum(count), .groups = "drop")

    # Create stacked bar chart with species on x-axis
    # Calculate total variants per species for proper stacking
    species_totals <- bar_chart_data %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(species_total = sum(total_count), .groups = "drop")

    # Define consistent colors for regions (same as pie chart)
    region_colors_bar <- region_colors

    # Create the stacked bar chart plot with species on x-axis
    bar_plot <- ggplot2::ggplot(bar_chart_data, ggplot2::aes(x = species, y = total_count, fill = region)) +
      ggplot2::geom_col(position = "stack") + # Create stacked bars
      # Add count labels in middle of each stacked segment
      ggplot2::geom_text(ggplot2::aes(label = ifelse(total_count > 0, total_count, "")), 
                         position = ggplot2::position_stack(vjust = 0.5), 
                         size = 2.5, color = "white", fontface = "bold") +
      # Add total count labels above bars
      ggplot2::geom_text(data = bar_chart_data %>% 
                           dplyr::group_by(species) %>% 
                           dplyr::summarise(total_count = sum(total_count)),
                         ggplot2::aes(x = species, y = total_count, label = total_count, fill = NULL), 
                         vjust = -0.5, size = 3, color = "black", fontface = "bold") +
      ggplot2::scale_fill_manual(values = region_colors_bar, name = "Region") +
      ggplot2::labs(
        y = "Total Variants", 
        x = "Species"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x =ggplot2::element_text(angle = 45, hjust = 1, size = 8), # Rotate species labels for readability
        axis.title.y =ggplot2::element_text(size = 10, angle = 90),
        axis.text.y =ggplot2::element_text(size = 9),
        legend.position = "right", # Show legend for region identification in stacked bars
        legend.title =ggplot2::element_text(size = 10, face = "bold"),
        legend.text =ggplot2::element_text(size = 9),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0, r = 10, b = 40, l = 10) # Extra bottom margin for rotated labels
      )

    # --- START: Advanced Grob Alignment for Perfect Width Matching ---

    # Convert plots to grob objects to manually align them
    pie_grob <- ggplot2::ggplotGrob(pie_matrix_plot)
    bar_grob <- ggplot2::ggplotGrob(bar_plot)

    # Identify the columns in each plot's grid layout that contain the main plotting panels
    pie_panel_cols <- pie_grob$layout[grepl("panel", pie_grob$layout$name), ]$l
    bar_panel_cols <- bar_grob$layout[grepl("panel", bar_grob$layout$name), ]$l
    
    # The core of the alignment: programmatically set the width of each pie chart panel
    # to be exactly equal to the width of the corresponding bar chart panel below it.
    # This forces the content areas to align perfectly.
    pie_grob$widths[pie_panel_cols] <- bar_grob$widths[bar_panel_cols]

    # Combine the now perfectly aligned grobs
    combined_plot <- gridExtra::grid.arrange(
      pie_grob, 
      bar_grob, 
      ncol = 1, 
      heights = c(0.75, 0.25) # Assign 75% height to pies, 25% to bars
    )
    
    # --- END: Advanced Grob Alignment ---
    
    # Prepare metadata
    metadata <- list(
      plot_dimensions = "4x17 matrix with summary bar chart",
      total_pie_charts = global_stats$total_combinations,
      color_scheme = region_colors,
      global_statistics = global_stats,
      matrix_statistics = matrix_stats,
      count_range = c(global_min_count, global_max_count),
      plot_type = "enhanced_pie_matrix_composite", 
      data_source = "detailed_normalized_frequencies",
      composite_design = TRUE
    )
    
    return(list(
      plot = combined_plot,
      data = pie_summary,
      statistics = matrix_stats,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Enhanced pie matrix generation failed: %s", e$message), level = "error")
    return(list(
      plot = ggplot2::ggplot() + 
        ggplot2::ggtitle("Enhanced Pie Matrix Generation Failed") +
        ggplot2::theme_void(),
      data = data,
      metadata = list(error = e$message)
    ))
  })
}
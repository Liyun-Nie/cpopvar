#' Generate region boxplot composite visualization
#' 
#' @title Generate region boxplot composite visualization
#' @description Creates composite boxplot with scatter overlay for region-based analysis
#' @param data Dataset for visualization
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @return List containing composite plot and data
#' @export
generate_region_boxplot_composite <- function(data, config, task_params) {
  
  log_message("Generating composite boxplot + scatter plot")
  
  tryCatch({
    # Get alpha scheme and color configuration
    alpha_scheme <- get_alpha_scheme()
    outlier_color <- get_task_parameter(task_params, config, "outlier_color", default = "red")
    normal_point_color <- get_task_parameter(task_params, config, "normal_point_color", default = "black")
    boxplot_border_color <- get_task_parameter(task_params, config, "boxplot_border_color", default = "black")
    
    # Validate input data
    if (is.null(data) || nrow(data) == 0) {
      stop("No data provided for composite boxplot")
    }
    
    # Check for required columns
    freq_col <- "frequency_per_kb"
    if (!freq_col %in% names(data)) {
      stop("Required frequency_per_kb column not found")
    }
    
    region_col <- "region_type"
    if (!region_col %in% names(data)) {
      stop("Required region_type column not found")
    }
    
    if (!"species" %in% names(data)) {
      stop("Species column required for boxplot")
    }
    
    # Trust upstream data preparation, no internal filtering
    boxplot_data <- data
    
    if (nrow(boxplot_data) == 0) {
      warning("No data available for composite boxplot")
      return(list(plot =ggplot2::ggplot() +ggplot2::ggtitle("No Data Available"), 
                 data = data, metadata = list()))
    }
    
    # Discover regions from the actual data (implicit ordering from data_scoping)
    actual_regions <- unique(boxplot_data[[region_col]])
    
    # Rename columns for consistency
    boxplot_data$frequency <- boxplot_data[[freq_col]]
    boxplot_data$region <- boxplot_data[[region_col]]
    
    # Apply species ordering using configuration
    species_order <- get_task_parameter(task_params, config, "species_order", NULL)
    if (!is.null(species_order)) {
      boxplot_data <- apply_species_ordering(boxplot_data, species_order_vector = species_order, species_col = "species")
    }
    
    # Calculate global y-axis limit for consistent comparison across panels
    global_y_max <- max(boxplot_data$frequency, na.rm = TRUE) * 1.05 # Add 5% buffer
    
    # Get color palette from the central manager using actual regions
    region_colors <- get_color_palette(config, "region_type", actual_regions)
    
    # Create list to store individual plots
    plot_list <- list()
    all_outliers <- data.frame()
    
    # Generate panel for each region (discovered from data)
    for (region in actual_regions) {
      region_data <- boxplot_data[boxplot_data$region == region, ]
      
      if (nrow(region_data) == 0) {
        # Create empty plot for missing regions
        plot_list[[region]] <- ggplot2::ggplot() + 
          ggplot2::ggtitle(paste("No", region, "Data")) + 
          theme_void()
        next
      }
      
      # Detect outliers using existing function
      outliers <- detect_outliers(
        data = region_data,
        value_col = "frequency",
        group_col = "species",
        method = "Z_score",
        #iqr_multiplier = get_task_parameter(task_params, config, "outlier_iqr_multiplier", default = 1.5),
        z_threshold = get_task_parameter(task_params, config, "outlier_z_threshold", default = 3),
        return_type = "outliers_only"
      )
      
      # Mark outliers in the data
      region_data$is_outlier <- region_data$frequency %in% outliers$frequency
      
      # Store outliers for metadata
      if (nrow(outliers) > 0) {
        outliers$region <- region
        all_outliers <- rbind(all_outliers, outliers)
      }
      
      # Get color for this region
      fill_color <- region_colors[[region]]
      
      # Create individual boxplot + scatter plot
      panel_plot <-ggplot2::ggplot(region_data,ggplot2::aes(x = species, y = frequency)) +
        # Boxplot without default outliers (we'll add custom ones)
        geom_boxplot(fill = fill_color, alpha = alpha_scheme$fill_main, outlier.shape = NA, 
                    color = boxplot_border_color, linewidth = 0.5) +
        # Add all points as scatter
        geom_point(data = region_data[!region_data$is_outlier, ],
                   ggplot2::aes(x = species, y = frequency), 
                   color = normal_point_color, size = 1.0, alpha = alpha_scheme$point_main, 
                   position = position_jitter(width = 0.25, height = 0)) +
        # Highlight outliers in configured color
        geom_point(data = region_data[region_data$is_outlier, ],
                   ggplot2::aes(x = species, y = frequency), 
                   color = outlier_color, size = 1.5, alpha = alpha_scheme$point_outlier, 
                   position = position_jitter(width = 0.25, height = 0)) +
        # Set a consistent y-axis limit across all panels for direct comparison
        ggplot2::coord_cartesian(ylim = c(0, global_y_max)) +
          ggplot2::labs(
          title = paste0(region, " Region"),
          subtitle = paste0("n=", nrow(region_data), ", outliers=", sum(region_data$is_outlier)),
          x = if(region == "intron") "Species" else "",  # Only show x-label on bottom plot
          y = "Frequency (per kb)"
        ) +
        { theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; 
          base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14;
          get_application_theme(theme_name, NULL, base_size) } +
        ggplot2::theme(
          axis.text.x =ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background =ggplot2::element_rect(fill = "white", color = NA)
        )
      
      # Add statistical annotations
      region_stats <- region_data %>%
        dplyr::summarise(
          mean_freq = mean(frequency, na.rm = TRUE),
          median_freq = median(frequency, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Add mean line
      panel_plot <- panel_plot +
        geom_hline(yintercept = region_stats$mean_freq, 
                   color = fill_color, linetype = "dashed", alpha = 0.7, linewidth = 0.5) +
        annotate("text", x = 1, y = region_stats$mean_freq,
                label = sprintf("mean=%.2f", region_stats$mean_freq),
                vjust = -0.5, size = 3, color = fill_color, fontface = "bold")
      
      plot_list[[region]] <- panel_plot
    }
    
    # Build a grob without drawing to the ambient device. save_plot() owns the
    # cross-platform PDF device and renders this object later.
    combined_plot <- do.call(gridExtra::arrangeGrob, c(plot_list, list(nrow = 1)))
    
    # Calculate summary statistics
    summary_stats <- boxplot_data %>%
      dplyr::group_by(region) %>%
      dplyr::summarise(
        count = dplyr::n(),
        species_count = dplyr::n_distinct(.data$species),
        mean_freq = mean(frequency, na.rm = TRUE),
        median_freq = median(frequency, na.rm = TRUE),
        sd_freq = sd(frequency, na.rm = TRUE),
        outliers = sum(frequency %in% all_outliers$frequency),
        .groups = "drop"
      )
    
    # Prepare metadata
    metadata <- list(
      total_points = nrow(boxplot_data),
      regions_analyzed = actual_regions,
      total_outliers = nrow(all_outliers),
      outliers_by_region = table(all_outliers$region),
      summary_statistics = summary_stats,
      color_scheme = region_colors,
      plot_type = "region_boxplot_composite",
      data_source = "gene_level_frequencies",
      outlier_method = "Z_score"
    )
    
    log_message(sprintf("Composite boxplot generated: %d regions, %d total outliers", 
                       length(actual_regions), nrow(all_outliers)))
    
    return(list(
      plot = combined_plot,
      individual_plots = plot_list,
      data = boxplot_data,
      outliers = all_outliers,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Composite boxplot generation failed: %s", e$message), level = "error")
    return(list(
      plot =ggplot2::ggplot() + 
        ggplot2::ggtitle("Composite Boxplot Generation Failed") +
        theme_void(),
      data = data,
      metadata = list(error = e$message)
    ))
  })
}
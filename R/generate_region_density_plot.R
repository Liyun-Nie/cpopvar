#' Generate region density plot
#' 
#' @title Generate region density plot
#' @description Creates density distribution plot for regional frequency analysis
#' @param data Dataset for visualization
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @return List containing density plot and analysis data
#' @export
generate_region_density_plot <- function(data, config, task_params) {
  
  log_message("Generating region-specific density distribution plot")
  
  tryCatch({
    # Get visual standards for consistent styling
    visual_standards <- get_visual_standards()
    geom_defaults <- visual_standards$ggplot2
    
    # Get alpha scheme
    alpha_scheme <- get_alpha_scheme()
    
    # Validate input data
    if (is.null(data) || nrow(data) == 0) {
      stop("No data provided for density plot")
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
    
    # Trust upstream data preparation, no internal filtering
    density_data <- data
    
    if (nrow(density_data) == 0) {
      warning("No data available for density plot")
      return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("No Data Available"), 
                 data = data, metadata = list()))
    }
    
    # Rename columns for consistency
    density_data$frequency <- density_data[[freq_col]]
    density_data$region <- density_data[[region_col]]
    
    # Define unique_regions for use in logging and metadata
    unique_regions <- unique(density_data$region)

    # Get color palette from the central manager
    region_colors <- get_color_palette(config, "region_type", unique_regions)
    
    # Calculate density statistics for each region, including peak coordinates
    density_stats <- density_data %>%
      dplyr::group_by(region) %>%
      dplyr::summarise(
        count = dplyr::n(),
        mean_freq = mean(frequency, na.rm = TRUE),
        median_freq = median(frequency, na.rm = TRUE),
        sd_freq = sd(frequency, na.rm = TRUE),
        .groups = "drop"
      )
    
    peak_coords <- density_data %>%
      dplyr::group_by(region) %>%
      dplyr::reframe({
        d <- density(frequency, na.rm = TRUE)
        tibble::tibble(
          peak_x = d$x[which.max(d$y)],
          peak_y = max(d$y)
        )
      })
    
    density_stats <- dplyr::left_join(density_stats, peak_coords, by = "region")
    
    log_message(sprintf("Density plot data: %d total points across %d regions", 
                       nrow(density_data), length(unique_regions)))
    
    # Create density plot
    density_plot <-ggplot2::ggplot(density_data, ggplot2::aes(x = frequency, fill = region, color = region)) +
      geom_density(alpha = alpha_scheme$fill_main, linewidth = geom_defaults$density_linewidth) +
      ggplot2::scale_fill_manual(values = region_colors, name = "Genomic Region") +
      ggplot2::scale_color_manual(values = region_colors, name = "Genomic Region") +
      ggplot2::labs(
        title = "Distribution of Variant Frequencies by Genomic Region",
        subtitle = "Density curves showing frequency patterns across CDS, IGS, and intron regions",
        x = "Variant Frequency (per kb)",
        y = "Density"
      ) +
      get_application_theme(config) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray40"),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        axis.text = ggplot2::element_text(size = 10),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_line(color = "gray90", size = 0.5),
        panel.grid.minor = ggplot2::element_line(color = "gray95", size = 0.3),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
    
    # Add statistics text annotation
    stats_text <- density_stats %>%
      dplyr::mutate(label = sprintf("%s: n=%d, mean=%.2f", region, count, mean_freq)) %>%
      dplyr::pull(label) %>%
      paste(collapse = "\n")
    
    # Get annotation text color from configuration
    annotation_text_color <- get_task_parameter(task_params, config, "annotation_text_color", default = "gray30")
    
    density_plot <- density_plot +
      annotate("text", x = Inf, y = Inf, label = stats_text, 
               hjust = 1.1, vjust = 1.1, size = geom_defaults$text_size, color = annotation_text_color,
               family = "sans") +
      # Add peak annotations using ggrepel for smart label placement
      ggrepel::geom_text_repel(
        data = density_stats,
        mapping = ggplot2::aes(
          x = peak_x, 
          y = peak_y, 
          label = sprintf("%.2f", peak_x),
          color = region
        ),
        nudge_y = 0.05,
        segment.color = "grey50",
        fontface = "bold",
        bg.color = "white",
        bg.r = 0.15,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    
    # Prepare metadata
    metadata <- list(
      total_points = nrow(density_data),
      regions_included = unique_regions,
      density_statistics = density_stats,
      color_scheme = region_colors,
      plot_type = "region_density",
      data_source = "gene_level_frequencies"
    )
    
    log_message("Region density plot generated successfully")
    
    return(list(
      plot = density_plot,
      data = density_data,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Region density plot generation failed: %s", e$message), level = "error")
    return(list(
      plot = ggplot2::ggplot() + 
        ggplot2::ggtitle("Region Density Plot Generation Failed") +
        ggplot2::theme_void(),
      data = data,
      metadata = list(error = e$message)
    ))
  })
}
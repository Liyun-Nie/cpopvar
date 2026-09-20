#' Generate distribution histogram
#' @param data Data frame with frequency data
#' @param config Configuration object
#' @param task_params Task parameters
#' @return List containing plot, data, and metadata
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_histogram labs scale_fill_manual theme_minimal theme element_text ggtitle expansion facet_wrap
#' @export
generate_distribution_histogram <- function(data, config, task_params) {
  log_message("Generating distribution histogram")
  
  tryCatch({
    # Get alpha scheme
    alpha_scheme <- get_alpha_scheme()
    
    # Use frequency data directly from normalization module
    if (!"frequency_per_kb" %in% names(data)) {
      stop("frequency_per_kb column not found in normalized data - check data source")
    }
    
    # Data is already filtered by orchestrator (NA and 0 values removed)
    hist_data <- data
    
    if (nrow(hist_data) == 0) {
      warning("No frequency data available for histogram")
      return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("No Data Available"), data = data, metadata = list()))
    }
    
    # Set appropriate labels for frequency data
    x_label <- "Frequency per Kilobase"
    plot_title <- "Variant Frequency Distribution"
    
    # Species ordering already applied by orchestrator
    # Data comes pre-ordered and ready for visualization
    
    # Get plot colors from configuration with fallback
    plot_colors <- tryCatch({
      get_color_palette(config, "distribution_histogram")
    }, error = function(e) {
      log_message(sprintf("Failed to get color palette: %s. Using default color.", e$message), level = "warning")
      c("#3498db")  # Default blue color
    })
    
    # Ensure plot_colors is valid
    if (is.null(plot_colors) || length(plot_colors) == 0) {
      plot_colors <- c("#3498db")  # Default blue color
    }
    
    # Create histogram plot for frequency data
    # Calculate reasonable x-axis limits based on data distribution
    freq_quantiles <- quantile(hist_data$frequency_per_kb, c(0.01, 0.99), na.rm = TRUE)
    x_max <- max(freq_quantiles[2], max(hist_data$frequency_per_kb, na.rm = TRUE))
    x_min <- min(0, freq_quantiles[1])
    
    # Get histogram bins from configuration
    histogram_bins <- get_task_parameter(task_params, config, "histogram_bins", default = 30)
    
    # Get histogram border color from configuration
    histogram_border_color <- get_task_parameter(task_params, config, "histogram_border_color", default = "white")
    
    # ============================================================
    # Data validation and NA cleanup
    # ============================================================
    # Ensure data passed to ggplot is clean and valid to prevent aesthetic mapping errors
    
    # Remove rows with NA or invalid frequency_per_kb values
    original_rows <- nrow(hist_data)
    hist_data <- hist_data[!is.na(hist_data$frequency_per_kb) & is.finite(hist_data$frequency_per_kb), ]
    cleaned_rows <- nrow(hist_data)
    
    if (cleaned_rows != original_rows) {
      log_message(sprintf("Data cleanup: Removed %d rows with NA/invalid frequency values", 
                         original_rows - cleaned_rows))
    }
    
    # Additional safety check for completely empty or problematic data
    if (nrow(hist_data) == 0 || all(hist_data$frequency_per_kb == 0)) {
      warning("No valid frequency data remaining after cleanup")
      return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("No Valid Data Available"), 
                  data = hist_data, metadata = list(error = "No valid data")))
    }
    # ============================================================
    
    # Enhanced data validation for ggplot aesthetics
    valid_data <- hist_data[!is.na(hist_data$frequency_per_kb) & 
                           is.finite(hist_data$frequency_per_kb) & 
                           hist_data$frequency_per_kb >= 0, ]
    
    if (nrow(valid_data) == 0) {
      return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("No Valid Data for Histogram"), 
                  data = hist_data, metadata = list(error = "No valid data for ggplot")))
    }
    
    # Simplify ggplot construction with explicit error handling
    hist_plot <- tryCatch({
      ggplot2::ggplot(valid_data, ggplot2::aes(x = frequency_per_kb)) +
        geom_histogram(bins = histogram_bins, 
                      fill = plot_colors[1], 
                      alpha = alpha_scheme$histogram_fill, 
                      color = histogram_border_color)
    }, error = function(e) {
      log_message(sprintf("Histogram construction failed: %s", e$message), level = "error")
      ggplot2::ggplot() + ggplot2::ggtitle("Histogram Construction Failed")
    })
    
    # Continue with the rest of the plot construction
    hist_plot <- hist_plot +
      scale_x_continuous(
        limits = c(x_min, x_max),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::labs(
        title = plot_title, 
        x = x_label,
        y = "Count"
      ) +
      # Get theme from configuration with error handling
      {
        theme_result <- tryCatch({
          theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"
          get_application_theme(theme_name)
        }, error = function(e) {
          log_message(sprintf("Failed to get application theme: %s. Using default theme.", e$message), level = "warning")
          theme_minimal()  # Fallback theme
        })
        theme_result
      } +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10)
      )
    
    # Add faceting if species grouping is requested
    if ("species" %in% names(valid_data) && length(unique(valid_data$species)) > 1) {
      # Get strip colors from configuration
      strip_fill_color <- get_task_parameter(task_params, config, "strip_fill_color", default = "gray95")
      strip_border_color <- get_task_parameter(task_params, config, "strip_border_color", default = "gray80")
      
      hist_plot <- hist_plot + 
        facet_wrap(~ species, scales = "free") +
        ggplot2::theme(
          strip.text = ggplot2::element_text(size = 8, face = "bold"),
          strip.background = ggplot2::element_rect(fill = strip_fill_color, color = strip_border_color)
        )
    }
    
    metadata <- list(
      total_variants = nrow(data),
      plotted_variants = nrow(valid_data),
      mean_frequency = mean(valid_data$frequency_per_kb, na.rm = TRUE),
      median_frequency = median(valid_data$frequency_per_kb, na.rm = TRUE)
    )
    
    return(list(plot = hist_plot, data = valid_data, metadata = metadata))
    
  }, error = function(e) {
    log_message(paste("Histogram generation failed:", e$message), level = "error")
    return(list(plot = ggplot2::ggplot() + ggplot2::ggtitle("Histogram Generation Failed"), data = data, metadata = list(error = e$message)))
  })
}
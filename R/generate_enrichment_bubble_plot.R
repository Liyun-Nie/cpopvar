#' Generate enrichment bubble plot
#' @param core_hotspots Core hotspots data
#' @param config Configuration object
#' @param task_params Task parameters
#' @return List containing plot, data, and metadata
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_gradient scale_x_continuous labs geom_vline annotate
#' @export
generate_enrichment_bubble_plot <- function(core_hotspots, config, task_params = NULL) {
  
  tryCatch({
    log_message("M03 Step 3: Generating functional enrichment bubble plot")
    
    # Get visual standards for consistent styling
    visual_standards <- get_visual_standards()
    geom_defaults <- visual_standards$ggplot2
    
    # Get structured color palette from configuration
    plot_colors <- get_color_palette(config, "dual_enrichment")
    
    # Get parameters
    max_categories <- get_task_parameter(task_params, config, "max_categories", 15)
    p_value_threshold <- get_task_parameter(task_params, config, "p_value_threshold", 0.05)
    min_enrichment_ratio <- get_task_parameter(task_params, config, "min_enrichment_ratio", 1.0)
    
    enrichment_data <- core_hotspots$enrichment
    
    if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
      log_message("No enrichment data available for bubble plot", level = "warning")
      return(NULL)
    }
    
    # Filter significant results and apply thresholds
    plot_data <- enrichment_data %>%
      dplyr::filter(
        p_adjusted <= p_value_threshold,
        enrichment_ratio >= min_enrichment_ratio,
        target_count >= 2  # Require at least 2 genes in category
      ) %>%
      head(max_categories) %>%
      dplyr::mutate(
        neg_log_p = -log10(p_adjusted),
        functional_category = factor(functional_category, levels = rev(functional_category))  # Reverse for correct ordering
      )
    
    if (nrow(plot_data) == 0) {
      log_message("No significant enrichment results meet the filtering criteria", level = "warning")
      # Create empty plot with message
      empty_plot <- ggplot2::ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = "No significant functional enrichment found", 
                size = 6, hjust = 0.5, vjust = 0.5) +
        # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
{ theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; get_application_theme(theme_name) } +
        ggplot2::labs(title = "Functional Enrichment Analysis",
             subtitle = "No significantly enriched functional categories")
      
      return(list(
        plot = empty_plot,
        data = data.frame(),
        metadata = list(message = "No significant results")
      ))
    }
    
    log_message(sprintf("Creating enrichment bubble plot with %d significant categories", nrow(plot_data)))
    
    # Create bubble plot
    p <-ggplot2::ggplot(plot_data,ggplot2::aes(x = enrichment_ratio, y = functional_category)) +
      geom_point(ggplot2::aes(size = target_count, color = neg_log_p), alpha = geom_defaults$point_alpha) +
      scale_size_continuous(
        name = "Gene Count",
        range = c(3, 12),
        breaks = c(2, 5, 10, 15),
        guide = guide_legend(override.aes = list(alpha = 1))
      ) +
      scale_color_gradient(
        name = "-log10(P adj)",
        low = plot_colors[["gradient"]][1],
        high = plot_colors[["gradient"]][2],
        guide = guide_colorbar(override.aes = list(alpha = 1))
      ) +
      scale_x_continuous(
        name = "Enrichment Ratio",
        limits = c(min_enrichment_ratio * 0.9, max(plot_data$enrichment_ratio) * 1.1)
      ) +
      ggplot2::labs(
        title = "Functional Enrichment of Core Shared Hotspot Genes",
        subtitle = "Bubble size = gene count, Color = significance",
        y = "Functional Group"
      ) +
      # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
{ theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; get_application_theme(theme_name) } +
      geom_vline(xintercept = 1, linetype = geom_defaults$line_ref_linetype, color = plot_colors[["ref_line_color"]], alpha = geom_defaults$line_ref_alpha)
    
    metadata <- list(
      significant_categories = nrow(plot_data),
      total_categories_tested = nrow(enrichment_data),
      max_enrichment_ratio = max(plot_data$enrichment_ratio),
      min_p_adjusted = min(plot_data$p_adjusted),
      parameters = list(
        max_categories = max_categories,
        p_value_threshold = p_value_threshold,
        min_enrichment_ratio = min_enrichment_ratio
      )
    )
    
    log_message("Functional enrichment bubble plot generated successfully")
    
    return(list(
      plot = p,
      data = plot_data,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate enrichment bubble plot: %s", e$message), level = "error")
    return(NULL)
  })
}

#' Generate dual enrichment plots
#' @param dual_enrichment_results Dual enrichment analysis results
#' @param config Configuration object
#' @param task_params Task parameters
#' @return List containing plots, data, and metadata
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_gradient labs guides guide_legend guide_colorbar
#' @export
generate_dual_enrichment_plots <- function(dual_enrichment_results,
                                          config,
                                          task_params = NULL) {
  
  tryCatch({
    log_message("M03 Step 3c: Generating dual enrichment analysis plots (core shared + private genes)")
    
    # Get structured color palette (V4 architecture)
    plot_colors <- get_color_palette(config, "dual_enrichment")
    
    # Get parameters
    max_categories <- get_task_parameter(task_params, config, "max_categories", 15)
    p_value_threshold <- get_task_parameter(task_params, config, "p_value_threshold", 0.05)
    min_enrichment_ratio <- get_task_parameter(task_params, config, "min_enrichment_ratio", 1.0)
    
    plots <- list()
    plot_data_list <- list()
    
    # Generate plot for core shared genes
    if ("core_shared_enrichment" %in% names(dual_enrichment_results)) {
      core_enrichment <- dual_enrichment_results$core_shared_enrichment
      
      if (nrow(core_enrichment) > 0) {
        # Filter and prepare data
        plot_data <- core_enrichment %>%
          dplyr::filter(
            p_adjusted <= p_value_threshold,
            enrichment_ratio >= min_enrichment_ratio
          ) %>%
          dplyr::arrange(p_adjusted, dplyr::desc(enrichment_ratio)) %>%
          head(max_categories)
        
        if (nrow(plot_data) > 0) {
          # Create bubble plot
          p_core <-ggplot2::ggplot(plot_data,ggplot2::aes(x = enrichment_ratio, y = reorder(functional_category, enrichment_ratio))) +
            geom_point(
            ggplot2::aes(size = target_count, color = -log10(p_adjusted)),
              alpha = 0.8
            ) +
            scale_size_continuous(
              name = "Gene Count",
              range = c(3, 15),
              breaks = pretty(plot_data$target_count, n = 4)
            ) +
            scale_color_gradient(
              name = "-log10(FDR)",
              low = plot_colors[["gradient"]][1], 
              high = plot_colors[["gradient"]][2],
              breaks = pretty(-log10(plot_data$p_adjusted), n = 4)
            ) +
            ggplot2::labs(
              title = "Functional Enrichment - Core Shared Hotspot Genes",
              subtitle = sprintf("Enriched in %d core shared genes (>=3 species)", 
                               dual_enrichment_results$metadata$core_shared_genes_count),
              x = "Enrichment Ratio",
              y = "Gene Category"
            ) +
            # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
{ theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; get_application_theme(theme_name) } +
            guides(
              size = guide_legend(override.aes = list(alpha = 1)),
              color = guide_colorbar(override.aes = list(alpha = 1))
            )
          
          plots$core_shared <- p_core
          plot_data_list$core_shared <- plot_data
          
          log_message(sprintf("Generated core shared enrichment plot with %d categories", nrow(plot_data)))
        } else {
          log_message("No significant enrichment found for core shared genes", level = "warning")
        }
      }
    }
    
    # Generate plot for private genes (independent analysis)
    if ("private_enrichment" %in% names(dual_enrichment_results)) {
      private_enrichment <- dual_enrichment_results$private_enrichment
      
      if (nrow(private_enrichment) > 0) {
        # Filter and prepare data
        plot_data <- private_enrichment %>%
          dplyr::filter(
            p_adjusted <= p_value_threshold,
            enrichment_ratio >= min_enrichment_ratio
          ) %>%
          dplyr::arrange(p_adjusted, dplyr::desc(enrichment_ratio)) %>%
          head(max_categories)
        
        if (nrow(plot_data) > 0) {
          # Create bubble plot for private genes
          p_private <- ggplot2::ggplot(plot_data, ggplot2::aes(x = enrichment_ratio, y = reorder(functional_category, enrichment_ratio))) +
            geom_point(
            ggplot2::aes(size = target_count, color = -log10(p_adjusted)),
              alpha = 0.8
            ) +
            scale_size_continuous(
              name = "Gene Count",
              range = c(3, 15),
              breaks = pretty(plot_data$target_count, n = 4)
            ) +
            scale_color_gradient(
              name = "-log10(FDR)",
              low = plot_colors[["gradient"]][1], 
              high = plot_colors[["gradient"]][2],
              breaks = pretty(-log10(plot_data$p_adjusted), n = 4)
            ) +
            ggplot2::labs(
              title = "Functional Enrichment - Private Hotspot Genes",
              subtitle = sprintf("Enriched in %d private genes (1 species only)", 
                               dual_enrichment_results$metadata$private_genes_count),
              x = "Enrichment Ratio",
              y = "Gene Category"
            ) +
            # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
{ theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; get_application_theme(theme_name) } +
            guides(
              size = guide_legend(override.aes = list(alpha = 1)),
              color = guide_colorbar(override.aes = list(alpha = 1))
            )
          
          plots$private <- p_private
          plot_data_list$private <- plot_data
          
          log_message(sprintf("Generated private enrichment plot with %d categories", nrow(plot_data)))
        } else {
          log_message("No significant enrichment found for private genes", level = "warning")
        }
      }
    }
    
    if (length(plots) == 0) {
      log_message("No significant functional enrichment found for either core shared or private genes", level = "warning")
      return(NULL)
    }
    
    metadata <- list(
      plots_generated = names(plots),
      total_plots = length(plots),
      parameters = list(
        max_categories = max_categories,
        p_value_threshold = p_value_threshold,
        min_enrichment_ratio = min_enrichment_ratio
      )
    )
    
    log_message(sprintf("Generated dual enrichment analysis plots: %s", paste(names(plots), collapse = ", ")))
    
    # Set the main plot as the core shared plot (or combined plot if gridExtra available)
    main_plot <- NULL
    if ("core_shared" %in% names(plots) && "private" %in% names(plots)) {
      # Try to create combined plot if both are available
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        main_plot <- gridExtra::arrangeGrob(
          plots$core_shared, plots$private,
          nrow = 1, ncol = 2,
          top = "Functional Enrichment Analysis: Core Shared vs Private Hotspot Genes"
        )
      } else {
        main_plot <- plots$core_shared  # Default to core shared if gridExtra not available
      }
    } else if ("core_shared" %in% names(plots)) {
      main_plot <- plots$core_shared
    } else if ("private" %in% names(plots)) {
      main_plot <- plots$private
    }
    
    return(list(
      plot = main_plot,  # Main plot for compatibility with calling code
      plots = plots,     # All plots (core_shared and/or private)
      data = plot_data_list,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate dual enrichment plots: %s", e$message), level = "error")
    return(NULL)
  })
}

# ============================================================
# M01 Relational Analysis Visualization Functions
# ============================================================
#
# This file contains visualization functions for M01-2: poiGS-Gene correlation.

#' Generate Species Correlation Heatmap
#' 
#' @description Creates heatmap showing correlation coefficients across species
#' @param correlation_results Data frame with correlation statistics
#' @param output_file Output file path
#' @param config Configuration object
#' @return Path to generated plot
#' @keywords internal
generate_correlation_heatmap <- function(correlation_results, output_file, config) {
  
  log_message("M01: Generating correlation coefficient heatmap")
  
  # Filter for species with sufficient data (aligned with M05)
  plot_data <- correlation_results %>%
    dplyr::filter(sufficient_data == TRUE)
  
  if (nrow(plot_data) == 0) {
    warning("generate_correlation_heatmap: No data available for plotting")
    return(NULL)
  }
  
  # Prepare data for heatmap (long format) - aligned with M05 data structure
  heatmap_data <- plot_data %>%
    dplyr::select(species, pearson_r, spearman_r, pearson_p, spearman_p) %>%
    tidyr::pivot_longer(
      cols = c(pearson_r, spearman_r),
      names_to = "method",
      values_to = "correlation"
    ) %>%
    dplyr::mutate(
      method = dplyr::case_when(
        method == "pearson_r" ~ "Pearson",
        method == "spearman_r" ~ "Spearman",
        TRUE ~ method
      )
    )
  
  # Add p-values for significance annotation
  heatmap_data <- heatmap_data %>%
    dplyr::left_join(
      plot_data %>%
        dplyr::select(species, pearson_p, spearman_p) %>%
        tidyr::pivot_longer(
          cols = c(pearson_p, spearman_p),
          names_to = "p_method",
          values_to = "p_value"
        ) %>%
        dplyr::mutate(
          method = dplyr::case_when(
            p_method == "pearson_p" ~ "Pearson",
            p_method == "spearman_p" ~ "Spearman",
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::select(species, method, p_value),
      by = c("species", "method")
    ) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      label_text = sprintf("%.2f%s", correlation, significance)
    )
  
  # Apply species ordering (use apply_species_ordering helper)
  heatmap_data_ordered <- apply_species_ordering(heatmap_data, config)
  
  # Reverse species order for heatmap (top-to-bottom display)
  if (is.factor(heatmap_data_ordered$species)) {
    heatmap_data_ordered$species <- factor(
      heatmap_data_ordered$species,
      levels = rev(levels(heatmap_data_ordered$species))
    )
  }
  
  # Create heatmap (aligned with M05 style)
  p <- ggplot2::ggplot(heatmap_data_ordered, ggplot2::aes(x = method, y = species, fill = correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label_text), size = 3, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#FFFFBF",
      high = "#D73027",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation\nCoefficient"
    ) +
    ggplot2::labs(
      title = "poiGS-Gene Correlation Coefficients Across Species (Neither-Zero Pairs Only)",
      subtitle = "Within-species correlation between IGS and flanking gene frequencies\n(*, **, *** indicate p < 0.05, 0.01, 0.001)",
      x = "Correlation Method",
      y = "Species"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "italic"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
  
  # Save plot
  ggplot2::ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = max(6, nrow(plot_data) * 0.3),
    dpi = 300,
    device = cpopvar_pdf_device
  )
  
  log_message(sprintf("Saved correlation heatmap: %s", basename(output_file)))
  
  return(output_file)
}


#' Generate poiGS-Gene Regression Plot
#' 
#' @description Creates global regression plot for poiGS vs gene frequencies
#' @param merged_data Merged poiGS-Gene data
#' @param correlation_results Correlation statistics
#' @param output_file Output file path
#' @param config Configuration object
#' @return List with plot file and regression stats
#' @keywords internal
generate_poigs_gene_regression_plot <- function(merged_data, correlation_results, output_file, config) {
  
  log_message("M01: Generating global poiGS-Gene regression plot")
  
  # Filter to neither_zero pattern
  plot_data <- merged_data %>%
    dplyr::filter(zero_pattern == "neither_zero")
  
  if (nrow(plot_data) == 0) {
    warning("No data available for regression plot")
    return(NULL)
  }
  
  # Calculate overall correlation (across all species) - aligned with M05
  overall_cor_pearson <- stats::cor.test(plot_data$poiGS_freq, plot_data$gene_freq_avg, 
                                          method = "pearson")
  overall_cor_spearman <- stats::cor.test(plot_data$poiGS_freq, plot_data$gene_freq_avg, 
                                           method = "spearman")
  
  # Calculate linear regression: poiGS_freq = slope * gene_freq_avg + intercept
  lm_fit <- stats::lm(poiGS_freq ~ gene_freq_avg, data = plot_data)
  slope <- stats::coef(lm_fit)[2]
  intercept <- stats::coef(lm_fit)[1]
  r_squared <- summary(lm_fit)$r.squared
  
  # Create annotation text with Overall statistics (aligned with M05)
  anno_text <- sprintf(
    "Overall Pearson r = %.3f (p %s)\nOverall Spearman rho = %.3f (p %s)\nRegression: poiGS = %.3f x gene_avg %s %.3f\nR^2 = %.3f\nn = %s poiGS-Gene pairs from %d species",
    overall_cor_pearson$estimate,
    ifelse(overall_cor_pearson$p.value < 0.001, "< 0.001", 
           sprintf("= %.3f", overall_cor_pearson$p.value)),
    overall_cor_spearman$estimate,
    ifelse(overall_cor_spearman$p.value < 0.001, "< 0.001", 
           sprintf("= %.3f", overall_cor_spearman$p.value)),
    slope,
    ifelse(intercept >= 0, "+", ""),
    intercept,
    r_squared,
    format(nrow(plot_data), big.mark = ","),
    length(unique(plot_data$species))
  )
  
  # Create scatter plot (aligned with M05: X=gene, Y=poiGS)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = gene_freq_avg, y = poiGS_freq)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.5, color = "#4DBBD5FF") +
    ggplot2::geom_smooth(method = "lm", color = "#E64B35FF", fill = "#E64B3550", 
                        linewidth = 1.2, alpha = 0.2) +
    ggplot2::annotate("text", x = Inf, y = Inf, 
                     label = anno_text,
                     hjust = 1.05, vjust = 1.1,
                     size = 3.5, color = "gray20",
                     fontface = "italic") +
    ggplot2::labs(
      title = "Global Correlation: poiGS vs Flanking Gene Variation Frequencies (Neither-Zero Pairs Only)",
      subtitle = "Linear regression fit with 95% confidence interval (all species pooled)",
      x = "Gene Frequency (variants/kb)",
      y = "poiGS Frequency (variants/kb)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  ggplot2::ggsave(
    output_file, p,
    width = 10, height = 8, dpi = 300,
    device = cpopvar_pdf_device
  )
  log_message(sprintf("Saved regression plot: %s", basename(output_file)))
  log_message(sprintf("  - Overall Pearson r = %.3f (p = %.2e)", 
                     overall_cor_pearson$estimate, overall_cor_pearson$p.value))
  log_message(sprintf("  - Overall Spearman rho = %.3f (p = %.2e)", 
                     overall_cor_spearman$estimate, overall_cor_spearman$p.value))
  log_message(sprintf("  - Regression equation: poiGS = %.3f x gene_avg %s %.3f (R^2 = %.3f)",
                     slope, ifelse(intercept >= 0, "+", ""), intercept, r_squared))
  
  return(list(
    plot_file = output_file,
    regression_params = list(
      slope = slope,
      intercept = intercept,
      r_squared = r_squared,
      pearson_r = overall_cor_pearson$estimate,
      pearson_p = overall_cor_pearson$p.value,
      spearman_r = overall_cor_spearman$estimate,
      spearman_p = overall_cor_spearman$p.value
    )
  ))
}


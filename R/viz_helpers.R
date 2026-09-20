############################################################
#### visualization_utils.R - Visualization Utilities ####
############################################################
#
# Provides themes, color schemes, and other helpers for
# creating plots. Migrated from src/utils and made compliant.
#
############################################################

# Load required libraries

#' Professional scientific publication theme (mytheme)
#' @param x_angle Angle for x-axis text labels
#' @param legend_position Position of the legend
#' @return ggplot2 theme object
#' @export
mytheme <- function(x_angle = 0, legend_position = "right") {
  theme_classic() +
    ggplot2::theme(
      text =ggplot2::element_text(colour = "black", size = 8),
      legend.text =ggplot2::element_text(colour = "black", size = 8),
      legend.title =ggplot2::element_text(colour = "black", size = 10),
      axis.line = element_line(size = 0.4, colour = "black"),
      axis.ticks = element_line(size = 0.4, colour = "black"),
      axis.text.x =ggplot2::element_text(angle = x_angle, hjust = if(x_angle > 0) 1 else 0.5, vjust = 0.5)
    )
}

#' Alternative theme for special cases (mytheme2)
mytheme2 <- function() {
  theme_minimal() +
    ggplot2::theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

#' Get standardized color palette for region types
get_region_type_colors <- function() {
  c("CDS" = "#1f77b4", "IGS" = "#ff7f0e", "intron" = "#2ca02c")
}

#' Get standardized color palette for variant types
get_variant_type_colors <- function() {
  c("SNP" = "#1f77b4", "INDEL" = "#ff7f0e", "complex" = "#2ca02c", "MNP" = "#9467bd")
}

#' Get frequency-based color palette
get_frequency_colors <- function() {
  c("CDS_frequency" = "#1f77b4", "intergenic_frequency" = "#ff7f0e", "intron_frequency" = "#2ca02c")
}

#' Create heatmap color gradient
#' @param type Type of heatmap colors to generate
#' @return Vector of color values
#' @export
get_heatmap_colors <- function(type = "hotspot") {
  if (type == "hotspot") return(c("white", "yellow", "orange", "red"))
  return(c("blue", "white", "red"))
}

#' Create explanatory subtitle for different plot types (Refactored for Compliance)
#' Rule #6 Compliant: No longer takes config object. Receives a simple list of parameters.
#' @param plot_type Type of plot for which to create subtitle
#' @param params_list List of parameters for subtitle generation
#' @return Character string with explanatory subtitle
#' @export
create_explanatory_subtitle <- function(plot_type, params_list = list()) {
    default_params <- list(
        frequency_unit = 1000,
        significance_level = 0.05,
        multiple_correction = "BH",
        correlation_method = "spearman",
        upset_min_degree = 3,
        outlier_method = "IQR",
        stat_method = "wilcox.test",
        version = "unknown"
    )
    params <- c(params_list, default_params[!names(default_params) %in% names(params_list)])

    switch(plot_type,
        "heatmap" = paste("Normalized frequency per", params$frequency_unit, "bp"),
        "hypergeometric" = paste("Significance level:", params$significance_level, "; Correction:", params$multiple_correction),
        "correlation" = paste("Correlation method:", params$correlation_method),
        "upset" = paste("Minimum degree:", params$upset_min_degree, "species"),
        "distribution" = paste("Outlier detection method:", params$outlier_method),
        "boxplot" = paste("Stat method:", params$stat_method, "; Sig. level:", params$significance_level),
        "stacked_bar" = "Variant count distribution by category and species",
        paste("Analysis based on", params$version, "workflow configuration")
  )
}

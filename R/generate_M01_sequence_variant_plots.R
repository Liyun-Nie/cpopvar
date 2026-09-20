# ============================================================
# M01 Sequence Variant Characterization Visualization Functions
# ============================================================
#
# This file contains visualization functions for M01-3.

# ============================================================
# Part 1: M05c - CDS SNP Annotation Plots
# ============================================================

#' Generate S/N Ratio Barplot
#' @keywords internal
generate_sn_ratio_barplot <- function(species_stats, output_file, config) {
  
  if (nrow(species_stats) == 0) {
    warning("No species statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering from file (same logic as M05)
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    species_stats <- apply_species_ordering(species_stats, 
                                             species_order_vector = species_order, 
                                             species_col = "species")
  }
  
  # Reverse factor levels for coord_flip() to maintain correct order
  if (is.factor(species_stats$species)) {
    species_stats$species <- factor(species_stats$species, 
                                     levels = rev(levels(species_stats$species)))
  }
  
  p <- ggplot2::ggplot(species_stats, ggplot2::aes(
    x = species, 
    y = sn_ratio,
    fill = sn_ratio
  )) +
    ggplot2::geom_col() +
    # Add significance stars next to the bars
    ggplot2::geom_text(
      ggplot2::aes(label = significance, y = sn_ratio), 
      hjust = -0.3, 
      vjust = 0.5, 
      color = "black", 
      size = 4
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 1,
      name = "S/N Ratio"
    ) +
    ggplot2::labs(
      title = "Synonymous/Nonsynonymous (S/N) Ratio by Species",
      subtitle = "Dashed line: neutral expectation (S/N = 1); Stars: binomial test significance",
      x = "Species",
      y = "S/N Ratio"
    ) +
    # ADDED: Significance legend annotation (aligned with M05)
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      label = "Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05",
      hjust = -0.05, vjust = 1.5,
      size = 3, color = "gray30"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 10, height = 8, device = "pdf")
  
  return(invisible(NULL))
}

#' Generate Syn vs Nonsyn Scatter
#' @keywords internal
generate_syn_vs_nonsyn_scatter <- function(species_stats, output_file, config) {
  
  if (nrow(species_stats) == 0) {
    warning("No species statistics to plot for scatter")
    return(invisible(NULL))
  }
  
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    species_stats <- apply_species_ordering(species_stats, species_order, "species")
  }
  
  # Use the M05 source color; a later color-system audit should route this through get_color_palette()
  p <- ggplot2::ggplot(species_stats, ggplot2::aes(
    x = synonymous, 
    y = nonsynonymous,
    label = species
  )) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = "#1F77B4") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    ggplot2::labs(
      title = "Synonymous vs Nonsynonymous Sites",
      subtitle = "Points above line: More nonsynonymous; Points below line: More synonymous",
      x = "Number of Synonymous Sites",
      y = "Number of Nonsynonymous Sites"
    ) +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 8, height = 7, device = "pdf")
  
  return(invisible(NULL))
}

#' Generate Site Type Composition Barplot
#' @keywords internal
generate_site_type_composition_barplot <- function(site_classification, output_file, config) {
  
  if (nrow(site_classification) == 0) {
    warning("No site classification data to plot")
    return(invisible(NULL))
  }
  
  # Aggregate by species and site_type
  plot_data <- site_classification %>%
    dplyr::group_by(species, site_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      site_type = factor(site_type, 
                        levels = c("synonymous", "nonsynonymous", "mixed"))
    )
  
  # Load species ordering from file (same logic as M05)
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    plot_data <- apply_species_ordering(plot_data, 
                                         species_order_vector = species_order, 
                                         species_col = "species")
  }
  
  # Reverse factor levels for coord_flip() to maintain correct order
  if (is.factor(plot_data$species)) {
    plot_data$species <- factor(plot_data$species, 
                                 levels = rev(levels(plot_data$species)))
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = species,
    y = count,
    fill = site_type
  )) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(
      values = c("synonymous" = "#4CAF50", 
                 "nonsynonymous" = "#F44336", 
                 "mixed" = "#FFC107"),
      labels = c("Synonymous", "Nonsynonymous", "Mixed"),
      name = "Site Type"
    ) +
    ggplot2::labs(
      title = "CDS SNP Site Type Composition Across Species",
      subtitle = "Stacked bar chart showing synonymous vs. nonsynonymous site distribution",
      x = "Species",
      y = "Number of Sites"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 10, height = 8, device = "pdf")
  
  return(invisible(NULL))
}

# ============================================================
# Part 2: M05d - IGS Variant Features Plots
# ============================================================

#' Generate SNP/INDEL Count Barplot
#' @keywords internal
generate_snp_indel_count_barplot <- function(species_stats, output_file, config) {
  
  plot_data <- species_stats %>%
    tidyr::pivot_longer(cols = c(snp, INDEL), names_to = "var_type", values_to = "count")
  
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    plot_data$species <- factor(plot_data$species, levels = species_order)
  }
  
  # Get unified colors for variant types
  var_type_colors <- get_color_palette(config, "var_type")
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = species, y = count, fill = var_type)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = var_type_colors) +
    ggplot2::labs(
      title = "IGS SNP and INDEL Site Counts",
      x = "Species",
      y = "Number of Sites",
      fill = "Variant Type"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  ggplot2::ggsave(output_file, p, width = 12, height = 6, dpi = 300)
  
  return(output_file)
}

#' Generate SNP/INDEL Ratio Lineplot
#' @keywords internal
generate_snp_indel_ratio_lineplot <- function(species_stats, output_file, config) {
  
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    species_stats$species <- factor(species_stats$species, levels = species_order)
  }
  
  # Use the M05 source colors; a later color-system audit should route these through get_color_palette()
  p <- ggplot2::ggplot(species_stats, ggplot2::aes(x = species, y = snp_indel_ratio, group = 1)) +
    ggplot2::geom_line(color = "#4DBBD5FF", size = 1) +
    ggplot2::geom_point(color = "#E64B35FF", size = 3) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = "IGS SNP/INDEL Ratio by Species",
      x = "Species",
      y = "SNP/INDEL Ratio"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  
  ggplot2::ggsave(output_file, p, width = 12, height = 6, dpi = 300)
  
  return(output_file)
}

#' Generate INDEL Length Boxplot
#' @keywords internal
generate_indel_length_boxplot <- function(indel_data, output_file, config, min_indel_length) {
  
  filtered_data <- indel_data %>% dplyr::filter(indel_length >= min_indel_length)
  
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    filtered_data$species <- factor(filtered_data$species, levels = species_order)
  }
  
  # Get unified colors
  boxplot_color <- get_color_palette(config, "primary")["fill"]
  
  p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = species, y = indel_length)) +
    ggplot2::geom_boxplot(fill = boxplot_color, alpha = 0.7) +
    ggplot2::labs(
      title = sprintf("INDEL Length Distribution (length >= %d)", min_indel_length),
      x = "Species",
      y = "INDEL Length (bp)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  ggplot2::ggsave(output_file, p, width = 12, height = 6, dpi = 300)
  
  return(output_file)
}

#' Generate INDEL Length Histogram
#' @keywords internal
generate_indel_length_histogram <- function(indel_data, output_file, config, min_indel_length) {
  
  filtered_data <- indel_data %>% dplyr::filter(indel_length >= min_indel_length)
  
  # Get unified colors
  hist_color <- get_color_palette(config, "primary")["fill"]
  
  p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = indel_length)) +
    ggplot2::geom_histogram(bins = 30, fill = hist_color, alpha = 0.7) +
    ggplot2::labs(
      title = sprintf("INDEL Length Distribution (length >= %d)", min_indel_length),
      x = "INDEL Length (bp)",
      y = "Count"
    ) +
    ggplot2::theme_classic()
  
  ggplot2::ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  
  return(output_file)
}

#' Generate INDEL Summary Comparison Plot
#' @keywords internal
generate_indel_summary_comparison_plot <- function(species_stats_filtered, output_file, config, min_indel_length) {
  
  if (missing(min_indel_length)) {
    min_indel_length <- 2
  }
  
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    species_stats_filtered <- apply_species_ordering(species_stats_filtered, species_order, "species")
  }
  
  plot_data <- species_stats_filtered %>%
    dplyr::select(species, n_indels, mean_length, median_length, max_length) %>%
    tidyr::pivot_longer(cols = c(n_indels, mean_length, median_length, max_length),
                       names_to = "metric", values_to = "value")
  
  metric_labels <- c(
    "n_indels" = sprintf("Number of INDELs (length >= %d)", min_indel_length),
    "mean_length" = "Mean Length (bp)",
    "median_length" = "Median Length (bp)",
    "max_length" = "Max Length (bp)"
  )
  
  plot_data$metric <- factor(plot_data$metric, levels = names(metric_labels))
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = species, y = value)) +
    ggplot2::geom_col(fill = "#4DBBD5FF", alpha = 0.7) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 2, 
                       labeller = ggplot2::as_labeller(metric_labels)) +
    ggplot2::labs(title = "INDEL Summary Statistics by Species", x = "Species", y = "Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  ggplot2::ggsave(output_file, p, width = 14, height = 10, dpi = 300)
  
  return(output_file)
}

#' Generate Shared INDEL Upset Plot
#' @keywords internal
generate_shared_indel_upset_plot <- function(upset_matrix, output_file, config) {
  
  if (nrow(upset_matrix) == 0 || ncol(upset_matrix) == 0) {
    warning("Empty UpSet matrix")
    return(NULL)
  }
  
  # Reorder columns by species order
  session_paths <- get_session_paths(config$session_info$session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    available_species <- intersect(species_order, colnames(upset_matrix))
    upset_matrix <- upset_matrix[, available_species, drop = FALSE]
  }
  
  upset_df <- as.data.frame(upset_matrix)
  
  pdf(output_file, width = 14, height = 8)
  print(UpSetR::upset(upset_df, 
                       nsets = ncol(upset_df),
                       order.by = "freq",
                       sets.bar.color = "#4575B4",
                       mainbar.y.label = "Shared INDEL Sequences",
                       sets.x.label = "Total Sequences per Species",
                       sets = rev(colnames(upset_df)),
                       keep.order = TRUE))
  dev.off()
  
  return(output_file)
}

#' Generate Species Clustering Dendrogram
#' @keywords internal
generate_species_clustering_dendrogram <- function(hc, output_file, config) {
  
  if (is.null(hc)) {
    warning("No clustering object provided")
    return(NULL)
  }
  
  pdf(output_file, width = 10, height = 8)
  plot(hc, 
       main = "Species Clustering Based on Shared INDEL Sequences (Jaccard Similarity)",
       sub = "UPGMA clustering using Jaccard distance",
       xlab = "Species",
       ylab = "Jaccard Distance (1 - Similarity)",
       cex = 0.8)
  dev.off()
  
  return(output_file)
}

#' Generate INDEL Length Histogram (Filtered)
#' @keywords internal
generate_indel_length_histogram_filtered <- function(indel_data, output_file, config, min_indel_length) {
  
  # Filter data
  indel_data_filtered <- indel_data %>%
    dplyr::filter(indel_length >= min_indel_length)
  
  if (nrow(indel_data_filtered) == 0) {
    warning("No INDEL data after filtering (length >= ", min_indel_length, ")")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    indel_data_filtered <- apply_species_ordering(indel_data_filtered, 
                                                   species_order_vector = species_order, 
                                                   species_col = "species")
  }
  
  p <- ggplot2::ggplot(indel_data_filtered, ggplot2::aes(x = indel_length, fill = indel_type)) +
    ggplot2::geom_histogram(binwidth = 1, alpha = 0.7, position = "stack") +
    ggplot2::scale_fill_manual(
      values = c("insertion" = "#4575B4", "deletion" = "#D73027"),
      name = "INDEL Type"
    ) +
    ggplot2::facet_wrap(~species, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = sprintf("INDEL Length Distribution (Filtered: length >= %d bp)", min_indel_length),
      subtitle = "Structural variants only - focusing on large INDELs",
      x = "INDEL Length (bp)",
      y = "Count"
    ) +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 14, height = 12, device = "pdf")
  
  return(invisible(NULL))
}


#' Generate INDEL Length Boxplot Comparison
#'
#' Creates side-by-side boxplots comparing INDEL length distributions for
#' all INDELs vs. filtered INDELs.
#'
#' @param species_stats_all Data frame with species-level statistics (all INDELs)
#' @param species_stats_filtered Data frame with species-level statistics (filtered INDELs)
#' @param output_file Path for output PDF file
#' @param config Configuration object
#' @param min_indel_length Minimum INDEL length threshold
#'
#' @return NULL (saves plot to file)
#'
#' @keywords internal
generate_indel_length_boxplot_comparison <- function(species_stats_all, species_stats_filtered, 
                                                      output_file, config, min_indel_length) {
  
  if (nrow(species_stats_all) == 0) {
    warning("No species statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
  }
  
  # Prepare data for comparison
  data_all <- species_stats_all %>%
    dplyr::select(species, mean_length, median_length) %>%
    dplyr::mutate(dataset = "All INDELs")
  
  if (nrow(species_stats_filtered) > 0) {
    data_filtered <- species_stats_filtered %>%
      dplyr::select(species, mean_length, median_length) %>%
      dplyr::mutate(dataset = sprintf("Filtered (>=%d bp)", min_indel_length))
    
    plot_data <- dplyr::bind_rows(data_all, data_filtered)
  } else {
    plot_data <- data_all
    warning("No filtered data available for comparison")
  }
  
  # Apply species ordering
  if (!is.null(species_order)) {
    plot_data <- apply_species_ordering(plot_data, 
                                        species_order_vector = species_order, 
                                        species_col = "species")
    
    # Reverse factor levels for coord_flip()
    if (is.factor(plot_data$species)) {
      plot_data$species <- factor(plot_data$species, 
                                   levels = rev(levels(plot_data$species)))
    }
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = species,
    y = mean_length,
    fill = dataset
  )) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(
      values = c("All INDELs" = "#ABDDA4", "Filtered (>=2 bp)" = "#2B83BA"),
      name = "Dataset"
    ) +
    ggplot2::labs(
      title = "Mean INDEL Length Comparison",
      subtitle = "All INDELs vs. Filtered structural variants",
      x = "Species",
      y = "Mean INDEL Length (bp)"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 10, height = 8, device = "pdf")
  
  return(invisible(NULL))
}

# ============================================================
# Part 3: M01-3e - Genome Region Statistics Plots (NEW)
# ============================================================

#' Generate S/N Ratio Barplot by Genome Region
#' @description Create grouped barplot showing S/N ratio across genome regions (LSC/SSC/IRA)
#' @param region_stats Data frame with columns: species, genome_region, sn_ratio
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_sn_ratio_barplot <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    region_stats <- apply_species_ordering(region_stats, 
                                           species_order_vector = species_order, 
                                           species_col = "species")
    
    # Reverse for coord_flip()
    if (is.factor(region_stats$species)) {
      region_stats$species <- factor(region_stats$species, 
                                     levels = rev(levels(region_stats$species)))
    }
  }
  
  # Ensure genome_region is a factor with correct order
  region_stats$genome_region <- factor(region_stats$genome_region, 
                                       levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Order species
  region_stats <- apply_species_ordering(region_stats, config$species_order)
  
  # Generate the bar plot with fixed bin width
  p <- ggplot2::ggplot(
    region_stats, 
    ggplot2::aes(
      x = species, 
      y = sn_ratio, 
      fill = genome_region
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8, preserve = "single"), 
      width = 0.7
    ) +
    ggplot2::scale_fill_manual(
      values = get_color_palette(config, "genome_region"),
      name = "Genome Region",
      labels = c("LSC" = "LSC", "IRB" = "IRB", "SSC" = "SSC", "IRA" = "IRA")
    ) +
    ggplot2::labs(
      title = "Synonymous/Nonsynonymous Ratio by Genome Region",
      subtitle = "CDS SNPs grouped by LSC/SSC/IR regions (* p<0.05, ** p<0.01, *** p<0.001)",
      x = "Species",
      y = "S/N Ratio (Synonymous / Nonsynonymous)"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.7) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  # Add significance markers
  if ("significance" %in% names(region_stats)) {
    # Calculate y positions for significance labels (slightly above the bars)
    max_ratio <- max(region_stats$sn_ratio, na.rm = TRUE)
    label_y <- max_ratio * 1.05
    
    # Filter for significant results
    sig_data <- region_stats %>%
      dplyr::filter(significance != "ns")
    
    if (nrow(sig_data) > 0) {
      p <- p + ggplot2::geom_text(
        data = sig_data,
        ggplot2::aes(x = species, y = label_y, label = significance, group = genome_region),
        size = 4,
        hjust = 0,
        position = ggplot2::position_dodge(width = 0.9),
        inherit.aes = FALSE
      )
    }
  }
  
  ggplot2::ggsave(output_file, p, width = 12, height = 10, device = "pdf")
  
  log_message(sprintf("Saved S/N ratio by genome region barplot to: %s", output_file))
  return(p)
}

#' Generate Site Type Composition Barplot by Genome Region
#' @description Create stacked barplot showing site type composition across genome regions
#' @param region_stats Data frame with columns: species, genome_region, synonymous, nonsynonymous, mixed
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_site_type_composition <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Reshape data for stacked bar plot
  plot_data <- region_stats %>%
    dplyr::select(species, genome_region, synonymous, nonsynonymous, mixed) %>%
    tidyr::pivot_longer(
      cols = c(synonymous, nonsynonymous, mixed),
      names_to = "site_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      site_type = factor(site_type, levels = c("synonymous", "nonsynonymous", "mixed"))
    )
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    plot_data <- apply_species_ordering(plot_data, 
                                        species_order_vector = species_order, 
                                        species_col = "species")
    
    # Reverse for coord_flip()
    if (is.factor(plot_data$species)) {
      plot_data$species <- factor(plot_data$species, 
                                  levels = rev(levels(plot_data$species)))
    }
  }
  
  # Ensure genome_region is a factor with correct order
  plot_data$genome_region <- factor(plot_data$genome_region, 
                                    levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Create faceted stacked bar plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = species,
    y = count,
    fill = site_type
  )) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::facet_wrap(~ genome_region, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c("synonymous" = "#4CAF50", 
                 "nonsynonymous" = "#F44336", 
                 "mixed" = "#FFC107"),
      labels = c("Synonymous", "Nonsynonymous", "Mixed"),
      name = "Site Type"
    ) +
    ggplot2::labs(
      title = "CDS SNP Site Type Composition by Genome Region",
      subtitle = "Stacked bar plot showing synonymous, nonsynonymous, and mixed sites",
      x = "Species",
      y = "Number of Sites"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 12, height = 14, device = "pdf")
  
  log_message(sprintf("Saved site type composition by genome region to: %s", output_file))
  return(p)
}

#' Generate Synonymous vs Nonsynonymous Scatter Plot by Genome Region
#' @description Create scatter plot showing syn vs nonsyn site counts, colored by genome region
#' @param region_stats Data frame with columns: species, genome_region, synonymous, nonsynonymous
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_syn_vs_nonsyn_scatter <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Ensure genome_region is a factor with correct order
  region_stats$genome_region <- factor(region_stats$genome_region, 
                                       levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Create scatter plot
  p <- ggplot2::ggplot(region_stats, ggplot2::aes(
    x = synonymous,
    y = nonsynonymous,
    color = genome_region,
    label = species
  )) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    ggplot2::scale_color_manual(
      values = get_color_palette(config, "genome_region"),
      name = "Genome Region",
      labels = c("LSC" = "LSC", "IRB" = "IRB", "SSC" = "SSC", "IRA" = "IRA")
    ) +
    ggplot2::labs(
      title = "Synonymous vs Nonsynonymous Sites by Genome Region",
      subtitle = "Each point represents a (species, genome_region) combination",
      x = "Number of Synonymous Sites",
      y = "Number of Nonsynonymous Sites"
    ) +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 12, height = 10, device = "pdf")
  
  log_message(sprintf("Saved syn vs nonsyn scatter by genome region to: %s", output_file))
  return(p)
}

#' Generate SNP/INDEL Ratio Plot by Genome Region
#' @description Create grouped barplot showing SNP/INDEL ratio across genome regions
#' @param region_stats Data frame with columns: species, genome_region, snp, INDEL, snp_indel_ratio
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_snp_indel_ratio_plot <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    region_stats <- apply_species_ordering(region_stats, 
                                           species_order_vector = species_order, 
                                           species_col = "species")
    
    # Reverse for coord_flip()
    if (is.factor(region_stats$species)) {
      region_stats$species <- factor(region_stats$species, 
                                     levels = rev(levels(region_stats$species)))
    }
  }
  
  # Ensure genome_region is a factor with correct order
  region_stats$genome_region <- factor(region_stats$genome_region, 
                                       levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Generate the bar plot with genome region colors and fixed bin width
  p <- ggplot2::ggplot(region_stats, ggplot2::aes(
    x = species,
    y = snp_indel_ratio,
    fill = genome_region
  )) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8, preserve = "single"), 
      width = 0.7
    ) +
    ggplot2::scale_fill_manual(
      values = get_color_palette(config, "genome_region"),
      name = "Genome Region",
      labels = c("LSC" = "LSC", "IRB" = "IRB", "SSC" = "SSC", "IRA" = "IRA")
    ) +
    ggplot2::labs(
      title = "SNP/INDEL Ratio by Genome Region",
      subtitle = "Comparison of variant ratios across genomic segments",
      x = "Species",
      y = "SNP/INDEL Ratio"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.7) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 12, height = 10, device = "pdf")
  
  log_message(sprintf("Saved SNP/INDEL ratio by genome region barplot to: %s", output_file))
  return(p)
}

#' Generate SNP/INDEL Count Stacked Barplot by Genome Region
#' @description Create stacked barplot showing SNP and INDEL counts across genome regions
#' @param region_stats Data frame with columns: species, genome_region, snp, INDEL
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_snp_indel_count_barplot <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Reshape data for stacked bar plot
  plot_data <- region_stats %>%
    dplyr::select(species, genome_region, snp, INDEL) %>%
    tidyr::pivot_longer(
      cols = c(snp, INDEL),
      names_to = "var_type",
      values_to = "count"
    )
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    plot_data <- apply_species_ordering(plot_data, 
                                        species_order_vector = species_order, 
                                        species_col = "species")
  }
  
  # Ensure genome_region is a factor with correct order
  plot_data$genome_region <- factor(plot_data$genome_region, 
                                    levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Get unified colors for variant types
  var_type_colors <- get_color_palette(config, "var_type")
  
  # Create faceted stacked bar plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = species,
    y = count,
    fill = var_type
  )) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::facet_wrap(~ genome_region, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = var_type_colors,
      labels = c("snp" = "SNP", "INDEL" = "INDEL"),
      name = "Variant Type"
    ) +
    ggplot2::labs(
      title = "SNP and INDEL Site Counts by Genome Region",
      subtitle = "Stacked bar plot showing variant counts across genomic segments",
      x = "Species",
      y = "Number of Sites"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  
  ggplot2::ggsave(output_file, p, width = 14, height = 12, device = "pdf")
  
  log_message(sprintf("Saved SNP/INDEL count barplot by genome region to: %s", output_file))
  return(p)
}

#' Generate Site Type Composition Plot by Genome Region
#' @description Create stacked barplot showing site type composition across genome regions
#' @param region_stats Data frame with columns: species, genome_region, synonymous, nonsynonymous, mixed
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_site_type_composition <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    region_stats <- apply_species_ordering(region_stats, 
                                           species_order_vector = species_order, 
                                           species_col = "species")
    
    # Reverse for coord_flip()
    if (is.factor(region_stats$species)) {
      region_stats$species <- factor(region_stats$species, 
                                     levels = rev(levels(region_stats$species)))
    }
  }
  
  # Reshape data for stacked bar plot
  plot_data <- region_stats %>%
    dplyr::select(species, genome_region, synonymous, nonsynonymous, mixed) %>%
    tidyr::pivot_longer(
      cols = c(synonymous, nonsynonymous, mixed),
      names_to = "site_type",
      values_to = "count"
    ) %>%
    dplyr::filter(count > 0) %>%  # Remove zero counts to avoid empty legend entries
    dplyr::mutate(
      site_type = factor(site_type, levels = c("synonymous", "nonsynonymous", "mixed")),
      genome_region = factor(genome_region, levels = c("LSC", "IRB", "SSC", "IRA"))
    )
  
  # Determine which site types are actually present in the data
  present_site_types <- unique(plot_data$site_type)
  
  # Get colors for site types (only for those present)
  all_site_type_colors <- c(
    "synonymous" = "#4CAF50",      # Green
    "nonsynonymous" = "#F44336",   # Red
    "mixed" = "#FFC107"            # Amber
  )
  site_type_colors <- all_site_type_colors[names(all_site_type_colors) %in% present_site_types]
  
  # Create labels only for present site types
  all_labels <- c("synonymous" = "Synonymous", 
                  "nonsynonymous" = "Nonsynonymous", 
                  "mixed" = "Mixed")
  site_type_labels <- all_labels[names(all_labels) %in% present_site_types]
  
  # Create faceted stacked barplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = species, y = count, fill = site_type)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::facet_wrap(~ genome_region, ncol = 4, scales = "free_x") +
    ggplot2::scale_fill_manual(
      values = site_type_colors,
      name = "Site Type",
      labels = site_type_labels,
      drop = TRUE  # Drop unused levels from legend
    ) +
    ggplot2::labs(
      title = "CDS SNP Site Type Composition by Genome Region",
      subtitle = "Distribution of synonymous, nonsynonymous, and mixed sites",
      x = "Species",
      y = "Number of Sites"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 16, height = 10, device = "pdf")
  
  log_message(sprintf("Saved site type composition by genome region to: %s", output_file))
  return(p)
}

#' Generate Synonymous vs Nonsynonymous Scatter Plot by Genome Region
#' @description Create scatter plot showing relationship between synonymous and nonsynonymous sites
#' @param region_stats Data frame with columns: species, genome_region, synonymous, nonsynonymous
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_syn_vs_nonsyn_scatter <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Ensure genome_region is a factor
  region_stats$genome_region <- factor(region_stats$genome_region, 
                                       levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Create scatter plot
  p <- ggplot2::ggplot(region_stats, ggplot2::aes(
    x = nonsynonymous, 
    y = synonymous, 
    color = genome_region,
    shape = genome_region
  )) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(
      values = get_color_palette(config, "genome_region"),
      name = "Genome Region"
    ) +
    ggplot2::scale_shape_manual(
      values = c("LSC" = 16, "IRB" = 17, "SSC" = 15, "IRA" = 18),
      name = "Genome Region"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = species),
      size = 2.5,
      max.overlaps = 20,
      box.padding = 0.5
    ) +
    ggplot2::labs(
      title = "Synonymous vs Nonsynonymous Sites by Genome Region",
      subtitle = "Dashed line represents S/N ratio = 1",
      x = "Nonsynonymous Sites",
      y = "Synonymous Sites"
    ) +
    get_application_theme(config) +
    ggplot2::theme(legend.position = "right")
  
  ggplot2::ggsave(output_file, p, width = 12, height = 10, device = "pdf")
  
  log_message(sprintf("Saved syn vs nonsyn scatter by genome region to: %s", output_file))
  return(p)
}

#' Generate SNP/INDEL Ratio Line Plot by Genome Region
#' @description Create line plot showing SNP/INDEL ratio trends across species
#' @param region_stats Data frame with columns: species, genome_region, snp_indel_ratio
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_snp_indel_ratio_lineplot <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    region_stats <- apply_species_ordering(region_stats, 
                                           species_order_vector = species_order, 
                                           species_col = "species")
  }
  
  # Ensure genome_region is a factor
  region_stats$genome_region <- factor(region_stats$genome_region, 
                                       levels = c("LSC", "IRB", "SSC", "IRA"))
  
  # Create line plot
  p <- ggplot2::ggplot(region_stats, ggplot2::aes(
    x = species, 
    y = snp_indel_ratio, 
    color = genome_region,
    group = genome_region
  )) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(
      values = get_color_palette(config, "genome_region"),
      name = "Genome Region"
    ) +
    ggplot2::labs(
      title = "SNP/INDEL Ratio Trends by Genome Region",
      subtitle = "Comparison of variant ratios across species",
      x = "Species",
      y = "SNP/INDEL Ratio"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 14, height = 8, device = "pdf")
  
  log_message(sprintf("Saved SNP/INDEL ratio lineplot by genome region to: %s", output_file))
  return(p)
}

#' Generate SNP/INDEL Count Stacked Barplot by Genome Region
#' @description Create stacked barplot showing SNP and INDEL counts across genome regions
#' @param region_stats Data frame with columns: species, genome_region, snp, INDEL
#' @param output_file Path to save the plot
#' @param config Configuration object
#' @keywords internal
generate_genome_region_snp_indel_count_barplot <- function(region_stats, output_file, config) {
  
  if (nrow(region_stats) == 0) {
    warning("No genome region statistics to plot")
    return(invisible(NULL))
  }
  
  # Load species ordering
  session_id <- config$session_info$session_id
  session_paths <- get_session_paths(session_id)
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    species_order <- species_order_data[[1]]
    
    # Apply species ordering
    region_stats <- apply_species_ordering(region_stats, 
                                           species_order_vector = species_order, 
                                           species_col = "species")
    
    # Reverse for coord_flip()
    if (is.factor(region_stats$species)) {
      region_stats$species <- factor(region_stats$species, 
                                     levels = rev(levels(region_stats$species)))
    }
  }
  
  # Reshape data for stacked bar plot
  plot_data <- region_stats %>%
    dplyr::select(species, genome_region, snp, INDEL) %>%
    tidyr::pivot_longer(
      cols = c(snp, INDEL),
      names_to = "var_type",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      var_type = factor(var_type, levels = c("snp", "INDEL")),
      genome_region = factor(genome_region, levels = c("LSC", "IRB", "SSC", "IRA"))
    )
  
  # Get colors from config (unified color system)
  var_type_colors <- get_color_palette(config, "var_type")
  
  # Create faceted stacked barplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = species, y = count, fill = var_type)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::facet_wrap(~ genome_region, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = var_type_colors,
      name = "Variant Type",
      labels = c("snp" = "SNP", "INDEL" = "INDEL")
    ) +
    ggplot2::labs(
      title = "SNP and INDEL Counts by Genome Region",
      subtitle = "Distribution of variant types across genomic segments",
      x = "Species",
      y = "Number of Variants"
    ) +
    ggplot2::coord_flip() +
    get_application_theme(config)
  
  ggplot2::ggsave(output_file, p, width = 16, height = 10, device = "pdf")
  
  log_message(sprintf("Saved SNP/INDEL count barplot by genome region to: %s", output_file))
  return(p)
}


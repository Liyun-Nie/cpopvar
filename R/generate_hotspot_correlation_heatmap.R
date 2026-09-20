#' Generate hotspot correlation heatmap
#' 
#' @title Generate hotspot correlation heatmap
#' @description Creates correlation heatmap showing relationships between species based on shared hotspots.
#' Supports both CDS-based (gene) and IGS-based (poiGS_ID) feature analysis.
#' @param candidate_hotspots Data containing identified hotspot candidates
#' @param normalized_data Normalized frequency dataset
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param feature_col Name of the feature column to use for grouping 
#'   (default: "gene" for M03, use "poiGS_ID" for M05)
#' @return List containing correlation heatmap and analysis results
#' @export
generate_hotspot_correlation_heatmap <- function(candidate_hotspots,
                                               normalized_data,
                                               config,
                                               task_params = NULL,
                                               feature_col = "gene") {
  
  tryCatch({
    log_message(sprintf("Hotspot Correlation: Analyzing species relationships based on core shared hotspots (%s)", feature_col))
    
    # Get visual standards for consistent styling
    visual_standards <- get_visual_standards()
    geom_defaults <- visual_standards$ggplot2
    
    # Get structured color palette from configuration
    plot_colors <- tryCatch({
      get_color_palette(config, "hotspot_correlation")
    }, error = function(e) {
      log_message(sprintf("Hotspot correlation colors not found, using defaults: %s", e$message), level = "warning")
      # Fallback colors for correlation heatmap
      list(
        gradient = c("#2166AC", "#FFFFBF", "#D73027")  # Blue-White-Red gradient for correlations
      )
    })
    
    if (is.null(candidate_hotspots) || nrow(candidate_hotspots) == 0) {
      log_message("No candidate hotspots available for correlation analysis", level = "warning")
      return(NULL)
    }
    
    # Get parameters from config, with a final fallback default.
    min_species_threshold <- get_task_parameter(task_params, config, "min_species_number", 3, log_default = TRUE)
    correlation_method <- get_task_parameter(task_params, config, "correlation_method", "spearman", log_default = TRUE)
    
    # Use parameterized feature column (gene for M03, poiGS_ID for M05)
    gene_col <- feature_col
    species_col <- "species"
    freq_col <- "frequency_per_kb"
    
    # Identify core shared genes (present in multiple species)
    gene_species_count <- candidate_hotspots %>%
      dplyr::group_by(.data[[gene_col]]) %>%
      dplyr::summarise(
        n_species = dplyr::n_distinct(.data[[species_col]]),
        species_list = paste(unique(.data[[species_col]]), collapse = ", "),
        avg_frequency = mean(.data[[freq_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_species),dplyr::desc(avg_frequency))
    
    core_genes <- gene_species_count %>%
      dplyr::filter(n_species >= min_species_threshold) %>%
      dplyr::pull(.data[[gene_col]])
    
    if (length(core_genes) == 0) {
      log_message(sprintf("No core %s found with threshold >= %d species", gsub("_ID$", "", feature_col), min_species_threshold), level = "warning")
      return(NULL)
    }
    
    log_message(sprintf("Identified %d core shared hotspot %s (>= %d species)", 
                       length(core_genes), gsub("_ID$", "", feature_col), min_species_threshold))
    
    # Filter to core genes and create gene-species matrix
    core_hotspots <- candidate_hotspots %>%
      dplyr::filter(.data[[gene_col]] %in% core_genes)
    
    # Create a complete matrix of core genes vs all species
    all_species <- unique(normalized_data[[species_col]])
    core_matrix <- matrix(NA, nrow = length(core_genes), ncol = length(all_species),
                         dimnames = list(core_genes, all_species))
    
    # Fill matrix with hotspot frequencies (NA if not a hotspot)
    for (i in 1:nrow(core_hotspots)) {
      gene <- core_hotspots[[gene_col]][i]
      species <- core_hotspots[[species_col]][i]
      freq <- core_hotspots[[freq_col]][i]
      core_matrix[gene, species] <- freq
    }
    
    # For species-species correlation, transpose matrix (species as rows, genes as columns)
    species_matrix <- t(core_matrix)
    
    # Remove species with no core hotspots
    species_with_data <- apply(species_matrix, 1, function(x) sum(!is.na(x)) > 0)
    species_matrix <- species_matrix[species_with_data, , drop = FALSE]
    
    if (nrow(species_matrix) < 3) {
      log_message("Insufficient species with core hotspots for correlation analysis", level = "warning")
      return(NULL)
    }
    
    log_message(sprintf("Core hotspot matrix: %d species x %d features (%s)", nrow(species_matrix), ncol(species_matrix), feature_col))
    
    # Calculate species-species correlation based on core hotspot patterns
    # Use presence/absence for correlation (convert NA to 0, others to 1)
    binary_matrix <- ifelse(is.na(species_matrix), 0, 1)
    
    # Validate binary matrix for correlation analysis
    if (any(apply(binary_matrix, 1, var) == 0)) {
      log_message("Some species have no variation in core hotspot presence, adjusting matrix", level = "warning")
      # Add small random noise to constant rows to allow correlation calculation
      constant_rows <- apply(binary_matrix, 1, var) == 0
      binary_matrix[constant_rows, ] <- binary_matrix[constant_rows, ] + 
        matrix(runif(sum(constant_rows) * ncol(binary_matrix), -1e-6, 1e-6), 
               nrow = sum(constant_rows))
    }
    
    # Calculate binary (presence/absence) correlation with robust error handling
    # For species-species correlation, we need to transpose (correlate rows, not columns)
    binary_correlation <- tryCatch({
      cor_matrix <- cor(t(binary_matrix), method = correlation_method)
      # Clean any remaining NaN/Inf values
      cor_matrix[is.na(cor_matrix) | is.infinite(cor_matrix)] <- 0
      cor_matrix
    }, error = function(e) {
      log_message(sprintf("Binary correlation calculation failed: %s", e$message), level = "warning")
      NULL
    })
    
    # REMOVED: Frequency-based correlation analysis (causes misleading results with sparse data)
    # As per user feedback: sparse hotspot gene matrices lead to zero variance and meaningless correlations
    # Only Binary (presence/absence) correlation analysis is retained for scientific accuracy
    
    # Calculate statistical significance for binary correlation only
    binary_pvalues <- NULL
    
    if (!is.null(binary_correlation)) {
      log_message("Calculating statistical significance for binary correlation (presence/absence only)")
      # Fix: Use binary_matrix directly (species as rows) to match correlation matrix dimensions
      binary_pvalues <- calculate_correlation_pvalues(binary_matrix, correlation_method)
    }
    
    # Create correlation heatmaps
    plots <- list()
    
    # 1. Binary presence/absence correlation heatmap
    if (!is.null(binary_correlation) && 
        !any(is.na(binary_correlation)) && 
        !any(is.infinite(binary_correlation)) &&
        nrow(binary_correlation) >= 2 && ncol(binary_correlation) >= 2) {
      
      # Create annotated correlation text with significance stars
      binary_annotations <- create_annotated_correlation_text(binary_correlation, binary_pvalues)
      
      binary_heatmap <- tryCatch({
        # Ensure plot_colors is a list and has gradient
        if (!is.list(plot_colors) || is.null(plot_colors[["gradient"]])) {
          log_message("Invalid plot_colors format, using fallback gradient", level = "warning")
          gradient_colors <- c("#2166AC", "#FFFFBF", "#D73027")
        } else {
          gradient_colors <- plot_colors[["gradient"]]
        }
        
        create_standardized_pheatmap(
          data_matrix = binary_correlation,
          title = sprintf("Species Correlation - Core Hotspot Presence/Absence\n(%s correlation, %d genes, *p<0.05, **p<0.01, ***p<0.001)", 
                         tools::toTitleCase(correlation_method), length(core_genes)),
          display_numbers = binary_annotations,
          color_palette = gradient_colors
        )
      }, error = function(e) {
        log_message(sprintf("pheatmap failed: %s. Creating ggplot2 fallback heatmap", e$message), level = "warning")
        # Create fallback ggplot2 heatmap with significance annotations
        heatmap_df <- expand.grid(Species1 = rownames(binary_correlation), 
                                 Species2 = colnames(binary_correlation))
        heatmap_df$Correlation <- as.vector(binary_correlation)
        heatmap_df$Annotation <- as.vector(binary_annotations)
        
        ggplot2::ggplot(heatmap_df,ggplot2::aes(x = Species1, y = Species2, fill = Correlation)) +
          geom_tile() +
          geom_text(ggplot2::aes(label = Annotation), size = geom_defaults$text_size, color = geom_defaults$text_color) +
          scale_fill_gradient2(low = plot_colors[["gradient"]][1], 
                              mid = plot_colors[["gradient"]][2], 
                              high = plot_colors[["gradient"]][3], 
                              midpoint = 0, limits = c(-1, 1)) +
          # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
{ theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; get_application_theme(theme_name) } +
          ggplot2::theme(axis.text.x =ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = sprintf("Species Correlation - Core Hotspot Presence/Absence\n(%s correlation, %d genes, *p<0.05, **p<0.01, ***p<0.001)", 
                              tools::toTitleCase(correlation_method), length(core_genes)),
               x = "Species", y = "Species")
      })
      plots[["binary_correlation"]] <- binary_heatmap
    } else {
      log_message("Binary correlation matrix not suitable for heatmap visualization", level = "warning")
    }
    
    # REMOVED: 2. Frequency-based correlation heatmap 
    # Deleted to avoid misleading results with sparse hotspot gene data
    # User feedback: sparse matrices → zero variance → meaningless frequency correlations
    log_message("Frequency-based correlation analysis removed to prevent misleading results with sparse data")
    
    # REMOVED: Gene distribution heatmap (identified as valueless graphic)
    # This visualization was removed per user feedback as it provides limited analytical value
    # compared to the species-species correlation analyses above
    
    # Prepare correlation summary data (simplified to binary-only)
    correlation_data <- list(
      binary_correlation = binary_correlation,
      binary_pvalues = binary_pvalues,
      core_genes = core_genes,
      gene_species_count = gene_species_count,
      core_matrix = core_matrix,
      analysis_mode = "binary_only"  # Indicates simplified analysis
    )
    
    # Create summary statistics
    metadata <- list(
      n_core_genes = length(core_genes),
      n_species_analyzed = nrow(species_matrix),
      min_species_threshold = min_species_threshold,
      correlation_method = correlation_method,
      total_core_hotspot_entries = nrow(core_hotspots),
      top_shared_genes = head(gene_species_count[[gene_col]], 10),
      analysis_date = Sys.time()
    )
    
    log_message(sprintf("Hotspot correlation analysis completed: %d core genes, %d species", 
                       length(core_genes), nrow(species_matrix)))
    
    return(list(
      plots = plots,
      correlation_data = correlation_data,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate hotspot correlation heatmap: %s", e$message), level = "error")
    return(NULL)
  })
}

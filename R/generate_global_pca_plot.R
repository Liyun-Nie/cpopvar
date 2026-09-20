#' Generate global PCA plot
#' 
#' @title Generate global PCA plot
#' @description Creates PCA visualization from global heatmap results.
#' Supports both CDS-based (gene) and IGS-based (poiGS_ID) feature analysis.
#' @param global_heatmap_result Results from global heatmap generation (optional)
#' @param normalized_data Normalized frequency dataset (optional)
#' @param config Configuration object
#' @param group_by_vars Variables to group by in PCA visualization
#' @param task_params Task-specific parameters
#' @param feature_col Name of the feature column to use for grouping 
#'   (default: "gene" for M03, use "poiGS_ID" for M05)
#' @return List containing PCA plot and analysis results
#' @import ggplot2
#' @importFrom rlang sym
#' @importFrom tibble column_to_rownames
#' @export
generate_global_pca_plot <- function(global_heatmap_result = NULL,
                                   normalized_data = NULL,
                                   config,
                                   group_by_vars = NULL,
                                   task_params = NULL,
                                   feature_col = "gene") {
  
  # Get parameters from config, with a final fallback default.
  n_components <- get_task_parameter(task_params, config, "pca_components", 4)
  min_variance_explained <- get_task_parameter(task_params, config, "min_variance_explained", 0.05)
  
  # Get visual standards for consistent styling
  visual_standards <- get_visual_standards()
  geom_defaults <- visual_standards$ggplot2
  
  tryCatch({
    log_message(sprintf("Global PCA: Performing principal component analysis on %s frequency patterns", feature_col))
    
    # Get structured color palette from configuration
    plot_colors <- get_color_palette(config, "global_pca")
    
    # Try to use global heatmap data first, fallback to creating matrix from normalized_data
    heatmap_matrix <- NULL
    
    if (!is.null(global_heatmap_result) && !is.null(global_heatmap_result$data)) {
      heatmap_matrix <- global_heatmap_result$data
      log_message("PCA: Using global heatmap matrix data")
    } else if (!is.null(normalized_data)) {
      log_message("PCA: Global heatmap data not available, creating matrix from normalized_data")
      
      # Create matrix from normalized_data (similar to global heatmap logic)
      species_col <- "species"
      gene_col <- feature_col  # Use parameterized feature column
      freq_col <- "frequency_per_kb"
      
      # Filter out NA frequency values
      filtered_data <- normalized_data %>%
        dplyr::filter(!is.na(.data[[freq_col]]))
      
      if (nrow(filtered_data) == 0) {
        log_message("No valid frequency data available for PCA", level = "warning")
        return(NULL)
      }
      
      # Create gene-species matrix
      gene_summary <- filtered_data %>%
      dplyr::group_by(.data[[gene_col]], .data[[species_col]]) %>%
      dplyr::summarise(!!sym(freq_col) := max(.data[[freq_col]], na.rm = TRUE), .groups = "drop")
      
      # Get all genes and species
      all_species <- unique(gene_summary[[species_col]])
      all_genes <- unique(gene_summary[[gene_col]])
      
      # Create complete combination matrix
      combination_data <- expand.grid(
        gene = factor(all_genes, levels = all_genes),
        species = factor(all_species, levels = all_species),
        stringsAsFactors = FALSE
      )
      names(combination_data) <- c(gene_col, species_col)
      
      # Merge and create matrix
      complete_data <- combination_data %>%
        dplyr::left_join(gene_summary, by = setNames(c(gene_col, species_col), c(gene_col, species_col))) %>%
        dplyr::mutate(frequency_display = ifelse(is.na(.data[[freq_col]]), NA_real_, .data[[freq_col]]))
      
      # Convert to matrix
      heatmap_matrix <- complete_data %>%
        dplyr::select(.data[[gene_col]], .data[[species_col]], frequency_display) %>%
        tidyr::spread(key = !!sym(species_col), value = frequency_display) %>%
        tibble::column_to_rownames(gene_col) %>%
        as.matrix()
        
      log_message(sprintf("PCA: Created matrix from normalized_data: %d features (%s) x %d species", 
                         nrow(heatmap_matrix), feature_col, ncol(heatmap_matrix)))
    } else {
      log_message("No data source available for PCA analysis (neither heatmap nor normalized_data)", level = "warning")
      return(NULL)
    }
    
    # Prepare data for PCA - transpose so species are observations and features are variables
    # Remove features with all missing values
    genes_with_data <- apply(heatmap_matrix, 1, function(x) sum(!is.na(x)) >= 2)
    species_with_data <- apply(heatmap_matrix, 2, function(x) sum(!is.na(x)) >= 2)
    
    if (sum(genes_with_data) < 3 || sum(species_with_data) < 3) {
      log_message(sprintf("Insufficient data for PCA analysis (need >=3 %s and >=3 species)", gsub("_ID$", "", feature_col)), level = "warning")
      return(NULL)
    }
    
    # Filter matrix and transpose for PCA (species as rows, genes as columns)
    pca_matrix <- t(heatmap_matrix[genes_with_data, species_with_data, drop = FALSE])
    
    # Handle missing values by imputation (use column means)
    for (i in 1:ncol(pca_matrix)) {
      pca_matrix[is.na(pca_matrix[, i]), i] <- mean(pca_matrix[, i], na.rm = TRUE)
    }
    
    # Remove any genes that still have no variation
    gene_variances <- apply(pca_matrix, 2, var, na.rm = TRUE)
    valid_genes <- !is.na(gene_variances) & gene_variances > 0
    
    if (sum(valid_genes) < 3) {
      log_message(sprintf("Insufficient %s variation for PCA analysis", gsub("_ID$", "", feature_col)), level = "warning")
      return(NULL)
    }
    
    pca_matrix <- pca_matrix[, valid_genes, drop = FALSE]
    
    log_message(sprintf("PCA input matrix: %d species x %d features (%s)", nrow(pca_matrix), ncol(pca_matrix), feature_col))
    
    # Perform PCA
    pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
    
    # Extract variance explained
    variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
    cumulative_variance <- cumsum(variance_explained)
    
    n_components <- min(n_components, length(variance_explained))
    
    # Create comprehensive PCA data including PC3 for 3D visualization
    pca_data <- data.frame(
      species = rownames(pca_result$x),
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      PC3 = if(ncol(pca_result$x) >= 3) pca_result$x[, 3] else 0,
      stringsAsFactors = FALSE
    )
    
    # Load species grouping information dynamically from group_info.csv
    # Extract session_id from global_heatmap_result metadata if available
    session_id <- if (!is.null(global_heatmap_result$metadata$session_id)) {
      global_heatmap_result$metadata$session_id
    } else {
      NULL
    }
    
    # Load group information using the robust, validated helper
    # Use config session_id as primary source, fallback to metadata session_id
    effective_session_id <- config$session_info$session_id %||% session_id
    session_paths <- get_session_paths(effective_session_id)
    pca_data_with_groups <- merge_with_group_info(
        main_data = pca_data,
        session_paths = session_paths
    )
    
    # Ensure we have valid data objects before proceeding
    if (is.null(pca_data)) {
      log_message("PCA data is NULL, cannot proceed with group mapping", level = "error")
      return(NULL)
    }
    
    # Use merged data directly (simplified logic - no redundant re-merging)
    if (!is.null(pca_data_with_groups) && is.data.frame(pca_data_with_groups) && 
        ncol(pca_data_with_groups) > ncol(pca_data)) {
      # Group info was successfully merged, use the complete merged dataset
      log_message("Group information loaded and merged for PCA analysis")
      pca_data <- pca_data_with_groups
      
      # V20 ARCHITECTURAL ALIGNMENT:
      # This function no longer assumes which columns to use for grouping.
      # It receives a 'group_by_vars' vector from the orchestrator, which is read from the config.
      # This enforces the "configuration-driven" principle.
      actual_group_cols <- character(0)
      if (!is.null(group_by_vars) && length(group_by_vars) > 0) {
          # Validate which of the requested group variables are actually in the merged data.
          actual_group_cols <- intersect(group_by_vars, names(pca_data))
          
          if (length(actual_group_cols) < length(group_by_vars)) {
              missing_vars <- setdiff(group_by_vars, actual_group_cols)
              log_message(sprintf("The following group variables requested from config were not found in the data and will be skipped: %s", 
                                  paste(missing_vars, collapse = ", ")), level = "warning")
          }
          
          log_message(sprintf("Will generate PCA plots grouped by: %s", paste(actual_group_cols, collapse = ", ")))
          
          # ============================================================
          # Comprehensive validation before processing
          # ============================================================
          # Add comprehensive NULL checks before any attribute operations in PCA
          
          # Validate pca_data is not NULL and is a valid data frame
          if (is.null(pca_data)) {
            log_message("PCA data is NULL before group column processing", level = "error")
            return(NULL)
          }
          
          if (!is.data.frame(pca_data)) {
            log_message("PCA data is not a valid data frame", level = "error")
            return(NULL)
          }
          
          if (nrow(pca_data) == 0) {
            log_message("PCA data frame is empty", level = "error")
            return(NULL)
          }
          
          # Clean up any empty factor levels in group columns to prevent NULL errors
          for (group_col in actual_group_cols) {
            # Validate each group column exists and is valid before processing
            if (!group_col %in% names(pca_data)) {
              log_message(sprintf("Group column '%s' not found in PCA data", group_col), level = "warning")
              next
            }
            
            if (is.null(pca_data[[group_col]])) {
              log_message(sprintf("Group column '%s' is NULL", group_col), level = "warning")
              next
            }
            
            # Safe attribute operations with error handling
            tryCatch({
              if (is.factor(pca_data[[group_col]])) {
                log_message(sprintf("Pre-sanitization: Group factor '%s' has %d levels.",
                                   group_col, nlevels(pca_data[[group_col]])))
                
                pca_data[[group_col]] <- droplevels(pca_data[[group_col]])
                
                log_message(sprintf("Post-sanitization: Group factor '%s' now has %d levels.",
                                   group_col, nlevels(pca_data[[group_col]])))
              } else if (is.character(pca_data[[group_col]])) {
                log_message(sprintf("Converting character column '%s' to factor for PCA", group_col))
                pca_data[[group_col]] <- factor(pca_data[[group_col]])
                pca_data[[group_col]] <- droplevels(pca_data[[group_col]])
                log_message(sprintf("Converted and sanitized: Group factor '%s' now has %d levels.",
                                   group_col, nlevels(pca_data[[group_col]])))
              }
            }, error = function(e) {
              log_message(sprintf("Failed to process group column '%s': %s", group_col, e$message), level = "error")
            })
          }
          # ============================================================
          
      } else {
          log_message("No grouping variables provided from config for PCA plot.", level = "info")
      }
    } else {
      # Fallback to generic grouping if no group info available
      log_message("No group information available for PCA analysis", level = "warning")
      
      # Safe fallback attribute setting
      if (!is.null(pca_data) && is.data.frame(pca_data)) {
        tryCatch({
          pca_data$group1 <- "All_Species"
          pca_data$group2 <- "All_Species"  
          pca_data$group3 <- "All_Species"
        }, error = function(e) {
          log_message(sprintf("Failed to set fallback group variables: %s", e$message), level = "error")
          return(NULL)
        })
      } else {
        log_message("PCA data is NULL or invalid, cannot set fallback groups", level = "error")
        return(NULL)
      }
    }
    
    # Create species ranking table based on PC coordinates
    species_ranking <- pca_data %>%
      dplyr::arrange(dplyr::desc(PC1)) %>%
      dplyr::mutate(
        PC1_rank = dplyr::row_number(),
        PC2_rank = rank(-PC2, ties.method = "min"),
        PC3_rank = rank(-PC3, ties.method = "min"),
        combined_rank = (PC1_rank + PC2_rank + PC3_rank) / 3
      ) %>%
      dplyr::arrange(combined_rank) %>%
      dplyr::select(species, PC1, PC2, PC3, PC1_rank, PC2_rank, PC3_rank, combined_rank)
    
    # Create comprehensive 3D PCA analysis with confidence ellipses for each group
    
    # Function to create individual PCA plot with confidence ellipse
    create_pca_plot <- function(data, x_var, y_var, group_var, title, x_label, y_label, colors) {
      p <- ggplot2::ggplot(data,ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]])) +
        geom_point(size = geom_defaults$point_size, alpha = geom_defaults$point_alpha) +
        stat_ellipse(level = 0.95, alpha = geom_defaults$violin_alpha, size = geom_defaults$line_size) +  # 95% confidence ellipse
        ggplot2::labs(title = title, x = x_label, y = y_label, color = "Group") +
        ggplot2::scale_color_manual(values = colors) +
        # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
        { theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; 
          base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14;
          get_application_theme(theme_name, NULL, base_size) }
      
      # Add text labels if ggrepel is available
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[["species"]]), 
                                         size = geom_defaults$text_size, max.overlaps = 15, alpha = geom_defaults$point_alpha)
      } else {
        p <- p + geom_text(ggplot2::aes(label = .data[["species"]]), 
                          size = geom_defaults$text_size, nudge_y = 0.1, check_overlap = TRUE, alpha = geom_defaults$point_alpha)
      }
      
      return(p)
    }
    
    # V17 CRITICAL FIX: Use actual column names instead of hardcoded group1/group2/group3
    # Generate individual PCA plots for each actual group variable
    group_plots <- list()
    
    # Use the actual group columns that exist in the data
    actual_group_names <- if (exists("actual_group_cols")) actual_group_cols else character(0)
    
    # Get descriptive names for groups based on actual available data
    get_group_description <- function(col_name, data) {
      if (!is.null(data) && col_name %in% names(data)) {
        unique_vals <- length(unique(data[[col_name]][!is.na(data[[col_name]])]))
        return(paste(tools::toTitleCase(gsub("_", " ", col_name)), "(", unique_vals, "categories)"))
      } else {
        return(tools::toTitleCase(gsub("_", " ", col_name)))
      }
    }
    
    # Define PC combinations for 1×3 layout
    pc_combinations <- list(
      list(x = "PC1", y = "PC2", label = "PC1 vs PC2"),
      list(x = "PC1", y = "PC3", label = "PC1 vs PC3"),
      list(x = "PC2", y = "PC3", label = "PC2 vs PC3")
    )
    
    for (group_var in actual_group_names) {
      # V21 FINAL POLISH: Add NA handling to match M02's robustness.
      
      # First, check if column exists (already implemented in V20)
      if (is.null(group_var) || is.na(group_var) || !group_var %in% names(pca_data)) {
        log_message(sprintf("Skipping PCA plot for group '%s': column not found in PCA data.", group_var), 
                   level = "warning")
        next # Skip to the next group
      }

      # Add the NA filter
      plot_data_for_group <- pca_data %>%
        dplyr::filter(!is.na(.data[[group_var]]))

      if(nrow(plot_data_for_group) == 0) {
        log_message(sprintf("Skipping PCA plot for group '%s': no valid data after removing NAs.", group_var), level = "warning")
        next
      }
      
      group_desc <- get_group_description(group_var, plot_data_for_group)
      
      # Check if group variable has variation (use filtered data)
      if (group_var %in% names(plot_data_for_group)) {
        group_values <- unique(plot_data_for_group[[group_var]])
        group_values <- group_values[!is.na(group_values)]
        
        if (length(group_values) > 1) {
          # Create color palette for this group using centralized management
          group_colors <- get_color_palette(config, group_var, group_values)
          
          # Generate three plots for this group (PC1/PC2, PC1/PC3, PC2/PC3)
          for (pc_idx in 1:3) {
            pc_combo <- pc_combinations[[pc_idx]]
            plot_name <- sprintf("pca_%d_%s", pc_idx, group_var)
            
            group_plots[[plot_name]] <- create_pca_plot(
              plot_data_for_group,
              x_var = pc_combo$x,
              y_var = pc_combo$y,
              group_var = group_var,
              title = sprintf("%s (%s)", pc_combo$label, group_desc),
              x_label = sprintf("%s (%.1f%% variance)", pc_combo$x, 
                               variance_explained[as.numeric(substr(pc_combo$x, 3, 3))] * 100),
              y_label = sprintf("%s (%.1f%% variance)", pc_combo$y,
                               if(pc_combo$y == "PC3" && length(variance_explained) >= 3) {
                                 variance_explained[3] * 100
                               } else if(pc_combo$y == "PC2") {
                                 variance_explained[2] * 100
                               } else {
                                 variance_explained[1] * 100
                               }),
              colors = group_colors
            )
          }
          
          log_message(sprintf("Generated PCA plots for %s with %d categories", group_var, length(group_values)))
        } else {
          log_message(sprintf("Skipping %s: insufficient variation (%d categories)", group_var, length(group_values)), level = "warning")
        }
      } else {
        log_message(sprintf("Group variable %s not found in PCA data", group_var), level = "warning")
      }
    }
    
    # Create scree plot
    scree_data <- data.frame(
      component = 1:length(variance_explained),
      variance = variance_explained,
      cumulative = cumulative_variance
    )
    
    scree_plot <-ggplot2::ggplot(scree_data[1:min(10, nrow(scree_data)), ], 
                                ggplot2::aes(x = component, y = variance)) +
      geom_col(fill = plot_colors[["scree_bar"]], alpha = geom_defaults$density_alpha) +
      geom_line(ggplot2::aes(y = cumulative), color = plot_colors[["scree_line"]], size = geom_defaults$line_size) +
      geom_point(ggplot2::aes(y = cumulative), color = plot_colors[["scree_line"]], size = geom_defaults$point_size) +
      scale_x_continuous(breaks = 1:min(10, nrow(scree_data))) +
      scale_y_continuous(
        name = "Proportion of Variance",
        sec.axis = sec_axis(~., name = "Cumulative Variance", 
                           labels = scales::percent_format())
      ) +
      ggplot2::labs(
        title = "PCA Scree Plot",
        subtitle = "Variance explained by each principal component",
        x = "Principal Component"
      ) +
      # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
      { theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; 
        base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14;
        get_application_theme(theme_name, NULL, base_size) } +
      ggplot2::theme(
        axis.title.y.right =ggplot2::element_text(color = plot_colors[["scree_line"]]),
        axis.text.y.right =ggplot2::element_text(color = plot_colors[["scree_line"]])
      )
    
    # Create comprehensive loadings plots for PC1/PC2, PC1/PC3, PC2/PC3
    loadings_data <- data.frame(
      gene = rownames(pca_result$rotation),
      PC1_loading = pca_result$rotation[, 1],
      PC2_loading = pca_result$rotation[, 2],
      PC3_loading = if(ncol(pca_result$rotation) >= 3) pca_result$rotation[, 3] else 0,
      stringsAsFactors = FALSE
    )
    
    # Get top contributing genes for each component combination
    loadings_data$contribution_12 <- sqrt(loadings_data$PC1_loading^2 + loadings_data$PC2_loading^2)
    loadings_data$contribution_13 <- sqrt(loadings_data$PC1_loading^2 + loadings_data$PC3_loading^2)
    loadings_data$contribution_23 <- sqrt(loadings_data$PC2_loading^2 + loadings_data$PC3_loading^2)
    
    top_genes_12 <- head(loadings_data[order(-loadings_data$contribution_12), ], 15)
    top_genes_13 <- head(loadings_data[order(-loadings_data$contribution_13), ], 15)
    top_genes_23 <- head(loadings_data[order(-loadings_data$contribution_23), ], 15)
    
    # Function to create loadings plot
    create_loadings_plot <- function(data, x_var, y_var, contrib_var, title, x_label, y_label) {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], size = .data[[contrib_var]])) +
        geom_point(alpha = geom_defaults$point_alpha, color = plot_colors[["loadings_points"]]) +
        geom_hline(yintercept = 0, linetype = geom_defaults$line_ref_linetype, alpha = geom_defaults$line_ref_alpha) +
        geom_vline(xintercept = 0, linetype = geom_defaults$line_ref_linetype, alpha = geom_defaults$line_ref_alpha) +
        ggplot2::labs(title = title, x = x_label, y = y_label, size = "Contribution") +
        # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
        { theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"; 
          base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14;
          get_application_theme(theme_name, NULL, base_size) }
      
      # Add text labels if ggrepel is available
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[["gene"]]), 
                                         size = geom_defaults$text_size, max.overlaps = 12, 
                                         box.padding = 0.3, point.padding = 0.3)
      } else {
        p <- p + geom_text(ggplot2::aes(label = .data[["gene"]]), 
                          size = geom_defaults$text_size, nudge_x = 0.01, nudge_y = 0.01, check_overlap = TRUE)
      }
      
      return(p)
    }
    
    # Create three loadings plots
    loadings_plot_12 <- create_loadings_plot(
      top_genes_12, "PC1_loading", "PC2_loading", "contribution_12",
      "PC1 vs PC2 Loadings", "PC1 Loading", "PC2 Loading"
    )
    
    loadings_plot_13 <- create_loadings_plot(
      top_genes_13, "PC1_loading", "PC3_loading", "contribution_13", 
      "PC1 vs PC3 Loadings", "PC1 Loading", "PC3 Loading"
    )
    
    loadings_plot_23 <- create_loadings_plot(
      top_genes_23, "PC2_loading", "PC3_loading", "contribution_23",
      "PC2 vs PC3 Loadings", "PC2 Loading", "PC3 Loading"
    )
    
    # Create group-specific 1×3 composite plots with confidence ellipses dynamically
    composite_plots <- list()
    
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      # Generate composite plots for each group that has plots
      for (group_idx in seq_along(actual_group_names)) {
        group_var <- actual_group_names[group_idx]
        group_desc <- get_group_description(group_var, pca_data)
        
        # Check if we have all three plots for this group
        plot1_name <- sprintf("pca_1_%s", group_var)
        plot2_name <- sprintf("pca_2_%s", group_var)
        plot3_name <- sprintf("pca_3_%s", group_var)
        
        if (all(c(plot1_name, plot2_name, plot3_name) %in% names(group_plots))) {
          composite_name <- sprintf("combined_%s", group_var)
          
          composite_plots[[composite_name]] <- gridExtra::grid.arrange(
            group_plots[[plot1_name]], 
            group_plots[[plot2_name]], 
            group_plots[[plot3_name]],
            nrow = 1, ncol = 3,
            top = sprintf("PCA Analysis: %s (with 95%% Confidence Ellipses)", group_desc)
          )
          
          log_message(sprintf("Created composite plot for %s", group_var))
        } else {
          log_message(sprintf("Skipping composite plot for %s: missing individual plots", group_var), level = "warning")
        }
      }
      
      # Loadings composite plot: PC1vsPC2, PC1vsPC3, PC2vsPC3
      composite_plots[["combined_loadings"]] <- gridExtra::grid.arrange(
        loadings_plot_12, loadings_plot_13, loadings_plot_23,
        nrow = 1, ncol = 3,
        top = "Gene Loadings Analysis: PC Component Contributions"
      )
      
      # Scree plot (individual)
      composite_plots[["scree"]] <- scree_plot
      
    } else {
      # Fallback: use individual plots if gridExtra not available
      log_message("gridExtra not available, using individual plots", level = "warning")
      for (group_idx in seq_along(actual_group_names)) {
        group_var <- actual_group_names[group_idx]
        plot1_name <- sprintf("pca_1_%s", group_var)
        if (plot1_name %in% names(group_plots)) {
          composite_plots[[sprintf("combined_%s", group_var)]] <- group_plots[[plot1_name]]
        }
      }
      composite_plots[["combined_loadings"]] <- loadings_plot_12
      composite_plots[["scree"]] <- scree_plot
    }
    
    # Prepare metadata
    metadata <- list(
      n_species = nrow(pca_data),
      n_genes = ncol(pca_matrix),
      variance_explained = variance_explained[1:n_components],
      cumulative_variance = cumulative_variance[1:n_components],
      total_variance_captured = cumulative_variance[n_components],
      top_contributing_genes = loadings_data$gene[1:min(10, nrow(loadings_data))],
      analysis_date = Sys.time()
    )
    
    log_message(sprintf("PCA analysis completed: %.1f%% variance explained by first %d components", 
                       cumulative_variance[n_components] * 100, n_components))
    
    # Combine all plots for return
    all_plots <- c(composite_plots, group_plots)
    all_plots[["loadings_12"]] <- loadings_plot_12
    all_plots[["loadings_13"]] <- loadings_plot_13  
    all_plots[["loadings_23"]] <- loadings_plot_23
    
    return(list(
      plots = all_plots,
      pca_data = pca_data,
      species_ranking = species_ranking,
      loadings_data = loadings_data,
      pca_result = pca_result,
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate global PCA plot: %s", e$message), level = "error")
    return(NULL)
  })
}

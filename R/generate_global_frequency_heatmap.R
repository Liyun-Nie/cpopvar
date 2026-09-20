#' Generate global frequency heatmap
#' 
#' @title Generate global frequency heatmap
#' @description Creates global frequency heatmap visualization for M03 hotspot analysis and M05 poiGS analysis
#' @param normalized_data Normalized frequency dataset
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param filename Output PDF filename (optional)
#' @param show_rownames Whether to show row names in heatmap (optional)
#' @param main Main title for the heatmap (optional)
#' @param width PDF width in inches (optional)
#' @param height PDF height in inches (optional)
#' @param feature_col Name of the feature column to use for grouping (default: "gene", M05 uses "poiGS_ID")
#' @return List containing heatmap plot and processed data
#' @importFrom rlang sym
#' @importFrom tibble column_to_rownames
#' @export
generate_global_frequency_heatmap <- function(normalized_data,
                                             config,
                                             task_params = NULL,
                                             filename = NULL,
                                             show_rownames = NULL,
                                             main = NULL,
                                             width = NULL,
                                             height = NULL,
                                             feature_col = "gene") {

  tryCatch({

    # Load required package for heatmap generation
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("Package 'pheatmap' is required for heatmap generation but not available")
    }

    # Default to M03 branding, but allow override for M05
    module_id <- if (!is.null(main) && grepl("IGS", main)) "M05" else "M03"
    data_type <- if (module_id == "M05") "IGS" else "gene"
    
    log_message(sprintf("%s Step 0a: Generating global %s frequency heatmap with functional ordering", module_id, data_type))

    # Get structured color palette from configuration
    plot_colors <- tryCatch({
      get_color_palette(config, "global_heatmap")
    }, error = function(e) {
      log_message(sprintf("Global heatmap colors not found, using defaults: %s", e$message), level = "warning")
      # Fallback colors for global heatmap
      list(
        gradient = c("#FFFFBF", "#FF7F0EFF", "#EE0000FF")  # Yellow-Orange-Red gradient for frequencies
      )
    })

    # Use internal standard column names (boundary principle)
    region_col <- "region_type"
    var_type_col <- "var_type"
    freq_col <- "frequency_per_kb"
    species_col <- "species"
    gene_col <- feature_col  # Use parameterized feature column (gene or poiGS_ID)

    # Trust the input data - it's already been scoped/filtered (Rule #11)
    # Only filter out NA values for data quality
    filtered_data <- normalized_data %>%
      dplyr::filter(!is.na(.data[[freq_col]]))

    if (nrow(filtered_data) == 0) {
      log_message("No valid data found for global heatmap", level = "warning")
      return(NULL)
    }

    log_message(sprintf("Processing data: %d rows for global frequency analysis", nrow(filtered_data)))

    # Load gene function mapping for ordering (using unified utils version)
    gene_function_map <- load_gene_function_mapping(config = config)

    # Get species ordering using unified helper (dependency injection)
    species_ordering_result <- load_and_apply_species_ordering(
      data = filtered_data,
      session_id = config$session_info$session_id,
      species_col = species_col
    )

    species_order <- species_ordering_result$species_order
    log_message(species_ordering_result$message)

    # Get all species and genes
    all_species <- unique(filtered_data[[species_col]])
    all_genes <- unique(filtered_data[[gene_col]])

    # Apply species ordering if available
    if (!is.null(species_order)) {
      ordered_species <- intersect(species_order, all_species)
      missing_species <- setdiff(all_species, species_order)
      all_species <- c(ordered_species, missing_species)
      log_message(sprintf("Applied phylogenetic ordering to %d species", length(ordered_species)))
    }

    # Apply gene functional ordering if gene function mapping is available (preliminary)
    if (!is.null(gene_function_map)) {
      log_message("Gene function mapping available for later ordering")
    } else {
      log_message("Gene function mapping not available, will use alphabetical ordering", level = "warning")
    }

    # Create complete gene-species combinations
    combination_data <- expand.grid(
      gene = factor(all_genes, levels = all_genes),
      species = factor(all_species, levels = all_species),
      stringsAsFactors = FALSE
    )
    names(combination_data) <- c(gene_col, species_col)
    complete_combinations <- combination_data

    # Merge with actual data, ensuring join columns are characters for robustness
    gene_summary <- filtered_data %>%
      dplyr::group_by(.data[[gene_col]], .data[[species_col]]) %>%
      dplyr::summarise(!!sym(freq_col) := max(.data[[freq_col]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        !!sym(gene_col) := as.character(.data[[gene_col]]),
        !!sym(species_col) := as.character(.data[[species_col]])
      )

    # Also convert combination_data columns to character for a robust join
    complete_combinations_char <- complete_combinations %>%
        dplyr::mutate(
            !!sym(gene_col) := as.character(.data[[gene_col]]),
            !!sym(species_col) := as.character(.data[[species_col]])
        )

    heatmap_data <- complete_combinations_char %>%
      dplyr::left_join(gene_summary, by = c(gene_col, species_col)) %>%
      dplyr::mutate(
        # Distinguish between true zero (gene present, no mutations) and NA (gene absent/no data)
        frequency_display = dplyr::case_when(
          is.na(.data[[freq_col]]) ~ NA_real_,  # Missing data - will be grey
          .data[[freq_col]] == 0 ~ 0,           # True zero - will be white
          TRUE ~ .data[[freq_col]]              # Actual frequency values
        )
      )

    # Convert to matrix format maintaining ordering
    
    # Validate heatmap_data structure before matrix creation

    heatmap_matrix <- heatmap_data %>%
      dplyr::select(.data[[gene_col]], .data[[species_col]], frequency_display) %>%
      tidyr::spread(key = !!sym(species_col), value = frequency_display) %>%
      tibble::column_to_rownames(var = gene_col) %>%
      as.matrix()

    # Matrix-derived ordering: get accurate gene and species lists from matrix dimensions
    genes_in_matrix <- rownames(heatmap_matrix)
    species_in_matrix <- colnames(heatmap_matrix)

    # Sort based on actual matrix dimensions
    if (!is.null(gene_function_map)) {
      # Filter map to only include genes present in our matrix
      gene_order_df <- gene_function_map %>%
        dplyr::filter(gene %in% genes_in_matrix) %>%
        dplyr::arrange(gene_category, gene)
      
      # Get the genes that are in the map and ordered
      ordered_genes <- gene_order_df$gene
      # Get the genes that are in the matrix but NOT in the map
      missing_genes <- setdiff(genes_in_matrix, ordered_genes)
      
      # Combine them: ordered genes first, then the rest alphabetically
      all_genes <- c(ordered_genes, sort(missing_genes))
      
      log_message(sprintf("Applied functional ordering to %d genes; %d genes without category appended.",
                         length(ordered_genes), length(missing_genes)))
    } else {
      all_genes <- sort(genes_in_matrix)
    }

    if (!is.null(species_order)) {
      ordered_species <- intersect(species_order, species_in_matrix)
      missing_species <- setdiff(species_in_matrix, ordered_species)
      all_species <- c(ordered_species, sort(missing_species))
    } else {
      all_species <- sort(species_in_matrix)
    }

    log_message(sprintf("Matrix reordering: %d genes, %d species based on matrix dimensions", 
                       length(all_genes), length(all_species)))

    # --- DEFENSIVE MATRIX SUBSETTING ---
    # Only use genes and species that actually exist in the matrix
    valid_genes <- intersect(all_genes, rownames(heatmap_matrix))
    valid_species <- intersect(all_species, colnames(heatmap_matrix))
    
    log_message(sprintf("Matrix validation: %d/%d genes valid, %d/%d species valid",
                       length(valid_genes), length(all_genes),
                       length(valid_species), length(all_species)))
    
    heatmap_matrix <- heatmap_matrix[valid_genes, valid_species, drop = FALSE]
    # 

    if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
      log_message("Empty matrix for global frequency heatmap", level = "warning")
      return(NULL)
    }

    log_message(sprintf("Global heatmap matrix: %d genes x %d species", nrow(heatmap_matrix), ncol(heatmap_matrix)))

    # Create color scheme for frequency values
    # Add comprehensive matrix diagnostics
    log_message(sprintf("Matrix diagnostics: class=%s, mode=%s, dim=%s", 
                       class(heatmap_matrix)[1], mode(heatmap_matrix), paste(dim(heatmap_matrix), collapse="x")))
    
    if (!is.matrix(heatmap_matrix)) {
      log_message(sprintf("CRITICAL: Not a matrix! class=%s", class(heatmap_matrix)[1]), level="error")
      return(NULL)
    }
    
    if (!is.numeric(heatmap_matrix)) {
      log_message(sprintf("CRITICAL: Matrix not numeric! mode=%s", mode(heatmap_matrix)), level="error")
      return(NULL)
    }
    
    # Add detailed matrix structure diagnostic
    tryCatch({
      log_message("Analyzing matrix structure...")
      log_message(sprintf("Matrix rownames length: %d", length(rownames(heatmap_matrix))))
      log_message(sprintf("Matrix colnames length: %d", length(colnames(heatmap_matrix))))
      log_message(sprintf("First few rownames: %s", paste(head(rownames(heatmap_matrix), 3), collapse=", ")))
      log_message(sprintf("First few colnames: %s", paste(head(colnames(heatmap_matrix), 3), collapse=", ")))
      
      # Test basic matrix operations step by step
      log_message("Testing is.na() operation...", level = "debug")
      test_na <- is.na(heatmap_matrix)
      log_message("is.na() passed", level = "debug")
      
      log_message("Testing is.finite() operation...", level = "debug")
      test_finite <- is.finite(heatmap_matrix)  
      log_message("is.finite() passed", level = "debug")
      
      log_message("Testing logical combination...", level = "debug")
      test_combined <- test_na | test_finite  # Use OR first to test
      log_message("Logical combination passed", level = "debug")
      
      log_message("Testing matrix subsetting...", level = "debug")
      # Try a simple subset first
      first_cell <- heatmap_matrix[1, 1]
      log_message(sprintf("First cell access successful: %s", first_cell), level = "debug")
      
      non_na_values <- heatmap_matrix[!test_na & test_finite]
      log_message(sprintf("Matrix indexing successful, extracted %d values", length(non_na_values)), level = "debug")
    }, error = function(e) {
      log_message(sprintf("MATRIX STRUCTURE ERROR: %s", e$message), level="error")
      # Emergency fallback - convert to simple numeric matrix
      tryCatch({
        log_message("Attempting emergency matrix reconstruction...")
        simple_matrix <- matrix(as.numeric(heatmap_matrix), 
                               nrow=nrow(heatmap_matrix), 
                               ncol=ncol(heatmap_matrix))
        rownames(simple_matrix) <- rownames(heatmap_matrix)
        colnames(simple_matrix) <- colnames(heatmap_matrix)
        heatmap_matrix <- simple_matrix
        non_na_values <- simple_matrix[!is.na(simple_matrix) & is.finite(simple_matrix)]
        log_message("Emergency matrix reconstruction successful", level = "debug")
      }, error = function(e2) {
        log_message(sprintf("Emergency reconstruction failed: %s", e2$message), level="error")
        return(NULL)
      })
    })
    if (length(non_na_values) == 0) {
      log_message("No valid values for color scaling", level = "warning")
      return(NULL)
    }

    # Add diagnostics for frequency value processing
    log_message(sprintf("Processing %d valid frequency values", length(non_na_values)))
    log_message(sprintf("Frequency value range: [%g, %g]", min(non_na_values, na.rm=TRUE), max(non_na_values, na.rm=TRUE)))
    
    tryCatch({
      max_freq <- max(non_na_values)
      log_message(sprintf("max_freq calculation successful: %g", max_freq), level = "debug")
      
      # Test the problematic indexing operation step by step
      log_message("Testing positive values filtering...", level = "debug")
      positive_condition <- non_na_values > 0
      log_message(sprintf("Positive condition created, length: %d", length(positive_condition)))
      positive_values <- non_na_values[positive_condition]
      log_message(sprintf("Positive values extracted: %d values", length(positive_values)))
      
      min_freq <- if(length(positive_values) > 0) min(positive_values) else 0
      log_message(sprintf("min_freq calculation successful: %g", min_freq), level = "debug")
    }, error = function(e) {
      log_message(sprintf("FREQUENCY PROCESSING ERROR: %s", e$message), level="error")
      return(NULL)
    })

    # Comprehensive solution for color and breaks generation edge cases
    
    # Step 1: Validate and prepare gradient colors
    heatmap_colors <- plot_colors[["gradient"]]
    if (is.null(heatmap_colors) || length(heatmap_colors) == 0) {
      log_message("WARNING: heatmap_colors configuration invalid, using fallback colors", level = "warning")
      heatmap_colors <- c("#3C5488", "#F39B7F", "#E64B35")  # NPG fallback gradient
    }
    
    # Ensure we have at least 2 colors for gradient generation
    if (length(heatmap_colors) < 2) {
      heatmap_colors <- c(heatmap_colors[1], "#FFFFFF")  # Add white as second color
    }
    
    # Step 2: Handle edge cases in frequency data
    freq_diff <- max_freq - min_freq
    
    # Case 1: No variation in data (all values the same)
    if (freq_diff == 0 || max_freq == 0) {
      colors <- c("white", heatmap_colors[1])
      color_breaks <- seq(0, max(max_freq, 0.001), length.out = length(colors) + 1)
      log_message(sprintf("Global heatmap: No frequency variation (max=%.3f). Using 2-color scale.", max_freq), level = "warning")
    }
    # Case 2: Very small variation that might cause numerical issues
    else if (freq_diff < 1e-6) {
      colors <- c("white", heatmap_colors[1], heatmap_colors[length(heatmap_colors)])
      color_breaks <- seq(0, max(max_freq, 0.001), length.out = length(colors) + 1)
      log_message(sprintf("Global heatmap: Very small frequency range (%.6f). Using 3-color scale.", freq_diff), level = "warning")
    }
    # Case 3: Normal case with meaningful range
    else {
      # Generate robust gradient with fixed 50-color resolution  
      n_breaks <- 50  # Total number of breaks (pheatmap requirement: colors = breaks - 1)
      n_colors <- n_breaks - 1  # 49 colors for 50 breaks
      
      # Generate breaks sequence: 0, min_freq, ..., max_freq
      break_sequence <- seq(min_freq, max_freq, length.out = n_breaks - 1)  # 49 breaks from min to max
      color_breaks <- c(0, break_sequence)  # 50 breaks total: 0 + 49 others
      
      # Generate corresponding colors: white for 0, then 48-color gradient  
      tryCatch({
        gradient_colors <- colorRampPalette(heatmap_colors)(n_colors - 1)  # 48 gradient colors
        colors <- c("white", gradient_colors)  # 49 colors total: white + 48 gradient
      }, error = function(e) {
        log_message(sprintf("Color gradient generation failed: %s. Using simple scale.", e$message), level = "error")
        colors <- c("white", heatmap_colors[1])
        color_breaks <- seq(0, max(max_freq, 0.001), length.out = length(colors) + 1)
      })
      
      log_message(sprintf("Global heatmap: Generated %d colors and %d breaks for range [%.3f, %.3f]",
                         length(colors), length(color_breaks), min_freq, max_freq))
    }
    
    # Step 3: Final validation - ensure pheatmap requirements
    if (length(color_breaks) != length(colors) + 1) {
      log_message(sprintf("Breaks/colors mismatch detected: %d breaks, %d colors. Using fallback scale.", 
                         length(color_breaks), length(colors)), level = "warning")
      # Emergency fallback to simple 2-color scale
      colors <- c("white", heatmap_colors[1])
      color_breaks <- seq(0, max(max_freq, 0.001), length.out = length(colors) + 1)
    }
    
    log_message(sprintf("Global heatmap: Final validation - %d colors, %d breaks (ratio: %.1f)",
                       length(colors), length(color_breaks), length(color_breaks)/length(colors)))
    # Color and breaks configuration completed

    # Generate annotations for gene categories with STRICT cleaning
    annotation_row <- NULL
    annotation_colors <- NULL
    if (!is.null(gene_function_map)) {
      # Expected 4 categories based on user feedback
      expected_categories <- c(
        "Other Functional Genes",
        "Unknown Function Genes", 
        "Photosynthesis Related",
        "Self-Replication Related"
      )

      # Apply STRICT gene category cleaning (same logic as enrichment analysis)
      gene_categories_raw <- gene_function_map %>%
        dplyr::filter(gene %in% rownames(heatmap_matrix)) %>%
        dplyr::select(gene, gene_category)

      log_message(sprintf("Global heatmap: Raw gene categories before cleaning: %d unique (%s)", 
                         length(unique(gene_categories_raw$gene_category)), 
                         paste(unique(gene_categories_raw$gene_category), collapse = ", ")))

      # Clean and standardize categories
      gene_categories_clean <- gene_categories_raw %>%
        dplyr::mutate(
          # Step 1: Basic cleaning
          gene_category_clean = trimws(gene_category),
          # Step 2: Remove empty/NA values
          gene_category_clean = ifelse(is.na(gene_category_clean) | gene_category_clean == "", NA, gene_category_clean)
        ) %>%
        dplyr::filter(!is.na(gene_category_clean)) %>%
        dplyr::mutate(
          # Step 3: Normalize and match with expected categories
          gene_category_final = sapply(gene_category_clean, function(cat) {
            cat_normalized <- tools::toTitleCase(tolower(trimws(cat)))

            # Try to match with expected categories
            for (expected in expected_categories) {
              if (grepl(gsub(" ", ".*", tolower(expected)), tolower(cat_normalized), ignore.case = TRUE) ||
                  grepl(gsub(" ", ".*", tolower(cat_normalized)), tolower(expected), ignore.case = TRUE)) {
                return(expected)
              }
            }
            return(cat_normalized)  # Keep original if no match
          })
        ) %>%
        # Step 4: Filter to only expected categories
        dplyr::filter(gene_category_final %in% expected_categories) %>%
        dplyr::select(gene, gene_category = gene_category_final) %>%
        tibble::column_to_rownames("gene")

      log_message(sprintf("Global heatmap: After cleaning - %d unique categories (%s)", 
                         length(unique(gene_categories_clean$gene_category)), 
                         paste(unique(gene_categories_clean$gene_category), collapse = ", ")))

      # Only include genes that are in our matrix and have valid categories
      genes_in_matrix_ann <- intersect(rownames(gene_categories_clean), rownames(heatmap_matrix))
      if (length(genes_in_matrix_ann) > 0) {
        annotation_row <- gene_categories_clean[genes_in_matrix_ann, , drop = FALSE]

        # Get colors from the central palette manager instead of hardcoding.
        unique_categories <- unique(annotation_row$gene_category)
        log_message(sprintf("Global heatmap: Requesting colors for categories: %s", paste(unique_categories, collapse = ", ")))

        category_colors <- tryCatch({
          get_color_palette(config, "gene_category", unique_categories)
        }, error = function(e) {
          log_message(sprintf("Color palette error: %s", e$message), level = "error")
          return(NULL)
        })

        if (is.null(category_colors) || length(category_colors) == 0) {
          log_message("No colors returned from get_color_palette, using fallback colors", level = "warning")
          # Create fallback colors
          fallback_colors <- rainbow(length(unique_categories))
          names(fallback_colors) <- unique_categories
          category_colors <- fallback_colors
        }

        annotation_colors <- list(gene_category = category_colors)
        log_message(sprintf("Global heatmap: Applied centrally managed colors to %d categories", length(unique_categories)))
      }
    }

    # Use fixed color palette (already properly configured above)
    # heatmap_color_palette configured via colors variable

    # Validate and set safe defaults for pheatmap colors
    na_color <- if (!is.null(plot_colors[["na_color"]]) && length(plot_colors[["na_color"]]) > 0) {
      plot_colors[["na_color"]]
    } else {
      "grey90"  # Safe fallback for missing data
    }

    border_color <- if (!is.null(plot_colors[["border_color"]]) && length(plot_colors[["border_color"]]) > 0) {
      plot_colors[["border_color"]]
    } else {
      "white"  # Safe fallback for borders
    }

    log_message(sprintf("Global heatmap colors: na_col=%s, border_color=%s", na_color, border_color))

    # Create pheatmap parameters with validated colors
    pheatmap_params <- list(
      mat = heatmap_matrix,
      cluster_rows = FALSE,  # Maintain functional ordering
      cluster_cols = FALSE,  # Maintain phylogenetic ordering
      color = colors,
      breaks = color_breaks,
      na_col = na_color,
      border_color = border_color,
      show_rownames = if (!is.null(show_rownames)) show_rownames else (nrow(heatmap_matrix) <= 100),  # Show gene names if not too many or if explicitly specified
      show_colnames = TRUE,
      fontsize = 8,
      fontsize_row = if(nrow(heatmap_matrix) > 50) 4 else 6,
      fontsize_col = 7,
      main = if (!is.null(main)) main else sprintf("Global Gene Frequency Distribution (%s %s)", 
                     tools::toTitleCase(get_task_parameter(task_params, config, "region_filter", "CDS")), 
                     toupper(get_task_parameter(task_params, config, "var_type_filter", "snp"))),
      angle_col = 45,
      silent = TRUE
    )
    
    # Add filename parameter if provided for automatic file saving
    if (!is.null(filename)) {
      pheatmap_params$filename <- filename
      
      # Add width and height parameters for PDF output (system architecture standard)
      # Default to M03-consistent dimensions if not specified
      pheatmap_params$width <- if (!is.null(width)) width else 14
      pheatmap_params$height <- if (!is.null(height)) height else 10
      
      log_message(sprintf("Heatmap will be saved to: %s (dimensions: %d x %d)", 
                         filename, pheatmap_params$width, pheatmap_params$height))
    }

    # Add row annotations if available with safety checks
    if (!is.null(annotation_row) && !is.null(annotation_colors)) {
      # Ensure gene_category is factor type for proper pheatmap levels matching
      if (!is.null(annotation_row) && "gene_category" %in% names(annotation_row)) {
        annotation_row$gene_category <- as.factor(annotation_row$gene_category)
      }
      pheatmap_params$annotation_row <- annotation_row
      pheatmap_params$annotation_colors <- annotation_colors
      log_message("Global heatmap: Added gene category annotations")

    } else if (!is.null(annotation_row)) {
      log_message("WARNING: annotation_row available but annotation_colors is NULL, skipping annotations", level = "warning")
    }

    # Generate the heatmap using namespace prefix
    p <- do.call(pheatmap::pheatmap, pheatmap_params)

    # Calculate statistics
    total_cells <- nrow(heatmap_matrix) * ncol(heatmap_matrix)
    na_count <- sum(is.na(heatmap_matrix))
    zero_count <- sum(heatmap_matrix == 0, na.rm = TRUE)
    positive_count <- sum(heatmap_matrix > 0, na.rm = TRUE)

    metadata <- list(
      total_genes = nrow(heatmap_matrix),
      total_species = ncol(heatmap_matrix),
      total_cells = total_cells,
      missing_data_cells = na_count,
      zero_frequency_cells = zero_count,
      positive_frequency_cells = positive_count,
      max_frequency = max_freq,
      mean_positive_frequency = ifelse(positive_count > 0, mean(heatmap_matrix[heatmap_matrix > 0], na.rm = TRUE), 0),
      functional_ordering_applied = !is.null(gene_function_map),
      phylogenetic_ordering_applied = !is.null(species_order),
      data_completeness = 1 - (na_count / total_cells),
      parameters = list(
        region_filter = get_task_parameter(task_params, config, "region_filter", "CDS"),
        var_type_filter = get_task_parameter(task_params, config, "var_type_filter", "snp")
      )
    )

    log_message(sprintf("Global heatmap completed: %d positive, %d zero, %d missing out of %d total cells", 
                       positive_count, zero_count, na_count, total_cells))

    return(list(
      plot = p,
      data = heatmap_matrix,
      metadata = metadata,
      gene_categories = if(!is.null(annotation_row)) annotation_row else NULL
    ))

  }, error = function(e) {
    # Error handling and environment diagnostics
    log_message("FATAL ERROR: Global heatmap generation failed", level = "error")
    log_message(paste("Error message:", e$message), level = "error")

    # Print detailed information about variables that might cause indexing issues
    if (exists("all_genes", inherits = FALSE)) {
      log_message("--- DUMPING all_genes VECTOR ---", level = "error")
      capture.output(str(all_genes), file = stderr())
      log_message(paste("Length:", length(all_genes), "Unique:", length(unique(all_genes))), level = "error")
    } else {
      log_message("all_genes variable not found in environment", level = "error")
    }

    if (exists("all_species", inherits = FALSE)) {
      log_message("--- DUMPING all_species VECTOR ---", level = "error")
      capture.output(str(all_species), file = stderr())
      log_message(paste("Length:", length(all_species), "Unique:", length(unique(all_species))), level = "error")
    } else {
      log_message("all_species variable not found in environment", level = "error")
    }

    if (exists("heatmap_matrix", inherits = FALSE)) {
      log_message("--- DUMPING heatmap_matrix DIMENSIONS ---", level = "error")
      log_message(paste("Dimensions:", paste(dim(heatmap_matrix), collapse = " x ")), level = "error")

      log_message("--- DUMPING heatmap_matrix ROWNAMES ---", level = "error")
      capture.output(str(rownames(heatmap_matrix)), file = stderr())
      log_message(paste("Length:", length(rownames(heatmap_matrix)), "Unique:", length(unique(rownames(heatmap_matrix)))), level = "error")

      log_message("--- DUMPING heatmap_matrix COLNAMES ---", level = "error")
      capture.output(str(colnames(heatmap_matrix)), file = stderr())
      log_message(paste("Length:", length(colnames(heatmap_matrix)), "Unique:", length(unique(colnames(heatmap_matrix)))), level = "error")

      # Check for mismatched genes between sorting vector and matrix
      if (exists("all_genes", inherits = FALSE)) {
        mismatched_genes <- setdiff(all_genes, rownames(heatmap_matrix))
        if (length(mismatched_genes) > 0) {
          log_message("--- MISMATCHED GENES DETECTED ---", level = "error")
          log_message(paste("The following", length(mismatched_genes), "gene(s) exist in the sorting vector but NOT in the matrix rownames:"), level = "error")
          log_message(paste(mismatched_genes, collapse = ", "), level = "error")
        } else {
          log_message("--- NO MISMATCHED GENES DETECTED ---", level = "error")
        }
      }

      # Check for mismatched species between sorting vector and matrix
      if (exists("all_species", inherits = FALSE)) {
        mismatched_species <- setdiff(all_species, colnames(heatmap_matrix))
        if (length(mismatched_species) > 0) {
          log_message("--- MISMATCHED SPECIES DETECTED ---", level = "error")
          log_message(paste("The following", length(mismatched_species), "species exist in the sorting vector but NOT in the matrix colnames:"), level = "error")
          log_message(paste(mismatched_species, collapse = ", "), level = "error")
        } else {
          log_message("--- NO MISMATCHED SPECIES DETECTED ---", level = "error")
        }
      }
    } else {
      log_message("heatmap_matrix variable not found in environment", level = "error")
    }

    # Re-throw the error to stop processing
    stop(e)
  })
}
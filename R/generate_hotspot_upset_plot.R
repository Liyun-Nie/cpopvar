#' Generate hotspot UpSet plot
#' 
#' @title Generate hotspot UpSet plot
#' @description Creates UpSet plot to visualize hotspot intersections across species.
#' Supports both CDS-based (gene) and IGS-based (poiGS_ID) feature analysis.
#' @param candidate_hotspots Data containing identified hotspot candidates
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param feature_col Name of the feature column to use for grouping 
#'   (default: "gene" for M03, use "poiGS_ID" for M05)
#' @return List containing UpSet plot and intersection data
#' @importFrom rlang sym
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr rowwise
#' @export
generate_hotspot_upset_plot <- function(candidate_hotspots, config, task_params = NULL, 
                                       feature_col = "gene") {
  
  # Load required package for UpSet plot generation
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' is required for UpSet plot generation but not available")
  }
  
  tryCatch({
    log_message(sprintf("Generating UpSet plot for hotspot sharing (feature: %s)", feature_col))
    
    # Declare identity, get structured color package (V4 architecture)
    plot_colors <- tryCatch({
      get_color_palette(config, "upset")
    }, error = function(e) {
      log_message(sprintf("UpSet plot colors not found, using defaults: %s", e$message), level = "warning")
      # Fallback colors for UpSet plot
      list(
        sets_bar_color = "#2E86AB",
        main_bar_color = "#F18F01", 
        matrix_dot_color = "#DC3545"
      )
    })
    
    # Get parameters
    top_genes_per_species <- get_task_parameter(task_params, config, "top_genes_per_species", 10)
    min_degree <- get_task_parameter(task_params, config, "min_degree", 3)
    
    # Use parameterized feature column (gene for M03, poiGS_ID for M05)
    gene_col <- feature_col
    species_col <- "species"
    
    # Create binary matrix for UpSet plot  
    # Rows = genes, Columns = species, Values = 1 if gene is hotspot in species, 0 otherwise
    
    # Add data validation before pivot_wider
    log_message(sprintf("Input data for UpSet matrix: %d rows", nrow(candidate_hotspots)))
    log_message(sprintf("Unique features (%s): %d, Unique species: %d", 
                       feature_col,
                       length(unique(candidate_hotspots[[gene_col]])),
                       length(unique(candidate_hotspots[[species_col]]))))
    
    upset_matrix <- candidate_hotspots %>%
      # Clean data before pivot_wider
      dplyr::filter(!is.na(.data[[species_col]]) & .data[[species_col]] != "") %>%
      dplyr::filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
      dplyr::select(.data[[gene_col]], .data[[species_col]]) %>%
      dplyr::distinct() %>%
      dplyr::mutate(is_hotspot = 1) %>%
      tidyr::pivot_wider(names_from = !!sym(species_col), values_from = is_hotspot, values_fill = 0)
    
    # Critical validation: ensure we have species columns after pivot_wider
    if (ncol(upset_matrix) <= 1) {
      log_message(sprintf("ERROR: UpSet matrix has insufficient columns (%d). Expected > 1.", ncol(upset_matrix)), level = "error")
      log_message("This suggests pivot_wider failed to create species columns", level = "error")
      return(NULL)
    }
    
    log_message(sprintf("UpSet matrix created: %d features (%s) x %d species", nrow(upset_matrix), feature_col, ncol(upset_matrix) - 1))
    
    # Load phylogenetic species ordering directly for matrix reordering
    species_order <- NULL
    species_ordering_message <- "Species ordering not attempted"
    
    tryCatch({
      # Get session paths for species ordering file
      session_paths <- get_session_paths(config$session_info$session_id)
      species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
      
      # Check if species ordering file exists
      if (file.exists(species_order_file)) {
        log_message(sprintf("Loading species ordering from: %s", species_order_file))
        species_order_data <- read.csv(species_order_file, stringsAsFactors = FALSE)
        
        # Extract species order from the first column
        if (ncol(species_order_data) >= 1) {
          species_order <- as.character(species_order_data[[1]])
          species_ordering_message <- sprintf("Loaded %d species for ordering", length(species_order))
          log_message(species_ordering_message)
        } else {
          log_message("Species ordering file is empty", level = "warning")
        }
      } else {
        log_message(sprintf("Species ordering file not found: %s", species_order_file), level = "warning")
      }
    }, error = function(e) {
      log_message(sprintf("WARNING: Failed to load species ordering: %s", e$message), level = "warning")
      species_order <- NULL
      species_ordering_message <- sprintf("Species ordering failed: %s", e$message)
    })
    
    if (!is.null(species_order) && length(species_order) > 0) {
      # Reorder columns according to phylogenetic order
      matrix_columns <- colnames(upset_matrix)[-1]  # Exclude 'gene' column
      available_species <- intersect(species_order, matrix_columns)
      
      log_message(sprintf("Species ordering validation: %d total species, %d available for ordering", 
                         length(matrix_columns), length(available_species)))
      
      if (length(available_species) > 1) {
        missing_species <- setdiff(matrix_columns, available_species)
        final_species_order <- c(available_species, missing_species)
        
        # Validate final order before selection
        if (length(final_species_order) > 0 && all(final_species_order %in% matrix_columns)) {
          # Reorder the upset_matrix columns
          upset_matrix <- upset_matrix %>%
            dplyr::select(!!sym(gene_col), dplyr::all_of(final_species_order))
          
          log_message(sprintf("Applied phylogenetic ordering to %d species", length(available_species)))
        } else {
          log_message("Species ordering validation failed, keeping original order", level = "warning")
        }
      } else {
        log_message("Insufficient species for ordering, keeping original order", level = "warning")
      }
    } else {
      log_message("Species ordering not available, using default order", level = "warning")
    }
    
    # Convert to data frame with gene_id as row names and ensure numeric values
    upset_df <- upset_matrix %>%
      tibble::column_to_rownames(gene_col) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(as.character(.))))  # Ensure all values are numeric
    
    if (nrow(upset_df) == 0) {
      log_message("No data available for UpSet plot", level = "warning")
      return(NULL)
    }
    
    # Ensure all values are 0 or 1 for binary matrix
    upset_df[is.na(upset_df)] <- 0
    upset_df[upset_df != 0 & upset_df != 1] <- 1
    
    log_message(sprintf("UpSet plot data: %d features (%s) across %d species", nrow(upset_df), feature_col, ncol(upset_df)))
    
    # Generate intersection details for CSV output
    intersection_details <- list()
    gene_names <- rownames(upset_df)
    species_names <- colnames(upset_df)
    
    # Calculate intersections for each possible combination
    for (i in 1:nrow(upset_df)) {
      gene <- gene_names[i]
      # Safely check for presence (value == 1) with type checking
      row_data <- as.numeric(as.character(upset_df[i, ]))
      # Additional safety check - replace any NAs with 0
      row_data[is.na(row_data)] <- 0
      species_with_gene <- species_names[row_data == 1]
      if (length(species_with_gene) > 0) {
        intersection_key <- paste(sort(species_with_gene), collapse = "_")
        if (is.null(intersection_details[[intersection_key]])) {
          intersection_details[[intersection_key]] <- list(
            species = species_with_gene,
            genes = character(),
            count = 0
          )
        }
        intersection_details[[intersection_key]]$genes <- c(intersection_details[[intersection_key]]$genes, gene)
        intersection_details[[intersection_key]]$count <- intersection_details[[intersection_key]]$count + 1
      }
    }
    
    # Create intersection CSV data (feature_count instead of gene_count for generality)
    feature_count_col <- paste0(gsub("_ID$", "", feature_col), "_count")
    features_col <- gsub("_ID$", "", feature_col)
    
    intersection_csv_data <- data.frame(
      intersection_id = character(),
      num_species = integer(),
      species_combination = character(),
      stringsAsFactors = FALSE
    )
    intersection_csv_data[[feature_count_col]] <- integer()
    intersection_csv_data[[features_col]] <- character()
    
    for (i in seq_along(intersection_details)) {
      detail <- intersection_details[[i]]
      row_data <- data.frame(
        intersection_id = names(intersection_details)[i],
        num_species = length(detail$species),
        species_combination = paste(detail$species, collapse = ", "),
        stringsAsFactors = FALSE
      )
      row_data[[feature_count_col]] <- detail$count
      row_data[[features_col]] <- paste(detail$genes, collapse = ", ")
      intersection_csv_data <- rbind(intersection_csv_data, row_data)
    }
    
    # Sort by feature count (descending) and number of species (descending)
    intersection_csv_data <- intersection_csv_data %>%
      dplyr::arrange(dplyr::desc(!!sym(feature_count_col)), dplyr::desc(num_species))
    
    # NOTE: File I/O operations moved to orchestrator layer for architecture compliance
    log_message(sprintf("Generated UpSet intersection analysis with %d intersections", nrow(intersection_csv_data)))
    
    # Final validation before UpSetR call
    if (ncol(upset_df) == 0) {
      log_message("ERROR: upset_df has no columns - cannot create UpSet plot", level = "error")
      return(NULL)
    }
    
    species_columns <- colnames(upset_df)
    log_message(sprintf("Final UpSet data validation: %d features (%s), %d species columns", nrow(upset_df), feature_col, length(species_columns)))
    log_message(sprintf("Species columns: %s", paste(species_columns, collapse = ", ")))
    
    # Additional validation for UpSetR parameters
    if (length(species_columns) == 0) {
      log_message("ERROR: No species columns available for UpSet plot", level = "error")
      return(NULL)
    }
    
    # Ensure species_columns is character vector and not empty
    species_columns <- as.character(species_columns)
    reversed_species <- rev(species_columns)
    
    if (length(reversed_species) == 0) {
      log_message("ERROR: Reversed species vector is empty", level = "error")
      return(NULL)
    }
    
    log_message(sprintf("UpSetR sets parameter will be: %s", paste(reversed_species, collapse = ", ")))
    
    # Create UpSet plot with forced species order
    # Use sets parameter to force the order of species on the left side
    upset_plot <- UpSetR::upset(
      upset_df,
      nsets = ncol(upset_df),
      nintersects = NA,  # Show ALL intersections (no limit) - ensures all patterns are visible
      order.by = "freq",  # Order by intersection size (number of features) - bars sorted by height for readability
      decreasing = TRUE,  # Show largest intersections first (tallest bars first)
      sets = reversed_species,  # Use verified, non-empty character vector
      sets.bar.color = plot_colors[["sets_bar_color"]],
      main.bar.color = plot_colors[["main_bar_color"]],
      matrix.color = plot_colors[["matrix_dot_color"]],
      sets.x.label = "Number of Hotspot Genes per Species",
      mainbar.y.label = "Number of Shared Hotspot Genes",
      text.scale = c(1.3, 1.3, 1, 1, 2, 1.3),
      keep.order = TRUE  # Try to maintain the order
    )
    
    # Calculate sharing statistics using robust base R approach
    # Using rowSums instead of c_across for better reliability with numeric matrices
    row_sums <- rowSums(upset_df, na.rm = TRUE)
    stats_table <- table(row_sums)
    sharing_stats <- data.frame(
      num_species = as.numeric(names(stats_table)),
      num_genes = as.numeric(stats_table),
      percentage = round(as.numeric(stats_table) / sum(stats_table) * 100, 1),
      stringsAsFactors = FALSE
    )
    
    metadata <- list(
      total_genes = nrow(upset_df),
      total_species = ncol(upset_df),
      sharing_stats = sharing_stats,
      intersection_details = intersection_csv_data,
      parameters = list(
        top_genes_per_species = top_genes_per_species,
        min_degree = min_degree
      )
    )
    
    log_message(sprintf("UpSet plot generated successfully with %d intersections", nrow(intersection_csv_data)))
    
    return(list(
      plot = upset_plot,
      data = upset_df,
      intersection_data = intersection_csv_data,  # Include intersection details for caller
      metadata = metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate UpSet plot: %s", e$message), level = "error")
    return(NULL)
  })
}

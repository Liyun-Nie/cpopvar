############################################################
#### load_gene_function_mapping.R - Gene Function Data Loader ####
############################################################
#
# Pure utility function to load gene function mapping data
# Follows Rule #6: Pure utility with dependency injection
# Follows Rule #10: Session ID from config object only
#
############################################################

#' Load Gene Function Mapping Data
#' 
#' @param config Configuration object containing session information
#' @return Gene function mapping data frame or NULL if not found
#' 
#' ARCHITECTURE: Pure utility function with dependency injection
#' Session ID obtained from config$session_info$session_id per Rule #10
#' load_gene_function_mapping
#' @export
load_gene_function_mapping <- function(config) {
  
  tryCatch({
    # Get session ID from config object (Rule #10 compliance)
    session_id <- config$session_info$session_id
    
    # V11 ARCHITECTURAL REFACTOR: Robust path finding for package data
    # Priority Order:
    # 1. User-uploaded file in the session directory.
    # 2. Path specified in the main config file.
    # 3. Default file included with the package.

    gene_function_file <- NULL

    # 1. Check for session-specific file
    if (!is.null(session_id)) {
      session_paths <- get_session_paths(session_id, create_dirs = FALSE)
      session_file <- file.path(session_paths$raw, "gene_function_map.csv")
      if (file.exists(session_file)) {
        gene_function_file <- session_file
        log_message("Using session-specific gene function map.", level = "info")
      }
    }

    # 2. Check for path in config if not already found
    # V3.16: Simplified configuration - direct file path
    if (is.null(gene_function_file) && !is.null(config$input_files$gene_function_mapping)) {
      config_path <- config$input_files$gene_function_mapping
      if (file.exists(config_path)) {
        gene_function_file <- config_path
        log_message(sprintf("Using gene function map from config: %s", config_path), level = "info")
      } else {
        log_message(sprintf("Path in config not found: %s", config_path), level = "warning")
      }
    }

    # 3. Fallback to default package file if still not found
    if (is.null(gene_function_file)) {
      package_file <- system.file("config", "gene_function_map.csv", package = "cpopvar")
      if (file.exists(package_file)) {
        gene_function_file <- package_file
        log_message("Using default gene function map from package.", level = "info")
      }
    }
    
    if (is.null(gene_function_file)) {
      log_message("Gene function mapping file not found", level = "warning")
      return(NULL)
    }
    
    # Read gene function mapping
    gene_function_data <- read.csv(gene_function_file, stringsAsFactors = FALSE)
    
    log_message(sprintf("Loaded gene function file with columns: %s", paste(names(gene_function_data), collapse = ", ")))
    
    # Check for gene_id column (original format) and rename if needed
    if ("gene_id" %in% names(gene_function_data) && !"gene" %in% names(gene_function_data)) {
      gene_function_data <- gene_function_data %>%
dplyr::rename(gene = gene_id)
      log_message("Renamed 'gene_id' column to 'gene' for consistency")
    }
    
    # Validate required columns
    required_cols <- c("gene", "functional_group", "gene_category")
    missing_cols <- setdiff(required_cols, names(gene_function_data))
    if (length(missing_cols) > 0) {
      log_message(sprintf("Gene function file missing columns: %s (available: %s)", 
                         paste(missing_cols, collapse = ", "), 
                         paste(names(gene_function_data), collapse = ", ")), level = "warning")
      return(NULL)
    }
    
    # Clean up gene names (remove any whitespace)
    gene_function_data$gene <- trimws(gene_function_data$gene)
    
    log_message(sprintf("Successfully loaded gene function mapping: %d genes with functional annotations", nrow(gene_function_data)))
    log_message(sprintf("Gene categories available: %s", paste(unique(gene_function_data$gene_category), collapse = ", ")))
    log_message(sprintf("Functional groups available (not used for enrichment): %s", paste(unique(gene_function_data$functional_group), collapse = ", ")))
    
    return(gene_function_data)
    
  }, error = function(e) {
    log_message(sprintf("Failed to load gene function mapping: %s", e$message), level = "error")
    return(NULL)
  })
}
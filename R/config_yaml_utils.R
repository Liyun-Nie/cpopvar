############################################################
#### config_yaml_utils.R - YAML Configuration Utilities ####
############################################################
#
# Utilities for working with YAML configuration files
# Includes validation, loading, and exporting functions
#
############################################################

# Load required libraries

#' Validate YAML configuration file
#' 
#' @param config_file Path to YAML configuration file
#' @param return_errors If TRUE, return detailed validation errors
#' @return TRUE if valid, FALSE if invalid, or list of errors if return_errors=TRUE
#' validate_yaml_config
#' @export
validate_yaml_config <- function(config_file, return_errors = FALSE) {
  
  validation_errors <- c()
  
  # Check if file exists
  if (!file.exists(config_file)) {
    error_msg <- paste("Configuration file not found:", config_file)
    if (return_errors) return(list(valid = FALSE, errors = error_msg))
    stop(error_msg, call. = FALSE)
  }
  
  # Try to parse YAML
  config_data <- tryCatch({
    yaml::read_yaml(config_file)
  }, error = function(e) {
    error_msg <- paste("Failed to parse YAML file:", e$message)
    if (return_errors) {
      validation_errors <<- c(validation_errors, error_msg)
      return(NULL)
    } else {
      stop(error_msg, call. = FALSE)
    }
  })
  
  if (is.null(config_data) && return_errors) {
    return(list(valid = FALSE, errors = validation_errors))
  }
  
  # Use the same top-level contract as the parser and dispatcher.
  structure_result <- validate_task_config_structure(config_data)
  validation_errors <- c(validation_errors, structure_result$errors)
  
  # Return results
  if (return_errors) {
    return(list(
      valid = length(validation_errors) == 0,
      errors = validation_errors,
      config_data = if (length(validation_errors) == 0) config_data else NULL
    ))
  } else {
    if (length(validation_errors) > 0) {
      stop(paste("Configuration validation failed:\n", 
                paste(validation_errors, collapse = "\n")), call. = FALSE)
    }
    return(TRUE)
  }
}

#' Load and validate YAML configuration
#' 
#' @param config_file Path to YAML configuration file
#' @return Validated configuration data
load_yaml_config <- function(config_file) {
  
  validation_result <- validate_yaml_config(config_file, return_errors = TRUE)
  
  if (!validation_result$valid) {
    stop(paste("Configuration validation failed:\n", 
              paste(validation_result$errors, collapse = "\n")), call. = FALSE)
  }
  
  return(validation_result$config_data)
}

#' Convert Shiny reactive values to YAML configuration
#' 
#' @param analysis_data Shiny reactiveValues object containing user selections
#' @param input Shiny input object containing UI parameters
#' @param session_paths Session directory paths
#' @return Configuration list ready for YAML export
shiny_to_yaml_config <- function(analysis_data, input, session_paths) {
  
  # Generate session info
  session_info <- list(
    session_id = basename(session_paths$base),
    run_name = input$run_name %||% "Shiny Frontend Export",
    description = input$run_description %||% "Configuration exported from Shiny frontend",
    timestamp = as.character(Sys.time())
  )
  
  # Collect input file information
  input_files <- list()
  
  # Main data file
  if (!is.null(analysis_data$raw_uploaded_data)) {
    input_files$main_data <- list(
      path = file.path(session_paths$raw, "all_combined_data.csv"),
      description = "Main variant data file uploaded via Shiny"
    )
  }
  
  # Other files
  file_types <- c("cds_lengths", "region_info", "genome_regions", "gene_function", "group_info")
  for (file_type in file_types) {
    if (!is.null(analysis_data$preview_data[[file_type]])) {
      input_files[[file_type]] <- list(
        path = file.path(session_paths$raw, paste0(file_type, ".csv")),
        description = paste("File uploaded via Shiny:", file_type),
        required = file_type %in% c("cds_lengths", "region_info", "genome_regions")
      )
    }
  }
  
  # Column mappings (if available)
  column_mappings <- list()
  if (!is.null(analysis_data$user_mapping)) {
    column_mappings$main_data <- analysis_data$user_mapping
  } else if (!is.null(analysis_data$suggested_mapping)) {
    column_mappings$main_data <- analysis_data$suggested_mapping
  } else {
    # Default mappings
    column_mappings$main_data <- list(
      species = "Species",
      variant_type = "Var_type", 
      position = "Position",
      gene = "Gene",
      region_type = "Region_type"
    )
  }
  
  # Default mappings for other file types
  column_mappings$cds_lengths <- list(
    gene_name = "Gene",
    gene_length = "Length"
  )
  
  column_mappings$region_info <- list(
    gene_name = "Gene",
    region_length = "Length",
    region_type = "Type"
  )
  
  # Extract parameters from Shiny inputs
  parameters <- list(
    filtering = list(
      snp_threshold = input$snp_threshold %||% 3,
      indel_threshold = input$indel_threshold %||% 3,
      complex_threshold = input$complex_threshold %||% 3,
      mnp_threshold = input$mnp_threshold %||% 3,
      remove_duplicates = input$remove_duplicates %||% TRUE,
      handle_missing_genes = input$handle_missing_genes %||% "keep_as_na",
      min_species_coverage = input$min_species_coverage %||% 1
    ),
    
    normalization = list(
      frequency_unit = input$frequency_unit %||% 1000,
      calculation_method = input$calculation_method %||% "simple",
      zero_frequency_as_na = input$zero_frequency_as_na %||% FALSE
    ),
    
    visualization = list(
      correlation_method = input$correlation_method %||% "spearman",
      hotspot_freq_percentile = input$hotspot_freq_percentile %||% 0.75,
      hotspot_min_species_ratio = input$hotspot_min_species_ratio %||% 0.5,
      
      module_02 = list(
        primary_factor = input$m02_primary_factor %||% "species",
        secondary_factor = input$m02_secondary_factor %||% "var_type",
        plot_type = input$m02_plot_type %||% "boxplot"
      ),
      
      module_03 = list(
        analysis_scope = input$m03_analysis_scope %||% "all",
        cluster_method = input$m03_cluster_method %||% "ward.D2"
      )
    )
  )
  
  # Module selection (check for module checkbox inputs)
  modules <- list(
    preprocessing = TRUE,
    filtering = TRUE,
    normalization = TRUE,
    module_01_distribution = input$enable_module_m01 %||% input$enable_module_01 %||% TRUE,
    module_02_comparative = input$enable_module_m02 %||% input$enable_module_02 %||% TRUE,
    module_03_hotspot = input$enable_module_m03 %||% input$enable_module_03 %||% TRUE,
    validation = input$enable_validation %||% FALSE
  )
  
  # Output configuration
  output <- list(
    save_csv = input$save_csv %||% TRUE,
    save_plots_png = input$save_plots_png %||% TRUE,
    add_timestamp = input$add_timestamp %||% TRUE,
    save_intermediate = input$save_intermediate %||% TRUE
  )
  
  # Combine all sections
  config_data <- list(
    session_info = session_info,
    input_files = input_files,
    column_mappings = column_mappings,
    parameters = parameters,
    modules = modules,
    output = output
  )
  
  return(config_data)
}

#' Export configuration to YAML file
#' 
#' @param config_data Configuration data list
#' @param output_file Path for output YAML file
#' @return TRUE if successful
export_yaml_config <- function(config_data, output_file) {
  
  tryCatch({
    # Ensure output directory exists
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Write YAML file
    yaml::write_yaml(config_data, output_file, 
                     handlers = list(
                       logical = function(x) if (x) "true" else "false"
                     ))
    
    return(TRUE)
    
  }, error = function(e) {
    warning(paste("Failed to export YAML configuration:", e$message))
    return(FALSE)
  })
}

#' Create a quick configuration file from minimal parameters
#' 
#' @param input_data_path Path to main data file
#' @param session_id Session identifier
#' @param output_file Path for output configuration file
#' @return Configuration file path
create_quick_config <- function(input_data_path, session_id = NULL, output_file = NULL) {
  
  if (is.null(session_id)) {
    session_id <- paste0("quick_run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  if (is.null(output_file)) {
    output_file <- paste0("config/quick_config_", session_id, ".yml")
  }
  
  # Basic configuration
  config_data <- list(
    session_info = list(
      session_id = session_id,
      run_name = "Quick Analysis Run",
      description = "Automatically generated configuration"
    ),
    
    input_files = list(
      main_data = list(path = input_data_path),
      cds_lengths = list(path = "data/raw/output_cds_lengths.csv"),
      region_info = list(path = "data/raw/region_length_statistics/output_gene_intergenic_intron_pos_length.csv"),
      genome_regions = list(path = system.file("config", "species_genome_regions.csv", package = "cpopvar"))
    ),
    
    column_mappings = list(
      main_data = list(
        species = "Species",
        variant_type = "Var_type",
        position = "Position",
        gene = "Gene",
        region_type = "Region_type"
      ),
      
      cds_lengths = list(
        gene_name = "Gene",
        gene_length = "Length"
      ),
      
      region_info = list(
        gene_name = "Gene",
        region_length = "Length",
        region_type = "Type"
      )
    ),
    
    parameters = list(
      filtering = list(
        snp_threshold = 3,
        indel_threshold = 3,
        complex_threshold = 3,
        mnp_threshold = 3
      ),
      normalization = list(
        frequency_unit = 1000,
        calculation_method = "simple"
      ),
      visualization = list(
        correlation_method = "spearman",
        hotspot_freq_percentile = 0.75
      )
    ),
    
    modules = list(
      preprocessing = TRUE,
      filtering = TRUE,
      normalization = TRUE,
      module_01_distribution = TRUE,
      module_02_comparative = TRUE,
      module_03_hotspot = TRUE
    ),
    
    output = list(
      save_csv = TRUE,
      save_plots_png = TRUE,
      save_intermediate = TRUE
    )
  )
  
  # Export configuration
  if (export_yaml_config(config_data, output_file)) {
    cat("Quick configuration created:", output_file, "\n")
    return(output_file)
  } else {
    stop("Failed to create quick configuration file", call. = FALSE)
  }
}

# cat("✓ YAML configuration utilities loaded successfully\n")
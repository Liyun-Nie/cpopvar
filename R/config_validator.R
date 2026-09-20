############################################################
#### config_validator.R - Configuration Validation System ####
############################################################
#
# Comprehensive configuration validation based on schema definitions
# Supports both CSV and R list configurations with detailed error reporting
#
############################################################

# Required imports handled by @importFrom directives

# Source utilities
# source("src/utils/common_utils.R")

#' Load and parse configuration schema
#' 
#' @param schema_file Path to YAML schema file
#' @return Parsed schema list
#' @importFrom yaml read_yaml
#' @export
load_config_schema <- function(schema_file = system.file("config", "species_genome_schema.yml", package = "cpopvar")) {
  
  if (!file.exists(schema_file)) {
    log_message(sprintf("WARNING: Schema file not found: %s", schema_file), level = "warning")
    return(NULL)
  }
  
  tryCatch({
    schema <- yaml::read_yaml(schema_file)
    log_message(sprintf("Successfully loaded configuration schema version %s", schema$schema_version))
    return(schema)
  }, error = function(e) {
    log_message(sprintf("ERROR: Failed to load schema file: %s", e$message), level = "error")
    return(NULL)
  })
}

#' Validate species configuration against schema
#' 
#' @param config_data Configuration data (data.frame or list)
#' @param schema_file Path to schema file
#' @param detailed_report Whether to return detailed validation report
#' @return List containing validation results
#' @importFrom magrittr %>%
#' @export
validate_species_configuration <- function(config_data, 
                                         schema_file = system.file("config", "species_genome_schema.yml", package = "cpopvar"),
                                         detailed_report = TRUE) {
  
  log_message("Starting comprehensive configuration validation")
  
  # Load schema
  schema <- load_config_schema(schema_file)
  if (is.null(schema)) {
    return(list(
      valid = FALSE,
      errors = "Failed to load validation schema",
      warnings = character(0),
      report = NULL
    ))
  }
  log_message(sprintf("Successfully loaded configuration schema version %s", schema$schema_version))
  
  # 2. Data Migration: Handle legacy 'special_handling' field
  # This step ensures backward compatibility by creating the new 'genome_structure_type' column
  # from the old 'special_handling' column if it doesn't already exist.
  if (!"genome_structure_type" %in% names(config_data) && "special_handling" %in% names(config_data)) {
    log_message("Performing legacy migration: 'special_handling' -> 'genome_structure_type'", level = "INFO")
    
    # Get migration mappings from schema, with a fallback for safety
    migration_map <- schema$legacy_migration$from_special_handling
    if (is.null(migration_map)) {
      errors <- c(errors, "FATAL: 'legacy_migration' rules not found in schema.")
      return(list(valid = FALSE, errors = errors, warnings = warnings))
    }
    
    # Replace NA in special_handling with "" to match the schema key
    config_data <- config_data %>%
      dplyr::mutate(special_handling = ifelse(is.na(special_handling), "", special_handling))

    # Apply the mapping
    config_data <- config_data %>%
      dplyr::mutate(genome_structure_type = purrr::map_chr(special_handling, ~migration_map[[.x]] %||% "unknown"))

    # Check if any unknown types were generated
    if (any(config_data$genome_structure_type == "unknown")) {
      warnings <- c(warnings, "Some 'special_handling' values had no corresponding migration rule and were mapped to 'unknown'.")
    }
  }
  
  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    report = list()
  )
  
  # Convert different input types to standardized format
  if (is.data.frame(config_data)) {
    config_df <- config_data
  } else if (is.list(config_data) && !is.null(names(config_data))) {
    # Convert from nested list (SPECIES_GENOME_REGIONS format) to data.frame
    config_df <- convert_nested_list_to_df(config_data)
  } else {
    validation_result$valid <- FALSE
    validation_result$errors <- c(validation_result$errors, 
                                 "Configuration data must be a data.frame or named list")
    return(validation_result)
  }
  
  # Validate column structure
  column_validation <- validate_columns(config_df, schema$columns)
  validation_result$valid <- validation_result$valid && column_validation$valid
  validation_result$errors <- c(validation_result$errors, column_validation$errors)
  validation_result$warnings <- c(validation_result$warnings, column_validation$warnings)
  
  # Validate data constraints
  constraint_validation <- validate_constraints(config_df, schema$columns)
  validation_result$valid <- validation_result$valid && constraint_validation$valid
  validation_result$errors <- c(validation_result$errors, constraint_validation$errors)
  validation_result$warnings <- c(validation_result$warnings, constraint_validation$warnings)
  
  # Validate business rules
  rules_validation <- validate_business_rules(config_df, schema$validation_rules)
  validation_result$valid <- validation_result$valid && rules_validation$valid
  validation_result$errors <- c(validation_result$errors, rules_validation$errors)
  validation_result$warnings <- c(validation_result$warnings, rules_validation$warnings)
  
  # Generate detailed report if requested
  if (detailed_report) {
    validation_result$report <- generate_validation_report(config_df, schema, validation_result)
  }
  
  # Log summary
  if (validation_result$valid) {
    log_message("Configuration validation PASSED")
  } else {
    log_message(sprintf("Configuration validation FAILED with %d errors", 
                       length(validation_result$errors)), level = "error")
  }
  
  if (length(validation_result$warnings) > 0) {
    log_message(sprintf("Configuration validation generated %d warnings", 
                       length(validation_result$warnings)), level = "warning")
  }
  
  return(validation_result)
}

#' Validate column structure against schema
#' 
#' @param config_df Configuration data frame
#' @param column_schema Schema column definitions
#' @return Validation result list
validate_columns <- function(config_df, column_schema) {
  
  result <- list(valid = TRUE, errors = character(0), warnings = character(0))
  
  # Check required columns
  required_columns <- names(column_schema)[sapply(column_schema, function(x) x$required)]
  missing_columns <- setdiff(required_columns, names(config_df))
  
  if (length(missing_columns) > 0) {
    result$valid <- FALSE
    result$errors <- c(result$errors, 
                      sprintf("Missing required columns: %s", 
                             paste(missing_columns, collapse = ", ")))
  }
  
  # Check for deprecated columns
  deprecated_columns <- names(column_schema)[sapply(column_schema, function(x) isTRUE(x$deprecated))]
  present_deprecated <- intersect(deprecated_columns, names(config_df))
  
  if (length(present_deprecated) > 0) {
    result$warnings <- c(result$warnings,
                        sprintf("Using deprecated columns: %s", 
                               paste(present_deprecated, collapse = ", ")))
  }
  
  # Check column types
  for (col_name in intersect(names(column_schema), names(config_df))) {
    expected_type <- column_schema[[col_name]]$type
    actual_data <- config_df[[col_name]]
    
    type_valid <- switch(expected_type,
      "character" = is.character(actual_data),
      "numeric" = is.numeric(actual_data),
      "logical" = is.logical(actual_data),
      TRUE  # Default to valid for unknown types
    )
    
    if (!type_valid) {
      result$valid <- FALSE
      result$errors <- c(result$errors,
                        sprintf("Column '%s' should be %s but is %s", 
                               col_name, expected_type, class(actual_data)[1]))
    }
  }
  
  return(result)
}

#' Validate data constraints
#' 
#' @param config_df Configuration data frame
#' @param column_schema Schema column definitions
#' @return Validation result list
validate_constraints <- function(config_df, column_schema) {
  
  result <- list(valid = TRUE, errors = character(0), warnings = character(0))
  
  for (col_name in names(column_schema)) {
    if (!col_name %in% names(config_df)) next
    
    col_data <- config_df[[col_name]]
    constraints <- column_schema[[col_name]]$constraints
    
    if (is.null(constraints)) next
    
    # Validate enum constraints
    if (!is.null(constraints$enum)) {
      invalid_values <- setdiff(col_data[!is.na(col_data)], constraints$enum)
      if (length(invalid_values) > 0) {
        result$valid <- FALSE
        result$errors <- c(result$errors,
                          sprintf("Column '%s' contains invalid values: %s. Allowed: %s",
                                 col_name, 
                                 paste(unique(invalid_values), collapse = ", "),
                                 paste(constraints$enum, collapse = ", ")))
      }
    }
    
    # Validate uniqueness
    if (isTRUE(constraints$unique)) {
      duplicates <- col_data[duplicated(col_data) & !is.na(col_data)]
      if (length(duplicates) > 0) {
        result$valid <- FALSE
        result$errors <- c(result$errors,
                          sprintf("Column '%s' contains duplicate values: %s",
                                 col_name, paste(unique(duplicates), collapse = ", ")))
      }
    }
    
    # Validate numeric constraints
    if (column_schema[[col_name]]$type == "numeric") {
      numeric_data <- col_data[!is.na(col_data)]
      
      if (!is.null(constraints$min) && any(numeric_data < constraints$min)) {
        result$valid <- FALSE
        result$errors <- c(result$errors,
                          sprintf("Column '%s' contains values below minimum %s",
                                 col_name, constraints$min))
      }
      
      if (!is.null(constraints$max) && any(numeric_data > constraints$max)) {
        result$valid <- FALSE
        result$errors <- c(result$errors,
                          sprintf("Column '%s' contains values above maximum %s",
                                 col_name, constraints$max))
      }
    }
  }
  
  return(result)
}

#' Validate business rules
#' 
#' @param config_df Configuration data frame
#' @param rules_schema Business rules from schema
#' @return Validation result list
validate_business_rules <- function(config_df, rules_schema) {
  
  result <- list(valid = TRUE, errors = character(0), warnings = character(0))
  
  if (is.null(rules_schema)) return(result)
  
  for (rule_name in names(rules_schema)) {
    rule <- rules_schema[[rule_name]]
    
    # Evaluate condition for each row
    for (i in 1:nrow(config_df)) {
      row_data <- config_df[i, ]
      
      # Simple condition evaluation (could be enhanced with proper expression parsing)
      condition_met <- evaluate_simple_condition(row_data, rule$condition)
      
      if (condition_met) {
        # Check requirements
        for (requirement in rule$requirements) {
          req_met <- evaluate_simple_condition(row_data, requirement)
          if (!req_met) {
            result$valid <- FALSE
            result$errors <- c(result$errors,
                              sprintf("Species '%s' violates rule '%s': %s",
                                     row_data$species, rule_name, requirement))
          }
        }
      }
    }
  }
  
  return(result)
}

#' Simple condition evaluator (basic implementation)
#' 
#' @param row_data Single row of configuration data
#' @param condition String condition to evaluate
#' @return Logical result
evaluate_simple_condition <- function(row_data, condition) {
  # A robust evaluator that uses an environment to avoid string substitution issues.
  tryCatch({
    # Create an environment from the row data, ensuring it inherits from the base
    # environment so that standard functions like `==` and `is.na` are available.
    # `list2env` is ideal for this. The previous `as.environment` method created
    # an isolated environment that could not find base R operators.
    eval_env <- list2env(as.list(row_data), parent = baseenv())
    
    # Evaluate the parsed condition expression within the correctly configured environment.
    result <- eval(parse(text = condition), envir = eval_env)
    
    # The result of a condition must be a single TRUE or FALSE.
    return(isTRUE(result))
    
  }, error = function(e) {
    # Log detailed error for easier debugging if evaluation fails.
    warning_message <- sprintf(
      "Rule evaluation failed for species '%s' on condition '%s'. Reason: %s",
      row_data$species, condition, e$message
    )
    log_message(warning_message, level = "warning")
    
    # For safety, if a rule cannot be evaluated, it is considered not met.
    return(FALSE)
  })
}

#' Convert nested list to data frame format
#' 
#' @param nested_list SPECIES_GENOME_REGIONS format list
#' @return Data frame representation
convert_nested_list_to_df <- function(nested_list) {
  
  species_names <- names(nested_list)
  df_rows <- list()
  
  for (species in species_names) {
    config <- nested_list[[species]]
    
    df_rows[[species]] <- data.frame(
      species = species,
      lsc_start = ifelse(is.null(config$LSC), NA, config$LSC["start"]),
      lsc_end = ifelse(is.null(config$LSC), NA, config$LSC["end"]),
      ir_start = ifelse(is.null(config$IR), NA, config$IR["start"]),
      ir_end = ifelse(is.null(config$IR), NA, config$IR["end"]),
      ssc_start = ifelse(is.null(config$SSC), NA, config$SSC["start"]),
      ssc_end = ifelse(is.null(config$SSC), NA, config$SSC["end"]),
      total_length = ifelse(is.null(config$total_length), NA, config$total_length),
      special_handling = ifelse(is.null(config$special_handling), NA, config$special_handling),
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, df_rows))
}

#' Generate detailed validation report
#' 
#' @param config_df Configuration data frame
#' @param schema Configuration schema
#' @param validation_result Validation results
#' @return Detailed report list
generate_validation_report <- function(config_df, schema, validation_result) {
  
  report <- list(
    summary = list(
      total_species = nrow(config_df),
      validation_status = ifelse(validation_result$valid, "PASSED", "FAILED"),
      error_count = length(validation_result$errors),
      warning_count = length(validation_result$warnings),
      schema_version = schema$schema_version
    ),
    species_breakdown = table(config_df$genome_structure_type),
    errors = validation_result$errors,
    warnings = validation_result$warnings,
    recommendations = generate_recommendations(config_df, validation_result)
  )
  
  return(report)
}

#' Generate recommendations based on validation results
#' 
#' @param config_df Configuration data frame
#' @param validation_result Validation results
#' @return Character vector of recommendations
generate_recommendations <- function(config_df, validation_result) {
  
  recommendations <- character(0)
  
  if (length(validation_result$errors) > 0) {
    recommendations <- c(recommendations,
                        "Fix all validation errors before proceeding with analysis")
  }
  
  if (length(validation_result$warnings) > 0) {
    recommendations <- c(recommendations,
                        "Consider addressing validation warnings for improved compatibility")
  }
  
  # Check for deprecated columns
  if ("special_handling" %in% names(config_df)) {
    recommendations <- c(recommendations,
                        "Consider migrating from 'special_handling' to 'genome_structure_type' field")
  }
  
  return(recommendations)
}

#' Quick validation function for interactive use
#' 
#' @param config_path Path to configuration file
#' @param session_id Session ID for session-specific configuration (default: NULL)
#' @return Simple validation result
quick_validate_config <- function(config_path = NULL, session_id = NULL) {
  # Auto-detect config path with session-aware priority
  if (is.null(config_path)) {
    if (!is.null(session_id)) {
      # Try session-specific path first using session paths
      if (exists("get_session_paths")) {
        tryCatch({
          session_paths <- get_session_paths(session_id, create_dirs = FALSE)
          session_path <- file.path(session_paths$raw, "species_genome_regions.csv")
        }, error = function(e) {
          session_path <<- file.path(getwd(), "app_data", "sessions", session_id, "raw", "species_genome_regions.csv")
        })
      } else {
        session_path <- file.path(getwd(), "app_data", "sessions", session_id, "raw", "species_genome_regions.csv")
      }
      if (file.exists(session_path)) {
        config_path <- session_path
      }
    }
    
    # Fallback to default global configuration in package
    if (is.null(config_path)) {
      config_path <- system.file("config", "species_genome_regions.csv", package = "cpopvar")
    }
  }
  
  if (!file.exists(config_path)) {
    cat("ERROR: Configuration file not found:", config_path, "\n")
    return(FALSE)
  }
  
  config_data <- read.csv(config_path, stringsAsFactors = FALSE)
  validation <- validate_species_configuration(config_data, detailed_report = FALSE)
  
  if (validation$valid) {
    cat("[SUCCESS] Configuration validation PASSED\n")
    if (length(validation$warnings) > 0) {
      cat("[WARNING] Warnings:", length(validation$warnings), "\n")
    }
  } else {
    cat("[ERROR] Configuration validation FAILED\n")
    cat("Errors:", length(validation$errors), "\n")
    for (error in validation$errors) {
      cat("  -", error, "\n")
    }
  }
  
  return(validation$valid)
}

# log_message("Configuration validation system loaded successfully")
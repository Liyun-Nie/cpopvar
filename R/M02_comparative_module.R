############################################################
#### M02_comparative_module.R - M02 Module Orchestrator ####
############################################################
#
# Thin orchestrator for M02 Comparative Analysis
# This is a Phase 2 "thin wrapper" that loads dependencies and delegates
# to the existing generate_M02_comparative_plots function
#
# Phase 2 Rule #3 Compliance: Does NOT modify existing function libraries,
# only provides a clean interface for task_dispatcher.R
#
############################################################

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

#' M02 Comparative Analysis Module Orchestrator
#' 
#' @param normalized_data Normalized frequency data
#' @param config Configuration object
#' @param output_dir Output directory for plots
#' @param task_name Task name for output organization
#' @param task_params Task-specific parameters
#' @return M02 analysis results
#' 
#' ARCHITECTURE: Implements M01-style standardized data cleaning and output generation
#' following the established golden standard for module orchestrators.
#' run_m02_comparative
#' @export
run_m02_comparative <- function(normalized_data,
                               config,
                               output_dir,
                               task_name,
                               task_params) {
  
  # Function dependencies are automatically loaded in R package context
  
  log_message(sprintf("=== M02 COMPARATIVE ANALYSIS: %s ===", task_name))
  
  # Initialize composite_plots to ensure variable exists in all code paths
  composite_plots <- list()
  
  # Detect the analysis mode from the configured factors
  # Determine analysis mode: explicit setting or auto-detect
  if (!is.null(task_params$analysis_mode)) {
    # Use explicit user setting
    analysis_mode <- task_params$analysis_mode
    log_message(sprintf("Analysis mode: %s (explicitly set by user)", analysis_mode))
  } else {
    # Auto-detect based on target_var_types
    has_var_types <- !is.null(task_params$data_scoping$target_var_types) && 
                     length(task_params$data_scoping$target_var_types) > 0
    analysis_mode <- if (has_var_types) "by_variant_type" else "by_all_types"
    log_message(sprintf("Analysis mode: %s (auto-detected)", analysis_mode))
    log_message(sprintf("  Reason: target_var_types %s in data_scoping", 
                       if (has_var_types) "specified" else "not specified"))
  }
  
  # Validate analysis_mode
  if (!analysis_mode %in% c("by_variant_type", "by_all_types")) {
    log_message(sprintf("ERROR: Invalid analysis_mode '%s'. Must be 'by_variant_type' or 'by_all_types'", 
                       analysis_mode), level = "error")
    stop("Invalid analysis_mode configuration")
  }
  
  # Filter to the factors requested for this task
  # Filter analysis factors based on mode
  original_single_count <- length(task_params$single_factor)
  original_dual_count <- length(task_params$dual_factor)
  
  if (analysis_mode == "by_all_types") {
    # Filter out var_type-related analyses
    task_params$single_factor <- task_params$single_factor[
      task_params$single_factor != "var_type"
    ]
    
    task_params$dual_factor <- task_params$dual_factor[
      !sapply(task_params$dual_factor, function(pair) "var_type" %in% pair)
    ]
    
    filtered_single_count <- original_single_count - length(task_params$single_factor)
    filtered_dual_count <- original_dual_count - length(task_params$dual_factor)
    
    if (filtered_single_count > 0 || filtered_dual_count > 0) {
      log_message("INFO: Auto-filtered var_type-related analyses (by_all_types mode):")
      if (filtered_single_count > 0) {
        log_message(sprintf("  - Removed %d single_factor analysis (var_type)", filtered_single_count))
      }
      if (filtered_dual_count > 0) {
        log_message(sprintf("  - Removed %d dual_factor combinations containing var_type", filtered_dual_count))
      }
    }
  } else {
    log_message("INFO: Retaining all analysis factors (by_variant_type mode)")
  }
  
  log_message(sprintf("Active analysis configuration:"))
  log_message(sprintf("  - single_factor: %d factors (%s)", 
                     length(task_params$single_factor),
                     paste(task_params$single_factor, collapse = ", ")))
  log_message(sprintf("  - dual_factor: %d combinations", length(task_params$dual_factor)))
  
  # Store analysis_mode in task_params for downstream use
  task_params$analysis_mode <- analysis_mode
  
  # Validate input data
  if (nrow(normalized_data) == 0) {
    log_message("ERROR: normalized_data is empty", level = "error")
    stop("M02 analysis cannot proceed with empty normalized_data")
  }
  
  # Apply data_scoping filters (PROJECT_CHARTER.md Principle 3)
  log_message("=== APPLYING DATA_SCOPING FILTERS ===")
  scoped_data <- apply_data_scoping(normalized_data, task_params)
  
  if (nrow(scoped_data) == 0) {
    log_message("ERROR: data_scoping resulted in empty dataset for M02", level = "error")
    stop("M02 analysis cannot proceed with empty data after data_scoping")
  }
  
  retention_rate <- ifelse(nrow(normalized_data) > 0, round(nrow(scoped_data)/nrow(normalized_data)*100, 2), 0)
  log_message(sprintf("M02 data scoping: %d -> %d rows (%.2f%% retained)", 
                      nrow(normalized_data), nrow(scoped_data), retention_rate))
  
  # Re-aggregate summary data when running in by_all_types mode
  # Architecture principle: M02 module creates its own summary data from filtered detailed data
  # This ensures upstream P03 normalization remains independent of downstream visualization needs
  if (analysis_mode == "by_all_types") {
    log_message("=== M02: RE-AGGREGATING SUMMARY DATA (by_all_types mode) ===")
    log_message(sprintf("Input: %d rows of detailed data (after data_scoping)", nrow(scoped_data)))
    
    # Re-aggregate: sum variant_count across all var_types for each region
    # Group by region_name rather than gene so IGS identities are preserved
    m02_summary_data <- scoped_data %>%
      dplyr::group_by(species, region_type, region_name, region_length, genome_region) %>%
      dplyr::summarise(
        total_variant_count = sum(variant_count, na.rm = TRUE),
        gene = dplyr::first(gene),  # Preserve gene column (CDS/intron have values, IGS is NA)
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        var_type = "all_types",
        variant_count = total_variant_count,
        frequency_per_kb = (total_variant_count / region_length) * 1000
      ) %>%
      dplyr::select(species, var_type, region_type, genome_region, region_name, gene,
                   variant_count, region_length, frequency_per_kb) %>%
      dplyr::filter(variant_count > 0)  # Only keep regions with variants
    
    log_message(sprintf("Output: %d rows of aggregated summary data", nrow(m02_summary_data)))
    log_message(sprintf("  Aggregation: %d detailed rows -> %d summary rows (%.1f%% reduction)", 
                       nrow(scoped_data), nrow(m02_summary_data), 
                       (1 - nrow(m02_summary_data)/nrow(scoped_data)) * 100))
    
    # Use the aggregated summary data for subsequent analysis
    scoped_data <- m02_summary_data
  } else {
    log_message("=== M02: USING DETAILED DATA (by_variant_type mode) ===")
    log_message(sprintf("Data: %d rows of detailed data (by variant type)", nrow(scoped_data)))
  }
  
  # Filter out NA frequencies and zero frequencies for comparative analysis
  cleaned_data <- scoped_data %>% 
    dplyr::filter(!is.na(frequency_per_kb) & frequency_per_kb > 0)
  
  log_message(sprintf("Data prepared: %d -> %d valid observations", 
                     nrow(normalized_data), nrow(cleaned_data)))
  
  # Confirm genome_region is present
  if (!"genome_region" %in% names(cleaned_data)) {
    log_message("ERROR: genome_region column missing after data scoping", level = "error")
    stop("genome_region field is required for M02 analysis. Please ensure preprocessing includes genome_region assignment.")
  }
  
  log_message(sprintf("genome_region distribution: %s", 
                     paste(names(table(cleaned_data$genome_region)), "=", table(cleaned_data$genome_region), collapse = ", ")))
  
  # Apply species ordering and create output directory
  species_ordering_result <- load_and_apply_species_ordering(
    data = cleaned_data,
    session_id = config$session_info$session_id,
    species_col = "species"
  )
  
  cleaned_data <- species_ordering_result$data
  
  # Create output directory
  data_output_dir <- file.path(output_dir, "M02_comparative", task_name)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Validate cleaned data
  if (nrow(cleaned_data) == 0) {
    log_message("ERROR: No valid data after cleaning", level = "error")
    stop("Cannot proceed with M02 analysis - no valid comparative data")
  }
  
  # Export an M02-specific data file for user inspection
  # Only export the actual data used for plotting (avoid redundancy)
  log_message("=== M02: EXPORTING DATA FILES ===")
  
  raw_data_file <- file.path(data_output_dir, "M02_raw_data_for_plotting.csv")
  tryCatch({
    utils::write.csv(cleaned_data, raw_data_file, row.names = FALSE)
    log_message(sprintf("Exported raw plotting data: %s (%d rows, %d columns)", 
                       basename(raw_data_file), nrow(cleaned_data), ncol(cleaned_data)))
    log_message(sprintf("  Data mode: %s", 
                       if (analysis_mode == "by_all_types") "Summary (aggregated)" else "Detailed (by variant type)"))
  }, error = function(e) {
    log_message(sprintf("WARNING: Failed to export raw plotting data: %s", e$message), level = "warning")
  })
  
  # Load group information using the robust, validated helper
  session_paths <- get_session_paths(config$session_info$session_id)
  final_plot_data <- merge_with_group_info(
    main_data = cleaned_data,
    session_paths = session_paths
  )
  
  # Execute universal analysis tasks
  all_plots <- list()
  all_results <- list()
  
  # Initialize dual-table collectors for "Collect-Enhance-Merge-Save" workflow
  single_factor_summaries <- list()
  single_factor_species_details <- list()
  dual_factor_summaries <- list()
  dual_factor_species_details <- list()
  
  # Single Factor Analysis Dispatch
  single_factor_list <- task_params$single_factor
  if (!is.null(single_factor_list) && length(single_factor_list) > 0) {
    log_message(sprintf("Analyzing %d single factors", length(single_factor_list)))
    
    for (factor_name in single_factor_list) {
      # Validate factor exists before calling analysis function
      if (!factor_name %in% names(final_plot_data)) {
        log_message(sprintf("Skipping factor '%s': not found in the data after merging group info.", 
                           factor_name), level = "warning")
        next # Skip to the next factor in the loop
      }
      
      single_result <- generate_universal_single_factor_comparison(
        plot_data = final_plot_data,
        factor_name = factor_name,
        config = config,
        task_params = task_params,
        output_dir = NULL
      )
      
      plot_key <- paste0("single_", factor_name)
      all_plots[[plot_key]] <- single_result$plot
      all_results[[plot_key]] <- single_result
      
      # Collect dual-table outputs for CSV export
      tryCatch({
        # Collect macro overview table (summary_table_gt)
        if (!is.null(single_result$summary_table_gt)) {
          if (requireNamespace("gt", quietly = TRUE)) {
            summary_df <- single_result$summary_table_gt$`_data` %>%
              dplyr::mutate(analysis_factor = factor_name, .before = 1)
            # Standardize column name: factor-specific column -> "level"
            names(summary_df)[names(summary_df) == factor_name] <- "level"
            single_factor_summaries[[factor_name]] <- summary_df
            log_message(sprintf("Collected summary table for single factor '%s' (%d rows)", factor_name, nrow(summary_df)))
          }
        }
        
        # Collect species detail table (species_summary_table_gt)
        if (!is.null(single_result$species_summary_table_gt)) {
          if (requireNamespace("gt", quietly = TRUE)) {
            species_df <- single_result$species_summary_table_gt$`_data` %>%
              dplyr::mutate(analysis_factor = factor_name, .before = 1)
            # Standardize column name: factor-specific column -> "level"
            names(species_df)[names(species_df) == factor_name] <- "level"
            single_factor_species_details[[factor_name]] <- species_df
            log_message(sprintf("Collected species detail table for single factor '%s' (%d rows)", factor_name, nrow(species_df)))
          }
        }
      }, error = function(e) {
        log_message(sprintf("Failed to collect tables for single factor '%s': %s", factor_name, e$message), level = "warning")
      })
    }
  } else {
    log_message("No single factor analyses configured - skipping")
  }
  
  # Dual Factor Analysis Dispatch
  dual_factor_pairs <- task_params$dual_factor
  if (!is.null(dual_factor_pairs) && length(dual_factor_pairs) > 0) {
    log_message(sprintf("Analyzing %d dual factor combinations", length(dual_factor_pairs)))
    
    for (factor_pair in dual_factor_pairs) {
      if (length(factor_pair) >= 2) {
        factor1 <- factor_pair[1]
        factor2 <- factor_pair[2]
        
        # Validate both factors exist in data
        if (factor1 %in% names(final_plot_data) && factor2 %in% names(final_plot_data)) {
          dual_result <- generate_universal_dual_factor_comparison(
            plot_data = final_plot_data,
            factor1 = factor1,
            factor2 = factor2,
            config = config,
            task_params = task_params,
            output_dir = NULL
          )
          
          plot_key <- paste0("dual_", factor1, "_", factor2)
          all_plots[[plot_key]] <- dual_result$plot
          all_results[[plot_key]] <- dual_result
          
          # Collect dual-table outputs for CSV export
          tryCatch({
            factor_combination <- paste(factor1, factor2, sep = "_x_")
            
            # Collect macro overview table (summary_table_gt)
            if (!is.null(dual_result$summary_table_gt)) {
              if (requireNamespace("gt", quietly = TRUE)) {
                summary_df <- dual_result$summary_table_gt$`_data` %>%
                  dplyr::mutate(analysis_factor = factor_combination, .before = 1)
                # Standardize column names: factor-specific columns -> "level1" and "level2"
                names(summary_df)[names(summary_df) == factor1] <- "level1"
                names(summary_df)[names(summary_df) == factor2] <- "level2"
                dual_factor_summaries[[factor_combination]] <- summary_df
                log_message(sprintf("Collected summary table for dual factor '%s' (%d rows)", factor_combination, nrow(summary_df)))
              }
            }
            
            # Collect species detail table (species_summary_table_gt)
            if (!is.null(dual_result$species_summary_table_gt)) {
              if (requireNamespace("gt", quietly = TRUE)) {
                species_df <- dual_result$species_summary_table_gt$`_data` %>%
                  dplyr::mutate(analysis_factor = factor_combination, .before = 1)
                # Standardize column names: factor-specific columns -> "level1" and "level2"
                names(species_df)[names(species_df) == factor1] <- "level1"
                names(species_df)[names(species_df) == factor2] <- "level2"
                dual_factor_species_details[[factor_combination]] <- species_df
                log_message(sprintf("Collected species detail table for dual factor '%s' (%d rows)", factor_combination, nrow(species_df)))
              }
            }
          }, error = function(e) {
            log_message(sprintf("Failed to collect tables for dual factor '%s x %s': %s", factor1, factor2, e$message), level = "warning")
          })
        } else {
          missing_factors <- setdiff(c(factor1, factor2), names(final_plot_data))
          log_message(sprintf("Factors '%s' not found", paste(missing_factors, collapse = "', '")), level = "warning")
        }
      }
    }
  }
  
  # ============================================================================
  # MERGE AND SAVE DUAL-TABLE OUTPUTS (Collect-Enhance-Merge-Save workflow)
  # ============================================================================
  log_message("=== MERGING AND SAVING DUAL-TABLE OUTPUTS ===")
  
  # Merge and save single factor summary tables
  if (length(single_factor_summaries) > 0) {
    tryCatch({
      final_single_summary_df <- dplyr::bind_rows(single_factor_summaries)
      single_summary_path <- file.path(data_output_dir, "M02_single_factor_summary.csv")
      success <- safe_write_file(final_single_summary_df, single_summary_path, row.names = FALSE)
      if (success) {
        log_message(sprintf("[SUCCESS] Saved single factor summary table: %s (%d rows)", basename(single_summary_path), nrow(final_single_summary_df)))
      }
    }, error = function(e) {
      log_message(sprintf("Failed to save single factor summary table: %s", e$message), level = "error")
    })
  } else {
    log_message("No single factor summary tables to save")
  }
  
  # Merge and save single factor species detail tables
  if (length(single_factor_species_details) > 0) {
    tryCatch({
      final_single_species_df <- dplyr::bind_rows(single_factor_species_details)
      single_species_path <- file.path(data_output_dir, "M02_single_factor_species_details.csv")
      success <- safe_write_file(final_single_species_df, single_species_path, row.names = FALSE)
      if (success) {
        log_message(sprintf("[SUCCESS] Saved single factor species detail table: %s (%d rows)", basename(single_species_path), nrow(final_single_species_df)))
      }
    }, error = function(e) {
      log_message(sprintf("Failed to save single factor species detail table: %s", e$message), level = "error")
    })
  } else {
    log_message("No single factor species detail tables to save")
  }
  
  # Merge and save dual factor summary tables
  if (length(dual_factor_summaries) > 0) {
    tryCatch({
      final_dual_summary_df <- dplyr::bind_rows(dual_factor_summaries)
      dual_summary_path <- file.path(data_output_dir, "M02_dual_factor_summary.csv")
      success <- safe_write_file(final_dual_summary_df, dual_summary_path, row.names = FALSE)
      if (success) {
        log_message(sprintf("[SUCCESS] Saved dual factor summary table: %s (%d rows)", basename(dual_summary_path), nrow(final_dual_summary_df)))
      }
    }, error = function(e) {
      log_message(sprintf("Failed to save dual factor summary table: %s", e$message), level = "error")
    })
  } else {
    log_message("No dual factor summary tables to save")
  }
  
  # Merge and save dual factor species detail tables
  if (length(dual_factor_species_details) > 0) {
    tryCatch({
      final_dual_species_df <- dplyr::bind_rows(dual_factor_species_details)
      dual_species_path <- file.path(data_output_dir, "M02_dual_factor_species_details.csv")
      success <- safe_write_file(final_dual_species_df, dual_species_path, row.names = FALSE)
      if (success) {
        log_message(sprintf("[SUCCESS] Saved dual factor species detail table: %s (%d rows)", basename(dual_species_path), nrow(final_dual_species_df)))
      }
    }, error = function(e) {
      log_message(sprintf("Failed to save dual factor species detail table: %s", e$message), level = "error")
    })
  } else {
    log_message("No dual factor species detail tables to save")
  }
  
  # FALLBACK: If no universal analysis configured, provide default analysis
  if (length(all_plots) == 0) {
    log_message("=== FALLBACK: No universal analysis configured - using default comprehensive analysis ===")
    log_message("WARNING: M02 task configuration is missing 'single_factor' and 'dual_factor' parameters")
    log_message("Providing basic region_type and var_type analysis as fallback")
    
    # Default single factor analysis
    default_single_factors <- c("region_type", "var_type")
    default_dual_factors <- list(c("region_type", "var_type"))
    
    # Check for group columns and add them to analysis if available
    available_group_cols <- intersect(c("Phylogeny", "Life_form", "Status"), names(final_plot_data))
    if (length(available_group_cols) > 0) {
      log_message(sprintf("Found group columns: %s", paste(available_group_cols, collapse = ", ")))
      default_single_factors <- c(default_single_factors, available_group_cols)
      
      # Add group-based dual factor combinations
      for (group_col in available_group_cols) {
        default_dual_factors <- append(default_dual_factors, list(c(group_col, "region_type")))
        default_dual_factors <- append(default_dual_factors, list(c(group_col, "var_type")))
      }
    }
    
    log_message(sprintf("Fallback analysis will include %d single factors and %d dual factor combinations", 
                       length(default_single_factors), length(default_dual_factors)))
    
    # Execute fallback analysis using the same logic as above
    fallback_plots <- list()
    fallback_results <- list()
    
    # Single factor analysis
    for (factor_name in default_single_factors) {
      if (factor_name %in% names(final_plot_data)) {
        tryCatch({
          single_result <- generate_universal_single_factor_comparison(
            plot_data = final_plot_data,
            factor_name = factor_name,
            config = config,
            task_params = task_params,
            output_dir = NULL
          )
          
          plot_key <- paste0("single_", factor_name)
          fallback_plots[[plot_key]] <- single_result$plot
          fallback_results[[plot_key]] <- single_result
          
          # Collect dual-table outputs for CSV export (fallback path)
          tryCatch({
            # Collect macro overview table (summary_table_gt)
            if (!is.null(single_result$summary_table_gt)) {
              if (requireNamespace("gt", quietly = TRUE)) {
                summary_df <- single_result$summary_table_gt$`_data` %>%
                  dplyr::mutate(analysis_factor = factor_name, .before = 1)
                # Standardize column name: factor-specific column -> "level"
                names(summary_df)[names(summary_df) == factor_name] <- "level"
                single_factor_summaries[[factor_name]] <- summary_df
                log_message(sprintf("Collected fallback summary table for single factor '%s' (%d rows)", factor_name, nrow(summary_df)))
              }
            }
            
            # Collect species detail table (species_summary_table_gt)
            if (!is.null(single_result$species_summary_table_gt)) {
              if (requireNamespace("gt", quietly = TRUE)) {
                species_df <- single_result$species_summary_table_gt$`_data` %>%
                  dplyr::mutate(analysis_factor = factor_name, .before = 1)
                # Standardize column name: factor-specific column -> "level"
                names(species_df)[names(species_df) == factor_name] <- "level"
                single_factor_species_details[[factor_name]] <- species_df
                log_message(sprintf("Collected fallback species detail table for single factor '%s' (%d rows)", factor_name, nrow(species_df)))
              }
            }
          }, error = function(e) {
            log_message(sprintf("Failed to collect fallback tables for single factor '%s': %s", factor_name, e$message), level = "warning")
          })
        }, error = function(e) {
          log_message(sprintf("Fallback single factor analysis failed for '%s': %s", factor_name, e$message), level = "warning")
        })
      }
    }
    
    # Dual factor analysis
    for (factor_pair in default_dual_factors) {
      if (length(factor_pair) >= 2) {
        factor1 <- factor_pair[1]
        factor2 <- factor_pair[2]
        
        if (factor1 %in% names(final_plot_data) && factor2 %in% names(final_plot_data)) {
          tryCatch({
            dual_result <- generate_universal_dual_factor_comparison(
              plot_data = final_plot_data,
              factor1 = factor1,
              factor2 = factor2,
              config = config,
              task_params = task_params,
              output_dir = NULL
            )
            
            plot_key <- paste0("dual_", factor1, "_", factor2)
            fallback_plots[[plot_key]] <- dual_result$plot
            fallback_results[[plot_key]] <- dual_result
            
            # Collect dual-table outputs for CSV export (fallback path)
            tryCatch({
              factor_combination <- paste(factor1, factor2, sep = "_x_")
              
              # Collect macro overview table (summary_table_gt)
              if (!is.null(dual_result$summary_table_gt)) {
                if (requireNamespace("gt", quietly = TRUE)) {
                  summary_df <- dual_result$summary_table_gt$`_data` %>%
                    dplyr::mutate(analysis_factor = factor_combination, .before = 1)
                  # Standardize column names: factor-specific columns -> "level1" and "level2"
                  names(summary_df)[names(summary_df) == factor1] <- "level1"
                  names(summary_df)[names(summary_df) == factor2] <- "level2"
                  dual_factor_summaries[[factor_combination]] <- summary_df
                  log_message(sprintf("Collected fallback summary table for dual factor '%s' (%d rows)", factor_combination, nrow(summary_df)))
                }
              }
              
              # Collect species detail table (species_summary_table_gt)
              if (!is.null(dual_result$species_summary_table_gt)) {
                if (requireNamespace("gt", quietly = TRUE)) {
                  species_df <- dual_result$species_summary_table_gt$`_data` %>%
                    dplyr::mutate(analysis_factor = factor_combination, .before = 1)
                  # Standardize column names: factor-specific columns -> "level1" and "level2"
                  names(species_df)[names(species_df) == factor1] <- "level1"
                  names(species_df)[names(species_df) == factor2] <- "level2"
                  dual_factor_species_details[[factor_combination]] <- species_df
                  log_message(sprintf("Collected fallback species detail table for dual factor '%s' (%d rows)", factor_combination, nrow(species_df)))
                }
              }
            }, error = function(e) {
              log_message(sprintf("Failed to collect fallback tables for dual factor '%s x %s': %s", factor1, factor2, e$message), level = "warning")
            })
          }, error = function(e) {
            log_message(sprintf("Fallback dual factor analysis failed for '%s x %s': %s", factor1, factor2, e$message), level = "warning")
          })
        }
      }
    }
    
    # Use fallback results
    all_plots <- fallback_plots
    all_results <- fallback_results
    
    # Initialize composite_plots for fallback branch (fix for "object not found" error)
    composite_plots <- list()
    composite_results <- list()
    
    log_message(sprintf("Fallback analysis completed: %d plots generated", length(all_plots)))
    
    # ============================================================================
    # MERGE AND SAVE DUAL-TABLE OUTPUTS (Fallback path)
    # ============================================================================
    log_message("=== MERGING AND SAVING DUAL-TABLE OUTPUTS (FALLBACK) ===")
    
    # Merge and save single factor summary tables
    if (length(single_factor_summaries) > 0) {
      tryCatch({
        final_single_summary_df <- dplyr::bind_rows(single_factor_summaries)
        single_summary_path <- file.path(data_output_dir, "M02_single_factor_summary.csv")
        success <- safe_write_file(final_single_summary_df, single_summary_path, row.names = FALSE)
        if (success) {
          log_message(sprintf("[SUCCESS] Saved fallback single factor summary table: %s (%d rows)", basename(single_summary_path), nrow(final_single_summary_df)))
        }
      }, error = function(e) {
        log_message(sprintf("Failed to save fallback single factor summary table: %s", e$message), level = "error")
      })
    } else {
      log_message("No fallback single factor summary tables to save")
    }
    
    # Merge and save single factor species detail tables
    if (length(single_factor_species_details) > 0) {
      tryCatch({
        final_single_species_df <- dplyr::bind_rows(single_factor_species_details)
        single_species_path <- file.path(data_output_dir, "M02_single_factor_species_details.csv")
        success <- safe_write_file(final_single_species_df, single_species_path, row.names = FALSE)
        if (success) {
          log_message(sprintf("[SUCCESS] Saved fallback single factor species detail table: %s (%d rows)", basename(single_species_path), nrow(final_single_species_df)))
        }
      }, error = function(e) {
        log_message(sprintf("Failed to save fallback single factor species detail table: %s", e$message), level = "error")
      })
    } else {
      log_message("No fallback single factor species detail tables to save")
    }
    
    # Merge and save dual factor summary tables
    if (length(dual_factor_summaries) > 0) {
      tryCatch({
        final_dual_summary_df <- dplyr::bind_rows(dual_factor_summaries)
        dual_summary_path <- file.path(data_output_dir, "M02_dual_factor_summary.csv")
        success <- safe_write_file(final_dual_summary_df, dual_summary_path, row.names = FALSE)
        if (success) {
          log_message(sprintf("[SUCCESS] Saved fallback dual factor summary table: %s (%d rows)", basename(dual_summary_path), nrow(final_dual_summary_df)))
        }
      }, error = function(e) {
        log_message(sprintf("Failed to save fallback dual factor summary table: %s", e$message), level = "error")
      })
    } else {
      log_message("No fallback dual factor summary tables to save")
    }
    
    # Merge and save dual factor species detail tables
    if (length(dual_factor_species_details) > 0) {
      tryCatch({
        final_dual_species_df <- dplyr::bind_rows(dual_factor_species_details)
        dual_species_path <- file.path(data_output_dir, "M02_dual_factor_species_details.csv")
        success <- safe_write_file(final_dual_species_df, dual_species_path, row.names = FALSE)
        if (success) {
          log_message(sprintf("[SUCCESS] Saved fallback dual factor species detail table: %s (%d rows)", basename(dual_species_path), nrow(final_dual_species_df)))
        }
      }, error = function(e) {
        log_message(sprintf("Failed to save fallback dual factor species detail table: %s", e$message), level = "error")
      })
    } else {
      log_message("No fallback dual factor species detail tables to save")
    }
    
  } else {
    single_factor_count <- length(single_factor_list %||% c())
    dual_factor_count <- length(dual_factor_pairs %||% c())
    
    # Generate composite plots
    composite_plots <- list()
    composite_results <- list()
    
    # Generate single factor composite plot if applicable
    single_factor_results <- all_results[grepl("^single_", names(all_results))]
    if (length(single_factor_results) > 0) {
      single_composite_result <- generate_single_factor_composite_plot(
        single_factor_results = single_factor_results,
        config = config,
        task_params = task_params
      )
      
      composite_plots[["single_factor_composite"]] <- single_composite_result$composite_plot
      composite_results[["single_factor_composite"]] <- single_composite_result
      
      # Save single factor composite plot
      single_composite_base <- file.path(data_output_dir, paste0(task_name, "_composite_single_factor"))
      tryCatch({
        saved_path <- save_plot(
          plot_object = single_composite_result$composite_plot,
          base_path = single_composite_base,
          config = config,
          width = 16,
          height = 12
        )
        log_message(sprintf("[SUCCESS] Single factor composite: composite_single_factor.png"))
      }, error = function(e) {
        log_message(sprintf("[ERROR] Failed composite plot: %s", e$message), level = "error")
      })
    }
    
    # Generate dual factor composite plot if applicable
    dual_factor_results <- all_results[grepl("^dual_", names(all_results))]
    if (length(dual_factor_results) > 0) {
      dual_composite_result <- generate_dual_factor_composite_plot(
        dual_factor_results = dual_factor_results,
        config = config,
        task_params = task_params
      )
      
      composite_plots[["dual_factor_composite"]] <- dual_composite_result$composite_plot
      composite_results[["dual_factor_composite"]] <- dual_composite_result
      
      # Save dual factor composite plot
      dual_composite_base <- file.path(data_output_dir, paste0(task_name, "_composite_dual_factor"))
      tryCatch({
        saved_path <- save_plot(
          plot_object = dual_composite_result$composite_plot,
          base_path = dual_composite_base,
          config = config,
          width = 16,
          height = 12
        )
        log_message(sprintf("[SUCCESS] Dual factor composite: composite_dual_factor.png"))
      }, error = function(e) {
        log_message(sprintf("[ERROR] Failed composite plot: %s", e$message), level = "error")
      })
    }
    
    # FLATTENED OUTPUT: Save individual plots with new naming convention
    single_factor_reports <- character(0)
    dual_factor_reports <- character(0)
    run_summary <- list(
      analyses = list(),
      session_id = config$session_info$session_id,
      timestamp = Sys.time(),
      total_analyses = length(all_results)
    )
    
    for (plot_name in names(all_results)) {
      plot_result <- all_results[[plot_name]]
      
      # STEP 1: Flattened file naming convention
      # Individual plots: M02_universal_analysis_plot_single_region_type.png, M02_universal_analysis_plot_dual_Phylogeny_var_type.png  
      plot_filename <- paste0(task_name, "_plot_", plot_name, ".png")
      plot_filepath <- file.path(data_output_dir, plot_filename)
      
      # Save the plot directly to main task directory
      if (!is.null(plot_result$plot)) {
        tryCatch({
          plot_base_path <- tools::file_path_sans_ext(plot_filepath)
          saved_path <- save_plot(
            plot_object = plot_result$plot,
            base_path = plot_base_path,
            config = config,
            width = 12,
            height = 8
          )
          # Simplified logging
        }, error = function(e) {
          log_message(sprintf("[ERROR] Failed to save %s: %s", plot_filename, e$message), level = "error")
        })
      }
      
      # STEP 2: Collect statistical reports for integration (implement in next step)
      if (!is.null(plot_result$statistical_report_text)) {
        if (grepl("^single_", plot_name)) {
          single_factor_reports <- c(single_factor_reports, plot_result$statistical_report_text)
        } else if (grepl("^dual_", plot_name)) {
          dual_factor_reports <- c(dual_factor_reports, plot_result$statistical_report_text)
        }
      }
      
      # STEP 3: Collect metadata for run summary  
      analysis_info <- list(
        name = plot_name,
        plot_file = plot_filename,
        factor_analyzed = if(!is.null(plot_result$metadata$factor_analyzed)) plot_result$metadata$factor_analyzed else plot_result$metadata$factors_analyzed,
        plot_type = plot_result$metadata$plot_type %||% "boxplot",
        total_observations = plot_result$metadata$total_observations %||% 0
      )
      run_summary$analyses[[plot_name]] <- analysis_info
    }
    
    # STEP 2: Write integrated statistical reports
    if (length(single_factor_reports) > 0) {
      single_report_text <- paste(single_factor_reports, collapse = "\n---\n")
      single_report_file <- file.path(data_output_dir, paste0(task_name, "_report_all_single_factor.txt"))
      writeLines(single_report_text, single_report_file)
      log_message("[SUCCESS] Single factor statistical report created")
    }
    
    if (length(dual_factor_reports) > 0) {
      dual_report_text <- paste(dual_factor_reports, collapse = "\n---\n")
      dual_report_file <- file.path(data_output_dir, paste0(task_name, "_report_all_dual_factor.txt"))
      writeLines(dual_report_text, dual_report_file)
      log_message("[SUCCESS] Dual factor statistical report created")
    }
    
    # STEP 3: Generate single run summary JSON
    run_summary$composite_plots <- list(
      single_factor_composite = if(exists("single_composite_file") && file.exists(single_composite_file)) basename(single_composite_file) else NULL,
      dual_factor_composite = if(exists("dual_composite_file") && file.exists(dual_composite_file)) basename(dual_composite_file) else NULL
    )
    run_summary$reports <- list(
      single_factor_report = if(length(single_factor_reports) > 0) "_report_all_single_factor.txt" else NULL,
      dual_factor_report = if(length(dual_factor_reports) > 0) "_report_all_dual_factor.txt" else NULL
    )
    run_summary$summary_counts <- list(
      single_factor_analyses = length(single_factor_results),
      dual_factor_analyses = length(dual_factor_results),
      composite_plots = length(composite_plots),
      total_plots = length(all_plots)
    )
    
    # Save as JSON
    run_summary_file <- file.path(data_output_dir, "_run_summary.json")
    tryCatch({
      if (requireNamespace("jsonlite", quietly = TRUE)) {
        jsonlite::write_json(run_summary, run_summary_file, pretty = TRUE, auto_unbox = TRUE)
        log_message("[SUCCESS] Run summary created")
      } else {
        # Fallback to basic text format if jsonlite not available
        run_summary_text <- capture.output(str(run_summary))
        writeLines(run_summary_text, gsub("\\.json$", ".txt", run_summary_file))
        log_message("[SUCCESS] Run summary created (text format)")
      }
    }, error = function(e) {
      log_message(sprintf("Failed to save run summary: %s", e$message), level = "warning")
    })
    
    # Create combined result structure
    result <- list(
      plots = all_plots,
      results = all_results,
      composite_plots = composite_plots,
      composite_results = composite_results,
      summary = list(
        total_analyses = length(all_plots),
        composite_plots_count = length(composite_plots),
        single_factor_count = length(single_factor_list %||% c()),
        dual_factor_count = length(dual_factor_pairs %||% c()),
        data_rows = nrow(final_plot_data),
        data_columns = ncol(final_plot_data)
      )
    )
  }
  
  # ============================================================================
  # GENERATE STATISTICAL DECISION FLOWCHART (TEXT VERSION)
  # ============================================================================
  tryCatch({
    # Generate text-based flowchart and save to file
    flowchart_text <- generate_statistics_flowchart_text()
    flowchart_file <- file.path(data_output_dir, "M02_workflow_flowchart.txt")
    
    writeLines(flowchart_text, flowchart_file)
    
    log_message(sprintf("[SUCCESS] Statistical decision flowchart (text) saved to: %s", flowchart_file))
  }, error = function(e) {
    log_message(sprintf("Warning: Failed to generate flowchart: %s", e$message), level = "warning")
  })
  
  # ================================================================
  # Phase 1 Pilot: Generate Manifest for Scientific Report System
  # ================================================================
  log_message("Generating output manifest for scientific report system")
  
  manifest_entries <- list()
  
  # Process composite plots (highest priority)
  if (!is.null(composite_plots) && length(composite_plots) > 0) {
    for (composite_name in names(composite_plots)) {
      
      if (composite_name == "single_factor_composite" && exists("single_composite_file") && file.exists(single_composite_file)) {
        title <- get_standard_title("M02", "composite_plot", "single_factor_composite")
        caption <- get_standard_caption("M02", "composite_plot")
        
        entry <- create_manifest_entry(
          path = single_composite_file,
          role = "composite_plot",
          title = title,
          caption = caption,
          module = "M02",
          file_type = "plot",
          priority = 1,
          parameters = task_params  # V7.1 Enhancement: Pass task parameters
        )
        manifest_entries[["single_factor_composite"]] <- entry
        log_message(sprintf("[SUCCESS] Created manifest entry for single factor composite plot"))
      }
      
      if (composite_name == "dual_factor_composite" && exists("dual_composite_file") && file.exists(dual_composite_file)) {
        title <- get_standard_title("M02", "composite_plot", "dual_factor_composite")
        caption <- get_standard_caption("M02", "composite_plot")
        
        entry <- create_manifest_entry(
          path = dual_composite_file,
          role = "composite_plot",
          title = title,
          caption = caption,
          module = "M02",
          file_type = "plot",
          priority = 2,
          parameters = task_params  # V7.1 Enhancement: Pass task parameters
        )
        manifest_entries[["dual_factor_composite"]] <- entry
        log_message(sprintf("[SUCCESS] Created manifest entry for dual factor composite plot"))
      }
    }
  }
  
  # Process individual analysis plots (medium priority)
  if (!is.null(all_results) && length(all_results) > 0) {
    for (plot_name in names(all_results)) {
      plot_result <- all_results[[plot_name]]
      
      # Reconstruct the file path from the saved naming convention
      plot_filename <- paste0(task_name, "_plot_", plot_name, ".png")
      plot_filepath <- file.path(data_output_dir, plot_filename)
      
      if (file.exists(plot_filepath)) {
        
        # Determine priority based on plot type
        priority <- if (grepl("^single_", plot_name)) 3 else 4
        
        # Generate academic title
        analysis_type <- if (grepl("^single_", plot_name)) "single factor" else "dual factor"
        title <- get_standard_title("M02", "individual_plot", analysis_type)
        caption <- get_standard_caption("M02", "individual_plot")
        
        entry <- create_manifest_entry(
          path = plot_filepath,
          role = "individual_plot",
          title = title,
          caption = caption,
          module = "M02",
          file_type = "plot",
          priority = priority,
          parameters = task_params  # V7.1 Enhancement: Pass task parameters
        )
        
        manifest_entries[[plot_name]] <- entry
        log_message(sprintf("[SUCCESS] Created manifest entry for %s analysis plot", analysis_type))
      }
    }
  }
  
  # Process statistical reports (low priority)
  report_files <- list()
  
  # Single factor report
  single_report_file <- file.path(data_output_dir, paste0(task_name, "_report_all_single_factor.txt"))
  if (file.exists(single_report_file)) {
    report_files[["single_factor_report"]] <- single_report_file
  }
  
  # Dual factor report
  dual_report_file <- file.path(data_output_dir, paste0(task_name, "_report_all_dual_factor.txt"))
  if (file.exists(dual_report_file)) {
    report_files[["dual_factor_report"]] <- dual_report_file
  }
  
  # Add statistical reports to manifest
  for (report_name in names(report_files)) {
    report_path <- report_files[[report_name]]
    
    title <- get_standard_title("M02", "statistical_report", report_name)
    caption <- get_standard_caption("M02", "statistical_report")
    
    entry <- create_manifest_entry(
      path = report_path,
      role = "statistical_report",
      title = title,
      caption = caption,
      module = "M02",
      file_type = "report",
      priority = 6,
      parameters = task_params  # V7.1 Enhancement: Pass task parameters
    )
    
    manifest_entries[[report_name]] <- entry
    log_message(sprintf("[SUCCESS] Created manifest entry for %s statistical report", gsub("_", " ", report_name)))
  }
  
  # Create manifest collection
  session_id <- config$session_info$session_id %||% "unknown"
  manifest <- create_manifest_collection(manifest_entries, "M02", session_id)
  
  # Add manifest to result
  if (exists("result") && !is.null(result)) {
    result$manifest <- manifest
  } else {
    # Create minimal result structure if it doesn't exist
    result <- list(
      manifest = manifest,
      summary = list(
        total_manifest_entries = length(manifest_entries)
      )
    )
  }
  
  log_message(sprintf("=== M02 COMPARATIVE ANALYSIS COMPLETED: %s ===", task_name))
  log_message(sprintf("Generated manifest with %d entries", length(manifest_entries)))
  
  return(result)
}
#' Generate M03 hotspot analysis plots
#' 
#' @title Generate M03 hotspot analysis plots
#' @description Comprehensive 8-step hotspot analysis pipeline including heatmaps, PCA, UpSet plots, and enrichment analysis
#' @param normalized_data Normalized frequency data from P03 stage
#' @param config Configuration list containing analysis parameters
#' @param output_dir Optional output directory for saving plots
#' @param task_name Optional task name for identification
#' @param task_params Task-specific parameters
#' @return List containing all generated plots and analysis results
#' @importFrom magrittr %>%
#' @export
generate_M03_hotspot_plots <- function(normalized_data, 
                                     config, 
                                     output_dir = NULL, 
                                     task_name = NULL, 
                                     task_params = NULL) {
  
  # Dependencies automatically loaded in R package context
  
  log_message("Starting M03 Enhanced 8-Step Hotspot Analysis Pipeline")
  
  tryCatch({
    # Comprehensive input validation
    if (is.null(normalized_data)) {
      stop("M03: normalized_data parameter is NULL")
    }
    
    if (!is.data.frame(normalized_data)) {
      stop("M03: normalized_data must be a data.frame")
    }
    
    if (nrow(normalized_data) == 0) {
      stop("M03: normalized_data is empty (0 rows)")
    }
    
    if (is.null(config)) {
      stop("M03: config parameter is NULL")
    }
    
    if (is.null(config$session_info$session_id)) {
      stop("M03: config$session_info$session_id is NULL")
    }
    
    # Use internal standard column names (boundary principle)
    species_col <- "species"
    gene_col <- "gene"
    freq_col <- "frequency_per_kb"
    region_col <- "region_type"
    var_type_col <- "var_type"
    
    # Validate required columns with detailed error reporting
    required_cols <- c(species_col, gene_col, freq_col, region_col, var_type_col)
    available_cols <- names(normalized_data)
    missing_cols <- setdiff(required_cols, available_cols)
    
    if (length(missing_cols) > 0) {
      log_message(sprintf("M03 ERROR: Missing required columns: %s", paste(missing_cols, collapse = ", ")), level = "error")
      log_message(sprintf("M03 ERROR: Available columns: %s", paste(available_cols, collapse = ", ")), level = "error")
      stop(sprintf("M03: Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    log_message(sprintf("M03: Input validation passed - %d rows, %d columns", nrow(normalized_data), ncol(normalized_data)))
    
    # Set up task parameters with defaults
    task_params <- task_params %||% list()
    
    # Extract M03-specific parameters with enhanced support for dual threshold modes
    frequency_percentile <- get_task_parameter(task_params, config, "frequency_percentile", 0.75, log_default = TRUE)
    min_species_ratio <- get_task_parameter(task_params, config, "min_species_ratio", 0.25, log_default = TRUE)
    min_species_number <- get_task_parameter(task_params, config, "min_species_number", NULL)
    threshold_mode <- get_task_parameter(task_params, config, "threshold_mode", "ratio")  # "ratio" or "absolute"
    region_filter <- get_task_parameter(task_params, config, "region_filter", "CDS")
    var_type_filter <- get_task_parameter(task_params, config, "var_type_filter", "snp")
    plot_types <- get_task_parameter(task_params, config, "plot_types", c("global_heatmap", "global_pca", "upset", "hotspot_correlation", "dual_enrichment", "statistical_test"))
    
    # Enhanced clustering parameters
    cluster_rows <- get_task_parameter(task_params, config, "cluster_rows", TRUE)
    cluster_cols <- get_task_parameter(task_params, config, "cluster_cols", TRUE)
    
    # Set up output directory - use output_dir provided by orchestrator 
    # ARCHITECTURE FIX: Orchestrator handles directory creation, we just use it
    if (is.null(output_dir)) {
      output_dir <- "plots"  # Simple fallback, caller should provide proper path
    }
    
    task_output_dir <- output_dir  # Directly use orchestrator-provided path
    
    log_message(sprintf("M03 Enhanced output directory: %s", task_output_dir))
    
    # Initialize results structure
    results <- list()
    plots_generated <- character()
    
    # =================================================================
    # PHASE 0: GLOBAL ANALYSIS - Understanding the Full Dataset
    # =================================================================
    
    # Step 0a: Generate Global Frequency Heatmap
    if ("global_heatmap" %in% plot_types) {
      log_message("M03 Step 0a: Generating global gene frequency heatmap")
      global_heatmap_result <- tryCatch({
        generate_global_frequency_heatmap(
          normalized_data = normalized_data,
          config = config,
          task_params = task_params
        )
      }, error = function(e) {
        log_message(sprintf("ERROR in global_heatmap generation: %s", e$message), level = "error")
        return(NULL)
      })
      
      if (!is.null(global_heatmap_result) && !is.null(global_heatmap_result$plot)) {
        plot_filename <- "M03_global_freq_heatmap.png"
        plot_path <- file.path(task_output_dir, plot_filename)
        
        tryCatch({
          plot_base_path <- tools::file_path_sans_ext(plot_path)
          saved_path <- save_plot(global_heatmap_result$plot, plot_base_path, config, width = 16, height = 12)
        }, error = function(e) {
          log_message(sprintf("ERROR saving global_heatmap plot: %s", e$message), level = "error")
        })
        
        results[["global_heatmap"]] <- list(
          plot = global_heatmap_result$plot,
          data = global_heatmap_result$data,
          file_path = plot_path,
          metadata = global_heatmap_result$metadata
        )
        plots_generated <- c(plots_generated, "global_heatmap")
        log_message(sprintf("M03: Saved global heatmap: %s", plot_path))
      }
    }
    
    # Note: category_correlation removed - not aligned with V3 "macro exploration + focus validation" logic
    
    # Step 0c: Generate Global PCA Analysis
    if ("global_pca" %in% plot_types) {
      log_message("M03 Step 0c: Generating global variation pattern PCA analysis")

      # V20 CONFIGURATION-DRIVEN LOGIC:
      # Read the list of variables to group by from the task parameters in config.yml.
      pca_grouping_vars <- task_params$pca_groups %||% character(0)

      pca_result <- generate_global_pca_plot(
        global_heatmap_result = global_heatmap_result,
        normalized_data = normalized_data,
        config = config,
        group_by_vars = pca_grouping_vars, # Pass the new parameter
        task_params = task_params
      )
      
      if (!is.null(pca_result)) {
        # Save individual PCA plots
        pca_files <- character()
        
        # Group-specific 1x3 composite plots with confidence ellipses (main outputs)
        # Dynamic detection of composite plots to avoid hardcoded naming issues
        composite_plot_names <- names(pca_result$plots)[grepl("^combined_", names(pca_result$plots))]
        
        if (length(composite_plot_names) > 0) {
          log_message(sprintf("M03: Found %d composite plots in PCA results: %s", 
                             length(composite_plot_names), 
                             paste(composite_plot_names, collapse = ", ")))
          
          for (plot_name in composite_plot_names) {
            if (!is.null(pca_result$plots[[plot_name]])) {
              # Extract group type from plot name (e.g., "combined_Phylogeny" -> "Phylogeny")
              group_type <- gsub("^combined_", "", plot_name)
              plot_filename <- sprintf("M03_global_pca_%s_composite.png", group_type)
              plot_base_path <- file.path(task_output_dir, sprintf("M03_global_pca_%s_composite", group_type))
              saved_path <- save_plot(pca_result$plots[[plot_name]], plot_base_path, config, width = 18, height = 6)
              plot_path <- saved_path %||% file.path(task_output_dir, plot_filename)
              pca_files <- c(pca_files, plot_path)
              log_message(sprintf("M03: [SUCCESS] Saved %s PCA composite plot", group_type))
            }
          }
        } else {
          log_message("M03: Warning: No composite plots found in PCA results")
        }
        
        # Comprehensive loadings composite plot
        if (!is.null(pca_result$plots$combined_loadings)) {
          plot_base_path <- file.path(task_output_dir, "M03_global_pca_loadings_composite")
          saved_path <- save_plot(pca_result$plots$combined_loadings, plot_base_path, config, width = 18, height = 6)
          plot_path <- saved_path %||% file.path(task_output_dir, "M03_global_pca_loadings_composite.png")
          pca_files <- c(pca_files, plot_path)
          log_message(sprintf("M03: Saved loadings composite plot: %s", plot_path))
        }
        
        # Individual component plots and scree plot
        for (plot_type in c("scree")) {
          if (!is.null(pca_result$plots[[plot_type]])) {
            plot_base_path <- file.path(task_output_dir, sprintf("M03_global_pca_%s", plot_type))
            saved_path <- save_plot(pca_result$plots[[plot_type]], plot_base_path, config, width = 12, height = 10)
            plot_path <- saved_path %||% file.path(task_output_dir, sprintf("M03_global_pca_%s.png", plot_type))
            pca_files <- c(pca_files, plot_path)
            log_message(sprintf("M03: Saved %s PCA plot: %s", plot_type, plot_path))
          }
        }
        
        # Save species ranking table
        if (!is.null(pca_result$species_ranking)) {
          ranking_filename <- "M03_global_pca_species_ranking.csv"
          ranking_path <- file.path(task_output_dir, ranking_filename)
          write.csv(pca_result$species_ranking, ranking_path, row.names = FALSE)
          log_message(sprintf("M03: Saved species PCA ranking table: %s", ranking_path))
        }
        
        results[["global_pca"]] <- list(
          plots = pca_result$plots,
          pca_data = pca_result$pca_data,
          species_ranking = pca_result$species_ranking,
          loadings_data = pca_result$loadings_data,
          file_paths = pca_files,
          metadata = pca_result$metadata
        )
        plots_generated <- c(plots_generated, "global_pca")
        log_message(sprintf("M03: Completed global PCA analysis with %.1f%% variance explained", 
                           pca_result$metadata$total_variance_captured * 100))
      }
    }
    
    # =================================================================
    # PHASE 1: HOTSPOT IDENTIFICATION - Enhanced Candidate Detection
    # =================================================================
    
    # Step 1: Get all ranked genes with percentile and threshold information
    log_message("M03 Step 1: Ranking all genes and calculating hotspot thresholds")
    all_genes_ranked <- identify_candidate_hotspots(
      normalized_data, 
      config = config,
      task_params = task_params
    )
    
    if (is.null(all_genes_ranked) || nrow(all_genes_ranked) == 0) {
      warning("M03: No genes available for ranking and hotspot identification")
      return(list(plots = list(), summary = list(error = "No genes available for ranking"), success = FALSE))
    }

    # Step 1b: Perform the definitive filtering to get the true candidate hotspots
    # CRITICAL FIX: This is the ONLY place where hotspot filtering should occur
    candidate_hotspots <- all_genes_ranked %>%
      dplyr::filter(is_potential_hotspot == TRUE) %>%
      dplyr::arrange(species, dplyr::desc(frequency_per_kb))
    
    if (is.null(candidate_hotspots) || nrow(candidate_hotspots) == 0) {
      warning("M03: No candidate hotspots identified after applying thresholds")
      return(list(plots = list(), summary = list(error = "No candidate hotspots found"), success = FALSE))
    }
    
    log_message(sprintf("M03: Identified %d candidate hotspot entries from %d unique genes across %d species",
                       nrow(candidate_hotspots), 
                       length(unique(candidate_hotspots$gene)),
                       length(unique(candidate_hotspots$species))))
    log_message(sprintf("M03: Unique genes identified as hotspots: %d", 
                       length(unique(candidate_hotspots$gene))))
    
    # =================================================================
    # PHASE 2: SHARING PATTERN ANALYSIS - Enhanced Visualization
    # =================================================================
    
    # Step 2a: Generate Enhanced UpSet Plot with proper species ordering
    if ("upset" %in% plot_types) {
      log_message("M03 Step 2a: Generating enhanced UpSet plot with phylogenetic ordering")
      plot_filename <- "M03_candidate_upset.pdf"
      plot_path <- file.path(task_output_dir, plot_filename)
      upset_result <- tryCatch({
        generate_hotspot_upset_plot(
          candidate_hotspots,
          config = config,
          task_params = task_params,
          output_path = plot_path
        )
      }, error = function(e) {
        log_message(sprintf("ERROR in upset plot generation: %s", e$message), level = "error")
        return(NULL)
      })
      
      if (!is.null(upset_result) && !is.null(upset_result$plot)) {
        if (file.exists(plot_path)) {
          log_message(sprintf("[SUCCESS] Saved UpSet plot using cross-platform PDF device: %s", basename(plot_path)))
        }
        
        # Save UpSet intersection data to CSV for user reference
        if (!is.null(upset_result$intersection_data)) {
          upset_csv_path <- file.path(task_output_dir, "M03_upset_intersections.csv")
          write.csv(upset_result$intersection_data, upset_csv_path, row.names = FALSE)
          log_message(sprintf("Saved UpSet intersection data: %s", basename(upset_csv_path)))
        }
        
        results[["upset"]] <- list(
          plot = upset_result$plot,
          data = upset_result$data,
          file_path = plot_path,
          metadata = upset_result$metadata
        )
        plots_generated <- c(plots_generated, "upset")
        log_message(sprintf("M03: Saved enhanced UpSet plot: %s", plot_path))
      }
    }
    
    # Step 2b: Hotspot Heatmap functionality removed per user request
    # The clustering heatmap functionality has been disabled due to fallback mechanism issues
    if ("heatmap" %in% plot_types) {
      log_message("M03 Step 2b: Hotspot clustering heatmap has been disabled per user request - skipping", level = "warning")
    }
    
    # Step 2c: Generate Hotspot Correlation Analysis
    if ("hotspot_correlation" %in% plot_types) {
      log_message("M03 Step 2c: Generating core hotspot correlation analysis")
      hotspot_correlation_result <- tryCatch({
        generate_hotspot_correlation_heatmap(
          candidate_hotspots = candidate_hotspots,
          normalized_data = normalized_data,
          config = config,
          task_params = task_params
        )
      }, error = function(e) {
        log_message(sprintf("ERROR in hotspot correlation generation: %s", e$message), level = "error")
        return(NULL)
      })
      
      if (!is.null(hotspot_correlation_result)) {
        # Save correlation heatmaps
        correlation_files <- character()
        
        for (plot_type in names(hotspot_correlation_result$plots)) {
          plot_filename <- sprintf("M03_core_shared_correlation_%s.png", plot_type)
          plot_path <- file.path(task_output_dir, plot_filename)
          
          # Save correlation plot
          tryCatch({
            plot_base_path <- tools::file_path_sans_ext(plot_path)
            saved_path <- save_plot(hotspot_correlation_result$plots[[plot_type]], plot_base_path, config, width = 12, height = 10)
            plot_path <- saved_path %||% plot_path
            correlation_files <- c(correlation_files, plot_path)
            log_message(sprintf("M03: Saved %s correlation heatmap: %s", plot_type, plot_path))
          }, error = function(e) {
            log_message(sprintf("Failed to save %s correlation heatmap: %s", plot_type, e$message), level = "error")
          })
        }
        
        # Save detailed correlation statistics with p-values
        correlation_csv_files <- character()
        
        # Save binary correlation matrix with significance
        if (!is.null(hotspot_correlation_result$correlation_data$binary_correlation)) {
          binary_corr_matrix <- hotspot_correlation_result$correlation_data$binary_correlation
          binary_corr_data <- expand.grid(
            Species1 = rownames(binary_corr_matrix),
            Species2 = colnames(binary_corr_matrix),
            stringsAsFactors = FALSE
          )
          binary_corr_data$Correlation <- as.vector(binary_corr_matrix)
          binary_corr_data$R_squared <- binary_corr_data$Correlation^2
          
          # Add p-values if available from correlation result
          if (!is.null(hotspot_correlation_result$correlation_data$binary_pvalues)) {
            binary_corr_data$P_value <- as.vector(hotspot_correlation_result$correlation_data$binary_pvalues)
          } else {
            binary_corr_data$P_value <- NA
          }
          
          binary_corr_data$Analysis_Type <- "Binary_Presence_Absence"
          binary_corr_data$Significant <- ifelse(is.na(binary_corr_data$P_value), FALSE, binary_corr_data$P_value < 0.05)
          
          binary_csv_filename <- "M03_core_shared_correlation_details.csv"
          binary_csv_path <- file.path(task_output_dir, binary_csv_filename)
          write.csv(binary_corr_data, binary_csv_path, row.names = FALSE)
          correlation_csv_files <- c(correlation_csv_files, binary_csv_path)
          log_message(sprintf("M03: Saved binary correlation details with %d significant pairs: %s", 
                             sum(binary_corr_data$Significant, na.rm = TRUE), binary_csv_path))
        }
        
        # REMOVED: Frequency correlation CSV export
        # Deleted to avoid misleading results with sparse hotspot gene data
        log_message("M03: Frequency correlation CSV export skipped (analysis mode: binary-only)")
        
        # Save core genes summary
        if (!is.null(hotspot_correlation_result$correlation_data$gene_species_count)) {
          core_genes_csv_filename <- "M03_core_shared_genes_summary.csv"
          core_genes_csv_path <- file.path(task_output_dir, core_genes_csv_filename)
          write.csv(hotspot_correlation_result$correlation_data$gene_species_count, core_genes_csv_path, row.names = FALSE)
          correlation_csv_files <- c(correlation_csv_files, core_genes_csv_path)
          log_message(sprintf("M03: Saved core genes summary: %s", core_genes_csv_path))
        }
        
        results[["hotspot_correlation"]] <- list(
          plots = hotspot_correlation_result$plots,
          correlation_data = hotspot_correlation_result$correlation_data,
          file_paths = correlation_files,
          csv_files = correlation_csv_files,
          metadata = hotspot_correlation_result$metadata
        )
        plots_generated <- c(plots_generated, "hotspot_correlation")
        log_message(sprintf("M03: Completed hotspot correlation analysis for %d core genes", 
                           hotspot_correlation_result$metadata$n_core_genes))
      }
    }
    
    # =================================================================
    # PHASE 3: DUAL THRESHOLD ANALYSIS - Core vs Private Genes
    # =================================================================
    
    # Step 3: Analyze Core Shared Hotspots with dual threshold support
    log_message("M03 Step 3: Analyzing core shared vs private hotspots with dual threshold logic")
    
    # Determine effective threshold based on mode
    effective_threshold <- if (threshold_mode == "absolute" && !is.null(min_species_number)) {
      min_species_number
    } else {
      min_species_ratio
    }
    
    core_hotspots <- analyze_core_shared_hotspots(
      candidate_hotspots,
      config = config,
      task_params = task_params
    )
    
    # =================================================================
    # PHASE 4: DUAL ENRICHMENT ANALYSIS - Separate Core/Private Analysis
    # =================================================================
    
    # Step 4: Generate Dual Enrichment Analysis (Core vs Private)
    if ("dual_enrichment" %in% plot_types && !is.null(core_hotspots) && !is.null(core_hotspots$genes) && nrow(core_hotspots$genes) > 0) {
      log_message("M03 Step 4: Generating dual functional enrichment analysis")
      # Load gene function data
      gene_function_data <- tryCatch({
        load_gene_function_mapping(config = config)
      }, error = function(e) {
        log_message(sprintf("Failed to load gene function data: %s", e$message), level = "warning")
        return(NULL)
      })
      
      dual_enrichment_result <- generate_dual_enrichment_analysis(
        hotspot_analysis_results = core_hotspots,
        gene_function_data = gene_function_data,
        background_genes = unique(normalized_data$gene),
        config = config,
        task_params = task_params
      )
      
      if (!is.null(dual_enrichment_result)) {
        # Generate visualization plots
        dual_enrichment_plots <- generate_dual_enrichment_plots(
          dual_enrichment_results = dual_enrichment_result,
          config = config,
          task_params = task_params
        )
        
        if (!is.null(dual_enrichment_plots)) {
          # Save dual enrichment plots
          if (!is.null(dual_enrichment_plots$plot)) {
            plot_filename <- "M03_private_dual_enrichment.png"
            plot_base_path <- file.path(task_output_dir, tools::file_path_sans_ext(plot_filename))
            
            saved_path <- save_plot(dual_enrichment_plots$plot, plot_base_path, config, width = 14, height = 10)
            
            results[["dual_enrichment"]] <- list(
              plot = dual_enrichment_plots$plot,
              data = dual_enrichment_result,
              file_path = plot_path,
              metadata = dual_enrichment_plots$metadata
            )
            plots_generated <- c(plots_generated, "dual_enrichment")
            log_message(sprintf("M03: Saved dual enrichment analysis plot: %s", plot_path))
          }
          
          # Save enrichment analysis CSV files
          # Remove timestamp variable - using standardized naming
          
          # Save core shared genes enrichment results
          if (!is.null(dual_enrichment_result$core_shared_enrichment)) {
            core_enrichment_data <- dual_enrichment_result$core_shared_enrichment
            log_message(sprintf("M03: Preparing to save core shared enrichment data (%d rows, %d cols)", 
                               nrow(core_enrichment_data), ncol(core_enrichment_data)))
            
            core_enrichment_csv <- "M03_core_shared_enrichment.csv"
            core_enrichment_path <- file.path(task_output_dir, core_enrichment_csv)
            
            tryCatch({
              write.csv(core_enrichment_data, core_enrichment_path, row.names = FALSE)
              log_message(sprintf("M03: [SUCCESS] Successfully saved core shared genes enrichment results: %s", core_enrichment_path))
            }, error = function(e) {
              log_message(sprintf("M03: [ERROR] Failed to save core shared enrichment: %s", e$message), level = "error")
            })
          } else {
            log_message("M03: Core shared enrichment data is NULL - skipping CSV export", level = "warning")
          }
          
          # Save private genes enrichment results
          if (!is.null(dual_enrichment_result$private_enrichment)) {
            private_enrichment_data <- dual_enrichment_result$private_enrichment
            log_message(sprintf("M03: Preparing to save private genes enrichment data (%d rows, %d cols)", 
                               nrow(private_enrichment_data), ncol(private_enrichment_data)))
                               
            private_enrichment_csv <- "M03_private_functional_enrichment.csv"
            private_enrichment_path <- file.path(task_output_dir, private_enrichment_csv)
            
            tryCatch({
              write.csv(private_enrichment_data, private_enrichment_path, row.names = FALSE)
              log_message(sprintf("M03: [SUCCESS] Successfully saved private genes enrichment results: %s", private_enrichment_path))
            }, error = function(e) {
              log_message(sprintf("M03: [ERROR] Failed to save private genes enrichment: %s", e$message), level = "error")
            })
          } else {
            log_message("M03: Private genes enrichment data is NULL - skipping CSV export", level = "warning")
          }
          
          # Save gene lists for both categories
          if (!is.null(dual_enrichment_result$gene_lists)) {
            gene_lists_csv <- "M03_general_enrichment_gene_lists.csv"
            gene_lists_path <- file.path(task_output_dir, gene_lists_csv)
            
            # Create comprehensive gene list data frame
            gene_lists_df <- data.frame(
              analysis_type = character(),
              gene_category = character(),
              gene_name = character(),
              p_value = numeric(),
              enrichment_ratio = numeric(),
              stringsAsFactors = FALSE
            )
            
            # Add core shared genes data
            if ("core_shared" %in% names(dual_enrichment_result$gene_lists)) {
              for (category in names(dual_enrichment_result$gene_lists$core_shared)) {
                genes <- dual_enrichment_result$gene_lists$core_shared[[category]]
                if (length(genes) > 0) {
                  gene_lists_df <- rbind(gene_lists_df, data.frame(
                    analysis_type = "Core_Shared_Genes",
                    gene_category = category,
                    gene_name = genes,
                    p_value = NA,  # Will be filled from enrichment results if available
                    enrichment_ratio = NA,
                    stringsAsFactors = FALSE
                  ))
                }
              }
            }
            
            # Add private genes data
            if ("private" %in% names(dual_enrichment_result$gene_lists)) {
              for (category in names(dual_enrichment_result$gene_lists$private)) {
                genes <- dual_enrichment_result$gene_lists$private[[category]]
                if (length(genes) > 0) {
                  gene_lists_df <- rbind(gene_lists_df, data.frame(
                    analysis_type = "Private_Genes",
                    gene_category = category,
                    gene_name = genes,
                    p_value = NA,
                    enrichment_ratio = NA,
                    stringsAsFactors = FALSE
                  ))
                }
              }
            }
            
            if (nrow(gene_lists_df) > 0) {
              write.csv(gene_lists_df, gene_lists_path, row.names = FALSE)
              log_message(sprintf("M03: Saved enrichment gene lists: %s", gene_lists_path))
            }
          }
        }
      }
    }
    
    # =================================================================
    # PHASE 5: STATISTICAL VALIDATION - Randomness Testing
    # =================================================================
    
    # Step 5: Perform Statistical Testing for Sharing Randomness
    if ("statistical_test" %in% plot_types) {
      log_message("M03 Step 5: Performing statistical tests for sharing pattern randomness")
      
      # CRITICAL FIX: Pass the complete core_hotspots structure which contains all_gene_sharing
      if (!is.null(core_hotspots) && !is.null(core_hotspots$all_gene_sharing)) {
        statistical_result <- perform_sharing_statistical_test(
          hotspot_analysis_results = core_hotspots,  # Pass the complete structure
          config = config,
          task_params = task_params
        )
      } else {
        log_message("Core hotspots analysis results missing - cannot perform statistical test", level = "error")
        statistical_result <- NULL
      }
      
      # Alternative approach if core_hotspots is missing - create minimal structure
      if (is.null(statistical_result) && !is.null(candidate_hotspots)) {
        log_message("Attempting statistical test with candidate hotspots data", level = "warning")
        
        # Create gene sharing summary from candidate hotspots
        temp_gene_sharing <- candidate_hotspots %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(
            num_species = dplyr::n_distinct(species),
            mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
            .groups = "drop"
          )
        
        statistical_result <- perform_sharing_statistical_test(
          hotspot_analysis_results = list(
            all_gene_sharing = temp_gene_sharing,
            parameters = list(
              total_species = length(unique(normalized_data$species))
            )
          ),
          config = config,
          task_params = task_params
        )
      }
      
      if (!is.null(statistical_result)) {
        log_message("M03: Statistical test completed successfully - preparing to save results")
        
        # Initialize stats_path to avoid undefined variable errors
        stats_filename <- "M03_core_shared_statistical_tests.txt"
        stats_path <- file.path(task_output_dir, stats_filename)
        
        # Save statistical test results - TXT report
        if (!is.null(statistical_result$report)) {
          tryCatch({
            writeLines(statistical_result$report, stats_path)
            log_message(sprintf("M03: [SUCCESS] Successfully saved statistical test report: %s", stats_path))
          }, error = function(e) {
            log_message(sprintf("M03: [ERROR] Failed to save statistical test report: %s", e$message), level = "error")
          })
        } else {
          log_message("M03: Statistical test report is NULL - skipping TXT export", level = "warning")
          # Create a placeholder report for debugging
          placeholder_report <- c(
            "=== STATISTICAL TEST REPORT ===",
            "Status: Statistical test completed but report data is NULL",
            sprintf("Analysis Date: %s", Sys.time()),
            "This indicates an issue with the statistical test function return values."
          )
          writeLines(placeholder_report, stats_path)
          log_message(sprintf("M03: Created placeholder statistical report: %s", stats_path))
        }
        
        # Save detailed CSV results - Chi-square test
        if (!is.null(statistical_result$comparison_table)) {
          comparison_data <- statistical_result$comparison_table
          log_message(sprintf("M03: Preparing to save chi-square test data (%d rows, %d cols)", 
                             nrow(comparison_data), ncol(comparison_data)))
                             
          chisq_csv_filename <- "M03_core_shared_chisquare_details.csv"
          chisq_csv_path <- file.path(task_output_dir, chisq_csv_filename)
          
          tryCatch({
            write.csv(comparison_data, chisq_csv_path, row.names = FALSE)
            log_message(sprintf("M03: [SUCCESS] Successfully saved chi-square test details: %s", chisq_csv_path))
          }, error = function(e) {
            log_message(sprintf("M03: [ERROR] Failed to save chi-square test details: %s", e$message), level = "error")
          })
        } else {
          log_message("M03: Chi-square test comparison table is NULL - skipping CSV export", level = "warning")
        }
        
        # Save detailed CSV results - Gene-by-gene hypergeometric test
        if (!is.null(statistical_result$tests$hypergeometric)) {
          log_message("M03: Preparing to save gene-by-gene hypergeometric test data")
          
          # Save individual gene-level test results
          if (!is.null(statistical_result$tests$hypergeometric$all_gene_results) && 
              nrow(statistical_result$tests$hypergeometric$all_gene_results) > 0) {
            
            gene_level_data <- statistical_result$tests$hypergeometric$all_gene_results
            
            log_message(sprintf("M03: Gene-level hypergeometric data ready (%d genes, %d cols)", 
                               nrow(gene_level_data), ncol(gene_level_data)))
            
            # Save complete gene-level results
            gene_level_csv_filename <- "M03_core_shared_hypergeometric_tests.csv"
            gene_level_csv_path <- file.path(task_output_dir, gene_level_csv_filename)
            
            tryCatch({
              write.csv(gene_level_data, gene_level_csv_path, row.names = FALSE)
              log_message(sprintf("M03: [SUCCESS] Successfully saved gene-level hypergeometric tests: %s", gene_level_csv_path))
            }, error = function(e) {
              log_message(sprintf("M03: [ERROR] Failed to save gene-level hypergeometric tests: %s", e$message), level = "error")
            })
            
            # Save driver genes summary
            driver_genes_summary <- data.frame(
              category = c(rep("Private", nrow(statistical_result$tests$hypergeometric$private_drivers %||% data.frame())),
                          rep("Core", nrow(statistical_result$tests$hypergeometric$core_drivers %||% data.frame()))),
              rbind(
                statistical_result$tests$hypergeometric$private_drivers %||% data.frame()[0,],
                statistical_result$tests$hypergeometric$core_drivers %||% data.frame()[0,]
              ),
              stringsAsFactors = FALSE
            )
            
            if (nrow(driver_genes_summary) > 0) {
              driver_csv_filename <- "M03_core_shared_driver_genes_summary.csv"
              driver_csv_path <- file.path(task_output_dir, driver_csv_filename)
              
              tryCatch({
                write.csv(driver_genes_summary, driver_csv_path, row.names = FALSE)
                log_message(sprintf("M03: [SUCCESS] Successfully saved driver genes summary: %s", driver_csv_path))
              }, error = function(e) {
                log_message(sprintf("M03: [ERROR] Failed to save driver genes summary: %s", e$message), level = "error")
              })
            }
            
            # Create legacy format summary for backward compatibility
            hypergeom_summary <- data.frame(
              test_method = "Gene-by-gene Hypergeometric Analysis",
              total_genes_tested = statistical_result$tests$hypergeometric$total_genes_tested %||% 0,
              private_genes_tested = statistical_result$tests$hypergeometric$private_genes_tested %||% 0,
              core_genes_tested = statistical_result$tests$hypergeometric$core_genes_tested %||% 0,
              significant_private_drivers = statistical_result$tests$hypergeometric$significant_private_drivers %||% 0,
              significant_core_drivers = statistical_result$tests$hypergeometric$significant_core_drivers %||% 0,
              alpha_level = statistical_result$tests$hypergeometric$alpha_level %||% 0.05,
              method_description = statistical_result$tests$hypergeometric$method_description %||% "Gene-by-gene analysis",
              stringsAsFactors = FALSE
            )
          } else {
            # Fallback summary if no detailed results
            hypergeom_summary <- data.frame(
              test_method = "Gene-by-gene Hypergeometric Analysis",
              total_genes_tested = 0,
              private_genes_tested = 0,
              core_genes_tested = 0,
              significant_private_drivers = 0,
              significant_core_drivers = 0,
              alpha_level = 0.05,
              method_description = "No gene-level results available",
              stringsAsFactors = FALSE
            )
          }
          
          # Use summary as the main hypergeometric data for legacy compatibility
          hypergeom_data <- hypergeom_summary
          
          log_message(sprintf("M03: Hypergeometric test data ready (%d rows, %d cols)", 
                             nrow(hypergeom_data), ncol(hypergeom_data)))
          
          hypergeom_csv_filename <- "M03_core_shared_hypergeometric_details.csv"
          hypergeom_csv_path <- file.path(task_output_dir, hypergeom_csv_filename)
          
          tryCatch({
            write.csv(hypergeom_data, hypergeom_csv_path, row.names = FALSE)
            log_message(sprintf("M03: [SUCCESS] Successfully saved hypergeometric test details: %s", hypergeom_csv_path))
          }, error = function(e) {
            log_message(sprintf("M03: [ERROR] Failed to save hypergeometric test details: %s", e$message), level = "error")
          })
        } else {
          log_message("M03: Hypergeometric test data is NULL - skipping CSV export", level = "warning")
        }
        
        # Save comprehensive statistics summary
        if (!is.null(statistical_result$sharing_summary)) {
          summary_data <- data.frame(
            statistic = c("total_genes", "total_species", "mean_sharing", "median_sharing", 
                         "max_sharing", "genes_in_multiple_species", "genes_in_all_species"),
            value = c(
              statistical_result$sharing_summary$total_genes,
              statistical_result$sharing_summary$total_species,
              statistical_result$sharing_summary$mean_sharing,
              statistical_result$sharing_summary$median_sharing,
              statistical_result$sharing_summary$max_sharing,
              statistical_result$sharing_summary$genes_in_multiple_species,
              statistical_result$sharing_summary$genes_in_all_species
            ),
            stringsAsFactors = FALSE
          )
          
          summary_csv_filename <- "M03_core_shared_sharing_statistics.csv"
          summary_csv_path <- file.path(task_output_dir, summary_csv_filename)
          write.csv(summary_data, summary_csv_path, row.names = FALSE)
          log_message(sprintf("M03: Saved sharing statistics summary: %s", summary_csv_path))
        }
        
        results[["statistical_test"]] <- list(
          data = statistical_result,
          file_path = stats_path,
          detailed_files = list(
            chisq_details = if(exists("chisq_csv_path")) chisq_csv_path else NULL,
            hypergeom_details = if(exists("hypergeom_csv_path")) hypergeom_csv_path else NULL,
            summary_stats = if(exists("summary_csv_path")) summary_csv_path else NULL
          ),
          metadata = list(
            tests_performed = names(statistical_result$tests),
            significant_results = statistical_result$significant_results
          )
        )
        plots_generated <- c(plots_generated, "statistical_test")
        log_message(sprintf("M03: Saved statistical test results: %s", stats_path))
        
        # =================================================================
        # PHASE 5b: POST-STATISTICAL FUNCTIONAL ENRICHMENT ANALYSIS
        # =================================================================
        
        # After statistical tests, perform separate functional enrichment for core and private genes
        log_message("M03 Step 5b: Performing functional enrichment analysis AFTER statistical testing")
        
        enrichment_results <- tryCatch({
          # Load gene function mapping
          gene_function_data <- load_gene_function_mapping(config = config)
          
          if (!is.null(gene_function_data) && nrow(gene_function_data) > 0 && 
              !is.null(core_hotspots) && !is.null(core_hotspots$all_gene_sharing)) {
            
            # Extract core shared and private genes from statistical results
            core_genes_list <- core_hotspots$core_shared_genes$gene
            private_genes_list <- core_hotspots$private_genes$gene
            all_candidate_genes <- unique(candidate_hotspots$gene)
            
            log_message(sprintf("Post-statistical enrichment: %d core genes, %d private genes, %d total candidates",
                               length(core_genes_list), length(private_genes_list), length(all_candidate_genes)))
            
            enrichment_results_combined <- list()
            
            # 1. Core shared genes functional enrichment
            if (length(core_genes_list) > 0) {
              core_genes_with_annotations <- intersect(core_genes_list, gene_function_data$gene)
              
              if (length(core_genes_with_annotations) > 0) {
                log_message(sprintf("Performing functional enrichment for %d core shared genes", length(core_genes_with_annotations)))
                
                # Get expected categories from configuration  
                expected_categories <- get_task_parameter(task_params, config, "expected_gene_categories", 
                                                         c("Other Functional Genes", "Unknown Function Genes", 
                                                           "Photosynthesis Related", "Self-Replication Related"))
                
                core_enrichment <- perform_functional_enrichment(
                  target_genes = core_genes_list,
                  gene_function_data = gene_function_data,
                  background_genes = all_candidate_genes,
                  expected_categories = expected_categories
                )
                
                enrichment_results_combined$core_shared <- core_enrichment
                
                # Save core enrichment results
                if (!is.null(core_enrichment) && is.data.frame(core_enrichment) && nrow(core_enrichment) > 0) {
                  core_enrichment_path <- file.path(task_output_dir, 
                    "M03_core_shared_functional_enrichment.csv")
                  write.csv(core_enrichment, core_enrichment_path, row.names = FALSE)
                  log_message(sprintf("M03: [SUCCESS] Saved core shared genes functional enrichment: %s", core_enrichment_path))
                } else {
                  log_message("M03: No significant enrichment found for core shared genes", level = "warning")
                }
              } else {
                log_message("M03: No core shared genes have functional annotations", level = "warning")
              }
            }
            
            # 2. Private genes functional enrichment  
            if (length(private_genes_list) > 0) {
              private_genes_with_annotations <- intersect(private_genes_list, gene_function_data$gene)
              
              if (length(private_genes_with_annotations) > 0) {
                log_message(sprintf("Performing functional enrichment for %d private genes", length(private_genes_with_annotations)))
                
                private_enrichment <- perform_functional_enrichment(
                  target_genes = private_genes_list,
                  gene_function_data = gene_function_data,
                  background_genes = all_candidate_genes,
                  expected_categories = expected_categories
                )
                
                enrichment_results_combined$private <- private_enrichment
                
                # Save private enrichment results
                if (!is.null(private_enrichment) && is.data.frame(private_enrichment) && nrow(private_enrichment) > 0) {
                  private_enrichment_path <- file.path(task_output_dir, 
                    "M03_private_functional_enrichment.csv")
                  write.csv(private_enrichment, private_enrichment_path, row.names = FALSE)
                  log_message(sprintf("M03: [SUCCESS] Saved private genes functional enrichment: %s", private_enrichment_path))
                } else {
                  log_message("M03: No significant enrichment found for private genes", level = "warning")
                }
              } else {
                log_message("M03: No private genes have functional annotations", level = "warning")
              }
            }
            
            # 3. Generate combined enrichment summary report with enhanced transparency
            # Always generate report to ensure consistency with logs
            enrichment_summary_lines <- c(
              "=== POST-STATISTICAL FUNCTIONAL ENRICHMENT ANALYSIS ===",
              sprintf("Analysis performed after statistical testing on %s", Sys.time()),
              ""
            )
            
            # Core shared genes enrichment reporting with detailed diagnostics
            core_shared_result <- enrichment_results_combined$core_shared
            if (!is.null(core_shared_result)) {
              if (is.data.frame(core_shared_result) && nrow(core_shared_result) > 0) {
                n_significant_core <- sum((core_shared_result$p_adjusted %||% c()) < 0.05, na.rm = TRUE)
                enrichment_summary_lines <- c(enrichment_summary_lines,
                  "=== CORE SHARED GENES FUNCTIONAL ENRICHMENT ===",
                  sprintf("Core shared genes enrichment: %d categories analyzed, %d significant", 
                         nrow(core_shared_result), n_significant_core),
                  sprintf("Categories tested: %s", 
                         paste(core_shared_result$functional_category %||% "Unknown", collapse = ", ")),
                  if (n_significant_core > 0) {
                    c("Significant categories:",
                      paste(sprintf("  - %s (p_adj=%.3f, ratio=%.2f)", 
                                  core_shared_result$functional_category[core_shared_result$p_adjusted < 0.05],
                                  core_shared_result$p_adjusted[core_shared_result$p_adjusted < 0.05],
                                  core_shared_result$enrichment_ratio[core_shared_result$p_adjusted < 0.05]), 
                           collapse = "\n"))
                  } else {
                    "No significantly enriched categories found"
                  })
              } else {
                enrichment_summary_lines <- c(enrichment_summary_lines,
                  "=== CORE SHARED GENES FUNCTIONAL ENRICHMENT ===",
                  sprintf("Core shared genes enrichment: Analysis attempted but returned empty results (data format: %s, length: %d)", 
                         class(core_shared_result)[1], length(core_shared_result)))
              }
            } else {
              enrichment_summary_lines <- c(enrichment_summary_lines,
                "=== CORE SHARED GENES FUNCTIONAL ENRICHMENT ===",
                "Core shared genes enrichment: 0 categories analyzed, 0 significant",
                "Reason: No core shared genes available for analysis")
            }
            
            enrichment_summary_lines <- c(enrichment_summary_lines, "")
            
            # Private genes enrichment reporting with detailed diagnostics
            private_result <- enrichment_results_combined$private
            if (!is.null(private_result)) {
              if (is.data.frame(private_result) && nrow(private_result) > 0) {
                n_significant_private <- sum((private_result$p_adjusted %||% c()) < 0.05, na.rm = TRUE)
                enrichment_summary_lines <- c(enrichment_summary_lines,
                  "=== PRIVATE GENES FUNCTIONAL ENRICHMENT ===",
                  sprintf("Private genes enrichment: %d categories analyzed, %d significant",
                         nrow(private_result), n_significant_private),
                  sprintf("Categories tested: %s", 
                         paste(private_result$functional_category %||% "Unknown", collapse = ", ")),
                  if (n_significant_private > 0) {
                    c("Significant categories:",
                      paste(sprintf("  - %s (p_adj=%.3f, ratio=%.2f)", 
                                  private_result$functional_category[private_result$p_adjusted < 0.05],
                                  private_result$p_adjusted[private_result$p_adjusted < 0.05],
                                  private_result$enrichment_ratio[private_result$p_adjusted < 0.05]), 
                           collapse = "\n"))
                  } else {
                    "No significantly enriched categories found"
                  })
              } else {
                enrichment_summary_lines <- c(enrichment_summary_lines,
                  "=== PRIVATE GENES FUNCTIONAL ENRICHMENT ===",
                  sprintf("Private genes enrichment: Analysis attempted but returned empty results (data format: %s, length: %d)", 
                         class(private_result)[1], length(private_result)))
              }
            } else {
              enrichment_summary_lines <- c(enrichment_summary_lines,
                "=== PRIVATE GENES FUNCTIONAL ENRICHMENT ===",
                "Private genes enrichment: 0 categories analyzed, 0 significant",
                "Reason: No private genes available for analysis")
            }
            
            enrichment_summary_lines <- c(enrichment_summary_lines,
              "",
              "=== FUNCTIONAL ENRICHMENT CONCLUSIONS ===",
              ifelse(!is.null(core_shared_result) && !is.null(private_result),
                "Both core shared and private genes functional enrichment analysis completed",
                ifelse(!is.null(core_shared_result),
                  "Only core shared genes functional enrichment analysis completed", 
                  ifelse(!is.null(private_result),
                    "Only private genes functional enrichment analysis completed",
                    "No functional enrichment analysis could be performed")))
            )
            
            # Save enrichment summary report
            enrichment_summary_path <- file.path(task_output_dir, 
              "M03_general_functional_enrichment_summary.txt")
            writeLines(enrichment_summary_lines, enrichment_summary_path)
            log_message(sprintf("M03: [SUCCESS] Saved functional enrichment summary: %s", enrichment_summary_path))
            
            enrichment_results_combined
          } else {
            log_message("M03: Gene function mapping or core hotspots data not available for post-statistical enrichment", level = "warning")
            NULL
          }
        }, error = function(e) {
          log_message(sprintf("M03: Failed to perform post-statistical functional enrichment: %s", e$message), level = "error")
          NULL
        })
        
        if (!is.null(enrichment_results)) {
          log_message("M03: Post-statistical functional enrichment analysis completed successfully")
          plots_generated <- c(plots_generated, "post_statistical_enrichment")
        }
      }
    }
    
    # =================================================================
    # PHASE 6: DATA PERSISTENCE - Save All Intermediate Results
    # =================================================================
    
    # CRITICAL FIX: Save the comprehensive ranked list of all genes for transparency
    log_message(sprintf("M03: Preparing to save full ranked gene data - all_genes_ranked: %s, nrow: %d", 
                       !is.null(all_genes_ranked), 
                       if (!is.null(all_genes_ranked)) nrow(all_genes_ranked) else 0))
    
    if (!is.null(all_genes_ranked) && nrow(all_genes_ranked) > 0) {
      # This file contains ALL genes from the initial filter, with ranks and thresholds.
      # This provides complete transparency for the hotspot identification process.
      ranked_data_path <- file.path(task_output_dir, "M03_all_genes_ranked_with_thresholds.csv")
      
      tryCatch({
        write.csv(all_genes_ranked, ranked_data_path, row.names = FALSE)
        log_message(sprintf("M03: [SUCCESS] Successfully saved full ranked gene data for transparency: %s", ranked_data_path))
        log_message(sprintf("M03: Full ranked genes file contains %d rows and %d columns", 
                           nrow(all_genes_ranked), ncol(all_genes_ranked)))
      }, error = function(e) {
        log_message(sprintf("M03: [ERROR] Failed to save full ranked genes CSV: %s", e$message), level = "error")
      })
    } else {
      log_message("M03: [ERROR] Cannot save full ranked genes CSV - data is NULL or empty", level = "error")
    }
    
    # Save the definitive list of true candidate hotspots used in all analyses
    log_message(sprintf("M03: Preparing to save definitive candidate hotspots - candidate_hotspots: %s, nrow: %d", 
                       !is.null(candidate_hotspots), 
                       if (!is.null(candidate_hotspots)) nrow(candidate_hotspots) else 0))
    
    if (!is.null(candidate_hotspots) && nrow(candidate_hotspots) > 0) {
      # This file contains ONLY the true hotspots used for all downstream analyses
      # This replaces the old, confusing file structure
      hotspots_path <- file.path(task_output_dir, "M03_candidate_hotspots.csv")
      
      tryCatch({
        write.csv(candidate_hotspots, hotspots_path, row.names = FALSE)
        log_message(sprintf("M03: [SUCCESS] Successfully saved definitive candidate hotspots used for analysis: %s", hotspots_path))
        log_message(sprintf("M03: Definitive hotspots file contains %d rows and %d columns", 
                           nrow(candidate_hotspots), ncol(candidate_hotspots)))
      }, error = function(e) {
        log_message(sprintf("M03: [ERROR] Failed to save definitive candidate hotspots CSV: %s", e$message), level = "error")
      })
      
      # Create and save a summary of the definitive hotspots
      hotspots_summary <- candidate_hotspots %>%
        dplyr::summarise(
          total_hotspot_entries = dplyr::n(),
          unique_genes = dplyr::n_distinct(gene),
          unique_species = dplyr::n_distinct(species),
          min_frequency = min(frequency_per_kb, na.rm = TRUE),
          max_frequency = max(frequency_per_kb, na.rm = TRUE),
          mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
          median_frequency = stats::median(frequency_per_kb, na.rm = TRUE),
          threshold_percentile = frequency_percentile,
          .groups = "drop"
        )
      
      summary_path <- file.path(task_output_dir, "M03_candidate_hotspots_summary.csv")
      write.csv(hotspots_summary, summary_path, row.names = FALSE)
      
      log_message(sprintf("M03: [SUCCESS] Definitive hotspots summary: %d entries from %d unique genes across %d species", 
                         hotspots_summary$total_hotspot_entries, hotspots_summary$unique_genes, hotspots_summary$unique_species))
      log_message(sprintf("M03: [SUCCESS] Saved definitive hotspots summary statistics: %s", summary_path))
    } else {
      log_message("M03: [ERROR] Cannot save definitive candidate hotspots CSV - data is NULL or empty", level = "error")
      log_message(sprintf("M03: candidate_hotspots diagnostic - is.null: %s, is.data.frame: %s", 
                         is.null(candidate_hotspots),
                         if (!is.null(candidate_hotspots)) is.data.frame(candidate_hotspots) else "N/A"))
    }
    
    if (!is.null(core_hotspots) && !is.null(core_hotspots$core_shared_genes) && nrow(core_hotspots$core_shared_genes) > 0) {
      core_hotspots_path <- file.path(task_output_dir, "M03_core_shared_hotspots.csv")
      write.csv(core_hotspots$core_shared_genes, core_hotspots_path, row.names = FALSE)
      log_message(sprintf("M03: Saved core hotspot data: %s", core_hotspots_path))
    }
    
    # Save private hotspots data if available
    if (!is.null(core_hotspots) && !is.null(core_hotspots$private_genes) && nrow(core_hotspots$private_genes) > 0) {
      private_hotspots_path <- file.path(task_output_dir, "M03_private_hotspots.csv")
      write.csv(core_hotspots$private_genes, private_hotspots_path, row.names = FALSE)
      log_message(sprintf("M03: Saved private hotspot data: %s", private_hotspots_path))
    }
    
    # =================================================================
    # PHASE 7: COMPREHENSIVE SUMMARY GENERATION
    # =================================================================
    
    # Get color palette for summary (following V4.1 architecture)
    color_palette <- get_color_palette(config, "M03_hotspot")
    
    # Generate enhanced task summary with all analysis components
    task_summary <- list(
      task_name = task_name,
      module = "M03_hotspot_enhanced",
      analysis_phases = c("global_analysis", "hotspot_identification", "sharing_patterns", "dual_threshold", "enrichment", "statistical_validation"),
      plots_generated = plots_generated,
      output_directory = task_output_dir,
      
      # Data metrics
      candidate_hotspots_count = nrow(candidate_hotspots),
      unique_genes_count = length(unique(candidate_hotspots$gene)),
      unique_species_count = length(unique(candidate_hotspots$species)),
      core_shared_genes_count = if (!is.null(core_hotspots) && !is.null(core_hotspots$genes)) nrow(core_hotspots$genes) else 0,
      
      # Enhanced parameters
      parameters_used = list(
        frequency_percentile = frequency_percentile,
        min_species_ratio = min_species_ratio,
        min_species_number = min_species_number,
        threshold_mode = threshold_mode,
        effective_threshold = effective_threshold,
        region_filter = region_filter,
        var_type_filter = var_type_filter,
        plot_types = plot_types,
        color_palette = color_palette,
        clustering_enabled = list(rows = cluster_rows, cols = cluster_cols)
      ),
      
      # Analysis workflow metadata
      workflow_version = "M03_Enhanced_8Step_v2.0",
      analysis_date = Sys.time(),
      total_analysis_steps = 8
    )
    
    log_message(sprintf("M03 Enhanced 8-Step Analysis completed: %d analysis components generated", length(plots_generated)))
    
    return(list(
      plots = results,
      summary = task_summary,
      success = length(plots_generated) > 0
    ))
    
  }, error = function(e) {
    log_message(sprintf("M03 Enhanced analysis failed: %s", e$message), level = "error")
    return(list(
      plots = list(),
      summary = list(error = e$message, workflow_version = "M03_Enhanced_8Step_v2.0"),
      success = FALSE
    ))
  })
}

############################################################
#### M01_distribution_module.R - M01 Module Orchestrator ####
############################################################
#
# Thin orchestrator for M01 Data Distribution Analysis
# This is a Phase 2 "thin wrapper" that loads dependencies and delegates
# to the existing generate_M01_data_distribution_plots function
#
# Phase 2 Rule #3 Compliance: Does NOT modify existing function libraries,
# only provides a clean interface for task_dispatcher.R
#
############################################################

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

#' Smart Inherit Data Scoping Configuration
#' @description Implements hierarchical inheritance with intelligent defaults based on submodule name
#' @param submodule_name Name of the submodule (e.g., "chromosome_distribution", "cds_snp_annotation")
#' @param submodule_config Submodule-specific configuration
#' @param global_config Global M01 configuration
#' @return Merged data_scoping configuration
#' @keywords internal
inherit_data_scoping <- function(submodule_name, submodule_config, global_config) {
  # Extract global data_scoping (M01 level)
  global_scoping <- global_config$data_scoping %||% list()
  
  # Extract explicit submodule data_scoping (if provided)
  submodule_scoping <- submodule_config$data_scoping %||% list()
  
  # Define intelligent defaults based on submodule name
  # CRITICAL: Use "all_types" string to indicate "no filtering" (vs NULL = "inherit global")
  smart_defaults <- list(
    # M01-2: Structural & Relational
    chromosome_distribution = list(
      target_region_types = "all_types",  # No region filtering (accept all)
      target_var_types = c("snp", "INDEL")
    ),
    poigs_gene_correlation = list(
      target_region_types = c("CDS", "IGS"),
      target_var_types = "all_types"  # No var_type filtering (accept all)
    ),
    chromosome_hotspot_map = list(
      target_region_types = c("CDS", "IGS"),  # CDS and IGS hotspots
      target_var_types = "all_types"  # No var_type filtering (accept all)
    ),
    
    # M01-3: Sequence Variant Characterization
    cds_snp_annotation = list(
      target_region_types = c("CDS"),
      target_var_types = c("snp")
    ),
    igs_variant_features = list(
      target_region_types = "IGS",
      target_var_types = c("snp", "INDEL")
    ),
    genome_region_stats = list(
      target_region_types = "all_types",  # No region filtering (accept all)
      target_var_types = c("snp", "INDEL")
    )
  )
  
  # Get smart default for this submodule
  default_for_submodule <- smart_defaults[[submodule_name]] %||% list()
  
  # Priority: explicit submodule > smart default > global
  # CRITICAL: Handle "all_types" as special value meaning "no filtering"
  merged_scoping <- list(
    target_region_types = {
      if (!is.null(submodule_scoping$target_region_types)) {
        submodule_scoping$target_region_types
      } else if (!is.null(default_for_submodule$target_region_types)) {
        default_for_submodule$target_region_types  # May be "all_types"
      } else {
        global_scoping$target_region_types
      }
    },
    target_var_types = {
      if (!is.null(submodule_scoping$target_var_types)) {
        submodule_scoping$target_var_types
      } else if (!is.null(default_for_submodule$target_var_types)) {
        default_for_submodule$target_var_types  # May be "all_types"
      } else {
        global_scoping$target_var_types
      }
    }
  )
  
  return(merged_scoping)
}

#' M01 Distribution Analysis Module Orchestrator
#' 
#' @title Run M01 distribution analysis
#' @description Implements differential data filtering for different plot types.
#' Applies biologically meaningful NA/0 distinction and species ordering.
#' @param pie_data Site-level data for pie charts (var_type separated)
#' @param freq_data Gene-level data for other plots (var_type="all_types")
#' @param config Configuration object
#' @param output_dir Output directory for plots
#' @param task_name Task name for output organization
#' @param task_params Task-specific parameters
#' @param normalized_data Optional complete normalized data for relational analyses
#' @param filtered_data Optional filtered variant data for sequence analyses
#' @param preprocessing_results Optional preprocessing outputs used by structural analyses
#' @return M01 analysis results
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
run_m01_distribution <- function(pie_data,
                                 freq_data,
                                 config,
                                 output_dir,
                                 task_name,
                                 task_params,
                                 normalized_data = NULL,
                                 filtered_data = NULL,
                                 preprocessing_results = NULL) {
  
  # Function dependencies are automatically loaded in R package context
  
  log_message(sprintf("=== M01 DATA DISTRIBUTION ANALYSIS: %s ===", task_name))
  log_message(sprintf("Raw pie data: %d rows, %d species", 
                     nrow(pie_data), 
                     length(unique(pie_data$species))))
  log_message(sprintf("Raw frequency data: %d rows, %d species", 
                     nrow(freq_data), 
                     length(unique(freq_data$species))))
  
  # Validate input data
  if (nrow(pie_data) == 0) {
    log_message("ERROR: pie_data is empty", level = "error")
    stop("M01 analysis cannot proceed with empty pie_data")
  }
  
  if (nrow(freq_data) == 0) {
    log_message("ERROR: freq_data is empty", level = "error")
    stop("M01 analysis cannot proceed with empty freq_data")
  }
  
  # Get session paths for species ordering file
  session_id <- config$session_info$session_id %||% "unknown"
  session_paths <- get_session_paths(session_id)
  
  # Load species ordering from user-provided file
  species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
  species_order <- NULL
  
  if (file.exists(species_order_file)) {
    log_message(sprintf("Loading species ordering from: %s", species_order_file))
    species_order_data <- utils::read.csv(species_order_file, stringsAsFactors = FALSE)
    # Assume first column contains species names, adjust column name as needed
    species_order <- species_order_data[[1]]  # or species_order_data$species
    log_message(sprintf("Loaded species ordering: %d species", length(species_order)))
  } else {
    log_message("Species ordering file not found, will use default alphabetical order", level = "warning")
  }
  
  # Fix data types to prevent "character vector argument expected" errors
  # Convert factor columns to character before calling visualization functions
  # Apply to pie_data
  if ("region_type" %in% names(pie_data) && is.factor(pie_data$region_type)) {
    pie_data$region_type <- as.character(pie_data$region_type)
  }
  if ("var_type" %in% names(pie_data) && is.factor(pie_data$var_type)) {
    pie_data$var_type <- as.character(pie_data$var_type)
  }
  if ("species" %in% names(pie_data) && is.factor(pie_data$species)) {
    pie_data$species <- as.character(pie_data$species)
  }
  
  # Apply to freq_data
  if ("region_type" %in% names(freq_data) && is.factor(freq_data$region_type)) {
    freq_data$region_type <- as.character(freq_data$region_type)
  }
  if ("var_type" %in% names(freq_data) && is.factor(freq_data$var_type)) {
    freq_data$var_type <- as.character(freq_data$var_type)
  }
  if ("species" %in% names(freq_data) && is.factor(freq_data$species)) {
    freq_data$species <- as.character(freq_data$species)
  }
  
  # Apply data_scoping filters (PROJECT_CHARTER.md Principle 3)
  log_message("=== APPLYING DATA_SCOPING FILTERS ===")
  
  # Validate data dimensions and quality before processing
  log_message(sprintf("Raw freq_data dimensions: %d rows", nrow(freq_data)))
  log_message(sprintf("Frequency range: %.3f to %.3f", 
                     min(freq_data$frequency_per_kb, na.rm=TRUE),
                     max(freq_data$frequency_per_kb, na.rm=TRUE)))
  
  # Unify M01 data aggregation so IGS identities are preserved
  # Architecture principle: M01 module creates its own summary data from detailed data
  # This ensures consistency with M02 and maintains upstream P03 independence
  
  # Step 1: Apply data_scoping to detailed data (pie_data)
  log_message("=== M01: APPLYING DATA SCOPING TO DETAILED DATA ===")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    log_message(jsonlite::toJSON(task_params, auto_unbox = TRUE, pretty = TRUE))
  } else {
    log_message(paste("task_params structure:", capture.output(str(task_params)), collapse = "\n"))
  }
  pie_data_scoped <- apply_data_scoping(pie_data, task_params)
  log_message(sprintf("Pie data scoping: %d -> %d rows", nrow(pie_data), nrow(pie_data_scoped)))
  
  # Step 2: Re-aggregate summary data from filtered detailed data
  log_message("=== M01: RE-AGGREGATING SUMMARY DATA FROM DETAILED DATA ===")
  log_message(sprintf("Input: %d rows of detailed data (after data_scoping)", nrow(pie_data_scoped)))
  
  # Re-aggregate: sum variant_count across all var_types for each region
  # Group by region_name rather than gene so IGS identities are preserved
  m01_summary_data <- pie_data_scoped %>%
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
  
  log_message(sprintf("Output: %d rows of aggregated summary data", nrow(m01_summary_data)))
  log_message(sprintf("  Aggregation: %d detailed rows -> %d summary rows (%.1f%% reduction)", 
                     nrow(pie_data_scoped), nrow(m01_summary_data), 
                     (1 - nrow(m01_summary_data)/nrow(pie_data_scoped)) * 100))
  
  # Update data references to use scoped and aggregated data
  pie_data <- pie_data_scoped
  freq_data <- m01_summary_data
  
  # DIFFERENTIAL DATA FILTERING: Apply biologically meaningful NA/0 distinction
  log_message("=== APPLYING DIFFERENTIAL DATA FILTERING ===")
  
  # V16 REGRESSION FIX: Revert to the simple, working logic from the backup.
  # All previous attempts to "fix" factor levels were the source of the error.
  
  # Frequency plots data: Remove NA frequencies and zero frequencies
  freq_plot_data <- freq_data %>% 
    dplyr::filter(!is.na(frequency_per_kb) & frequency_per_kb > 0)
  
  # Prepare pie-chart data
  # Pie charts only display variant COUNTS, not frequencies
  # Remove region_length and frequency_per_kb to avoid NA-related filtering issues
  pie_plot_data <- pie_data %>% 
    dplyr::filter(variant_count > 0) %>%
    dplyr::select(species, var_type, region_type, genome_region, gene, variant_count)
  
  log_message(sprintf("Data filtering results:"))
  log_message(sprintf("  Frequency plots: %d rows (after removing NA/0 frequencies)", nrow(freq_plot_data)))
  log_message(sprintf("  Pie chart: %d rows (after removing zero variants, frequency columns removed)", nrow(pie_plot_data)))
  
  # Validate RNA regions handling
  pie_rna_count <- sum(pie_plot_data$region_type == "RNA", na.rm = TRUE)
  freq_rna_count <- sum(freq_plot_data$region_type == "RNA", na.rm = TRUE)
  
  log_message(sprintf("RNA regions: Pie chart=%d rows, Frequency plots=%d rows", 
                     pie_rna_count, freq_rna_count))
  
  # Apply species ordering if available
  if (!is.null(species_order)) {
    log_message("Applying species ordering to both datasets")
    freq_plot_data <- apply_species_ordering(freq_plot_data, species_order_vector = species_order, species_col = "species")
    pie_plot_data <- apply_species_ordering(pie_plot_data, species_order_vector = species_order, species_col = "species")
  }
  
  # V3.16.14: Unified directory structure - data and plots in same directory
  # Save to M01_distribution_standard (flat structure, not nested)
  data_output_dir <- file.path(output_dir, "M01_distribution", "M01_distribution_standard")
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save cleaned datasets for verification
  pie_data_file <- file.path(data_output_dir, "M01_pie_data_cleaned.csv")
  freq_data_file <- file.path(data_output_dir, "M01_freq_data_cleaned.csv")
  utils::write.csv(pie_plot_data, pie_data_file, row.names = FALSE)
  utils::write.csv(freq_plot_data, freq_data_file, row.names = FALSE)
  
  # Generate data quality report
  quality_report <- data.frame(
    data_type = c("pie_chart", "frequency_plots"),
    total_rows = c(nrow(pie_plot_data), nrow(freq_plot_data)),
    rna_regions = c(pie_rna_count, freq_rna_count),
    total_variants = c(
      sum(pie_plot_data$variant_count, na.rm = TRUE),
      sum(freq_plot_data$variant_count, na.rm = TRUE)
    ),
    species_count = c(
      length(unique(pie_plot_data$species)),
      length(unique(freq_plot_data$species))
    )
  )
  utils::write.csv(quality_report, file.path(data_output_dir, "data_quality_report.csv"), row.names = FALSE)
  
  log_message(sprintf("Cleaned data saved to: %s", data_output_dir))
  
  # Validate filtered data
  if (nrow(freq_plot_data) == 0) {
    log_message("ERROR: No valid frequency data after filtering", level = "error")
    stop("Cannot proceed with frequency-based plots - no valid data")
  }
  
  if (nrow(pie_plot_data) == 0) {
    log_message("ERROR: No variants found for pie chart", level = "error")
    stop("Cannot proceed with pie chart - no variant data")
  }
  
  log_message(sprintf("Final data validation - Pie data var_types: %s", 
                     paste(unique(pie_plot_data$var_type), collapse = ", ")))
  log_message(sprintf("Final data validation - Freq data var_types: %s", 
                     paste(unique(freq_plot_data$var_type), collapse = ", ")))

  # Pass cleaned and ordered data to plotting functions
  result <- generate_M01_data_distribution_plots(
    pie_data = pie_plot_data,
    freq_data = freq_plot_data,
    config = config,
    output_dir = output_dir,
    task_name = task_name,
    task_params = task_params
  )
  
  # ================================================================
  # Phase 1 Pilot: Generate Manifest for Scientific Report System
  # ================================================================
  log_message("Generating output manifest for scientific report system")
  
  manifest_entries <- list()
  
  if (!is.null(result$plots) && length(result$plots) > 0) {
    
    for (plot_type in names(result$plots)) {
      plot_info <- result$plots[[plot_type]]
      
      if (!is.null(plot_info$file_path) && file.exists(plot_info$file_path)) {
        
        # Determine role based on plot type
        role <- if (plot_type %in% c("region_boxplot_composite")) {
          "composite_plot"  # High priority for final report
        } else {
          "individual_plot"  # Medium priority
        }
        
        # Generate academic title
        title <- get_standard_title("M01", role, plot_type)
        
        # Generate standard caption
        caption <- get_standard_caption("M01", role, plot_type)
        
        # Create manifest entry
        entry <- create_manifest_entry(
          path = plot_info$file_path,
          role = role,
          title = title,
          caption = caption,
          module = "M01",
          file_type = "plot",
          priority = if (role == "composite_plot") 1 else 3,
          parameters = task_params  # V7.1 Enhancement: Pass task parameters
        )
        
        manifest_entries[[plot_type]] <- entry
        log_message(sprintf("[SUCCESS] Created manifest entry for %s (%s)", plot_type, role))
      }
    }
  }
  
  # Create manifest collection
  manifest <- create_manifest_collection(manifest_entries, "M01", session_id)
  
  # Add manifest to result
  result$manifest <- manifest
  
  # ================================================================
  # Part B: M01-2 Structural & Relational Analysis (NEW)
  # ================================================================
  
  # Check if structural/relational analyses are enabled
  enabled_submodules <- get_task_parameter(
    task_params,
    config,
    "enabled_submodules",
    default = c("foundational"),
    log_default = TRUE
  )
  
  if ("structural_relational" %in% enabled_submodules) {
    log_message("=== M01 PART B: STRUCTURAL & RELATIONAL ANALYSIS ===")
    
    # Part B-1: Chromosome Distribution
    structural_params <- get_task_parameter(task_params, config, "structural_relational", default = list(), log_default = FALSE)
    chr_params <- get_task_parameter(structural_params, config, "chromosome_distribution", default = list(), log_default = FALSE)
    
    if (!is.null(chr_params$enabled) && chr_params$enabled) {
      log_message("M01: Executing chromosome distribution analysis")
      
      if (is.null(filtered_data)) {
        warning("M01: filtered_data is NULL, skipping chromosome distribution analysis")
      } else {
        # Create output directory (NEW: three-tier architecture)
        chr_output_dir <- file.path(output_dir, "M01_distribution", "structural_relational", "chromosome_distribution")
        dir.create(chr_output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Apply smart data_scoping inheritance
        chr_params$data_scoping <- inherit_data_scoping("chromosome_distribution", chr_params, task_params)
        
        # Call helper function (migrated from M05)
        chr_dist_results <- generate_chromosome_distribution_analysis(
          filtered_data = filtered_data,
          config = config,
          output_dir = chr_output_dir,
          task_name = paste0(task_name, "_chromosome"),
          task_params = chr_params
        )
        
        if (chr_dist_results$success) {
          log_message("M01: Chromosome distribution completed successfully")
          log_message(sprintf("  - Plots generated: %d", chr_dist_results$plots_generated))
          result$chromosome_distribution_results <- chr_dist_results
        } else {
          warning("M01: Chromosome distribution encountered errors")
        }
      }
    }
    
    # Part B-2: poiGS-Gene Correlation
    poigs_params <- get_task_parameter(structural_params, config, "poigs_gene_correlation", default = list(), log_default = FALSE)
    
    if (!is.null(poigs_params$enabled) && poigs_params$enabled) {
      log_message("M01: Executing poiGS-Gene correlation analysis")
      
      if (is.null(normalized_data)) {
        warning("M01: normalized_data is NULL, skipping poiGS correlation analysis")
      } else {
        # CRITICAL: Wrap Part B-2 in tryCatch to ensure independence
        tryCatch({
          # Create output directory (NEW: three-tier architecture)
          corr_output_dir <- file.path(output_dir, "M01_distribution", "structural_relational", "poigs_gene_correlation")
          dir.create(corr_output_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Apply smart data_scoping inheritance
          poigs_params$data_scoping <- inherit_data_scoping("poigs_gene_correlation", poigs_params, task_params)
        
        # Step 1: Apply M01-specific data_scoping
        log_message("M01: Applying M01-specific data_scoping for poiGS correlation")
        
        m01_scoping_params <- get_task_parameter(
          poigs_params,
          config,
          "data_scoping",
          default = list(
            target_region_types = c("CDS", "IGS"),
            target_var_types = c("snp", "INDEL")
          ),
          log_default = TRUE
        )
        
        m01_scoped_data <- apply_data_scoping(
          data = normalized_data,
          task_params = list(data_scoping = m01_scoping_params)
        )
        
        log_message(sprintf("  - M01 scoped data: %d records from %d species",
                           nrow(m01_scoped_data),
                           length(unique(m01_scoped_data$species))))
        
        # Step 2: Build poiGS mapping
        log_message("M01: Building poiGS mapping from IGS region statistics")
        
        poigs_mapping <- build_poigs_mapping(
          session_id = config$session_info$session_id,
          config = config
        )
        
        # Step 3: Build poiGS frequencies
        log_message("M01: Building poiGS frequencies from scoped data")
        
        poigs_result <- build_poigs_frequencies(
          normalized_data = m01_scoped_data,
          poigs_mapping = poigs_mapping,
          config = config
        )
        
        poigs_freq_data <- poigs_result$frequencies
        
        if (nrow(poigs_freq_data) > 0) {
          # Step 4: Extract gene frequencies
          log_message("M01: Extracting gene frequencies from scoped data")
          
          # Calculate gene length by deduplicating exons
          # S8 FIX introduced a bug: sum(region_length) triple-counted lengths because
          # normalized_data contains multiple rows per exon (snp, INDEL, complex)
          # Correct logic: Sum unique exon lengths
          gene_freq_data <- m01_scoped_data %>%
            dplyr::filter(region_type == "CDS", !is.na(gene)) %>%
            dplyr::group_by(species, gene) %>%
            dplyr::summarise(
              # Step 1: Aggregate raw counts across all exons of the gene
              total_variant_count = sum(variant_count, na.rm = TRUE),
              # Deduplicate regions to avoid triple-counting length (snp+indel+complex rows)
              total_region_length = sum(region_length[!duplicated(region_name)], na.rm = TRUE),
              # Step 2: Recalculate frequency from aggregated totals (correct density)
              frequency_per_kb = dplyr::if_else(
                total_region_length > 0,
                (total_variant_count / total_region_length) * 1000,
                0
              ),
              .groups = "drop"
            ) %>%
            dplyr::select(species, gene, frequency_per_kb)
          
          log_message(sprintf("  - Extracted %d gene frequency records (density uses deduplicated gene length)", nrow(gene_freq_data)))
          
          # Step 5: Merge poiGS and gene frequencies
          log_message("M01: Merging poiGS and gene frequency data")
          
          merged_data <- merge_poigs_gene_frequencies(
            poigs_freq = poigs_freq_data,
            gene_freq = gene_freq_data,
            poigs_mapping = poigs_mapping,
            keep_all_zeros = TRUE
          )
          
          if (nrow(merged_data) > 0) {
            # Save merged data
            merged_data_file <- file.path(corr_output_dir, "M01_poigs_gene_merged_data.csv")
            readr::write_csv(merged_data, merged_data_file)
            log_message(sprintf("  - Saved merged data: %s", merged_data_file))
            
            # Step 6: Calculate correlation (filter [0,0] pairs)
            log_message("M01: Calculating within-species correlation")
            
            correlation_data <- merged_data %>%
              dplyr::filter(zero_pattern == "neither_zero")
            
            if (nrow(correlation_data) > 0) {
              min_pairs_per_species <- get_task_parameter(poigs_params, config, "min_pairs_per_species", default = 5, log_default = TRUE)
              
              # Calculate correlations on neither_zero data
              correlation_result <- calculate_within_species_correlation(
                merged_data = correlation_data,
                min_pairs_per_species = min_pairs_per_species
              )
              
              # Apply FDR correction (aligned with the M05 option A-2 path)
              log_message("M01: Applying FDR correction to correlation p-values")
              correlation_result$correlation_results <- apply_fdr_correction(
                correlation_result$correlation_results,
                fdr_method = "BH"
              )
              
              # Save additional correlation outputs (aligned with M05)
              log_message("M01: Saving correlation analysis results")
              
              # Save correlation results CSV
              correlation_results_file <- file.path(corr_output_dir, "M01_within_species_correlation.csv")
              readr::write_csv(correlation_result$correlation_results, correlation_results_file)
              log_message(sprintf("  - Saved: %s", basename(correlation_results_file)))
              
              # Save zero pattern summary
              zero_pattern_summary <- merged_data %>%
                dplyr::group_by(species, zero_pattern) %>%
                dplyr::summarise(count = dplyr::n(), .groups = "drop")
              
              zero_summary_file <- file.path(corr_output_dir, "M01_zero_pattern_summary.csv")
              readr::write_csv(zero_pattern_summary, zero_summary_file)
              log_message(sprintf("  - Saved: %s", basename(zero_summary_file)))
              
              # Save summary statistics
              summary_stats_file <- file.path(corr_output_dir, "M01_correlation_summary_stats.txt")
              sink(summary_stats_file)
              cat("=== M01 poiGS-Gene Correlation Analysis Summary ===\n\n")
              cat(sprintf("Total species analyzed: %d\n", correlation_result$summary_stats$total_species))
              cat(sprintf("Species with sufficient data: %d\n", correlation_result$summary_stats$species_with_sufficient_data))
              cat(sprintf("Species excluded (insufficient pairs): %d\n", correlation_result$summary_stats$species_excluded))
              cat(sprintf("Total poiGS-Gene pairs: %d\n", correlation_result$summary_stats$total_pairs))
              cat(sprintf("Mean pairs per species: %.1f\n", correlation_result$summary_stats$mean_pairs_per_species))
              cat(sprintf("\nMean Pearson correlation: %.3f\n", correlation_result$summary_stats$mean_pearson_r))
              cat(sprintf("Mean Spearman correlation: %.3f\n", correlation_result$summary_stats$mean_spearman_r))
              cat(sprintf("\nSignificant correlations (p < 0.05):\n"))
              cat(sprintf("  - Pearson: %d species\n", correlation_result$summary_stats$significant_pearson_count))
              cat(sprintf("  - Spearman: %d species\n", correlation_result$summary_stats$significant_spearman_count))
              sink()
              log_message(sprintf("  - Saved: %s", basename(summary_stats_file)))
              
              # Step 7: Generate visualizations
              log_message("M01: Generating poiGS-Gene correlation visualizations")
              
              # Heatmap (FIXED: correct filename to match content)
              heatmap_file <- file.path(corr_output_dir, "M01_species_correlation_heatmap.pdf")
              tryCatch({
                generate_correlation_heatmap(
                  correlation_result$correlation_results,  # correct field name
                  heatmap_file,
                  config
                )
                log_message(sprintf("  - Saved heatmap: %s", basename(heatmap_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate correlation heatmap: %s", e$message))
              })
              
              # Regression plot
              regression_file <- file.path(corr_output_dir, "M01_poigs_gene_regression_plot.pdf")
              tryCatch({
                generate_poigs_gene_regression_plot(
                  merged_data = correlation_data,
                  correlation_results = correlation_result$correlation_results,
                  output_file = regression_file,
                  config = config
                )
                log_message(sprintf("  - Saved regression plot: %s", basename(regression_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate regression plot: %s", e$message))
              })
              
              # Zero pattern distribution plot (FIXED: apply species ordering)
              zero_plot_file <- file.path(corr_output_dir, "M01_zero_distribution.pdf")
              tryCatch({
                # Apply species ordering (use species_order from module initialization)
                zero_pattern_summary_ordered <- apply_species_ordering(zero_pattern_summary, species_order)
                
                p_zero <- zero_pattern_summary_ordered %>%
                  ggplot2::ggplot(ggplot2::aes(x = species, y = count, fill = zero_pattern)) +
                  ggplot2::geom_bar(stat = "identity", position = "fill") +
                  ggplot2::scale_fill_manual(
                    values = c(
                      "neither_zero" = "#00A087FF",
                      "both_zero" = "#DC0000FF",
                      "poiGS_only" = "#4DBBD5FF",
                      "gene_only" = "#F39B7FFF"
                    ),
                    labels = c(
                      "neither_zero" = "Both non-zero",
                      "both_zero" = "Both zero",
                      "poiGS_only" = "poiGS only",
                      "gene_only" = "Gene only"
                    )
                  ) +
                  ggplot2::labs(
                    title = "M01: Zero Value Pattern Distribution (Percentage)",
                    x = "Species",
                    y = "Percentage (%)",
                    fill = "Zero Pattern"
                  ) +
                  ggplot2::scale_y_continuous(labels = scales::percent_format()) +
                  ggplot2::theme_classic() +
                  ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                    legend.position = "bottom"
                  )
                
                ggplot2::ggsave(zero_plot_file, p_zero, width = 12, height = 6, device = "pdf")
                log_message(sprintf("  - Saved: %s", basename(zero_plot_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate zero pattern plot: %s", e$message))
              })
              
              # STRATEGY CHANGE: Generate 3 INDEPENDENT histogram density plots
              # Each plot will have its own PDF file for complete independence from facet constraints
              tryCatch({
                # Prepare data for 3 independent plots
                hist_data_neither <- merged_data %>%
                  dplyr::filter(zero_pattern == "neither_zero") %>%
                  tidyr::pivot_longer(
                    cols = c(poiGS_freq, gene_freq_avg),
                    names_to = "variable",
                    values_to = "frequency"
                  ) %>%
                  dplyr::mutate(
                    variable = dplyr::recode(variable,
                      "poiGS_freq" = "poiGS Frequency",
                      "gene_freq_avg" = "Gene Frequency (Avg)"
                    )
                  )
                
                hist_data_gene_only <- merged_data %>%
                  dplyr::filter(zero_pattern == "gene_only") %>%
                  dplyr::select(species, gene_freq_avg) %>%
                  dplyr::rename(frequency = gene_freq_avg) %>%
                  dplyr::mutate(variable = "Gene Frequency (Avg)")
                
                hist_data_poigs_only <- merged_data %>%
                  dplyr::filter(zero_pattern == "poiGS_only") %>%
                  dplyr::select(species, poiGS_freq) %>%
                  dplyr::rename(frequency = poiGS_freq) %>%
                  dplyr::mutate(variable = "poiGS Frequency")
                
                # Calculate sample sizes
                n_neither <- nrow(merged_data %>% dplyr::filter(zero_pattern == "neither_zero"))
                n_gene_only <- nrow(merged_data %>% dplyr::filter(zero_pattern == "gene_only"))
                n_poigs_only <- nrow(merged_data %>% dplyr::filter(zero_pattern == "poiGS_only"))
                
                # Plot 1: Neither Zero (Both non-zero, used for correlation)
                # CRITICAL FIX: Density must be grouped by 'fill' aesthetic to match histogram
                # Using after_stat(count) WITHOUT group causes density to use ALL data
                p1 <- ggplot2::ggplot(hist_data_neither, ggplot2::aes(x = frequency, fill = variable)) +
                  ggplot2::geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
                  # FIXED: Add 'group = variable' to ensure density calculated per variable
                  ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count), group = variable), 
                                       alpha = 0.3) +
                  ggplot2::annotate("text", x = Inf, y = Inf, 
                                   label = sprintf("N = %d", n_neither),
                                   hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
                  ggplot2::scale_fill_manual(values = c(
                    "poiGS Frequency" = "#E64B35FF",
                    "Gene Frequency (Avg)" = "#4DBBD5FF"
                  )) +
                  ggplot2::labs(
                    title = "Neither-zero pairs (used for correlation)",
                    x = "Frequency (per kb)",
                    y = "Count",
                    fill = "Variable"
                  ) +
                  ggplot2::theme_classic() +
                  ggplot2::theme(legend.position = "bottom")
                
                # Plot 2: Gene Only (poiGS = 0)
                # CRITICAL FIX: Manually scale density to histogram bin width
                # after_stat(count) still uses density-based calculation, needs explicit scaling
                bin_width_2 <- diff(range(hist_data_gene_only$frequency, na.rm = TRUE)) / 30
                p2 <- ggplot2::ggplot(hist_data_gene_only, ggplot2::aes(x = frequency)) +
                  ggplot2::geom_histogram(bins = 30, alpha = 0.6, fill = "#4DBBD5FF") +
                  # Scale density to match histogram: density * N * bin_width
                  ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density) * n_gene_only * bin_width_2), 
                                       alpha = 0.3, fill = "#4DBBD5FF") +
                  ggplot2::annotate("text", x = Inf, y = Inf, 
                                   label = sprintf("N = %d", n_gene_only),
                                   hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
                  ggplot2::labs(
                    title = "Gene Only (poiGS = 0)",
                    x = "Gene Frequency (per kb)",
                    y = "Count"
                  ) +
                  ggplot2::theme_classic()
                
                # Plot 3: poiGS Only (Gene = 0)
                # CRITICAL FIX: Manually scale density to histogram bin width
                bin_width_3 <- diff(range(hist_data_poigs_only$frequency, na.rm = TRUE)) / 30
                p3 <- ggplot2::ggplot(hist_data_poigs_only, ggplot2::aes(x = frequency)) +
                  ggplot2::geom_histogram(bins = 30, alpha = 0.6, fill = "#E64B35FF") +
                  # Scale density to match histogram: density * N * bin_width
                  ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(density) * n_poigs_only * bin_width_3), 
                                       alpha = 0.3, fill = "#E64B35FF") +
                  ggplot2::annotate("text", x = Inf, y = Inf, 
                                   label = sprintf("N = %d", n_poigs_only),
                                   hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
                  ggplot2::labs(
                    title = "poiGS Only (Gene = 0)",
                    x = "poiGS Frequency (per kb)",
                    y = "Count"
                  ) +
                  ggplot2::theme_classic()
                
                # Save 3 independent PDF files
                ggplot2::ggsave(file.path(corr_output_dir, "M01_histogram_density_1_neither_zero.pdf"), 
                               p1, width = 10, height = 6, device = "pdf")
                log_message("  - Saved: M01_histogram_density_1_neither_zero.pdf")
                
                ggplot2::ggsave(file.path(corr_output_dir, "M01_histogram_density_2_gene_only.pdf"), 
                               p2, width = 10, height = 6, device = "pdf")
                log_message("  - Saved: M01_histogram_density_2_gene_only.pdf")
                
                ggplot2::ggsave(file.path(corr_output_dir, "M01_histogram_density_3_poigs_only.pdf"), 
                               p3, width = 10, height = 6, device = "pdf")
                log_message("  - Saved: M01_histogram_density_3_poigs_only.pdf")
              }, error = function(e) {
                warning(sprintf("Failed to generate histogram density plots: %s", e$message))
              })
              
              result$poigs_gene_correlation_results <- list(
                merged_data = merged_data,
                correlation_results = correlation_result,
                output_dir = corr_output_dir  # correct variable name
              )
              
              log_message("M01: poiGS-Gene correlation completed successfully")
            } else {
              warning("M01: No valid pairs for correlation (all zeros filtered out)")
            }
          } else {
            warning("M01: No poiGS-Gene pairs could be matched")
          }
        } else {
          warning("M01: No poiGS data could be built")
        }
        }, error = function(e) {
          # CRITICAL ERROR HANDLER: Isolate Part B-2 failures
          error_msg <- sprintf("M01 Part B-2 (poiGS-Gene Correlation) encountered an error: %s", e$message)
          log_message(error_msg, level = "error")
          warning(error_msg)
          # Store error for diagnostic purposes
          result$poigs_gene_correlation_error <- e$message
        })
      }
    }
    
    # Part B-3: Chromosome Hotspot Map (MIGRATED TO M05)
    # Chromosome hotspot maps are produced by the M05_hotspot_ideogram module
    # for better architectural separation and clearer dependency management.
    # 
    # Migration Notes:
    # - M05 is now a standalone module that depends on M03 and M04 outputs
    # - M05 should be executed AFTER M03 and M04 in the task sequence
    # - Configuration: Add M05_hotspot_ideogram task to analysis_tasks in YAML
    # 
    # The chromosome_hotspot_map parameter in M01 configuration is now deprecated.
    # Please use M05_hotspot_ideogram instead.
    
    hotspot_params <- get_task_parameter(structural_params, config, "chromosome_hotspot_map", default = list(), log_default = FALSE)
    
    if (!is.null(hotspot_params$enabled) && hotspot_params$enabled) {
      warning_msg <- paste(
        "M01 Part B-3 (Chromosome Hotspot Map) has been migrated to M05_hotspot_ideogram module.",
        "Please update your configuration:",
        "  1. Disable chromosome_hotspot_map in M01 configuration",
        "  2. Add M05_hotspot_ideogram task to analysis_tasks",
        "  3. Ensure M05 is executed AFTER M03 and M04",
        "See full_refactor_test.yml for the complete configuration example.",
        sep = "\n"
      )
      log_message(warning_msg, level = "warning")
      warning(warning_msg)
    }
  }
  
  # ================================================================
  # Part C: M01-3 Sequence Variant Characterization (NEW)
  # ================================================================
  
  if ("sequence_variants" %in% enabled_submodules) {
    log_message("=== M01 PART C: SEQUENCE VARIANT CHARACTERIZATION ===")
    
    # Extract enabled sub-modules for sequence variants
    seq_params <- get_task_parameter(task_params, config, "sequence_variants", default = list(), log_default = FALSE)
    enabled_submodules_seq_variants <- get_task_parameter(seq_params, config, "enabled_submodules", default = c(), log_default = FALSE)
    
    # Part C-1: CDS SNP Annotation
    cds_params <- get_task_parameter(seq_params, config, "cds_snp_annotation", default = list(), log_default = FALSE)
    
    if (!is.null(cds_params$enabled) && cds_params$enabled) {
      # CRITICAL: Wrap entire Part C-1 in tryCatch to ensure independence from Part C-2
      tryCatch({
        log_message("M01: Executing CDS SNP annotation analysis")
        
        # Apply smart data_scoping inheritance
        cds_params$data_scoping <- inherit_data_scoping("cds_snp_annotation", cds_params, task_params)
        
        # Create output directory (NEW: three-tier architecture)
        cds_output_dir <- file.path(output_dir, "M01_distribution", "sequence_variants", "cds_snp_annotation")
        dir.create(cds_output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Get data file paths from session
      filtered_data_file <- file.path(
        dirname(output_dir),
        "processed_data", "P02_filtered", "variants_filtered.csv"
      )
      
      preprocessed_data_file <- file.path(
        dirname(output_dir),
        "processed_data", "P01_preprocessed", "annotated_variants_ira_only.csv"
      )
      
      if (!file.exists(filtered_data_file)) {
        warning(sprintf("M01: Filtered data file not found: %s", filtered_data_file))
      } else if (!file.exists(preprocessed_data_file)) {
        warning(sprintf("M01: Preprocessed data file not found: %s", preprocessed_data_file))
      } else {
        # Step 1: Extract CDS SNP positions (whitelist)
        log_message("M01: Extracting filtered CDS SNP positions...")
        
        cds_snp_whitelist <- extract_filtered_cds_snp_positions(filtered_data_file)
        log_message(sprintf("  - Extracted %d unique CDS SNP positions", nrow(cds_snp_whitelist)))
        
        if (nrow(cds_snp_whitelist) > 0) {
          # Step 2: Extract sample-level annotations
          log_message("M01: Extracting sample-level annotations from P01...")
          
          sample_annotations <- extract_sample_level_annotations(
            preprocessed_data_file,
            cds_snp_whitelist
          )
          log_message(sprintf("  - Extracted %d sample-level annotations", nrow(sample_annotations)))
          
          if (nrow(sample_annotations) > 0) {
            # Step 3: Parse variant types from V12
            log_message("M01: Parsing variant types from V12 column...")
            
            annotated_variants <- parse_variant_type_from_v12(sample_annotations)
            log_message(sprintf("  - Classified %d variants", nrow(annotated_variants)))
            
            if (nrow(annotated_variants) > 0) {
              # Step 4: Classify sites using Majority Rule
              log_message("M01: Classifying sites using Majority Rule...")
              
              majority_threshold <- get_task_parameter(cds_params, config, "majority_rule_threshold", default = 0.65, log_default = TRUE)
              
              site_classification <- classify_snp_sites_majority_rule(
                annotated_variants,
                majority_threshold
              )
              
              log_message(sprintf("  - Classified %d sites", nrow(site_classification)))
              
              # Step 5: Calculate S/N ratios
              log_message("M01: Calculating S/N ratios...")
              
              pseudocount <- get_task_parameter(cds_params, config, "pseudocount", default = 1, log_default = TRUE)
              min_sites_per_gene <- get_task_parameter(cds_params, config, "min_sites_per_gene", default = 5, log_default = TRUE)
              
              ratio_results <- calculate_sn_ratio(
                site_classification,
                pseudocount,
                min_sites_per_gene
              )
              
              # Save results
              species_stats_file <- file.path(cds_output_dir, "M01_cds_snp_species_stats.csv")
              readr::write_csv(ratio_results$species_stats, species_stats_file)
              
              gene_stats_file <- file.path(cds_output_dir, "M01_cds_snp_gene_stats.csv")
              readr::write_csv(ratio_results$gene_stats, gene_stats_file)
              
              log_message(sprintf("  - Saved species stats: %s", basename(species_stats_file)))
              log_message(sprintf("  - Saved gene stats: %s", basename(gene_stats_file)))
              
            # Save sample-level annotations (intermediate data for diagnostic purposes)
            # CRITICAL: Save annotated_variants (with variant_class) instead of sample_annotations
            sample_annotations_file <- file.path(cds_output_dir, "M01_cds_snp_sample_annotations.csv")
            readr::write_csv(annotated_variants, sample_annotations_file)
            log_message(sprintf("  - Saved sample annotations: %s", basename(sample_annotations_file)))
              
              # Save site classification (intermediate data for diagnostic purposes)
              site_classification_file <- file.path(cds_output_dir, "M01_cds_snp_site_classification.csv")
              readr::write_csv(site_classification, site_classification_file)
              log_message(sprintf("  - Saved site classification: %s", basename(site_classification_file)))
              
              # Step 6: Generate visualizations
              log_message("M01: Generating CDS SNP annotation visualizations...")
              
              # Plot 1: S/N ratio barplot
              sn_barplot_file <- file.path(cds_output_dir, "M01_cds_sn_ratio_barplot.pdf")
              tryCatch({
                generate_sn_ratio_barplot(
                  ratio_results$species_stats,
                  sn_barplot_file,
                  config
                )
                log_message(sprintf("  - Saved: %s", basename(sn_barplot_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate S/N ratio barplot: %s", e$message))
              })
              
              # Plot 2: Scatter plot
              scatter_file <- file.path(cds_output_dir, "M01_syn_vs_nonsyn_scatter.pdf")
              tryCatch({
                generate_syn_vs_nonsyn_scatter(
                  ratio_results$species_stats,
                  scatter_file,
                  config
                )
                log_message(sprintf("  - Saved: %s", basename(scatter_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate scatter plot: %s", e$message))
              })
              
              # Plot 3: Site type composition
              composition_file <- file.path(cds_output_dir, "M01_site_type_composition.pdf")
              tryCatch({
                generate_site_type_composition_barplot(
                  site_classification,
                  composition_file,
                  config
                )
                log_message(sprintf("  - Saved: %s", basename(composition_file)))
              }, error = function(e) {
                warning(sprintf("Failed to generate composition barplot: %s", e$message))
              })
              
              result$cds_snp_annotation_results <- list(
                species_stats = ratio_results$species_stats,
                gene_stats = ratio_results$gene_stats,
                output_dir = cds_output_dir
              )
              
              log_message("M01: CDS SNP annotation completed successfully")
            }
          }
        }
      }
      }, error = function(e) {
        # CRITICAL ERROR HANDLER: Isolate Part C-1 failures to prevent blocking Part C-2
        error_msg <- sprintf("M01 Part C-1 (CDS SNP Annotation) encountered an error: %s", e$message)
        log_message(error_msg, level = "error")
        warning(error_msg)
        # Store error but continue to Part C-2
        result$cds_snp_annotation_error <- e$message
      })
    }
    
    # Part C-2: IGS Variant Features
    igs_params <- get_task_parameter(seq_params, config, "igs_variant_features", default = list(), log_default = FALSE)
    
    if (!is.null(igs_params$enabled) && igs_params$enabled) {
      # CRITICAL: Wrap entire Part C-2 in tryCatch to ensure independence from Part C-1
      tryCatch({
        log_message("M01: Executing IGS variant features analysis")
        
        # Apply smart data_scoping inheritance
        igs_params$data_scoping <- inherit_data_scoping("igs_variant_features", igs_params, task_params)
        
        # Create output directory (NEW: three-tier architecture)
        igs_output_dir <- file.path(output_dir, "M01_distribution", "sequence_variants", "igs_variant_features")
        dir.create(igs_output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Get task parameters
      pseudocount <- get_task_parameter(igs_params, config, "pseudocount", default = 1, log_default = TRUE)
      min_indel_length <- get_task_parameter(igs_params, config, "min_indel_length", default = 2, log_default = TRUE)
      min_species_for_sharing <- get_task_parameter(igs_params, config, "min_species_for_sharing", default = 3, log_default = TRUE)
      
      # Define data file paths
      filtered_data_file <- file.path(
        dirname(output_dir),
        "processed_data", "P02_filtered", "variants_filtered.csv"
      )
      
      preprocessed_data_file <- file.path(
        dirname(output_dir),
        "processed_data", "P01_preprocessed", "annotated_variants_ira_only.csv"
      )
      
      if (!file.exists(filtered_data_file)) {
        warning(sprintf("M01: Filtered data file not found: %s", filtered_data_file))
      } else {
        # Part 1: SNP/INDEL Ratio Analysis
        log_message("M01: Part 1 - Analyzing SNP/INDEL ratio for IGS regions")
        
        snp_indel_stats <- calculate_igs_snp_indel_ratio(filtered_data_file, pseudocount)
        
        if (nrow(snp_indel_stats) > 0) {
          # Save results
          snp_indel_file <- file.path(igs_output_dir, "M01_igs_snp_indel_ratio.csv")
          readr::write_csv(snp_indel_stats, snp_indel_file)
          log_message(sprintf("  - Saved SNP/INDEL ratio statistics: %s", basename(snp_indel_file)))
          
          # Generate visualizations
          log_message("  - Generating SNP/INDEL ratio visualizations...")
          
          # Plot 1: Count barplot
          count_barplot_file <- file.path(igs_output_dir, "M01_igs_snp_indel_count_barplot.pdf")
          tryCatch({
            generate_snp_indel_count_barplot(snp_indel_stats, count_barplot_file, config)
            log_message(sprintf("    Saved: %s", basename(count_barplot_file)))
          }, error = function(e) {
            warning(sprintf("Failed to generate SNP/INDEL count barplot: %s", e$message))
          })
          
          # Plot 2: Ratio lineplot
          ratio_lineplot_file <- file.path(igs_output_dir, "M01_igs_snp_indel_ratio_lineplot.pdf")
          tryCatch({
            generate_snp_indel_ratio_lineplot(snp_indel_stats, ratio_lineplot_file, config)
            log_message(sprintf("    Saved: %s", basename(ratio_lineplot_file)))
          }, error = function(e) {
            warning(sprintf("Failed to generate SNP/INDEL ratio lineplot: %s", e$message))
          })
          
          log_message("M01: Part 1 completed - SNP/INDEL ratio analysis")
        }
        
        # Part 2-4: INDEL Features Analysis (requires preprocessed data)
        if (!file.exists(preprocessed_data_file)) {
          warning(sprintf("M01: Preprocessed data file not found: %s", preprocessed_data_file))
        } else {
          # Extract IGS INDEL positions
          igs_indel_whitelist <- extract_filtered_igs_indel_positions(filtered_data_file)
          log_message(sprintf("  - Extracted %d unique IGS INDEL positions", nrow(igs_indel_whitelist)))
          
          if (nrow(igs_indel_whitelist) > 0) {
            # Extract INDEL length information
            indel_data <- extract_indel_length_info(preprocessed_data_file, igs_indel_whitelist)
            log_message(sprintf("  - Extracted %d INDEL records from samples", nrow(indel_data)))
            
            if (nrow(indel_data) > 0) {
              # Part 2: INDEL Length Features
              log_message("M01: Part 2 - Analyzing INDEL length features")
              
              indel_stats <- calculate_indel_length_stats(indel_data, min_indel_length)
              
              # Save results
              indel_data_file <- file.path(igs_output_dir, "M01_igs_indel_sample_data.csv")
              readr::write_csv(indel_data, indel_data_file)
              log_message(sprintf("  - Saved INDEL sample data: %s", basename(indel_data_file)))
              
              # Save species stats (all INDELs) for diagnostic purposes
              species_all_file <- file.path(igs_output_dir, "M01_igs_indel_species_stats_all.csv")
              readr::write_csv(indel_stats$species_stats_all, species_all_file)
              log_message(sprintf("  - Saved species stats (all INDELs): %s", basename(species_all_file)))
              
              if (nrow(indel_stats$species_stats_filtered) > 0) {
                species_filtered_file <- file.path(igs_output_dir, "M01_igs_indel_species_stats_filtered.csv")
                readr::write_csv(indel_stats$species_stats_filtered, species_filtered_file)
                log_message(sprintf("  - Saved species stats (filtered): %s", basename(species_filtered_file)))
              }
              
              # Save position-level statistics (ADDED: align with M05)
              if (!is.null(indel_stats$position_stats) && nrow(indel_stats$position_stats) > 0) {
                position_stats_file <- file.path(igs_output_dir, "M01_igs_indel_position_stats.csv")
                readr::write_csv(indel_stats$position_stats, position_stats_file)
                log_message(sprintf("  - Saved position stats: %s", basename(position_stats_file)))
              }
              
              if (nrow(indel_stats$species_stats_filtered) > 0) {
                
                # Generate visualizations
                log_message("  - Generating INDEL length visualizations...")
                
                # Histogram (filtered) - only if filtered data exists
                if (nrow(indel_stats$species_stats_filtered) > 0) {
                  hist_filtered_file <- file.path(igs_output_dir, "M01_igs_indel_length_histogram_filtered.pdf")
                  tryCatch({
                    generate_indel_length_histogram_filtered(indel_data, hist_filtered_file, config, min_indel_length)
                    log_message(sprintf("    Saved: %s", basename(hist_filtered_file)))
                  }, error = function(e) {
                    warning(sprintf("Failed to generate histogram (filtered): %s", e$message))
                  })
                }
                
                # Boxplot comparison
                if (nrow(indel_stats$species_stats_all) > 0) {
                  boxplot_comp_file <- file.path(igs_output_dir, "M01_igs_indel_length_boxplot_comparison.pdf")
                  tryCatch({
                    generate_indel_length_boxplot_comparison(
                      indel_stats$species_stats_all,
                      indel_stats$species_stats_filtered,
                      boxplot_comp_file,
                      config,
                      min_indel_length
                    )
                    log_message(sprintf("    Saved: %s", basename(boxplot_comp_file)))
                  }, error = function(e) {
                    warning(sprintf("Failed to generate boxplot comparison: %s", e$message))
                  })
                }
                
                # Summary statistics plot
                summary_comp_file <- file.path(igs_output_dir, "M01_igs_indel_summary_comparison.pdf")
                tryCatch({
                  generate_indel_summary_comparison_plot(
                    indel_stats$species_stats_filtered,
                    summary_comp_file,
                    config,
                    min_indel_length
                  )
                  log_message(sprintf("    Saved: %s", basename(summary_comp_file)))
                }, error = function(e) {
                  warning(sprintf("Failed to generate summary statistics plot: %s", e$message))
                })
                
                log_message("M01: Part 2 completed - INDEL length features analysis")
              }
              
              # Part 3: Shared INDEL Sequences
              log_message("M01: Part 3 - Identifying shared INDEL sequences")
              
              shared_sequences <- identify_shared_indel_sequences(indel_data, min_species_for_sharing)
              
              if (nrow(shared_sequences) > 0) {
                # Save results
                shared_seqs_file <- file.path(igs_output_dir, "M01_shared_indel_sequences.csv")
                readr::write_csv(shared_sequences, shared_seqs_file)
                log_message(sprintf("  - Identified %d shared INDEL sequences", nrow(shared_sequences)))
                
                # Generate UpSet plot
                log_message("  - Generating UpSet plot for shared sequences...")
                
                filtered_indel_data <- indel_data %>%
                  dplyr::semi_join(shared_sequences, by = "indel_sequence")
                
                upset_matrix <- prepare_upset_matrix_by_sequence(filtered_indel_data)
                
                if (nrow(upset_matrix) > 0) {
                  upset_file <- file.path(igs_output_dir, "M01_shared_indel_upset.pdf")
                  tryCatch({
                    generate_shared_indel_upset_plot(upset_matrix, upset_file, config)
                    log_message(sprintf("    Saved: %s", basename(upset_file)))
                  }, error = function(e) {
                    warning(sprintf("Failed to generate UpSet plot: %s", e$message))
                  })
                }
                
                log_message("M01: Part 3 completed - Shared INDEL sequences analysis")
                
                # Part 4: Species Clustering
                log_message("M01: Part 4 - Performing species clustering (Jaccard similarity)")
                
                if (exists("upset_matrix") && !is.null(upset_matrix) && nrow(upset_matrix) > 0) {
                  log_message("  - Calculating Jaccard distances from shared sequences...")
                  clustering_result <- cluster_species_by_shared_indels(upset_matrix)
                  
                  if (!is.null(clustering_result)) {
                    # Save Jaccard distance matrix
                    jaccard_dist_file <- file.path(igs_output_dir, "M01_jaccard_distance_matrix.csv")
                    utils::write.csv(clustering_result$dist_matrix, jaccard_dist_file, row.names = TRUE)
                    log_message(sprintf("  - Saved Jaccard distance matrix: %s", basename(jaccard_dist_file)))
                    
                    # Save pairwise sharing statistics
                    pairwise_stats_file <- file.path(igs_output_dir, "M01_species_pairwise_sharing.csv")
                    readr::write_csv(clustering_result$pairwise_stats, pairwise_stats_file)
                    log_message(sprintf("  - Saved pairwise sharing stats: %s", basename(pairwise_stats_file)))
                    
                    # Generate dendrogram
                    dendrogram_file <- file.path(igs_output_dir, "M01_species_clustering_dendrogram.pdf")
                    tryCatch({
                      generate_species_clustering_dendrogram(clustering_result$hc, dendrogram_file, config)
                      log_message(sprintf("  - Saved dendrogram: %s", basename(dendrogram_file)))
                    }, error = function(e) {
                      warning(sprintf("Failed to generate dendrogram: %s", e$message))
                    })
                    
                    log_message("M01: Part 4 completed - Jaccard-based species clustering")
                  }
                }
              }
              
              result$igs_variant_features_results <- list(
                snp_indel_stats = snp_indel_stats,
                indel_stats = indel_stats,
                shared_sequences = shared_sequences,
                output_dir = igs_output_dir
              )
              
              log_message("M01: IGS variant features completed successfully")
            }
          }
        }
      }
      }, error = function(e) {
        # CRITICAL ERROR HANDLER: Isolate Part C-2 failures
        error_msg <- sprintf("M01 Part C-2 (IGS Variant Features) encountered an error: %s", e$message)
        log_message(error_msg, level = "error")
        warning(error_msg)
        # Store error for diagnostic purposes
        result$igs_variant_features_error <- e$message
      })
    }
    
    # ============================================================
    # M01-3e: Genome Region Statistics (NEW)
    # ============================================================
    # CRITICAL FIX: Check genome_region_config$enabled directly instead of relying on enabled_submodules
    genome_region_config <- get_task_parameter(seq_params, config, "genome_region_stats", default = list(), log_default = FALSE)
    
    if (!is.null(genome_region_config$enabled) && genome_region_config$enabled) {
      
      tryCatch({
        log_message("M01: Executing genome region statistics analysis")
        
        # Apply smart data_scoping inheritance
        genome_region_config$data_scoping <- inherit_data_scoping("genome_region_stats", genome_region_config, task_params)
        
        grstat_output_dir <- file.path(output_dir, "M01_distribution", "sequence_variants", "genome_region_stats")
        dir.create(grstat_output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Part 1: S/N Ratio by Genome Region
        log_message("M01: Part 1 - Calculating S/N ratio by genome region")
        
        sn_region_results <- calculate_sn_ratio_by_genome_region(
          preprocessed_data_file = preprocessed_data_file,
          filtered_data_file = filtered_data_file,
          majority_rule_threshold = get_task_parameter(genome_region_config, config, "majority_rule_threshold"),
          pseudocount = get_task_parameter(genome_region_config, config, "pseudocount")
        )
        
        sn_region_stats <- sn_region_results$region_stats
        sn_site_classification <- sn_region_results$site_classification
        
        if (nrow(sn_region_stats) > 0) {
          # Save main data table with all statistics
          sn_file <- file.path(grstat_output_dir, "M01_genome_region_sn_ratio.csv")
          readr::write_csv(sn_region_stats, sn_file)
          log_message("  - Saved: M01_genome_region_sn_ratio.csv")
          
          # Save detailed site classification (for verification)
          sn_site_file <- file.path(grstat_output_dir, "M01_genome_region_cds_snp_site_classification.csv")
          readr::write_csv(sn_site_classification, sn_site_file)
          log_message("  - Saved: M01_genome_region_cds_snp_site_classification.csv")
          
          # Generate visualization 1: S/N Ratio Barplot with significance
          sn_plot_file <- file.path(grstat_output_dir, "M01_genome_region_sn_ratio_barplot.pdf")
          tryCatch({
            generate_genome_region_sn_ratio_barplot(sn_region_stats, sn_plot_file, config)
            log_message("  - Saved: M01_genome_region_sn_ratio_barplot.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate S/N ratio barplot: %s", e$message))
          })
          
          # Generate visualization 2: Site Type Composition
          composition_plot_file <- file.path(grstat_output_dir, "M01_genome_region_site_type_composition.pdf")
          tryCatch({
            generate_genome_region_site_type_composition(sn_region_stats, composition_plot_file, config)
            log_message("  - Saved: M01_genome_region_site_type_composition.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate site type composition plot: %s", e$message))
          })
          
          # Generate visualization 3: Syn vs Nonsyn Scatter
          scatter_plot_file <- file.path(grstat_output_dir, "M01_genome_region_syn_vs_nonsyn_scatter.pdf")
          tryCatch({
            generate_genome_region_syn_vs_nonsyn_scatter(sn_region_stats, scatter_plot_file, config)
            log_message("  - Saved: M01_genome_region_syn_vs_nonsyn_scatter.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate syn vs nonsyn scatter plot: %s", e$message))
          })
        } else {
          warning("M01: No S/N ratio data by genome region (all species may lack valid IR regions)")
        }
        
        # Part 2: SNP/INDEL Ratio by Genome Region
        log_message("M01: Part 2 - Calculating SNP/INDEL ratio by genome region")
        
        snp_indel_region_stats <- calculate_snp_indel_ratio_by_genome_region(
          filtered_data_file = filtered_data_file,
          pseudocount = get_task_parameter(genome_region_config, config, "pseudocount")
        )
        
        if (nrow(snp_indel_region_stats) > 0) {
          # Save main data table with all statistics
          snp_indel_file <- file.path(grstat_output_dir, "M01_genome_region_snp_indel_ratio.csv")
          readr::write_csv(snp_indel_region_stats, snp_indel_file)
          log_message("  - Saved: M01_genome_region_snp_indel_ratio.csv")
          
          # Save detailed counts for transparency
          snp_indel_detail_file <- file.path(grstat_output_dir, "M01_genome_region_variant_counts.csv")
          snp_indel_detail_data <- snp_indel_region_stats %>%
            dplyr::select(species, genome_region, snp, INDEL, total_sites, snp_indel_ratio)
          readr::write_csv(snp_indel_detail_data, snp_indel_detail_file)
          log_message("  - Saved: M01_genome_region_variant_counts.csv")
          
          # Generate visualization 1: SNP/INDEL Ratio Barplot (grouped by region)
          snp_indel_plot_file <- file.path(grstat_output_dir, "M01_genome_region_snp_indel_ratio_barplot.pdf")
          tryCatch({
            generate_genome_region_snp_indel_ratio_plot(snp_indel_region_stats, snp_indel_plot_file, config)
            log_message("  - Saved: M01_genome_region_snp_indel_ratio_barplot.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate SNP/INDEL ratio barplot: %s", e$message))
          })
          
          # Generate visualization 2: SNP/INDEL Ratio Lineplot (trends across species)
          lineplot_file <- file.path(grstat_output_dir, "M01_genome_region_snp_indel_ratio_lineplot.pdf")
          tryCatch({
            generate_genome_region_snp_indel_ratio_lineplot(snp_indel_region_stats, lineplot_file, config)
            log_message("  - Saved: M01_genome_region_snp_indel_ratio_lineplot.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate SNP/INDEL ratio lineplot: %s", e$message))
          })
          
          # Generate visualization 3: SNP/INDEL Count Stacked Barplot
          count_plot_file <- file.path(grstat_output_dir, "M01_genome_region_snp_indel_count_barplot.pdf")
          tryCatch({
            generate_genome_region_snp_indel_count_barplot(snp_indel_region_stats, count_plot_file, config)
            log_message("  - Saved: M01_genome_region_snp_indel_count_barplot.pdf")
          }, error = function(e) {
            warning(sprintf("Failed to generate SNP/INDEL count barplot: %s", e$message))
          })
        } else {
          warning("M01: No SNP/INDEL ratio data by genome region (all species may lack valid IR regions)")
        }
        
        log_message("M01: Genome region statistics completed successfully")
        
      }, error = function(e) {
        error_msg <- sprintf("M01 Part C-3 (Genome Region Statistics) encountered an error: %s", e$message)
        log_message(error_msg, level = "error")
        warning(error_msg)
        result$genome_region_stats_error <- e$message
      })
    }
  }
  
  log_message(sprintf("=== M01 DATA DISTRIBUTION ANALYSIS COMPLETED: %s ===", task_name))
  log_message(sprintf("Generated manifest with %d entries", length(manifest_entries)))
  
  return(result)
}
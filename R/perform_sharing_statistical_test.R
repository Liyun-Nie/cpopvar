#' Perform sharing statistical test for Two-Speed Evolution Model
#' 
#' @title Perform sharing statistical test for Two-Speed Evolution Model
#' @description Executes stepwise statistical tests to evaluate the Two-Speed Evolution Model.
#' Supports both CDS-based (gene) and IGS-based (poiGS_ID) feature analysis.
#' @param hotspot_analysis_results Results from hotspot analysis containing feature sharing data
#' @param config Configuration list containing analysis parameters
#' @param task_params Task-specific parameters including alpha level
#' @param feature_col Name of the feature column used in sharing data 
#'   (default: "gene" for M03, use "poiGS_ID" for M05)
#' @return Statistical test results with significance values
#' @importFrom magrittr %>%
#' @importFrom dplyr rowwise
#' @export
perform_sharing_statistical_test <- function(hotspot_analysis_results,
                                            config,
                                            task_params = NULL,
                                            feature_col = "gene") {
  
  tryCatch({
    log_message("M03 Step 3d: Performing STEPWISE statistical test for Two-Speed Evolution Model")
    
    # Get parameters (simplified for stepwise approach)
    alpha_level <- get_task_parameter(task_params, config, "alpha_level", 0.05, log_default = TRUE)
    
    gene_sharing <- hotspot_analysis_results$all_gene_sharing
    total_species <- hotspot_analysis_results$parameters$total_species
    
    # Validate input data
    if (is.null(gene_sharing) || nrow(gene_sharing) == 0) {
      log_message("No gene sharing data available for statistical testing", level = "warning")
      return(NULL)
    }
    
    if (is.null(total_species) || length(total_species) == 0 || total_species <= 0) {
      log_message("Invalid total_species parameter for statistical testing", level = "warning")
      return(NULL)
    }
    
    if (!"num_species" %in% names(gene_sharing)) {
      log_message("Missing 'num_species' column in gene sharing data", level = "warning")
      return(NULL)
    }
    
    # Create observed frequency distribution
    sharing_counts <- gene_sharing %>%
      dplyr::count(num_species, name = "observed_count")
    
    # Calculate expected frequencies under random model
    # Under null hypothesis, each gene has equal probability of being a hotspot in each species
    total_genes <- sum(sharing_counts$observed_count)
    
    # Estimate probability that a gene is a hotspot in any given species
    total_gene_species_pairs <- sum(gene_sharing$num_species)
    total_possible_pairs <- total_genes * total_species
    
    # Validate calculation parameters
    if (total_possible_pairs == 0) {
      log_message("Total possible gene-species pairs is zero - cannot calculate hotspot probability", level = "warning")
      return(NULL)
    }
    
    p_hotspot <- total_gene_species_pairs / total_possible_pairs
    
    # Validate probability
    if (is.na(p_hotspot) || is.infinite(p_hotspot) || p_hotspot <= 0) {
      log_message("Invalid hotspot probability calculated - cannot proceed with statistical test", level = "warning")
      return(NULL)
    }
    
    log_message(sprintf("Estimated hotspot probability per gene-species pair: %.4f", p_hotspot))
    
    # Calculate expected distribution using binomial model
    expected_dist <- data.frame(
      num_species = 1:total_species,
      expected_count = total_genes * dbinom(1:total_species, total_species, p_hotspot)
    )
    
    # Merge observed and expected
    test_data <- sharing_counts %>%
      dplyr::full_join(expected_dist, by = "num_species") %>%
      tidyr::replace_na(list(observed_count = 0)) %>%
      dplyr::filter(expected_count > 0)
    
    # Calculate residuals for all categories
    test_data <- test_data %>%
      dplyr::mutate(
        residual = observed_count - expected_count,
        std_residual = residual / sqrt(expected_count),
        chi_contribution = residual^2 / expected_count,
        abs_std_residual = abs(std_residual)
      )
    
    # =================================================================
    # STEP 1: PRIVATE HOTSPOT GENES ENRICHMENT TEST (num_species = 1)
    # =================================================================
    log_message("Step 1: Testing private hotspot gene enrichment (num_species = 1)")
    
    private_genes_data <- test_data %>% dplyr::filter(num_species == 1)
    
    if (nrow(private_genes_data) > 0) {
      private_observed <- private_genes_data$observed_count
      private_expected <- private_genes_data$expected_count
      private_std_residual <- private_genes_data$std_residual
      
      # Z-test for this specific category
      private_z_score <- private_std_residual
      private_p_value <- 2 * (1 - stats::pnorm(abs(private_z_score)))  # Two-tailed test
      
      private_test_result <- list(
        category = "Private Hotspot Genes (num_species = 1)",
        observed = private_observed,
        expected = private_expected,
        standardized_residual = private_std_residual,
        z_score = private_z_score,
        p_value = private_p_value,
        significant = private_p_value < alpha_level,
        chi_contribution = private_genes_data$chi_contribution,
        interpretation = ifelse(private_p_value < alpha_level,
                              sprintf("Private genes SIGNIFICANTLY ENRICHED (obs=%d vs exp=%.1f, z=%.2f, p=%.4f)", 
                                     private_observed, private_expected, private_z_score, private_p_value),
                              sprintf("Private genes show no significant enrichment (obs=%d vs exp=%.1f, z=%.2f, p=%.4f)", 
                                     private_observed, private_expected, private_z_score, private_p_value))
      )
      
      log_message(sprintf("Private genes test: %s", private_test_result$interpretation))
    } else {
      private_test_result <- list(category = "Private Hotspot Genes", interpretation = "No private genes found")
    }
    
    # =================================================================  
    # STEP 2: CORE SHARED GENES DISTRIBUTION TEST (num_species >= 3)
    # =================================================================
    log_message("Step 2: Testing core shared gene distribution (num_species >= 3)")
    
    core_genes_data <- test_data %>% dplyr::filter(num_species >= 3)
    
    if (nrow(core_genes_data) > 0) {
      # Test if core shared genes follow random distribution
      core_chi_stat <- sum(core_genes_data$chi_contribution)
      core_df <- nrow(core_genes_data) - 1
      core_p_value <- ifelse(core_df > 0, 1 - stats::pchisq(core_chi_stat, core_df), 1)
      
      # Check if most residuals are small (< 1.96)
      small_residuals <- sum(core_genes_data$abs_std_residual < 1.96)
      total_core_categories <- nrow(core_genes_data)
      
      core_test_result <- list(
        category = "Core Shared Genes (num_species >= 3)",
        n_categories = total_core_categories,
        chi_square = core_chi_stat,
        df = core_df,
        p_value = core_p_value,
        significant = core_p_value < alpha_level,
        small_residuals = small_residuals,
        pct_small_residuals = small_residuals / total_core_categories * 100,
        interpretation = ifelse(core_p_value < alpha_level,
                              sprintf("Core genes show SIGNIFICANT deviation from random (chi-square=%.2f, df=%d, p=%.4f)", 
                                     core_chi_stat, core_df, core_p_value),
                              sprintf("Core genes follow RANDOM distribution (chi-square=%.2f, df=%d, p=%.4f, %.0f%% small residuals)", 
                                     core_chi_stat, core_df, core_p_value, small_residuals / total_core_categories * 100))
      )
      
      log_message(sprintf("Core genes test: %s", core_test_result$interpretation))
    } else {
      core_test_result <- list(category = "Core Shared Genes", interpretation = "No core shared genes found")
    }
    
    # =================================================================
    # STEP 3: TWO-SPEED EVOLUTION MODEL SYNTHESIS  
    # =================================================================
    log_message("Step 3: Synthesizing Two-Speed Evolution Model")
    
    # Determine the evolutionary pattern
    private_enriched <- !is.null(private_test_result$significant) && private_test_result$significant
    core_random <- !is.null(core_test_result$significant) && !core_test_result$significant
    
    two_speed_model <- list(
      private_genes_enriched = private_enriched,
      core_genes_random = core_random,
      model_supported = private_enriched && core_random
    )
    
    if (two_speed_model$model_supported) {
      model_interpretation <- "TWO-SPEED EVOLUTION MODEL SUPPORTED: Strong non-random forces drive private gene formation, while core shared genes disperse randomly"
    } else if (private_enriched && !core_random) {
      model_interpretation <- "MODIFIED TWO-SPEED MODEL: Both private and core gene formation show non-random patterns"
    } else if (!private_enriched && core_random) {
      model_interpretation <- "WEAK EVOLUTIONARY SIGNAL: Neither private nor core genes show strong non-random patterns"
    } else {
      model_interpretation <- "UNIFORM NON-RANDOM MODEL: Both private and core genes show non-random patterns"
    }
    
    log_message(sprintf("Two-Speed Model: %s", model_interpretation))
    
    # =================================================================
    # STEP 4: GENE-BY-GENE HYPERGEOMETRIC DRIVER ANALYSIS
    # =================================================================
    log_message("Step 4: Performing gene-by-gene hypergeometric tests to identify driver genes")
    
    # Perform correct hypergeometric tests: test each gene individually to identify "driver genes"
    hypergeometric_results <- tryCatch({
      if (is.null(gene_sharing) || nrow(gene_sharing) == 0) {
        log_message("No gene sharing data available for gene-by-gene hypergeometric analysis", level = "warning")
        return(NULL)
      }
      
      # Step 1: Prepare data for feature-by-feature analysis
      # Get all hotspot events (feature-species pairs) from the original data
      all_hotspot_events <- gene_sharing %>%
        # Create one row for each feature-species occurrence
        dplyr::group_by(!!sym(feature_col)) %>%
        dplyr::reframe(
          species_count = num_species,
          # Each feature appears once per species it's found in
          hotspot_events = pmax(1, num_species)  # Ensure at least 1 event per feature
        ) %>%
        # Total hotspot events in the system
        dplyr::ungroup()
      
      total_hotspot_events <- sum(all_hotspot_events$hotspot_events)
      
      # Step 2: Define categories for analysis
      private_genes <- gene_sharing %>% dplyr::filter(num_species == 1)
      core_genes <- gene_sharing %>% dplyr::filter(num_species >= 3)
      
      log_message(sprintf("Gene-by-gene analysis setup: %d total hotspot events, %d private genes, %d core genes",
                         total_hotspot_events, nrow(private_genes), nrow(core_genes)))
      
      if (nrow(private_genes) == 0 && nrow(core_genes) == 0) {
        log_message("No genes in either private or core categories for gene-by-gene analysis", level = "warning")
        return(NULL)
      }
      
      # Step 3: Perform hypergeometric test for each gene in each category
      gene_level_results <- list()
      
      # Test private genes (num_species = 1)
      if (nrow(private_genes) > 0) {
        log_message(sprintf("Testing %d private genes for significant contribution", nrow(private_genes)))
        
        private_category_events <- sum(private_genes$num_species)  # Total events in private category
        
        private_gene_tests <- private_genes %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            # Hypergeometric test parameters:
            # Population (N): Total hotspot events across all genes and species
            # Success states in population (K): Events for this specific gene
            # Sample size (n): Total events in private category
            # Observed successes (x): Events for this gene in private category (always 1 for private)
            
            gene_events = num_species,  # Events for this gene (1 for private genes)
            population_size = total_hotspot_events,
            success_states_in_pop = gene_events,
            sample_size = private_category_events,
            observed_successes = 1,  # This gene appears once in private category
            
            # Expected number of times this gene should appear in private category under random model
            expected_in_sample = (gene_events / population_size) * sample_size,
            
            # Hypergeometric p-value (upper tail test)
            p_value = stats::phyper(observed_successes - 1, success_states_in_pop, 
                                   population_size - success_states_in_pop, sample_size, 
                                   lower.tail = FALSE),
            
            # Enrichment metrics
            enrichment_ratio = observed_successes / pmax(expected_in_sample, 0.001),
            significant = p_value < alpha_level,
            category = "Private"
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(!!sym(feature_col), category, gene_events, population_size, success_states_in_pop, sample_size,
                observed_successes, expected_in_sample, enrichment_ratio, p_value, significant)
        
        gene_level_results$private <- private_gene_tests
        
        # Summary for private genes
        n_significant_private <- sum(private_gene_tests$significant, na.rm = TRUE)
        log_message(sprintf("Private genes: %d/%d genes show significant enrichment (p < %.3f)",
                           n_significant_private, nrow(private_gene_tests), alpha_level))
      }
      
      # Test core shared genes (num_species >= 3)
      if (nrow(core_genes) > 0) {
        log_message(sprintf("Testing %d core shared genes for significant contribution", nrow(core_genes)))
        
        core_category_events <- sum(core_genes$num_species)  # Total events in core category
        
        core_gene_tests <- core_genes %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            # For core genes, they appear multiple times (num_species times)
            gene_events = num_species,
            population_size = total_hotspot_events,
            success_states_in_pop = gene_events,
            sample_size = core_category_events,
            observed_successes = num_species,  # This gene appears num_species times in core
            
            # Expected number of times this gene should appear in core category
            expected_in_sample = (gene_events / population_size) * sample_size,
            
            # Hypergeometric p-value (upper tail test)
            p_value = phyper(observed_successes - 1, success_states_in_pop,
                            population_size - success_states_in_pop, sample_size,
                            lower.tail = FALSE),
            
            # Enrichment metrics
            enrichment_ratio = observed_successes / pmax(expected_in_sample, 0.001),
            significant = p_value < alpha_level,
            category = "Core"
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(!!sym(feature_col), category, gene_events, population_size, success_states_in_pop, sample_size,
                observed_successes, expected_in_sample, enrichment_ratio, p_value, significant)
        
        gene_level_results$core <- core_gene_tests
        
        # Summary for core genes
        n_significant_core <- sum(core_gene_tests$significant, na.rm = TRUE)
        log_message(sprintf("Core genes: %d/%d genes show significant enrichment (p < %.3f)",
                           n_significant_core, nrow(core_gene_tests), alpha_level))
      }
      
      # Step 4: Combine results and create summary
      all_gene_tests <- dplyr::bind_rows(gene_level_results)
      
      if (nrow(all_gene_tests) > 0) {
        # Identify top driver genes in each category
        top_private_drivers <- if (!is.null(gene_level_results$private)) {
          gene_level_results$private %>%
            dplyr::filter(significant) %>%
            dplyr::arrange(p_value) %>%
            utils::head(5)
        } else { data.frame() }
        
        top_core_drivers <- if (!is.null(gene_level_results$core)) {
          gene_level_results$core %>%
            dplyr::filter(significant) %>%
            dplyr::arrange(p_value) %>%
            utils::head(5)
        } else { data.frame() }
        
        # Create summary results
        summary_results <- list(
          test_method = "gene_by_gene_hypergeometric",
          total_genes_tested = nrow(all_gene_tests),
          private_genes_tested = nrow(gene_level_results$private %||% data.frame()),
          core_genes_tested = nrow(gene_level_results$core %||% data.frame()),
          significant_private_drivers = nrow(top_private_drivers),
          significant_core_drivers = nrow(top_core_drivers),
          alpha_level = alpha_level,
          
          # Detailed results
          all_gene_results = all_gene_tests,
          private_drivers = top_private_drivers,
          core_drivers = top_core_drivers,
          
          # Method description
          method_description = "Gene-by-gene hypergeometric test identifies individual genes that contribute disproportionately to private or core shared categories"
        )
        
        log_message(sprintf("Gene-by-gene hypergeometric analysis completed: %d genes tested, %d significant private drivers, %d significant core drivers",
                           summary_results$total_genes_tested, summary_results$significant_private_drivers, summary_results$significant_core_drivers))
        
        summary_results
      } else {
        NULL
      }
    }, error = function(e) {
      log_message(sprintf("Failed to perform gene-by-gene hypergeometric tests: %s", e$message), level = "error")
      NULL
    })
    
    # =================================================================
    # STEP 5: BIPOLARITY INDEX CALCULATION
    # =================================================================
    log_message("Step 5: Calculating Bipolarity Index for two-pole gene differentiation")
    
    # Calculate bipolarity index based on chi-square test results
    bipolarity_results <- tryCatch({
      # Define poles: private genes (num_species = 1) and core shared genes (num_species >= 3)
      private_data <- test_data %>% dplyr::filter(num_species == 1)
      core_data <- test_data %>% dplyr::filter(num_species >= 3)
      
      if (nrow(private_data) == 0 || nrow(core_data) == 0) {
        log_message("Cannot calculate Bipolarity Index - insufficient data for both poles", level = "warning")
        NULL
      } else {
        # Step 1: Extract observed and expected counts for both poles
        private_observed <- private_data$observed_count[1]
        private_expected <- private_data$expected_count[1]
        core_observed <- sum(core_data$observed_count)
        core_expected <- sum(core_data$expected_count)
        
        # Step 2: Calculate observed index (exclude intermediate genes)
        total_bipolar_observed <- private_observed + core_observed
        observed_index <- private_observed / total_bipolar_observed
        
        # Step 3: Calculate expected index  
        total_bipolar_expected <- private_expected + core_expected
        expected_index <- private_expected / total_bipolar_expected
        
        # Step 4: Calculate bipolarity metrics
        bipolarity_deviation <- observed_index - expected_index
        bipolarity_ratio <- observed_index / expected_index
        bipolarity_strength <- abs(bipolarity_deviation) / expected_index  # Relative deviation
        
        # Interpret bipolarity
        bipolarity_interpretation <- if (observed_index > expected_index + 0.1) {
          sprintf("STRONG PRIVATE GENE BIAS: %.1f%% observed vs %.1f%% expected private genes", 
                 observed_index * 100, expected_index * 100)
        } else if (observed_index < expected_index - 0.1) {
          sprintf("STRONG CORE SHARED BIAS: %.1f%% observed vs %.1f%% expected private genes", 
                 observed_index * 100, expected_index * 100)
        } else {
          sprintf("BALANCED BIPOLARITY: %.1f%% observed vs %.1f%% expected private genes", 
                 observed_index * 100, expected_index * 100)
        }
        
        bipolarity_result <- list(
          private_genes = list(observed = private_observed, expected = private_expected),
          core_genes = list(observed = core_observed, expected = core_expected),
          observed_index = observed_index,
          expected_index = expected_index,
          bipolarity_deviation = bipolarity_deviation,
          bipolarity_ratio = bipolarity_ratio,
          bipolarity_strength = bipolarity_strength,
          interpretation = bipolarity_interpretation,
          calculation_details = sprintf(
            "Observed Index = %d / (%d + %d) = %.3f; Expected Index = %.1f / (%.1f + %.1f) = %.3f",
            private_observed, private_observed, core_observed, observed_index,
            private_expected, private_expected, core_expected, expected_index
          )
        )
        
        log_message(sprintf("Bipolarity Index calculated: %s", bipolarity_interpretation))
        log_message(sprintf("Bipolarity metrics - Deviation: %.3f, Ratio: %.3f, Strength: %.3f", 
                           bipolarity_deviation, bipolarity_ratio, bipolarity_strength))
        
        bipolarity_result
      }
    }, error = function(e) {
      log_message(sprintf("Failed to calculate Bipolarity Index: %s", e$message), level = "error")
      NULL
    })
    
    # Generate comprehensive report
    report_lines <- c(
      "=== STEPWISE HOTSPOT GENE SHARING ANALYSIS ===",
      "=== Two-Speed Evolution Model Testing ===",
      "",
      sprintf("Analysis Date: %s", Sys.time()),
      sprintf("Total Genes: %d", total_genes),
      sprintf("Total Species: %d", total_species),
      sprintf("Estimated Hotspot Probability: %.4f per gene-species pair", p_hotspot),
      "",
      "=== STEP 1: PRIVATE HOTSPOT GENE ENRICHMENT ===",
      private_test_result$interpretation,
      ifelse(!is.null(private_test_result$significant) && private_test_result$significant,
             "CONCLUSION: Private hotspot genes show SIGNIFICANT enrichment - stronger selection/constraint forces",
             "CONCLUSION: Private hotspot genes show no significant enrichment"),
      "",
      "=== STEP 2: CORE SHARED GENE DISTRIBUTION ===", 
      core_test_result$interpretation,
      ifelse(!is.null(core_test_result$significant) && !core_test_result$significant,
             "CONCLUSION: Core shared genes follow RANDOM distribution - neutral dispersal process",
             "CONCLUSION: Core shared genes show non-random distribution patterns"),
      "",
      "=== STEP 3: TWO-SPEED EVOLUTION MODEL SYNTHESIS ===",
      model_interpretation,
      ""
    )
    
    # Add evolutionary interpretation based on model support
    if (two_speed_model$model_supported) {
      report_lines <- c(report_lines,
        "EVOLUTIONARY INTERPRETATION:",
        "- FIRST SPEED: Strong evolutionary forces (selection/constraint) create species-specific hotspots.",
        "- SECOND SPEED: Once genes become 'shareable', further dispersal follows random/neutral processes.", 
        "- This suggests two distinct evolutionary mechanisms operating at different scales.",
        ""
      )
    } else {
      report_lines <- c(report_lines,
        "ALTERNATIVE EVOLUTIONARY PATTERN:",
        "- The data does not support the classic two-speed model.",
        "- May indicate uniform evolutionary processes or more complex patterns.",
        ""
      )
    }
    
    # === STEP 4: GENE-BY-GENE HYPERGEOMETRIC DRIVER ANALYSIS ===
    report_lines <- c(report_lines, "=== STEP 4: GENE-BY-GENE HYPERGEOMETRIC DRIVER ANALYSIS ===")
    
    if (!is.null(hypergeometric_results)) {
      report_lines <- c(report_lines,
        "GENE-BY-GENE HYPERGEOMETRIC TESTS:",
        sprintf("- Test Method: %s", hypergeometric_results$method_description %||% "Gene-by-gene hypergeometric analysis"),
        sprintf("- Total Genes Tested: %d", hypergeometric_results$total_genes_tested %||% 0),
        sprintf("- Private Genes Tested: %d", hypergeometric_results$private_genes_tested %||% 0),
        sprintf("- Core Genes Tested: %d", hypergeometric_results$core_genes_tested %||% 0),
        sprintf("- Alpha Level: %.3f", hypergeometric_results$alpha_level %||% 0.05),
        "",
        "DRIVER GENE IDENTIFICATION:",
        sprintf("- Significant Private Driver Genes: %d", hypergeometric_results$significant_private_drivers %||% 0),
        sprintf("- Significant Core Driver Genes: %d", hypergeometric_results$significant_core_drivers %||% 0),
        ""
      )
      
      if ((hypergeometric_results$significant_private_drivers %||% 0) > 0) {
        top_drivers <- hypergeometric_results$private_drivers
        feature_name <- toupper(gsub("_ID$", "", feature_col))
        driver_lines <- paste(sprintf("  - %s (p=%.4f, ratio=%.2f)", 
                              top_drivers[[feature_col]], top_drivers$p_value, top_drivers$enrichment_ratio))
        report_lines <- c(report_lines, sprintf("TOP PRIVATE DRIVER %sS:", feature_name), driver_lines, "")
      } else {
        feature_name <- tolower(gsub("_ID$", "", feature_col))
        report_lines <- c(report_lines, sprintf("- No significant private driver %ss identified.", feature_name), "")
      }
      
      if ((hypergeometric_results$significant_core_drivers %||% 0) > 0) {
        top_drivers <- hypergeometric_results$core_drivers
        feature_name <- toupper(gsub("_ID$", "", feature_col))
        driver_lines <- paste(sprintf("  - %s (p=%.4f, ratio=%.2f)", 
                              top_drivers[[feature_col]], top_drivers$p_value, top_drivers$enrichment_ratio))
        report_lines <- c(report_lines, sprintf("TOP CORE DRIVER %sS:", feature_name), driver_lines, "")
      } else {
        feature_name <- tolower(gsub("_ID$", "", feature_col))
        report_lines <- c(report_lines, sprintf("- No significant core driver %ss identified.", feature_name), "")
      }
      
      interpretation <- "BIOLOGICAL INTERPRETATION: No significant driver genes identified - suggests random distribution patterns."
      if ((hypergeometric_results$significant_private_drivers %||% 0) > 0 && (hypergeometric_results$significant_core_drivers %||% 0) > 0) {
        interpretation <- "BIOLOGICAL INTERPRETATION: Both private and core categories have significant driver genes."
      } else if ((hypergeometric_results$significant_private_drivers %||% 0) > 0) {
        interpretation <- "BIOLOGICAL INTERPRETATION: Only private category shows significant driver genes - supports species-specific evolution."
      } else if ((hypergeometric_results$significant_core_drivers %||% 0) > 0) {
        interpretation <- "BIOLOGICAL INTERPRETATION: Only core category shows significant driver genes - supports shared evolutionary pressures."
      }
      report_lines <- c(report_lines, interpretation, "")
      
    } else {
      report_lines <- c(report_lines, 
        "GENE-BY-GENE HYPERGEOMETRIC TESTS: Cannot be performed due to insufficient data.",
        "- Requires gene sharing data with both private and core shared categories.",
        ""
      )
    }
    
    # === STEP 5: BIPOLARITY INDEX ANALYSIS ===
    report_lines <- c(report_lines, "=== STEP 5: BIPOLARITY INDEX ANALYSIS ===")
    
    if (!is.null(bipolarity_results)) {
      report_lines <- c(report_lines,
        "BIPOLARITY INDEX CALCULATION:",
        sprintf("- %s", bipolarity_results$interpretation),
        sprintf("- Calculation Details: %s", bipolarity_results$calculation_details),
        sprintf("- Bipolarity Deviation: %.3f (%.1f%% points difference)", 
               bipolarity_results$bipolarity_deviation, 
               bipolarity_results$bipolarity_deviation * 100),
        sprintf("- Bipolarity Ratio: %.3f (%.1fx stronger than expected)", 
               bipolarity_results$bipolarity_ratio, bipolarity_results$bipolarity_ratio),
        sprintf("- Bipolarity Strength: %.3f (relative deviation from expectation)", 
               bipolarity_results$bipolarity_strength),
        ""
      )
      
      interpretation <- "BIOLOGICAL INTERPRETATION: Balanced bipolarity indicates neutral evolutionary processes."
      if (bipolarity_results$bipolarity_deviation > 0.1) {
        interpretation <- c(
          "BIOLOGICAL INTERPRETATION:",
          "- Strong systematic force pushes genes toward species-specific (private) direction.",
          "- Evidence of non-random evolutionary constraint favoring private gene formation.",
          "- Suggests selective pressure for genomic differentiation between species."
        )
      } else if (bipolarity_results$bipolarity_deviation < -0.1) {
        interpretation <- c(
          "BIOLOGICAL INTERPRETATION:",
          "- Strong systematic force pushes genes toward shared (core) direction.", 
          "- Evidence of non-random evolutionary constraint favoring core gene formation.",
          "- Suggests selective pressure for genomic conservation across species."
        )
      }
      report_lines <- c(report_lines, interpretation, "")
      
    } else {
      report_lines <- c(report_lines, 
        "BIPOLARITY INDEX: Cannot be calculated due to insufficient data.",
        "- Requires both private genes (num_species=1) and core genes (num_species>=3).",
        ""
      )
    }
    
    sharing_summary <- list(
      total_genes = total_genes,
      total_species = total_species,
      mean_sharing = mean(gene_sharing$num_species, na.rm = TRUE),
      median_sharing = stats::median(gene_sharing$num_species, na.rm = TRUE),
      max_sharing = max(gene_sharing$num_species, na.rm = TRUE),
      genes_in_multiple_species = sum(gene_sharing$num_species > 1),
      genes_in_all_species = sum(gene_sharing$num_species == total_species),
      hotspot_probability = p_hotspot
    )
    
    test_results <- list(
      private_genes = private_test_result,
      core_genes = core_test_result,
      two_speed_model = two_speed_model,
      hypergeometric = hypergeometric_results,
      bipolarity_index = bipolarity_results
    )
    
    results <- list(
      tests = test_results,
      comparison_table = test_data,
      sharing_summary = sharing_summary,
      report = report_lines,
      test_metadata = list(
        test_method = "stepwise_two_speed",
        alpha_level = alpha_level,
        analysis_date = Sys.time(),
        model_interpretation = model_interpretation
      )
    )
    
    # Debug logging to identify NULL components
    log_message(sprintf("Statistical test results validation - report: %s, comparison_table: %s, hypergeometric: %s",
                       !is.null(results$report), !is.null(results$comparison_table), 
                       !is.null(results$tests$hypergeometric)))
    
    if (!is.null(results$report)) {
      log_message(sprintf("Report contains %d lines", length(results$report)))
    }
    if (!is.null(results$comparison_table)) {
      log_message(sprintf("Comparison table contains %d rows", nrow(results$comparison_table)))
    }
    if (!is.null(results$tests$hypergeometric)) {
      log_message("Hypergeometric test results available")
    }
    
    log_message("Stepwise statistical testing completed successfully")
    return(results)
    
  }, error = function(e) {
    log_message(sprintf("Failed to perform sharing statistical test: %s", e$message), level = "error")
    return(NULL)
  })
}

#' Generate universal single factor comparison analysis
#' 
#' @title Generate universal single factor comparison analysis
#' @description An intelligent statistical engine that automatically selects the most appropriate
#' statistical approach for single factor analysis based on data characteristics 
#' (normality, homogeneity of variance). Implements both parametric (ANOVA + TukeyHSD) 
#' and non-parametric (Kruskal-Wallis + Dunn's test) pathways with complete diagnostic transparency.
#' @importFrom magrittr %>%
#' @param plot_data Pre-processed plot data (already cleaned, ordered, and merged with group_info by orchestrator)
#' @param factor_name Name of the factor to analyze (e.g., "region_type", "var_type", "Phylogeny")
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param output_dir Optional output directory for saving statistical report TXT file
#' @param generate_plots Whether to generate plot objects in addition to statistics
#' @return List containing plot, data, summary statistics, statistical results (including diagnostic info), statistical report text, and metadata
#' @export
generate_universal_single_factor_comparison <- function(plot_data,
                                                      factor_name,
                                                      config,
                                                      task_params = NULL,
                                                      output_dir = NULL,
                                                      generate_plots = TRUE) {
  
  tryCatch({
    log_message(sprintf("Generating intelligent single-factor comparison for factor: %s", factor_name))
    
    # Get visual standards for consistent styling
    visual_standards <- get_visual_standards()
    geom_defaults <- visual_standards$ggplot2
    
    # V19 FINAL FIX: The defensive checks and NA filtering from V16 are correct.
    # The root cause is further downstream. We will add more robust data handling throughout.
    
    # V16 ENHANCED DEFENSIVE CHECK: Validate that factor exists in the data before proceeding
    if (!factor_name %in% names(plot_data)) {
      log_message(sprintf("Skipping factor '%s': column not found in the provided data.", factor_name), level = "warning")
      # Return an empty/fallback result to allow pipeline to continue (fixed function call)
      return(list(
        plot = ggplot2::ggplot() + 
          ggplot2::ggtitle("Factor Analysis Skipped") + 
          ggplot2::labs(subtitle = sprintf("Factor '%s' not found in data", factor_name)) +
          ggplot2::theme_minimal(),
        data = plot_data,
        summary_statistics = data.frame(),
        summary_table_gt = NULL,
        species_summary_table_gt = NULL,
        statistical_results = list(error = sprintf("Factor '%s' not found in data", factor_name)),
        statistical_report_text = NULL,
        metadata = list(
          factor_analyzed = factor_name,
          error = sprintf("Factor '%s' not found in data", factor_name),
          n_observations = 0
        )
      ))
    }
    
    # V16 NA Handling: Filter out NA values for the specific factor being analyzed
    analysis_data <- plot_data %>% 
      dplyr::filter(!is.na(.data[[factor_name]]))
    
    if(nrow(analysis_data) == 0) {
      log_message(sprintf("Skipping factor '%s': no valid data after removing NAs.", factor_name), level = "warning")
      return(list(
        plot = ggplot2::ggplot() + 
          ggplot2::ggtitle("Factor Analysis Skipped") + 
          ggplot2::labs(subtitle = sprintf("No data for factor '%s'", factor_name)) +
          ggplot2::theme_minimal(),
        data = plot_data,
        summary_statistics = data.frame(),
        summary_table_gt = NULL,
        species_summary_table_gt = NULL,
        statistical_results = list(error = sprintf("No data for factor '%s'", factor_name)),
        statistical_report_text = NULL,
        metadata = list(
          factor_analyzed = factor_name,
          error = sprintf("No data for factor '%s'", factor_name),
          n_observations = 0
        )
      ))
    }
    
    # Get plot configuration - statistical_test no longer needed (intelligent auto-selection)
    plot_type <- get_task_parameter(task_params, config, "plot_type", "boxplot")
    significance_level <- get_task_parameter(task_params, config, "significance_level", 0.05)
    
    # Get alpha scheme for consistent transparency
    alpha_scheme <- get_alpha_scheme()
    
    # Use internal standard column names (Rule #9 - Data Boundary Principle)
    freq_col <- "frequency_per_kb"
    species_col <- "species"
    
    # Validate required columns (check against original data structure)
    required_cols <- c(freq_col, species_col, factor_name)
    if (!all(required_cols %in% names(analysis_data))) {
      missing_cols <- setdiff(required_cols, names(analysis_data))
      stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Filter out zero frequencies for meaningful comparisons
    filtered_data <- analysis_data[analysis_data[[freq_col]] > 0, ]
    
    if (nrow(filtered_data) == 0) {
      warning("No non-zero frequency data available for single-factor comparison")
      return(list(
        plot = ggplot2::ggplot() + 
          ggplot2::ggtitle("No Data Available") + 
          ggplot2::labs(subtitle = sprintf("No valid data for factor: %s", factor_name)),
        data = plot_data, 
        summary_statistics = data.frame(),
        summary_table_gt = NULL,
        species_summary_table_gt = NULL,
        statistical_results = list(),
        statistical_report_text = NULL,
        metadata = list(factor_analyzed = factor_name, n_observations = 0)
      ))
    }
    
    # Auto-detect factor type
    factor_type <- if (factor_name %in% c("region_type", "var_type", "genome_region")) "data_factor" else "group_factor"
    
    # Get factor-level order from data_scoping
    # The order specified in data_scoping configuration determines the display order
    preferred_order <- NULL
    if (factor_name == "region_type" && !is.null(task_params$data_scoping$target_region_types)) {
      preferred_order <- task_params$data_scoping$target_region_types
      log_message(sprintf("Using data_scoping order for region_type: %s", 
                         paste(preferred_order, collapse = ", ")))
    } else if (factor_name == "var_type" && !is.null(task_params$data_scoping$target_var_types)) {
      preferred_order <- task_params$data_scoping$target_var_types
      log_message(sprintf("Using data_scoping order for var_type: %s", 
                         paste(preferred_order, collapse = ", ")))
    } else if (factor_name == "genome_region" && !is.null(task_params$data_scoping$target_genome_regions)) {
      preferred_order <- task_params$data_scoping$target_genome_regions
      log_message(sprintf("Using data_scoping order for genome_region: %s", 
                         paste(preferred_order, collapse = ", ")))
    }
    
    # ============================================================
    # Remove empty factor levels
    # ============================================================
    # After all filtering, some factor levels might have 0 rows.
    # We must drop these empty levels to prevent downstream errors.
    
    # Convert character to factor FIRST, before any factor operations
    if (is.character(filtered_data[[factor_name]])) {
      log_message(sprintf("Converting character column '%s' to factor before sanitization", factor_name))
      filtered_data[[factor_name]] <- factor(filtered_data[[factor_name]])
    }
    
    # NOW it's safe to check levels and apply droplevels
    if (is.factor(filtered_data[[factor_name]])) {
      log_message(sprintf("Pre-sanitization: Factor '%s' has %d levels.",
                         factor_name, nlevels(filtered_data[[factor_name]])))
      
      # This line is the core of the fix:
      filtered_data[[factor_name]] <- droplevels(filtered_data[[factor_name]])
      
      log_message(sprintf("Post-sanitization: Factor '%s' now has %d levels.",
                         factor_name, nlevels(filtered_data[[factor_name]])))
    } else {
      log_message(sprintf("Column '%s' is not a factor after conversion, skipping droplevels", factor_name), level = "warning")
    }
    # ============================================================
    
    # Apply ordering logic
    factor_values <- unique(filtered_data[[factor_name]])
    if (!is.null(preferred_order)) {
      ordered_values <- intersect(preferred_order, factor_values)
      remaining_values <- setdiff(factor_values, ordered_values)
      factor_levels <- c(ordered_values, sort(remaining_values))
    } else {
      factor_levels <- sort(factor_values)
    }
    
    # Apply factor levels to data
    filtered_data[[factor_name]] <- factor(filtered_data[[factor_name]], levels = factor_levels)
    
    # Get colors using V4.1 unified color system (Rule #12) with error handling
    colors <- tryCatch({
      get_color_palette(config, factor_name, factor_levels)
    }, error = function(e) {
      log_message(sprintf("Color palette retrieval failed for factor '%s': %s", factor_name, e$message), level = "error")
      NULL
    })
    
    # Ensure colors is not NULL - provide fallback
    if (is.null(colors) || length(colors) == 0) {
      log_message(sprintf("Color palette returned NULL for factor '%s', using default colors", factor_name), level = "warning")
      colors <- grDevices::rainbow(length(factor_levels))
      names(colors) <- factor_levels
    }
    
    # Create base plot
    p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = .data[[factor_name]], y = .data[[freq_col]], 
                                                     fill = .data[[factor_name]], color = .data[[factor_name]]))
    
    # Add main plot layer based on plot_type
    if (plot_type == "violin") {
      p <- p + 
        ggplot2::geom_violin(alpha = alpha_scheme$fill_main, scale = "width") +
        ggplot2::geom_boxplot(width = 0.1, alpha = alpha_scheme$fill_main, outlier.shape = NA)
    } else {
      # Default to boxplot
      p <- p + ggplot2::geom_boxplot(alpha = alpha_scheme$fill_main, outlier.shape = NA)
    }
    
    # Add scatter points with jitter for data visualization
    p <- p + 
      ggplot2::geom_jitter(width = geom_defaults$jitter_width, alpha = alpha_scheme$point_main, size = geom_defaults$point_size) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors)
    
    # ============================================
    # INTELLIGENT ADAPTIVE STATISTICAL ENGINE
    # ============================================
    statistical_results <- list()
    
    if (length(factor_levels) > 1) {
      # Step 1: Prerequisite Testing (Diagnostic Phase)
      log_message("=== INTELLIGENT SINGLE-FACTOR ENGINE: DIAGNOSTIC PHASE ===")
      
      diagnostic_results <- list()
      use_parametric <- FALSE
      
      tryCatch({
        # Test 1: Homogeneity of Variance (Levene's Test)
        if (requireNamespace("car", quietly = TRUE)) {
          levene_formula <- as.formula(paste(freq_col, "~", factor_name))
          levene_result <- car::leveneTest(levene_formula, data = filtered_data)
          diagnostic_results$levene_test <- levene_result
          homogeneity_ok <- levene_result$`Pr(>F)`[1] >= significance_level
          diagnostic_results$homogeneity_ok <- homogeneity_ok
          log_message(sprintf("Levene's Test for Homogeneity: p = %.6f (%s)", 
                             levene_result$`Pr(>F)`[1], 
                             ifelse(homogeneity_ok, "PASSED", "FAILED")))
        } else {
          log_message("car package not available, assuming homogeneity violated", level = "warning")
          homogeneity_ok <- FALSE
          diagnostic_results$homogeneity_ok <- FALSE
        }
        
        # Test 2: Normality Testing (Shapiro-Wilk for each factor level)
        unique_levels <- unique(filtered_data[[factor_name]])
        normality_results <- list()
        normality_ok <- TRUE
        
        log_message(sprintf("Testing normality for %d factor levels...", length(unique_levels)))
        
        for (level in unique_levels) {
          level_data <- filtered_data[filtered_data[[factor_name]] == level, ]
          if (nrow(level_data) >= 3) {  # Minimum sample size for Shapiro-Wilk
            shapiro_result <- stats::shapiro.test(level_data[[freq_col]])
            normality_results[[as.character(level)]] <- shapiro_result
            if (shapiro_result$p.value < significance_level) {
              normality_ok <- FALSE
              log_message(sprintf("Level %s: Shapiro-Wilk p = %.6f (FAILED)", level, shapiro_result$p.value))
            }
          } else {
            log_message(sprintf("Level %s: Insufficient data (n=%d) for normality test", level, nrow(level_data)))
            normality_ok <- FALSE  # Conservative approach
          }
        }
        
        diagnostic_results$normality_tests <- normality_results
        diagnostic_results$normality_ok <- normality_ok
        
        # Decision Logic: IF (Levene's p < 0.05) OR (Any Shapiro-Wilk p < 0.05) THEN Non-parametric
        use_parametric <- homogeneity_ok && normality_ok
        diagnostic_results$statistical_pathway <- ifelse(use_parametric, "parametric", "non-parametric")
        
        log_message(sprintf("=== DECISION: %s pathway selected ===", 
                           ifelse(use_parametric, "PARAMETRIC (ANOVA)", "NON-PARAMETRIC (Kruskal-Wallis)")))
        
      }, error = function(e) {
        log_message(sprintf("Diagnostic testing failed: %s. Defaulting to non-parametric.", e$message), level = "warning")
        use_parametric <- FALSE
        diagnostic_results$error <- e$message
        diagnostic_results$statistical_pathway <- "non-parametric"
      })
      
      # Step 2: Execute Selected Statistical Pathway
      statistical_results$diagnostic <- diagnostic_results
      
      if (use_parametric) {
        # ========== PARAMETRIC PATHWAY: ANOVA + TukeyHSD ==========
        log_message("Executing PARAMETRIC pathway: One-way ANOVA + TukeyHSD")
        tryCatch({
          aov_formula <- as.formula(paste(freq_col, "~", factor_name))
          aov_result <- stats::aov(aov_formula, data = filtered_data)
          statistical_results$anova_test <- summary(aov_result)
          statistical_results$pathway <- "parametric"
          
          # Check if ANOVA is significant for post-hoc
          anova_summary <- statistical_results$anova_test
          if (is.list(anova_summary) && length(anova_summary) > 0) {
            anova_table <- anova_summary[[1]]
            if (is.data.frame(anova_table) && nrow(anova_table) > 0) {
              main_pvalue <- anova_table[1, "Pr(>F)"]
              statistical_results$main_effect_pvalue <- main_pvalue
              
              # Perform TukeyHSD if main effect is significant
              if (!is.na(main_pvalue) && main_pvalue < significance_level) {
                log_message(sprintf("Significant main effect detected (p = %.6f). Performing TukeyHSD post-hoc.", main_pvalue))
                tukey_result <- stats::TukeyHSD(aov_result)
                statistical_results$pairwise_tests <- tukey_result
                statistical_results$post_hoc_significant <- TRUE
              } else {
                statistical_results$post_hoc_significant <- FALSE
                log_message(sprintf("Main effect not significant (p = %.6f). No post-hoc needed.", main_pvalue))
              }
            }
          }
        }, error = function(e) {
          log_message(sprintf("Parametric analysis failed: %s", e$message), level = "error")
          statistical_results$parametric_error <- e$message
        })
        
      } else {
        # ========== NON-PARAMETRIC PATHWAY: Kruskal-Wallis + Dunn's Test ==========
        log_message("Executing NON-PARAMETRIC pathway: Kruskal-Wallis + Dunn's test")
        tryCatch({
          # Kruskal-Wallis test
          kruskal_result <- stats::kruskal.test(filtered_data[[freq_col]], filtered_data[[factor_name]])
          statistical_results$overall_test <- kruskal_result
          statistical_results$pathway <- "non-parametric"
          
          # Dunn's Test for post-hoc pairwise comparisons if Kruskal-Wallis is significant
          if (kruskal_result$p.value < significance_level) {
            log_message(sprintf("Kruskal-Wallis significant (p = %.6f).", kruskal_result$p.value))
            
            # Only perform Dunn's test if there are more than 2 levels
            if (length(factor_levels) > 2) {
              log_message("Performing Dunn's test for post-hoc pairwise comparisons.")
              
              if (requireNamespace("FSA", quietly = TRUE)) {
                dunn_result <- FSA::dunnTest(filtered_data[[freq_col]], 
                                           filtered_data[[factor_name]],
                                           method = "bh")  # Benjamini-Hochberg correction
                statistical_results$pairwise_tests <- dunn_result
                statistical_results$post_hoc_significant <- TRUE
                log_message(sprintf("Dunn's test completed with %d pairwise comparisons", 
                                   ifelse(is.data.frame(dunn_result$res), nrow(dunn_result$res), 0)))
              } else {
                log_message("FSA package not available for Dunn's test", level = "warning")
                statistical_results$dunn_unavailable <- TRUE
              }
            } else {
              log_message("Factor has only 2 levels. Kruskal-Wallis result equivalent to Wilcoxon test. No post-hoc comparison needed.")
              statistical_results$pairwise_tests <- list(message = "Dunn's test not applicable for 2 groups (equivalent to Wilcoxon test).")
              statistical_results$post_hoc_significant <- FALSE
            }
          } else {
            statistical_results$post_hoc_significant <- FALSE
            log_message(sprintf("Kruskal-Wallis not significant (p = %.6f). No post-hoc needed.", kruskal_result$p.value))
          }
        }, error = function(e) {
          log_message(sprintf("Non-parametric analysis failed: %s", e$message), level = "error")
          statistical_results$nonparametric_error <- e$message
        })
      }
    } else {
      log_message("Insufficient factor levels for statistical comparison", level = "warning")
      statistical_results$insufficient_levels <- TRUE
    }
    
    # Add statistical comparison annotations to plot.
    # ggplot2 4.x cannot draw ggpubr::stat_compare_means() overall tests
    # (method = "kruskal.test"/"anova" without comparisons); that layer
    # produces a blank PDF. Pairwise `comparisons =` still renders.
    if (!is.null(statistical_results) && length(statistical_results) > 0) {
      tryCatch({
        overall_p <- if (identical(statistical_results$pathway, "parametric")) {
          statistical_results$main_effect_pvalue
        } else if (!is.null(statistical_results$overall_test$p.value)) {
          statistical_results$overall_test$p.value
        } else {
          NA_real_
        }

        y_max <- max(filtered_data[[freq_col]], na.rm = TRUE)
        if (length(factor_levels) >= 2 && is.finite(overall_p)) {
          p <- p +
            ggplot2::annotate(
              "text",
              x = (length(factor_levels) + 1) / 2,
              y = y_max * 1.12,
              label = sprintf("p = %s", format.pval(overall_p, digits = 2, eps = 1e-16)),
              size = 3.5
            ) +
            ggplot2::scale_y_continuous(
              expand = ggplot2::expansion(mult = c(0.05, 0.28))
            )
        }

        if (requireNamespace("ggpubr", quietly = TRUE) &&
            length(factor_levels) > 2 &&
            length(factor_levels) <= 6) {
          comparison_list <- utils::combn(as.character(factor_levels), 2, simplify = FALSE)

          if (length(comparison_list) <= 15) {
            pairwise_method <- "wilcox.test"
            if (!is.null(statistical_results$pathway) &&
                statistical_results$pathway == "parametric") {
              pairwise_method <- "t.test"
            }

            p <- p + ggpubr::stat_compare_means(
              comparisons = comparison_list,
              method = pairwise_method,
              step.increase = 0.12,
              label = "p.signif",
              hide.ns = FALSE,
              tip.length = 0.02,
              vjust = 0.5,
              bracket.size = 0.6,
              size = 3.5
            )
          }
        }
      }, error = function(e) {
        log_message(sprintf("Plot annotation failed: %s", e$message), level = "warning")
      })
    }
    
    # Apply unified theme (Rule #12 - Unified Color System)
    # V24 FINAL FIX: Correctly call get_application_theme by extracting the theme name string.
    theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"
    base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14
    theme_settings <- get_application_theme(theme_name, NULL, base_size)
    
    # Ensure theme_settings is not NULL - provide fallback
    if (is.null(theme_settings)) {
      log_message("Application theme returned NULL, using default theme", level = "warning")
      theme_settings <- ggplot2::theme_minimal()  # Fallback to default theme
    }
    
    # Create intelligent caption and labels
    factor_display <- tools::toTitleCase(gsub("_", " ", factor_name))
    
    # Build a title that includes the statistical method
    statistical_method <- "Intelligent Adaptive Engine"
    if (!is.null(statistical_results$pathway)) {
      if (statistical_results$pathway == "parametric") {
        statistical_method <- "ANOVA (parametric pathway)"
      } else if (statistical_results$pathway == "non-parametric") {
        statistical_method <- "Kruskal-Wallis (non-parametric pathway)"
      }
    }
    
    # Calculate sample size per factor level for x-axis labels
    level_counts <- filtered_data %>%
      dplyr::group_by(.data[[factor_name]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    
    # Add sample sizes to x-axis labels
    enhanced_labels <- sapply(factor_levels, function(level) {
      count <- level_counts$n[level_counts[[factor_name]] == level]
      if (length(count) == 0) count <- 0
      sprintf("%s\n(n=%d)", level, count)
    })
    names(enhanced_labels) <- factor_levels
    
    # Title format: Factor | Statistical Method (N=total)
    plot_title <- sprintf("%s | %s (N=%d)", 
                         factor_display, 
                         statistical_method, 
                         nrow(filtered_data))
    
    p <- p + 
      theme_settings +
      ggplot2::labs(
        title = plot_title,
        subtitle = NULL,         # Remove subtitle to avoid redundancy
        x = NULL,                # Remove redundant x-axis label
        y = "Frequency per Kilobase",
        fill = factor_display,
        color = factor_display,
        caption = NULL           # Remove caption (now in title)
      ) +
      ggplot2::scale_x_discrete(labels = enhanced_labels) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.caption = ggplot2::element_text(size = geom_defaults$text_size, color = geom_defaults$text_color),
        legend.position = "bottom"
      )
    
    # Calculate comprehensive summary statistics
    summary_stats <- tryCatch({
      filtered_data %>%
        dplyr::group_by(.data[[factor_name]]) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean_freq = mean(.data[[freq_col]], na.rm = TRUE),
          median_freq = stats::median(.data[[freq_col]], na.rm = TRUE),
          sd_freq = stats::sd(.data[[freq_col]], na.rm = TRUE),
          q25 = stats::quantile(.data[[freq_col]], 0.25, na.rm = TRUE),
          q75 = stats::quantile(.data[[freq_col]], 0.75, na.rm = TRUE),
          min_freq = min(.data[[freq_col]], na.rm = TRUE),
          max_freq = max(.data[[freq_col]], na.rm = TRUE),
          .groups = "drop"
        )
    }, error = function(e) {
      log_message(sprintf("Failed to generate summary statistics: %s", e$message), level = "warning")
      data.frame()  # Return empty data frame as fallback
    })
    
    # Ensure summary_stats is not NULL
    if (is.null(summary_stats)) {
      log_message("Summary statistics is NULL, using empty data frame", level = "warning")
      summary_stats <- data.frame()
    }
    
    # Generate professional summary table using gt package if available
    summary_table_gt <- NULL
    if (requireNamespace("gt", quietly = TRUE) && nrow(summary_stats) > 0) {
      tryCatch({
        factor_display <- tools::toTitleCase(gsub("_", " ", factor_name))
        summary_table_gt <- summary_stats %>%
          gt::gt() %>%
          gt::tab_header(
            title = paste("Single Factor Analysis Summary:", factor_display)
          ) %>%
          gt::cols_label(
            !!factor_name := factor_display,
            n = "Sample Size (n)",
            mean_freq = "Mean Frequency",
            median_freq = "Median Frequency", 
            sd_freq = "Standard Deviation",
            q25 = "Q1 (25th percentile)",
            q75 = "Q3 (75th percentile)",
            min_freq = "Minimum",
            max_freq = "Maximum"
          ) %>%
          gt::fmt_number(
            columns = c("mean_freq", "median_freq", "sd_freq", "q25", "q75", "min_freq", "max_freq"),
            decimals = 2
          )
        
        log_message(sprintf("Generated professional summary table for factor '%s' with %d rows", factor_name, nrow(summary_stats)))
      }, error = function(e) {
        log_message(sprintf("Failed to generate gt summary table for factor '%s': %s", factor_name, e$message), level = "warning")
        summary_table_gt <- NULL
      })
    } else if (!requireNamespace("gt", quietly = TRUE)) {
      log_message("gt package not available. Summary table generation skipped.", level = "warning")
    }
    
    # Calculate species-level detailed summary statistics
    species_summary_stats <- tryCatch({
      filtered_data %>%
        dplyr::group_by(.data[[factor_name]], species) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean_freq = mean(.data[[freq_col]], na.rm = TRUE),
          median_freq = stats::median(.data[[freq_col]], na.rm = TRUE),
          sd_freq = stats::sd(.data[[freq_col]], na.rm = TRUE),
          q25 = stats::quantile(.data[[freq_col]], 0.25, na.rm = TRUE),
          q75 = stats::quantile(.data[[freq_col]], 0.75, na.rm = TRUE),
          min_freq = min(.data[[freq_col]], na.rm = TRUE),
          max_freq = max(.data[[freq_col]], na.rm = TRUE),
          .groups = "drop"
        )
    }, error = function(e) {
      log_message(sprintf("Failed to generate species-level summary statistics: %s", e$message), level = "warning")
      data.frame()  # Return empty data frame as fallback
    })
    
    # Ensure species_summary_stats is not NULL
    if (is.null(species_summary_stats)) {
      log_message("Species-level summary statistics is NULL, using empty data frame", level = "warning")
      species_summary_stats <- data.frame()
    }
    
    # Generate professional species-level summary table using gt package if available
    species_summary_table_gt <- NULL
    if (requireNamespace("gt", quietly = TRUE) && nrow(species_summary_stats) > 0) {
      tryCatch({
        factor_display <- tools::toTitleCase(gsub("_", " ", factor_name))
        species_summary_table_gt <- species_summary_stats %>%
          gt::gt() %>%
          gt::tab_header(
            title = paste("Single Factor Species-Level Details:", factor_display)
          ) %>%
          gt::cols_label(
            !!factor_name := factor_display,
            species = "Species",
            n = "Sample Size (n)",
            mean_freq = "Mean Frequency",
            median_freq = "Median Frequency", 
            sd_freq = "Standard Deviation",
            q25 = "Q1 (25th percentile)",
            q75 = "Q3 (75th percentile)",
            min_freq = "Minimum",
            max_freq = "Maximum"
          ) %>%
          gt::fmt_number(
            columns = c("mean_freq", "median_freq", "sd_freq", "q25", "q75", "min_freq", "max_freq"),
            decimals = 2
          )
        
        log_message(sprintf("Generated professional species-level summary table for factor '%s' with %d rows", factor_name, nrow(species_summary_stats)))
      }, error = function(e) {
        log_message(sprintf("Failed to generate gt species-level summary table for factor '%s': %s", factor_name, e$message), level = "warning")
        species_summary_table_gt <- NULL
      })
    } else if (!requireNamespace("gt", quietly = TRUE)) {
      log_message("gt package not available. Species-level summary table generation skipped.", level = "warning")
    }
    
    # Prepare metadata for report generation (including diagnostic information)
    analysis_metadata <- list(
      factor_analyzed = factor_name,
      factor_type = factor_type,
      plot_type = plot_type,
      statistical_pathway = if (!is.null(statistical_results$pathway)) statistical_results$pathway else "unknown",
      intelligent_engine = TRUE,  # Flag indicating this uses the new adaptive engine
      diagnostic_results = if (!is.null(statistical_results$diagnostic)) statistical_results$diagnostic else NULL,
      n_categories = length(factor_levels),
      category_names = factor_levels,
      total_observations = nrow(filtered_data),
      original_observations = nrow(plot_data),
      significance_level = significance_level
    )
    
    # Generate statistical report text for aggregation (no file output here)
    statistical_report_text <- NULL
    if (!is.null(statistical_results) && length(statistical_results) > 0) {
      tryCatch({
        # Generate the statistical report using unified helper (no file output)
        statistical_report_text <- generate_statistical_report_text(
          analysis_results = statistical_results,
          analysis_metadata = analysis_metadata,
          output_file = NULL  # Return text only, no file output
        )
        
      }, error = function(e) {
        log_message(sprintf("Failed to generate statistical report for %s: %s", 
                           factor_name, e$message), level = "warning")
      })
    }
    
    # Return comprehensive results
    return(list(
      plot = p,
      data = filtered_data,
      summary_statistics = summary_stats,
      summary_table_gt = summary_table_gt,
      species_summary_table_gt = species_summary_table_gt,
      statistical_results = statistical_results,
      statistical_report_text = statistical_report_text,
      metadata = analysis_metadata
    ))
    
  }, error = function(e) {
    log_message(sprintf("Universal single-factor comparison failed for factor '%s': %s", 
                       factor_name, e$message), level = "error")
    return(list(
      plot = ggplot2::ggplot() + 
        ggplot2::ggtitle("Universal Single-Factor Analysis Failed") + 
        ggplot2::labs(subtitle = sprintf("Factor: %s | Error: %s", factor_name, e$message)),
      data = if(exists("plot_data")) plot_data else data.frame(),
      summary_statistics = data.frame(),
      summary_table_gt = NULL,
      species_summary_table_gt = NULL,
      statistical_results = list(error = e$message),
      statistical_report_text = NULL,
      metadata = list(
        factor_analyzed = factor_name,
        error = e$message,
        n_observations = 0
      )
    ))
  })
}
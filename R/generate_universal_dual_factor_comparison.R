#' Intelligent Adaptive Dual Factor Statistical Engine
#' 
#' @title Generate universal dual factor comparison
#' @description An intelligent statistical engine that automatically selects the most appropriate
#' statistical approach based on data characteristics (normality, homogeneity of variance).
#' Implements both parametric (ANOVA + TukeyHSD) and non-parametric (Kruskal-Wallis + Dunn's test)
#' pathways with complete diagnostic transparency.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_point scale_fill_manual scale_color_manual labs theme element_text ggtitle position_dodge position_jitterdodge
#' @param plot_data Pre-processed plot data (already cleaned, ordered, and merged with group_info by orchestrator)
#' @param factor1 Name of the first factor (X-axis grouping, e.g., "region_type", "var_type", "Phylogeny")
#' @param factor2 Name of the second factor (fill/color grouping, e.g., "var_type", "Phylogeny", "Life_form")
#' @param config Configuration object 
#' @param task_params Task-specific parameters
#' @param output_dir Optional output directory for saving statistical report TXT file
#' @param generate_plots Whether to generate plot objects in addition to statistics
#' @return List containing plot, data, summary statistics, statistical results (including diagnostic info), statistical report text, and metadata
#' generate_universal_dual_factor_comparison
#' @export
generate_universal_dual_factor_comparison <- function(plot_data,
                                                    factor1,
                                                    factor2,
                                                    config,
                                                    task_params = NULL,
                                                    output_dir = NULL,
                                                    generate_plots = TRUE) {
  
  tryCatch({
    log_message(sprintf("Generating universal dual-factor comparison: %s x %s", factor1, factor2))
    
    # Validate that both factors exist in the data
    if (!factor1 %in% names(plot_data)) {
      available_factors <- names(plot_data)
      stop(sprintf("Factor1 '%s' not found in plot data. Available factors: %s", 
                  factor1, paste(available_factors, collapse = ", ")))
    }
    
    if (!factor2 %in% names(plot_data)) {
      available_factors <- names(plot_data)
      stop(sprintf("Factor2 '%s' not found in plot data. Available factors: %s", 
                  factor2, paste(available_factors, collapse = ", ")))
    }
    
    # Get plot configuration - statistical_test no longer needed (intelligent auto-selection)
    plot_type <- get_task_parameter(task_params, config, "plot_type", "boxplot")
    significance_level <- get_task_parameter(task_params, config, "significance_level", 0.05)
    
    # Get alpha scheme for consistent transparency
    alpha_scheme <- get_alpha_scheme()
    
    # Use internal standard column names (Rule #9 - Data Boundary Principle)
    freq_col <- "frequency_per_kb"
    species_col <- "species"
    
    # Validate required columns
    required_cols <- c(freq_col, species_col, factor1, factor2)
    if (!all(required_cols %in% names(plot_data))) {
      missing_cols <- setdiff(required_cols, names(plot_data))
      stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Drop NA values for both factors being analysed
    # Architecture principle: NA filtering should be per-analysis, not global
    # This ensures species with NA in one factor can still be used in other analyses
    filtered_data <- plot_data %>%
      dplyr::filter(!is.na(.data[[factor1]]) & !is.na(.data[[factor2]]))
    
    na_count <- nrow(plot_data) - nrow(filtered_data)
    if (na_count > 0) {
      log_message(sprintf("Filtered out %d rows with NA values in factors '%s' or '%s'", 
                         na_count, factor1, factor2))
    }
    
    # Filter out zero frequencies for meaningful comparisons
    filtered_data <- filtered_data[filtered_data[[freq_col]] > 0, ]
    
    if (nrow(filtered_data) == 0) {
      warning("No non-zero frequency data available for dual-factor comparison")
      return(list(
        plot = ggplot2::ggplot() + 
          ggplot2::ggtitle("No Data Available") +
          ggplot2::labs(subtitle = sprintf("No valid data for factors: %s x %s", factor1, factor2)),
        data = plot_data, 
        summary_statistics = data.frame(),
        summary_table_gt = NULL,
        species_summary_table_gt = NULL,
        statistical_results = list(),
        statistical_report_text = NULL,
        metadata = list(factors_analyzed = c(factor1, factor2), n_observations = 0)
      ))
    }
    
    # Auto-detect factor types for appropriate ordering
    factor1_type <- if (factor1 %in% c("region_type", "var_type", "genome_region")) "data_factor" else "group_factor"
    factor2_type <- if (factor2 %in% c("region_type", "var_type", "genome_region")) "data_factor" else "group_factor"
    
    # Get factor1 (x-axis) level order from data_scoping
    factor1_preferred_order <- NULL
    if (factor1 == "region_type" && !is.null(task_params$data_scoping$target_region_types)) {
      factor1_preferred_order <- task_params$data_scoping$target_region_types
    } else if (factor1 == "var_type" && !is.null(task_params$data_scoping$target_var_types)) {
      factor1_preferred_order <- task_params$data_scoping$target_var_types
    } else if (factor1 == "genome_region" && !is.null(task_params$data_scoping$target_genome_regions)) {
      factor1_preferred_order <- task_params$data_scoping$target_genome_regions
    }
    
    # Apply ordering logic for factor1
    factor1_values <- unique(filtered_data[[factor1]])
    if (!is.null(factor1_preferred_order)) {
      ordered_values1 <- intersect(factor1_preferred_order, factor1_values)
      remaining_values1 <- setdiff(factor1_values, ordered_values1)
      factor1_levels <- c(ordered_values1, sort(remaining_values1))
    } else {
      factor1_levels <- sort(factor1_values)
    }
    
    # Get factor2 (legend/fill) level order from data_scoping
    factor2_values <- unique(filtered_data[[factor2]])
    factor2_preferred_order <- NULL
    if (factor2 == "region_type" && !is.null(task_params$data_scoping$target_region_types)) {
      factor2_preferred_order <- task_params$data_scoping$target_region_types
    } else if (factor2 == "var_type" && !is.null(task_params$data_scoping$target_var_types)) {
      factor2_preferred_order <- task_params$data_scoping$target_var_types
    } else if (factor2 == "genome_region" && !is.null(task_params$data_scoping$target_genome_regions)) {
      factor2_preferred_order <- task_params$data_scoping$target_genome_regions
    }
    
    # Apply ordering logic for factor2
    if (!is.null(factor2_preferred_order)) {
      ordered_values2 <- intersect(factor2_preferred_order, factor2_values)
      remaining_values2 <- setdiff(factor2_values, ordered_values2)
      factor2_levels <- c(ordered_values2, sort(remaining_values2))
    } else {
      factor2_levels <- sort(factor2_values)
    }
    
    # Apply factor levels to data
    filtered_data[[factor1]] <- factor(filtered_data[[factor1]], levels = factor1_levels)
    filtered_data[[factor2]] <- factor(filtered_data[[factor2]], levels = factor2_levels)
    
    # Get colors using V4.1 unified color system (Rule #12)
    # Priority: factor2 (fill/color aesthetic) gets the primary color mapping
    colors <- get_color_palette(config, factor2, factor2_levels)
    
    # Create base plot with architect-approved strategy: position_dodge instead of facet_wrap
    p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = .data[[factor1]], y = .data[[freq_col]], 
                                                     fill = .data[[factor2]], color = .data[[factor2]]))
    
    # Add main plot layer based on plot_type with position_dodge
    if (plot_type == "violin") {
      p <- p + 
        ggplot2::geom_violin(alpha = alpha_scheme$fill_main, position = ggplot2::position_dodge(0.8), scale = "width") +
        ggplot2::geom_boxplot(width = 0.3, alpha = alpha_scheme$fill_main, 
                             outlier.shape = NA, position = ggplot2::position_dodge(0.8))
    } else {
      # Default to boxplot with architect's recommended position_dodge strategy
      p <- p + ggplot2::geom_boxplot(alpha = alpha_scheme$fill_main, 
                                     outlier.shape = NA, 
                                     position = ggplot2::position_dodge(0.8))
    }
    
    # Add scatter points with position_jitterdodge for better visibility
    p <- p + 
      ggplot2::geom_point(alpha = alpha_scheme$point_main, 
                          size = 1.0, 
                          position = ggplot2::position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors)
    
    # ============================================
    # INTELLIGENT ADAPTIVE STATISTICAL ENGINE
    # ============================================
    statistical_results <- list()
    
    # Step 1: Prerequisite Testing (Diagnostic Phase)
    log_message("=== INTELLIGENT STATISTICAL ENGINE: DIAGNOSTIC PHASE ===")
    
    diagnostic_results <- list()
    use_parametric <- FALSE
    
    if (length(factor1_levels) >= 2 && length(factor2_levels) >= 2) {
      tryCatch({
        # Test 1: Homogeneity of Variance (Levene's Test)
        if (requireNamespace("car", quietly = TRUE)) {
          levene_formula <- as.formula(paste(freq_col, "~", factor1, "*", factor2))
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
        
        # Test 2: Normality Testing (Shapiro-Wilk for each interaction group)
        filtered_data$interaction_group <- interaction(filtered_data[[factor1]], filtered_data[[factor2]])
        unique_groups <- unique(filtered_data$interaction_group)
        normality_results <- list()
        normality_ok <- TRUE
        
        log_message(sprintf("Testing normality for %d interaction groups...", length(unique_groups)))
        
        for (group in unique_groups) {
          group_data <- filtered_data[filtered_data$interaction_group == group, ]
          if (nrow(group_data) >= 3) {  # Minimum sample size for Shapiro-Wilk
            shapiro_result <- stats::shapiro.test(group_data[[freq_col]])
            normality_results[[as.character(group)]] <- shapiro_result
            if (shapiro_result$p.value < significance_level) {
              normality_ok <- FALSE
              log_message(sprintf("Group %s: Shapiro-Wilk p = %.6f (FAILED)", group, shapiro_result$p.value))
            }
          } else {
            log_message(sprintf("Group %s: Insufficient data (n=%d) for normality test", group, nrow(group_data)))
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
        log_message("Executing PARAMETRIC pathway: Two-way ANOVA + TukeyHSD")
        tryCatch({
          aov_formula <- as.formula(paste(freq_col, "~", factor1, "*", factor2))
          aov_result <- stats::aov(aov_formula, data = filtered_data)
          statistical_results$anova_test <- summary(aov_result)
          statistical_results$pathway <- "parametric"
          
          # Check interaction effect significance
          anova_summary <- statistical_results$anova_test
          if (is.list(anova_summary) && length(anova_summary) > 0) {
            anova_table <- anova_summary[[1]]
            if (is.data.frame(anova_table) && !is.null(rownames(anova_table))) {
              interaction_term <- paste(factor1, factor2, sep = ":")
              if (interaction_term %in% rownames(anova_table)) {
                interaction_pvalue <- anova_table[interaction_term, "Pr(>F)"]
                statistical_results$interaction_pvalue <- interaction_pvalue
                
                # Perform TukeyHSD if interaction is significant
                if (!is.na(interaction_pvalue) && interaction_pvalue < significance_level) {
                  log_message(sprintf("Significant interaction detected (p = %.6f). Performing TukeyHSD post-hoc.", interaction_pvalue))
                  tukey_result <- stats::TukeyHSD(aov_result)
                  statistical_results$interaction_pairwise <- tukey_result
                  statistical_results$interaction_significant <- TRUE
                } else {
                  statistical_results$interaction_significant <- FALSE
                  log_message(sprintf("Interaction not significant (p = %.6f). No post-hoc needed.", interaction_pvalue))
                }
              }
            }
          }
        }, error = function(e) {
          log_message(sprintf("Parametric analysis failed: %s", e$message), level = "error")
          statistical_results$parametric_error <- e$message
        })
        
      } else {
        # ========== NON-PARAMETRIC PATHWAY: Proper Two-Way Factorial Analysis ==========
        log_message("Executing NON-PARAMETRIC pathway: Two-way factorial analysis with separate main effects")
        tryCatch({
          # Test main effect of factor1
          factor1_kruskal <- stats::kruskal.test(filtered_data[[freq_col]], filtered_data[[factor1]])
          statistical_results$factor1_main_effect <- factor1_kruskal
          log_message(sprintf("Main effect %s: Kruskal-Wallis p = %.6f", factor1, factor1_kruskal$p.value))
          
          # Test main effect of factor2
          factor2_kruskal <- stats::kruskal.test(filtered_data[[freq_col]], filtered_data[[factor2]])
          statistical_results$factor2_main_effect <- factor2_kruskal
          log_message(sprintf("Main effect %s: Kruskal-Wallis p = %.6f", factor2, factor2_kruskal$p.value))
          
          # For interaction testing in non-parametric context, we use the interaction groups approach
          # but acknowledge this is a limitation of non-parametric factorial analysis
          filtered_data$interaction_group <- interaction(filtered_data[[factor1]], filtered_data[[factor2]])
          interaction_kruskal <- stats::kruskal.test(filtered_data[[freq_col]], filtered_data$interaction_group)
          statistical_results$overall_test <- interaction_kruskal  # This represents the omnibus test
          statistical_results$pathway <- "non-parametric"
          
          log_message(sprintf("Interaction analysis (omnibus): Kruskal-Wallis p = %.6f", interaction_kruskal$p.value))
          
          # Determine which effects are significant and perform post-hoc tests
          significant_effects <- list()
          
          # Post-hoc for factor1 main effect
          if (factor1_kruskal$p.value < significance_level) {
            log_message(sprintf("Significant main effect for %s (p = %.6f).", factor1, factor1_kruskal$p.value))
            
            # Only perform Dunn's test if factor1 has more than 2 levels
            if (length(factor1_levels) > 2) {
              log_message(sprintf("Performing Dunn's test for %s main effect.", factor1))
              if (requireNamespace("FSA", quietly = TRUE)) {
                factor1_dunn <- FSA::dunnTest(filtered_data[[freq_col]], 
                                            filtered_data[[factor1]],
                                            method = "bh")
                statistical_results$factor1_pairwise <- factor1_dunn
                significant_effects$factor1 <- TRUE
                log_message(sprintf("Factor1 Dunn's test: %d pairwise comparisons", 
                                   ifelse(is.data.frame(factor1_dunn$res), nrow(factor1_dunn$res), 0)))
              }
            } else {
              log_message(sprintf("Factor %s has only 2 levels. Kruskal-Wallis result equivalent to Wilcoxon test. No post-hoc comparison needed.", factor1))
              statistical_results$factor1_pairwise <- list(message = "Dunn's test not applicable for 2 groups (equivalent to Wilcoxon test).")
              significant_effects$factor1 <- FALSE
            }
          } else {
            significant_effects$factor1 <- FALSE
          }
          
          # Post-hoc for factor2 main effect
          if (factor2_kruskal$p.value < significance_level) {
            log_message(sprintf("Significant main effect for %s (p = %.6f).", factor2, factor2_kruskal$p.value))
            
            # Only perform Dunn's test if factor2 has more than 2 levels
            if (length(factor2_levels) > 2) {
              log_message(sprintf("Performing Dunn's test for %s main effect.", factor2))
              if (requireNamespace("FSA", quietly = TRUE)) {
                factor2_dunn <- FSA::dunnTest(filtered_data[[freq_col]], 
                                            filtered_data[[factor2]],
                                            method = "bh")
                statistical_results$factor2_pairwise <- factor2_dunn
                significant_effects$factor2 <- TRUE
                log_message(sprintf("Factor2 Dunn's test: %d pairwise comparisons", 
                                   ifelse(is.data.frame(factor2_dunn$res), nrow(factor2_dunn$res), 0)))
              }
            } else {
              log_message(sprintf("Factor %s has only 2 levels. Kruskal-Wallis result equivalent to Wilcoxon test. No post-hoc comparison needed.", factor2))
              statistical_results$factor2_pairwise <- list(message = "Dunn's test not applicable for 2 groups (equivalent to Wilcoxon test).")
              significant_effects$factor2 <- FALSE
            }
          } else {
            significant_effects$factor2 <- FALSE
          }
          
          # Post-hoc for interaction (if overall test is significant)
          if (interaction_kruskal$p.value < significance_level) {
            log_message(sprintf("Significant overall effect (p = %.6f). Performing interaction Dunn's test.", interaction_kruskal$p.value))
            if (requireNamespace("FSA", quietly = TRUE)) {
              interaction_dunn <- FSA::dunnTest(filtered_data[[freq_col]], 
                                              filtered_data$interaction_group,
                                              method = "bh")
              statistical_results$interaction_pairwise <- interaction_dunn
              statistical_results$interaction_significant <- TRUE
              log_message(sprintf("Interaction Dunn's test: %d pairwise comparisons", 
                                 ifelse(is.data.frame(interaction_dunn$res), nrow(interaction_dunn$res), 0)))
            }
          } else {
            statistical_results$interaction_significant <- FALSE
          }
          
          statistical_results$significant_effects <- significant_effects
          
        }, error = function(e) {
          log_message(sprintf("Non-parametric analysis failed: %s", e$message), level = "error")
          statistical_results$nonparametric_error <- e$message
        })
      }
    } else {
      log_message("Insufficient factor levels for dual-factor analysis", level = "warning")
      statistical_results$insufficient_factors <- TRUE
    }
        
    # Add statistical comparison annotations to plot
    if (requireNamespace("ggpubr", quietly = TRUE) && 
        !is.null(statistical_results) && 
        length(statistical_results) > 0) {
      
      tryCatch({
        # Add overall test result based on pathway used
        if (!is.null(statistical_results$pathway)) {
          test_method <- ifelse(statistical_results$pathway == "parametric", "anova", "kruskal.test")
          
          p <- p + ggpubr::stat_compare_means(
            ggplot2::aes(group = .data[[factor1]]), 
            method = test_method,
            label = "p.signif",
            label.y = max(filtered_data[[freq_col]], na.rm = TRUE) * 1.1,
            size = 3,
            hide.ns = FALSE
          )
          
          # Add pairwise comparisons for factor1 if reasonable number of levels
          if (length(factor1_levels) >= 2 && length(factor1_levels) <= 4) {
            comparison_list <- utils::combn(factor1_levels, 2, simplify = FALSE)
            
            if (length(comparison_list) <= 6) {
              pairwise_method <- ifelse(statistical_results$pathway == "parametric", "t.test", "wilcox.test")
              
              p <- p + ggpubr::stat_compare_means(
                comparisons = comparison_list,
                method = pairwise_method,
                step.increase = 0.1,
                label = "p.signif",
                hide.ns = FALSE,
                tip.length = 0.01
              )
            }
          }
        }
      }, error = function(e) {
        log_message(sprintf("Plot annotation failed: %s", e$message), level = "warning")
      })
    }
    
    # Apply unified theme (Rule #12 - Unified Color System)
    theme_name <- config$visualization_settings$theme_and_sizing$theme %||% "professional_light"
    base_size <- config$visualization_settings$theme_and_sizing$base_size %||% 14
    theme_settings <- get_application_theme(theme_name, NULL, base_size)
    
    # Create informative labels (using base R)
    factor1_display <- tools::toTitleCase(gsub("_", " ", factor1))
    factor2_display <- tools::toTitleCase(gsub("_", " ", factor2))
    
    # Build a title that includes the statistical method
    statistical_method <- "Intelligent Adaptive Engine"
    if (!is.null(statistical_results$pathway)) {
      if (statistical_results$pathway == "parametric") {
        statistical_method <- "ANOVA (parametric pathway)"
      } else if (statistical_results$pathway == "non-parametric") {
        statistical_method <- "Kruskal-Wallis (non-parametric pathway)"
      }
    }
    
    # Sample sizes for factor1 levels (x-axis labels)
    factor1_counts <- filtered_data %>%
      dplyr::group_by(.data[[factor1]]) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    
    factor1_levels <- unique(filtered_data[[factor1]])
    enhanced_x_labels <- sapply(factor1_levels, function(level) {
      count <- factor1_counts$n[factor1_counts[[factor1]] == level]
      if (length(count) == 0) count <- 0
      sprintf("%s\n(n=%d)", level, count)
    })
    names(enhanced_x_labels) <- factor1_levels
    
    # Title format: Factor1 x Factor2 | Statistical Method (N=total)
    plot_title <- sprintf("%s x %s | %s (N=%d)", 
                         factor1_display, 
                         factor2_display, 
                         statistical_method, 
                         nrow(filtered_data))
    
    p <- p + 
      theme_settings +
      ggplot2::labs(
        title = plot_title,
        subtitle = NULL,         # Remove subtitle to avoid redundancy
        x = factor1_display,
        y = "Frequency per Kilobase", 
        fill = factor2_display,
        color = factor2_display,
        caption = NULL           # Remove caption (now in title)
      ) +
      ggplot2::scale_x_discrete(labels = enhanced_x_labels) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    # Calculate comprehensive summary statistics by both factors
    summary_stats <- filtered_data %>%
      dplyr::group_by(.data[[factor1]], .data[[factor2]]) %>%
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
    
    # Generate professional summary table using gt package if available
    summary_table_gt <- NULL
    if (requireNamespace("gt", quietly = TRUE) && nrow(summary_stats) > 0) {
      tryCatch({
        factor1_display <- tools::toTitleCase(gsub("_", " ", factor1))
        factor2_display <- tools::toTitleCase(gsub("_", " ", factor2))
        summary_table_gt <- summary_stats %>%
          gt::gt() %>%
          gt::tab_header(
            title = paste("Dual Factor Analysis Summary:", factor1_display, "x", factor2_display)
          ) %>%
          gt::cols_label(
            !!factor1 := factor1_display,
            !!factor2 := factor2_display,
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
        
        log_message(sprintf("Generated professional summary table for factors '%s' x '%s' with %d rows", factor1, factor2, nrow(summary_stats)))
      }, error = function(e) {
        log_message(sprintf("Failed to generate gt summary table for factors '%s' x '%s': %s", factor1, factor2, e$message), level = "warning")
        summary_table_gt <- NULL
      })
    } else if (!requireNamespace("gt", quietly = TRUE)) {
      log_message("gt package not available. Summary table generation skipped.", level = "warning")
    }
    
    # Calculate species-level detailed summary statistics by both factors
    species_summary_stats <- filtered_data %>%
      dplyr::group_by(.data[[factor1]], .data[[factor2]], species) %>%
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
    
    # Generate professional species-level summary table using gt package if available
    species_summary_table_gt <- NULL
    if (requireNamespace("gt", quietly = TRUE) && nrow(species_summary_stats) > 0) {
      tryCatch({
        factor1_display <- tools::toTitleCase(gsub("_", " ", factor1))
        factor2_display <- tools::toTitleCase(gsub("_", " ", factor2))
        species_summary_table_gt <- species_summary_stats %>%
          gt::gt() %>%
          gt::tab_header(
            title = paste("Dual Factor Species-Level Details:", factor1_display, "x", factor2_display)
          ) %>%
          gt::cols_label(
            !!factor1 := factor1_display,
            !!factor2 := factor2_display,
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
        
        log_message(sprintf("Generated professional species-level summary table for factors '%s' x '%s' with %d rows", factor1, factor2, nrow(species_summary_stats)))
      }, error = function(e) {
        log_message(sprintf("Failed to generate gt species-level summary table for factors '%s' x '%s': %s", factor1, factor2, e$message), level = "warning")
        species_summary_table_gt <- NULL
      })
    } else if (!requireNamespace("gt", quietly = TRUE)) {
      log_message("gt package not available. Species-level summary table generation skipped.", level = "warning")
    }
    
    # Prepare metadata for statistical report (including diagnostic information)
    analysis_metadata <- list(
      factors_analyzed = c(factor1, factor2),
      factor1_type = factor1_type,
      factor2_type = factor2_type,
      plot_type = plot_type,
      statistical_pathway = if (!is.null(statistical_results$pathway)) statistical_results$pathway else "unknown",
      intelligent_engine = TRUE,  # Flag indicating this uses the new adaptive engine
      diagnostic_results = if (!is.null(statistical_results$diagnostic)) statistical_results$diagnostic else NULL,
      n_factor1_categories = length(factor1_levels),
      n_factor2_categories = length(factor2_levels),
      factor1_levels = factor1_levels,
      factor2_levels = factor2_levels,
      total_observations = nrow(filtered_data),
      original_observations = nrow(plot_data),
      factorial_design = paste0(length(factor1_levels), "x", length(factor2_levels)),
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
        log_message(sprintf("Failed to generate statistical report for %s x %s: %s", 
                           factor1, factor2, e$message), level = "warning")
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
    log_message(sprintf("Universal dual-factor comparison failed for factors '%s' x '%s': %s", 
                       factor1, factor2, e$message), level = "error")
    return(list(
      plot = ggplot2::ggplot() + 
        ggplot2::ggtitle("Universal Dual-Factor Analysis Failed") + 
        ggplot2::labs(subtitle = sprintf("Factors: %s x %s | Error: %s", factor1, factor2, e$message)),
      data = if(exists("plot_data")) plot_data else data.frame(),
      summary_statistics = data.frame(),
      summary_table_gt = NULL,
      species_summary_table_gt = NULL,
      statistical_results = list(error = e$message),
      statistical_report_text = NULL,
      metadata = list(
        factors_analyzed = c(factor1, factor2),
        error = e$message,
        n_observations = 0
      )
    ))
  })
}
# Helper operators for string formatting (defined outside function to ensure availability)
`%rep%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

`%&%` <- function(x, y) {
  paste0(x, y)
}

#' Smart P-value Formatter
#' 
#' Automatically selects the best display format based on P-value magnitude:
#' - P < 0.001: Scientific notation (e.g., 1.234E-05)
#' - P >= 0.001: Fixed-point decimal (e.g., 0.045000)
#' - Non-numeric/NA: Returns "N/A"
#'
#' @param p Numeric P-value or vector of P-values
#' @return Formatted string representation
#' @noRd
format_p_value <- function(p) {
  if (is.null(p) || (length(p) == 1 && is.na(p))) return("N/A")
  if (!is.numeric(p)) return(as.character(p))
  
  sapply(p, function(val) {
    if (is.na(val)) return("N/A")
    if (val < 0.001) {
      return(sprintf("%.3E", val))
    } else {
      return(sprintf("%.6f", val))
    }
  })
}

#' Generate Statistical Report Text
#' 
#' Creates a formatted text report from statistical analysis results.
#' This unified helper function handles both single-factor and dual-factor
#' statistical analyses, providing consistent formatting across all M02 analyses.
#'
#' @param analysis_results List containing statistical results from analysis functions
#' @param analysis_metadata Metadata about the analysis (factor names, test methods, etc.)
#' @param output_file Path where the text report should be saved
#' @return Formatted text content (also saved to file if output_file provided)
#'
#' @examples
#' \dontrun{
#' # For single factor analysis
#' generate_statistical_report_text(
#'   single_factor_results$statistical_results,
#'   single_factor_results$metadata,
#'   "single_factor_stats.txt"
#' )
#' 
#' # For dual factor analysis  
#' generate_statistical_report_text(
#'   dual_factor_results$statistical_results,
#'   dual_factor_results$metadata,
#'   "dual_factor_stats.txt"
#' )
#' }
#' @export
generate_statistical_report_text <- function(analysis_results, 
                                           analysis_metadata, 
                                           output_file = NULL) {
  
  tryCatch({
    
    # Initialize report content
    report_lines <- c()
    
    # Determine analysis type and extract key information for executive summary
    is_dual_factor <- length(analysis_metadata$factors_analyzed %||% analysis_metadata$factor_analyzed) > 1
    analysis_type <- ifelse(is_dual_factor, "Dual-Factor Analysis", "Single-Factor Analysis")
    
    # Header section
    report_lines <- c(report_lines,
      paste(rep("=", 80), collapse = ""),
      "STATISTICAL ANALYSIS REPORT",
      paste(rep("=", 80), collapse = ""),
      paste("Generated:", Sys.time()),
      paste("Analysis Type:", analysis_type),
      ""
    )
    
    # ============================================================================
    # SECTION 1: EXECUTIVE SUMMARY (New pyramid structure)
    # ============================================================================
    report_lines <- c(report_lines,
      "EXECUTIVE SUMMARY",
      paste(rep("-", 40), collapse = "")
    )
    
    # Generate dynamic summary based on results
    summary_lines <- c()
    
    if (is_dual_factor) {
      # For dual factor analyses
      factors <- analysis_metadata$factors_analyzed
      factor1_significant <- !is.null(analysis_results$factor1_main_effect) && 
                           (analysis_results$factor1_main_effect$p.value %||% 1) < 0.05
      factor2_significant <- !is.null(analysis_results$factor2_main_effect) && 
                           (analysis_results$factor2_main_effect$p.value %||% 1) < 0.05
      interaction_significant <- !is.null(analysis_results$interaction_pvalue) && 
                               (analysis_results$interaction_pvalue %||% 1) < 0.05
      
      if (interaction_significant) {
        summary_lines <- c(summary_lines, 
          sprintf("Key Finding: Significant interaction between %s and %s detected.", factors[1], factors[2]),
          "Interpretation: The effect of one factor depends on the level of the other factor."
        )
      } else if (factor1_significant && factor2_significant) {
        summary_lines <- c(summary_lines,
          sprintf("Key Finding: Both %s and %s show significant main effects.", factors[1], factors[2]),
          "Interpretation: Both factors independently influence the outcome."
        )
      } else if (factor1_significant) {
        summary_lines <- c(summary_lines,
          sprintf("Key Finding: %s has a significant main effect.", factors[1]),
          "Interpretation: The primary factor shows meaningful differences between groups."
        )
      } else if (factor2_significant) {
        summary_lines <- c(summary_lines,
          sprintf("Key Finding: %s has a significant main effect.", factors[2]),
          "Interpretation: The secondary factor shows meaningful differences between groups."
        )
      } else {
        summary_lines <- c(summary_lines,
          "Key Finding: No significant differences detected between any factor combinations.",
          "Interpretation: Current data shows no evidence of systematic group differences."
        )
      }
    } else {
      # For single factor analyses
      factor_name <- analysis_metadata$factor_analyzed
      overall_significant <- !is.null(analysis_results$overall_test) && 
                           (analysis_results$overall_test$p.value %||% 1) < 0.05
      
      if (overall_significant) {
        summary_lines <- c(summary_lines,
          sprintf("Key Finding: %s shows significant differences between groups.", factor_name),
          "Interpretation: Post-hoc comparisons reveal which specific groups differ."
        )
      } else {
        summary_lines <- c(summary_lines,
          sprintf("Key Finding: No significant differences found between %s groups.", factor_name),
          "Interpretation: Current data shows no evidence of systematic group differences."
        )
      }
    }
    
    report_lines <- c(report_lines, summary_lines, "")
    
    # ============================================================================
    # SECTION 2: KEY RESULTS TABLE (Simplified overview)
    # ============================================================================
    report_lines <- c(report_lines,
      "KEY RESULTS TABLE",
      "-" %rep% 40
    )
    
    if (is_dual_factor) {
      # Simplified dual-factor results table
      factors <- analysis_metadata$factors_analyzed
      report_lines <- c(report_lines,
        sprintf("%-25s | %12s | %15s", "Effect", "P-value", "Significance"),
        "-" %rep% 60
      )
      
      # Main effect 1
      if (!is.null(analysis_results$factor1_main_effect)) {
        p_val <- analysis_results$factor1_main_effect$p.value %||% "N/A"
        sig_status <- if (is.numeric(p_val) && p_val < 0.05) "Significant" else "Not Significant"
        p_display <- format_p_value(p_val)
        report_lines <- c(report_lines,
          sprintf("%-25s | %12s | %15s", paste("Main Effect:", factors[1]), p_display, sig_status)
        )
      }
      
      # Main effect 2
      if (!is.null(analysis_results$factor2_main_effect)) {
        p_val <- analysis_results$factor2_main_effect$p.value %||% "N/A"
        sig_status <- if (is.numeric(p_val) && p_val < 0.05) "Significant" else "Not Significant"
        p_display <- format_p_value(p_val)
        report_lines <- c(report_lines,
          sprintf("%-25s | %12s | %15s", paste("Main Effect:", factors[2]), p_display, sig_status)
        )
      }
      
      # Interaction effect
      if (!is.null(analysis_results$interaction_pvalue)) {
        p_val <- analysis_results$interaction_pvalue
        sig_status <- if (is.numeric(p_val) && p_val < 0.05) "Significant" else "Not Significant"
        p_display <- format_p_value(p_val)
        report_lines <- c(report_lines,
          sprintf("%-25s | %12s | %15s", "Interaction Effect", p_display, sig_status)
        )
      }
      
    } else {
      # Simplified single-factor results table
      factor_name <- analysis_metadata$factor_analyzed
      report_lines <- c(report_lines,
        sprintf("%-25s | %12s | %15s", "Test", "P-value", "Significance"),
        "-" %rep% 60
      )
      
      if (!is.null(analysis_results$overall_test)) {
        p_val <- analysis_results$overall_test$p.value %||% "N/A"
        test_method <- analysis_results$overall_test$method %||% "Unknown Test"
        sig_status <- if (is.numeric(p_val) && p_val < 0.05) "Significant" else "Not Significant"
        p_display <- format_p_value(p_val)
        report_lines <- c(report_lines,
          sprintf("%-25s | %12s | %15s", test_method, p_display, sig_status)
        )
      }
    }
    
    report_lines <- c(report_lines, "", "")
    
    # ============================================================================
    # SECTION 3: ANALYSIS OVERVIEW (Moved from original position)
    # ============================================================================
    report_lines <- c(report_lines, 
      "ANALYSIS OVERVIEW",
      "-" %rep% 40
    )
    
    if (!is.null(analysis_metadata$factor_analyzed)) {
      # Single factor analysis
      report_lines <- c(report_lines,
        paste("Factor Analyzed:", analysis_metadata$factor_analyzed),
        paste("Factor Type:", analysis_metadata$factor_type %||% "Unknown"),
        paste("Number of Categories:", analysis_metadata$n_categories %||% "Unknown"),
        paste("Category Names:", paste(analysis_metadata$category_names %||% "Unknown", collapse = ", "))
      )
    } else if (!is.null(analysis_metadata$factors_analyzed)) {
      # Dual factor analysis
      factors <- analysis_metadata$factors_analyzed
      report_lines <- c(report_lines,
        paste("Primary Factor (X-axis):", factors[1]),
        paste("Secondary Factor (Fill/Color):", factors[2]),
        paste("Primary Factor Type:", analysis_metadata$factor1_type %||% "Unknown"),
        paste("Secondary Factor Type:", analysis_metadata$factor2_type %||% "Unknown"),
        paste("Factorial Design:", analysis_metadata$factorial_design %||% "Unknown"),
        paste("Primary Factor Levels:", paste(analysis_metadata$factor1_levels %||% "Unknown", collapse = ", ")),
        paste("Secondary Factor Levels:", paste(analysis_metadata$factor2_levels %||% "Unknown", collapse = ", "))
      )
    }
    
    report_lines <- c(report_lines,
      paste("Plot Type:", analysis_metadata$plot_type %||% "Unknown"),
      paste("Statistical Pathway:", 
            if (!is.null(analysis_metadata$intelligent_engine) && analysis_metadata$intelligent_engine) {
              paste(analysis_metadata$statistical_pathway %||% "Unknown", "(intelligent auto-selection)")
            } else {
              analysis_metadata$statistical_test %||% "Legacy manual selection"
            }),
      paste("Total Observations:", analysis_metadata$total_observations %||% "Unknown"),
      paste("Original Observations:", analysis_metadata$original_observations %||% "Unknown"),
      ""
    )
    
    # Store diagnostic information for later use in methodology appendix
    diagnostic_section <- c()
    if (!is.null(analysis_metadata$intelligent_engine) && analysis_metadata$intelligent_engine &&
        !is.null(analysis_metadata$diagnostic_results)) {
      
      diagnostic <- analysis_metadata$diagnostic_results
      
      diagnostic_section <- c(diagnostic_section,
        "AUTOMATED STATISTICAL STRATEGY DIAGNOSIS",
        "-" %rep% 65,
        paste("Significance Level Used:", analysis_metadata$significance_level %||% "0.05")
      )
      
      # Homogeneity Test Results
      if (!is.null(diagnostic$levene_test)) {
        homogeneity_status <- ifelse(diagnostic$homogeneity_ok, "PASSED", "FAILED")
        levene_p <- diagnostic$levene_test$`Pr(>F)`[1]
        diagnostic_section <- c(diagnostic_section,
          sprintf("Homogeneity of Variance (Levene's Test): p = %s (%s)", format_p_value(levene_p), homogeneity_status)
        )
      } else {
        diagnostic_section <- c(diagnostic_section, "Homogeneity Test: Not available (car package missing)")
      }
      
      # Normality Test Results
      if (!is.null(diagnostic$normality_tests)) {
        normality_status <- ifelse(diagnostic$normality_ok, "PASSED", "FAILED")
        n_groups <- length(diagnostic$normality_tests)
        failed_groups <- sum(sapply(diagnostic$normality_tests, function(x) x$p.value < (analysis_metadata$significance_level %||% 0.05)))
        
        diagnostic_section <- c(diagnostic_section,
          sprintf("Normality Tests (Shapiro-Wilk): %s", normality_status),
          sprintf("  - Tested %d groups/levels", n_groups),
          sprintf("  - Failed normality: %d groups", failed_groups)
        )
        
        # Show details for failed groups
        if (failed_groups > 0 && failed_groups <= 5) {  # Limit detail output
          failed_details <- sapply(names(diagnostic$normality_tests), function(group_name) {
            test <- diagnostic$normality_tests[[group_name]]
            if (test$p.value < (analysis_metadata$significance_level %||% 0.05)) {
              sprintf("    %s: p = %s", group_name, format_p_value(test$p.value))
            } else NULL
          })
          failed_details <- failed_details[!sapply(failed_details, is.null)]
          diagnostic_section <- c(diagnostic_section, failed_details)
        }
      } else {
        diagnostic_section <- c(diagnostic_section, "Normality Tests: Not performed")
      }
      
      # Decision Logic Explanation
      decision_explanation <- paste(
        "Decision Logic: IF (Levene p < 0.05) OR (Any Shapiro-Wilk p < 0.05)",
        "THEN Non-parametric ELSE Parametric"
      )
      
      final_decision <- sprintf("System Decision: %s pathway automatically selected", 
                               diagnostic$statistical_pathway %||% analysis_metadata$statistical_pathway %||% "Unknown")
      
      diagnostic_section <- c(diagnostic_section,
        "",
        decision_explanation,
        final_decision,
        ""
      )
    }
    
    # ============================================================================
    # SECTION 4: DETAILED COMPARISONS (Complete statistical results)
    # ============================================================================
    if (!is.null(analysis_results) && length(analysis_results) > 0) {
      report_lines <- c(report_lines,
        "DETAILED COMPARISONS",
        "-" %rep% 40
      )
      
      # Handle overall test results
      if (!is.null(analysis_results$overall_test)) {
        overall_test <- analysis_results$overall_test
        report_lines <- c(report_lines,
          "Overall Test (Multiple Group Comparison):",
          paste("  Test Method:", overall_test$method %||% "Kruskal-Wallis"),
          paste("  Test Statistic:", 
                if (!is.null(overall_test$statistic)) {
                  sprintf("%.4f", overall_test$statistic)
                } else {
                  "Not available"
                }),
          paste("  P-value:", 
                if (!is.null(overall_test$p.value)) {
                  format_p_value(overall_test$p.value)
                } else {
                  "Not available"
                }),
          paste("  Significance:", 
                if (!is.null(overall_test$p.value)) {
                  if (overall_test$p.value < 0.001) "*** (p < 0.001)"
                  else if (overall_test$p.value < 0.01) "** (p < 0.01)" 
                  else if (overall_test$p.value < 0.05) "* (p < 0.05)"
                  else if (overall_test$p.value < 0.1) ". (p < 0.1)"
                  else "ns (not significant)"
                } else {
                  "Cannot determine"
                }),
          ""
        )
      }
      
      # Handle ANOVA test results (for dual factor analysis)
      if (!is.null(analysis_results$anova_test)) {
        anova_result <- analysis_results$anova_test
        report_lines <- c(report_lines,
          "Two-Way ANOVA Results:",
          "  Source of Variation | Sum of Squares | Mean Square | F-value | P-value",
          "  " %&% ("-" %rep% 70)
        )
        
        # Extract ANOVA table if available
        if (is.list(anova_result) && length(anova_result) > 0) {
          anova_table <- anova_result[[1]]  # First element is usually the ANOVA table
          if (!is.null(anova_table) && is.data.frame(anova_table)) {
            for (i in 1:nrow(anova_table)) {
              row_name <- rownames(anova_table)[i]
              if (row_name != "Residuals") {
                report_lines <- c(report_lines,
                  sprintf("  %-18s | %13.4f | %11.4f | %7.4f | %s", 
                          row_name,
                          anova_table[i, "Sum Sq"] %||% 0,
                          anova_table[i, "Mean Sq"] %||% 0,
                          anova_table[i, "F value"] %||% 0,
                          if (!is.na(anova_table[i, "Pr(>F)"])) {
                            format_p_value(anova_table[i, "Pr(>F)"])
                          } else {
                            "Not available"
                          }
                  )
                )
              }
            }
          }
        }
        report_lines <- c(report_lines, "")
      }
      
      # Handle non-parametric dual-factor main effects analysis
      if (!is.null(analysis_results$factor1_main_effect) && !is.null(analysis_results$factor2_main_effect)) {
        report_lines <- c(report_lines,
          "NON-PARAMETRIC TWO-WAY FACTORIAL ANALYSIS:",
          "----------------------------------------"
        )
        
        # Factor 1 main effect
        factor1_test <- analysis_results$factor1_main_effect
        report_lines <- c(report_lines,
          sprintf("Main Effect - %s:", names(analysis_results)[names(analysis_results) == "factor1_main_effect"] %||% "Factor 1"),
          sprintf("  Test Method: %s", factor1_test$method %||% "Kruskal-Wallis"),
          sprintf("  Test Statistic: %.4f", factor1_test$statistic %||% 0),
          sprintf("  P-value: %s", format_p_value(factor1_test$p.value %||% 1)),
          sprintf("  Significance: %s", 
                 if (!is.null(factor1_test$p.value)) {
                   if (factor1_test$p.value < 0.001) "*** (p < 0.001)"
                   else if (factor1_test$p.value < 0.01) "** (p < 0.01)" 
                   else if (factor1_test$p.value < 0.05) "* (p < 0.05)"
                   else if (factor1_test$p.value < 0.1) ". (p < 0.1)"
                   else "ns (not significant)"
                 } else "Cannot determine"),
          ""
        )
        
        # Factor 2 main effect
        factor2_test <- analysis_results$factor2_main_effect
        report_lines <- c(report_lines,
          sprintf("Main Effect - %s:", names(analysis_results)[names(analysis_results) == "factor2_main_effect"] %||% "Factor 2"),
          sprintf("  Test Method: %s", factor2_test$method %||% "Kruskal-Wallis"),
          sprintf("  Test Statistic: %.4f", factor2_test$statistic %||% 0),
          sprintf("  P-value: %s", format_p_value(factor2_test$p.value %||% 1)),
          sprintf("  Significance: %s", 
                 if (!is.null(factor2_test$p.value)) {
                   if (factor2_test$p.value < 0.001) "*** (p < 0.001)"
                   else if (factor2_test$p.value < 0.01) "** (p < 0.01)" 
                   else if (factor2_test$p.value < 0.05) "* (p < 0.05)"
                   else if (factor2_test$p.value < 0.1) ". (p < 0.1)"
                   else "ns (not significant)"
                 } else "Cannot determine"),
          ""
        )
      }
      
      # Handle main effects post-hoc tests (Dunn's test for each factor)
      if (!is.null(analysis_results$factor1_pairwise)) {
        dunn_result <- analysis_results$factor1_pairwise
        report_lines <- c(report_lines,
          sprintf("MAIN EFFECT POST-HOC ANALYSIS - Factor 1 (Dunn's Test):"),
          sprintf("  P-value Adjustment Method: %s", dunn_result$method %||% "Benjamini-Hochberg"),
          "",
          "  Pairwise Comparisons:",
          "  " %&% ("-" %rep% 70)
        )
        
        # Format Dunn's test results for factor 1
        if (!is.null(dunn_result$res) && is.data.frame(dunn_result$res) && nrow(dunn_result$res) > 0) {
          report_lines <- c(report_lines,
            sprintf("  %-25s | %10s | %10s | %12s", 
                   "Comparison", "Z.statistic", "P.unadj", "P.adj"),
            "  " %&% ("-" %rep% 65)
          )
          
          # +++ START OF NEW, CORRECTED CODE +++
          for (i in 1:nrow(dunn_result$res)) {
            comparison <- dunn_result$res[i, "Comparison"] %||% paste("Row", i)
            z_stat <- dunn_result$res[i, "Z"] %||% 0
            p_unadj_val <- dunn_result$res[i, "P.unadj"] %||% "N/A"
            p_adj_val <- dunn_result$res[i, "P.adj"] %||% "N/A"

            # Format p-values using format_p_value function
            formatted_p_unadj <- format_p_value(p_unadj_val)
            formatted_p_adj <- format_p_value(p_adj_val)
            
            # Determine significance stars
            sig_stars <- ""
            if (is.numeric(p_adj_val) && !is.na(p_adj_val)) {
              if (p_adj_val < 0.001) sig_stars <- "***"
              else if (p_adj_val < 0.01) sig_stars <- "**"
              else if (p_adj_val < 0.05) sig_stars <- "*"
              else if (p_adj_val < 0.1) sig_stars <- "."
            }

            # Final sprintf uses string formats (%s) for p-values to handle both numbers and text
            # and adjusts spacing for proper alignment.
            report_lines <- c(report_lines,
              sprintf("  %-25s | %11.4f | %11s | %-12s %s",
                      comparison, z_stat, formatted_p_unadj, formatted_p_adj, sig_stars)
            )
          }
          # +++ END OF NEW, CORRECTED CODE +++
          report_lines <- c(report_lines, "")
        }
      }
      
      if (!is.null(analysis_results$factor2_pairwise)) {
        dunn_result <- analysis_results$factor2_pairwise
        report_lines <- c(report_lines,
          sprintf("MAIN EFFECT POST-HOC ANALYSIS - Factor 2 (Dunn's Test):"),
          sprintf("  P-value Adjustment Method: %s", dunn_result$method %||% "Benjamini-Hochberg"),
          "",
          "  Pairwise Comparisons:",
          "  " %&% ("-" %rep% 70)
        )
        
        # Format Dunn's test results for factor 2
        if (!is.null(dunn_result$res) && is.data.frame(dunn_result$res) && nrow(dunn_result$res) > 0) {
          report_lines <- c(report_lines,
            sprintf("  %-25s | %10s | %10s | %12s", 
                   "Comparison", "Z.statistic", "P.unadj", "P.adj"),
            "  " %&% ("-" %rep% 65)
          )
          
          # +++ START OF NEW, CORRECTED CODE +++
          for (i in 1:nrow(dunn_result$res)) {
            comparison <- dunn_result$res[i, "Comparison"] %||% paste("Row", i)
            z_stat <- dunn_result$res[i, "Z"] %||% 0
            p_unadj_val <- dunn_result$res[i, "P.unadj"] %||% "N/A"
            p_adj_val <- dunn_result$res[i, "P.adj"] %||% "N/A"

            # Format p-values using format_p_value function
            formatted_p_unadj <- format_p_value(p_unadj_val)
            formatted_p_adj <- format_p_value(p_adj_val)
            
            # Determine significance stars
            sig_stars <- ""
            if (is.numeric(p_adj_val) && !is.na(p_adj_val)) {
              if (p_adj_val < 0.001) sig_stars <- "***"
              else if (p_adj_val < 0.01) sig_stars <- "**"
              else if (p_adj_val < 0.05) sig_stars <- "*"
              else if (p_adj_val < 0.1) sig_stars <- "."
            }

            # Final sprintf uses string formats (%s) for p-values to handle both numbers and text
            # and adjusts spacing for proper alignment.
            report_lines <- c(report_lines,
              sprintf("  %-25s | %11.4f | %11s | %-12s %s",
                      comparison, z_stat, formatted_p_unadj, formatted_p_adj, sig_stars)
            )
          }
          # +++ END OF NEW, CORRECTED CODE +++
          report_lines <- c(report_lines, "")
        }
      }
      
      # Handle interaction post-hoc tests (Dunn's Test or TukeyHSD results)
      if (!is.null(analysis_results$interaction_pairwise)) {
        interaction_result <- analysis_results$interaction_pairwise
        
        # Check if this is a Dunn's test result (FSA package) or TukeyHSD result
        if (!is.null(interaction_result$res) && is.data.frame(interaction_result$res)) {
          # This is a Dunn's test result from FSA package
          dunn_result <- interaction_result
          
          report_lines <- c(report_lines,
            "INTERACTION EFFECT POST-HOC ANALYSIS (Dunn's Test):",
            sprintf("  Interaction Test Method: %s", dunn_result$method %||% "Benjamini-Hochberg"),
            "",
            "  Pairwise Comparisons Between All Factor Level Combinations:",
            "  " %&% ("-" %rep% 70)
          )
          
          # Format Dunn's test results using our robust method
          dunn_table <- dunn_result$res
          if (is.data.frame(dunn_table) && nrow(dunn_table) > 0) {
            
            # Add header
            report_lines <- c(report_lines,
              sprintf("  %-35s | %10s | %10s | %12s", 
                     "Comparison", "Z.statistic", "P.unadj", "P.adj"),
              "  " %&% ("-" %rep% 75)
            )
            
            # +++ START OF ROBUST DUNN'S TEST FORMATTING +++
            for (i in 1:nrow(dunn_table)) {
              comparison <- dunn_table[i, "Comparison"] %||% paste("Row", i)
              z_stat <- dunn_table[i, "Z"] %||% 0
              p_unadj_val <- dunn_table[i, "P.unadj"] %||% "N/A"
              p_adj_val <- dunn_table[i, "P.adj"] %||% "N/A"

              # Initialize formatted values
              formatted_p_unadj <- ""
              formatted_p_adj <- ""
              sig_stars <- ""

              # Format unadjusted p-value robustly
              if (is.numeric(p_unadj_val)) {
                formatted_p_unadj <- format_p_value(p_unadj_val)
              } else {
                formatted_p_unadj <- as.character(p_unadj_val)
              }

              # Format adjusted p-value robustly and determine significance stars
              if (is.numeric(p_adj_val)) {
                formatted_p_adj <- format_p_value(p_adj_val)
                if (!is.na(p_adj_val)) {
                  if (p_adj_val < 0.001) sig_stars <- "***"
                  else if (p_adj_val < 0.01) sig_stars <- "**"
                  else if (p_adj_val < 0.05) sig_stars <- "*"
                  else if (p_adj_val < 0.1) sig_stars <- "."
                  else sig_stars <- ""
                }
              } else {
                # If not numeric, format as string and ensure no significance stars
                formatted_p_adj <- as.character(p_adj_val)
                sig_stars <- ""
              }

              # Final sprintf uses string formats (%s) for p-values to handle both numbers and text
              report_lines <- c(report_lines,
                sprintf("  %-35s | %11.4f | %11s | %-12s %s",
                        comparison, z_stat, formatted_p_unadj, formatted_p_adj, sig_stars)
              )
            }
            # +++ END OF ROBUST DUNN'S TEST FORMATTING +++
            
            report_lines <- c(report_lines, "")
          }
          
          # Add interpretation note for Dunn's test
          report_lines <- c(report_lines,
            "  Interaction Dunn's Test Notes:",
            "  - Non-parametric post-hoc test for interaction effects",
            "  - Z-statistics show the magnitude of rank differences between factor combinations",
            "  - P-values are adjusted using Benjamini-Hochberg method (controls FDR)",
            "  - Appropriate when interaction violates normality or homogeneity assumptions",
            ""
          )
          
        } else if (is.list(interaction_result) && length(interaction_result) > 0) {
          # This is a TukeyHSD result (original logic)
          tukey_result <- interaction_result
        
        report_lines <- c(report_lines,
          "INTERACTION EFFECT POST-HOC ANALYSIS (Tukey's HSD):",
          "",
          "  Pairwise Comparisons Between All Factor Level Combinations:",
          "  " %&% ("-" %rep% 70)
        )
        
          # Format TukeyHSD results (original logic preserved)
          for (term_name in names(tukey_result)) {
            term_results <- tukey_result[[term_name]]
            if (is.matrix(term_results) || is.data.frame(term_results)) {
              
              # Add header for this comparison term
              report_lines <- c(report_lines,
                sprintf("  %s:", term_name),
                sprintf("  %-35s | %10s | %10s | %10s | %12s", 
                       "Comparison", "Diff", "Lower CI", "Upper CI", "Adj. P-value"),
                "  " %&% ("-" %rep% 85)
              )
              
              # Add each pairwise comparison
              comparison_names <- rownames(term_results)
              if (!is.null(comparison_names)) {
                for (i in 1:nrow(term_results)) {
                  comparison <- comparison_names[i]
                  diff_val <- term_results[i, "diff"] %||% 0
                  lwr_val <- term_results[i, "lwr"] %||% 0
                  upr_val <- term_results[i, "upr"] %||% 0
                  adj_pval <- term_results[i, "p adj"] %||% 1
                  
                  # Add significance stars with type checking for safety
                  sig_stars <- ""
                  if (is.numeric(adj_pval) && !is.na(adj_pval)) {
                    if (adj_pval < 0.001) sig_stars <- "***"
                    else if (adj_pval < 0.01) sig_stars <- "**" 
                    else if (adj_pval < 0.05) sig_stars <- "*"
                    else if (adj_pval < 0.1) sig_stars <- "."
                    else sig_stars <- ""
                  }
                  
                  # Safe sprintf with type checking
                  if (is.numeric(adj_pval)) {
                  report_lines <- c(report_lines,
                    sprintf("  %-35s | %10.4f | %10.4f | %10.4f | %8.6f %s", 
                           comparison, diff_val, lwr_val, upr_val, adj_pval, sig_stars)
                  )
                  } else {
                    report_lines <- c(report_lines,
                      sprintf("  %-35s | %10.4f | %10.4f | %10.4f | %-8s %s", 
                             comparison, diff_val, lwr_val, upr_val, as.character(adj_pval), sig_stars)
                    )
                  }
                }
              }
              report_lines <- c(report_lines, "")
          }
        }
        
          # Add interpretation note for TukeyHSD
        report_lines <- c(report_lines,
          "  Post-hoc Analysis Notes:",
          "  - Tukey's HSD provides pairwise comparisons between all factor combinations",
          "  - P-values are adjusted for multiple comparisons (family-wise error control)",
          "  - Confidence intervals show the range of plausible differences",
          "  - Positive differences indicate higher values in the first group",
          ""
        )
        }
      }
      
      # Handle non-parametric post-hoc results (Dunn's Test from FSA package)
      if (!is.null(analysis_results$pairwise_tests) && 
          !is.null(analysis_results$pathway) && 
          analysis_results$pathway == "non-parametric" &&
          is.list(analysis_results$pairwise_tests) &&
          !is.null(analysis_results$pairwise_tests$res)) {
        
        dunn_result <- analysis_results$pairwise_tests
        
        report_lines <- c(report_lines,
          "NON-PARAMETRIC POST-HOC ANALYSIS (Dunn's Test):",
          sprintf("  Overall Test P-value: %.6f (significant)", 
                 analysis_results$overall_test$p.value %||% 
                 analysis_results$main_effect_pvalue %||% "Unknown"),
          sprintf("  P-value Adjustment Method: %s", dunn_result$method %||% "Benjamini-Hochberg"),
          "",
          "  Pairwise Comparisons:",
          "  " %&% ("-" %rep% 70)
        )
        
        # Format Dunn's test results
        dunn_table <- dunn_result$res
        if (is.data.frame(dunn_table) && nrow(dunn_table) > 0) {
          
          # Add header
          report_lines <- c(report_lines,
            sprintf("  %-25s | %10s | %10s | %12s", 
                   "Comparison", "Z.statistic", "P.unadj", "P.adj"),
            "  " %&% ("-" %rep% 65)
          )
          
          # +++ START OF NEW, CORRECTED CODE +++
          for (i in 1:nrow(dunn_table)) {
            comparison <- dunn_table[i, "Comparison"] %||% paste("Row", i)
            z_stat <- dunn_table[i, "Z"] %||% 0
            p_unadj_val <- dunn_table[i, "P.unadj"] %||% "N/A"
            p_adj_val <- dunn_table[i, "P.adj"] %||% "N/A"

            # Format p-values using format_p_value function
            formatted_p_unadj <- format_p_value(p_unadj_val)
            formatted_p_adj <- format_p_value(p_adj_val)
            
            # Determine significance stars
            sig_stars <- ""
            if (is.numeric(p_adj_val) && !is.na(p_adj_val)) {
              if (p_adj_val < 0.001) sig_stars <- "***"
              else if (p_adj_val < 0.01) sig_stars <- "**"
              else if (p_adj_val < 0.05) sig_stars <- "*"
              else if (p_adj_val < 0.1) sig_stars <- "."
            }

            # Final sprintf uses string formats (%s) for p-values to handle both numbers and text
            # and adjusts spacing for proper alignment.
            report_lines <- c(report_lines,
              sprintf("  %-25s | %11.4f | %11s | %-12s %s",
                      comparison, z_stat, formatted_p_unadj, formatted_p_adj, sig_stars)
            )
          }
          # +++ END OF NEW, CORRECTED CODE +++
          
          report_lines <- c(report_lines, "")
        }
        
        # Add interpretation note for Dunn's test
        report_lines <- c(report_lines,
          "  Dunn's Test Notes:",
          "  - Non-parametric post-hoc test following significant Kruskal-Wallis",
          "  - Z-statistics show the magnitude of rank differences between groups",
          "  - P-values are adjusted using Benjamini-Hochberg method (controls FDR)",
          "  - Appropriate when data violate normality or homogeneity assumptions",
          ""
        )
      }
      
      # Handle parametric pairwise test results (traditional format)
      if (!is.null(analysis_results$pairwise_tests)) {
        pairwise <- analysis_results$pairwise_tests
        report_lines <- c(report_lines,
          "Pairwise Comparisons:",
          paste("  Method:", pairwise$method %||% "Unknown"),
          paste("  P-value Adjustment:", pairwise$p.adjust.method %||% "Unknown"),
          ""
        )
        
        # Format pairwise comparison matrix
        if (!is.null(pairwise$p.value) && is.matrix(pairwise$p.value)) {
          p_matrix <- pairwise$p.value
          report_lines <- c(report_lines, "  P-value Matrix:")
          
          # Add column headers
          col_names <- colnames(p_matrix)
          if (!is.null(col_names)) {
            header_line <- "    " %&% paste(sprintf("%12s", col_names), collapse = " ")
            report_lines <- c(report_lines, header_line)
          }
          
          # Add matrix rows
          for (i in 1:nrow(p_matrix)) {
            row_name <- rownames(p_matrix)[i]
            if (!is.null(row_name)) {
              row_values <- p_matrix[i, ]
              formatted_values <- sapply(row_values, function(x) {
                if (is.na(x)) {
                  sprintf("%12s", "-")
                } else {
                  sprintf("%12.6f", x)
                }
              })
              row_line <- sprintf("%-12s", row_name) %&% paste(formatted_values, collapse = " ")
              report_lines <- c(report_lines, row_line)
            }
          }
          report_lines <- c(report_lines, "")
        }
      }
      
      # Handle any error messages
      if (!is.null(analysis_results$error)) {
        report_lines <- c(report_lines,
          "ERROR ENCOUNTERED:",
          paste("  ", analysis_results$error),
          ""
        )
      }
      
    } else {
      report_lines <- c(report_lines,
        "DETAILED COMPARISONS",
        "-" %rep% 40,
        "No statistical test results available.",
        ""
      )
    }
    
    # ============================================================================
    # SECTION 5: METHODOLOGY APPENDIX (Technical details and diagnostics)
    # ============================================================================
    report_lines <- c(report_lines,
      "METHODOLOGY APPENDIX",
      "=" %rep% 40
    )
    
    # Add the diagnostic section that was stored earlier
    if (length(diagnostic_section) > 0) {
      report_lines <- c(report_lines, diagnostic_section)
    }
    
    # Statistical interpretation and methodology details
    report_lines <- c(report_lines,
      "STATISTICAL INTERPRETATION GUIDE",
      "-" %rep% 40,
      "Significance Levels:",
      "  *** p < 0.001 (highly significant)",
      "  **  p < 0.01  (very significant)", 
      "  *   p < 0.05  (significant)",
      "  .   p < 0.1   (marginally significant)",
      "  ns  p >= 0.1  (not significant)",
      "",
      "Methodological Notes:",
      "- Statistical tests assess whether observed differences are",
      "  likely due to chance or represent real biological variation.",
      "- Multiple comparison corrections (e.g., Benjamini-Hochberg) are applied",
      "  to control family-wise error rates in pairwise tests.",
      "- Non-parametric tests (Wilcoxon, Kruskal-Wallis, Dunn's) are preferred",
      "  when data do not meet normality or homogeneity assumptions.",
      "- Intelligent statistical engine automatically selects appropriate tests",
      "  based on diagnostic results from assumption testing.",
      "",
      "Test Selection Criteria:",
      "- Levene's Test: Assesses homogeneity of variance across groups",
      "- Shapiro-Wilk Test: Evaluates normality within each group",
      "- Kruskal-Wallis Test: Non-parametric alternative to ANOVA",
      "- Dunn's Test: Post-hoc comparisons following significant Kruskal-Wallis",
      ""
    )
    
    # Footer
    report_lines <- c(report_lines,
      "=" %rep% 80,
      "End of Statistical Report",
      "=" %rep% 80
    )
    
    # Combine all lines into single text
    report_text <- paste(report_lines, collapse = "\n")
    
    # Save to file if output_file provided
    if (!is.null(output_file)) {
      tryCatch({
        writeLines(report_text, output_file)
        log_message(sprintf("Statistical report saved to: %s", output_file))
      }, error = function(e) {
        log_message(sprintf("Failed to save statistical report to %s: %s", output_file, e$message), level = "error")
      })
    }
    
    return(report_text)
    
  }, error = function(e) {
    error_msg <- sprintf("Failed to generate statistical report: %s", e$message)
    log_message(error_msg, level = "error")
    
    # Return basic error report
    error_report <- paste(
      "STATISTICAL ANALYSIS REPORT - ERROR",
      paste("Generated:", Sys.time()),
      paste("Error:", e$message),
      "",
      "The statistical analysis encountered an error and could not generate a complete report.",
      "Please check the input data and analysis parameters.",
      sep = "\n"
    )
    
    if (!is.null(output_file)) {
      writeLines(error_report, output_file)
    }
    
    return(error_report)
  })
}


# Helper function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || (is.atomic(x) && length(x) == 1 && is.na(x))) y else x
}
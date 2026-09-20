############################################################
#### task_specific_visualizations.R - Core Visualization Engine ####
############################################################
#
# This script is the central engine for generating all task-driven visualizations.
# It implements the core logic for the following analysis modules:
# - M01: Data Distribution Overview
# - M02: Comparative Statistical Analysis
# - M03: Hotspot Gene Identification and Analysis
#
# Key Features:
# - Task-specific parameter application from YAML configuration
# - Multiple plot variants per module (e.g., boxplot, violin, heatmap)
# - Unified theme and color management for visual consistency
# - Advanced, publication-quality plots using ggplot2 and other libraries
# - Session-aware data loading and output organization
#
############################################################

# Load required libraries

# Source base visualization functions
# source("src/functions/generate_visualizations.R")

# Dependencies automatically loaded in R package context

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

# Visualization modules automatically loaded in R package context

# PATH_CONFIG initialization removed for R package compatibility

#' Load existing processed data for visualization-only mode
#' @param session_id Session identifier
#' @param data_type Data type ("filtered" or "unfiltered") 
#' @return Normalized frequency data or NULL if not found
# load_existing_normalized_data function removed - using version from data_loading_helpers.R

#' Aggregate mutation sites by gene to match original script approach
#' 
#' This function transforms mutation-site level data into gene-level frequency data,
#' exactly matching the approach used in the original 03_filtered_analysis.R script:
#' 1. Groups mutations by Species + Gene + Region_type
#' 2. Counts mutations per gene 
#' 3. Calculates normalized frequency per gene: (count / gene_length) * 1000
#' 4. Creates gene × species frequency matrix for proper statistical analysis
#'
#' @param mutation_data Raw mutation site data (each row = one mutation)
#' @param gene_lengths_data CDS lengths data with columns: Species, Gene, CDS_length
#' @param session_id Session ID for data validation
#' @return Data frame with gene-level aggregated frequencies
# Note: aggregate_mutations_by_gene function removed to follow "single data source" principle
# Use normalized frequency data directly from normalization module
  
# =================================================================
# M01 DATA DISTRIBUTION ANALYSIS MODULE
# =================================================================
# 
# NOTICE: M01 functions have been successfully refactored into:
# src/functions/visualization_modules/viz_module_M01.R
# 
# The legacy generate_M01_data_distribution_plots function has been
# removed to prevent function overriding that caused "unused argument (config = config)" error.
# The task dispatcher now correctly calls the function loaded from the module.
# =================================================================

# =================================================================
# M02 COMPARATIVE ANALYSIS MODULE
# =================================================================
# 
# NOTICE: M02 functions have been successfully refactored into:
# src/functions/visualization_modules/viz_module_M02.R
# 
# All M02 functions including generate_M02_comparative_plots, 
# generate_M02_single_factor_comparison, generate_M02_dual_factor_comparison,
# generate_M02_custom_dual_factor_comparison, and generate_M02_comprehensive_analysis
# have been removed to prevent function overriding conflicts.
# The task dispatcher now correctly calls functions loaded from the module.
# =================================================================

# =================================================================
# M03 HOTSPOT ANALYSIS MODULE
# =================================================================
# 
# NOTICE: M03 functions have been successfully refactored into:
# src/functions/visualization_modules/viz_module_M03.R
# 
# All M03-related functions have been moved to prevent function overriding.
# Task-specific visualization functions loaded successfully.
# ================================================================= 

# log_message("Task-specific visualization functions loaded successfully")
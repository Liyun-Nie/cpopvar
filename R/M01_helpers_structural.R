# ============================================================
# M01 Structural Analysis Helper Functions
# ============================================================
#
# This file contains helper functions for M01-2: Structural & Relational Analysis.

#' M01 Chromosome Distribution Analysis Helper Function
#'
#' @title Generate chromosome distribution analysis for M01 module
#' @description Analyzes and visualizes variant spatial distribution across 
#' chromosomes using windowed density calculations and ideogram plots.
#'
#' @param filtered_data Data frame containing filtered variant data
#' @param config Configuration object
#' @param output_dir Output directory for generated plots
#' @param task_name Task identifier
#' @param task_params Task-specific parameters
#'
#' @return List containing analysis results and file manifest
#'
#' @keywords internal
generate_chromosome_distribution_analysis <- function(filtered_data,
                                                     config,
                                                     output_dir,
                                                     task_name,
                                                     task_params) {

  log_message(sprintf("=== M01 CHROMOSOME DISTRIBUTION ANALYSIS: %s ===", task_name))
  log_message(sprintf("Input data: %d rows across %d species", 
                     nrow(filtered_data), 
                     length(unique(filtered_data$species))))
  
  # Initialize results structure
  analysis_results <- list(
    success = FALSE,
    plots_generated = 0,
    species_processed = character(0),
    manifest = list(),
    error_message = NULL
  )
  
  tryCatch({
    
    # =================================================================
    # Part A: Data Loading and Integration
    # =================================================================

    log_message("STEP 1: Applying data_scoping filters")
    
    # Apply data scoping
    scoped_data <- apply_data_scoping(filtered_data, task_params)
    
    if (is.null(scoped_data) || nrow(scoped_data) == 0) {
      stop("No data remaining after applying data_scoping filters")
    }
    
    log_message(sprintf("Data after scoping: %d rows across %d species", 
                       nrow(scoped_data), 
                       length(unique(scoped_data$species))))

    # Load genome regions data
    log_message("Loading genome regions data for _vs_regions modes")
    session_paths <- get_session_paths(config$session_info$session_id)
    regions_file_path <- file.path(session_paths$processed_data, "P01_preprocessed", "processed_genome_regions.csv")
    
    genome_regions_data <- NULL
    if (file.exists(regions_file_path)) {
      genome_regions_data <- readr::read_csv(regions_file_path, show_col_types = FALSE)
      log_message(sprintf("Successfully loaded genome regions data: %d rows", nrow(genome_regions_data)))
    } else {
      log_message(sprintf("Warning: Genome regions file not found: %s", regions_file_path), level = "warning")
    }

    # Get species list and apply ordering
    species_list <- unique(scoped_data$species)
    
    ordering_result <- load_and_apply_species_ordering(
      data = data.frame(species = species_list),
      session_id = config$session_info$session_id,
      species_col = "species"
    )
    
    if (ordering_result$success) {
      species_list <- ordering_result$species_order
      log_message(sprintf("Applied user-defined species ordering: %s", 
                         paste(species_list, collapse = ", ")))
    } else {
      log_message(sprintf("Using default species order: %s", ordering_result$message))
    }

    # V20 ENHANCEMENT: Extract custom_label_mapping
    custom_label_mapping <- ordering_result$custom_label_mapping
    
    # V20 FIX: Check if custom labels should be used (task-level parameter only)
    use_custom_labels <- FALSE
    if (!is.null(task_params$use_custom_labels)) {
      use_custom_labels <- task_params$use_custom_labels
      log_message(sprintf("M01: use_custom_labels = %s (from task parameters)", use_custom_labels))
    }
    
    # Create species-to-Chr mapping
    # If custom labels are enabled and available, use them; otherwise use numeric IDs
    if (use_custom_labels && !is.null(custom_label_mapping)) {
      species_id_mapping <- custom_label_mapping[species_list]
      log_message(sprintf("M01: Using custom labels for Chr column: %s", 
                         paste(head(species_id_mapping, 3), collapse = ", ")))
    } else {
      species_id_mapping <- setNames(as.character(seq_along(species_list)), species_list)
      log_message("M01: Using numeric IDs for Chr column")
    }
    log_message(sprintf("Created species mapping for %d species", length(species_id_mapping)))

    # Extract display mode
    display_mode <- get_task_parameter(task_params, config, "chromosome_display_mode", "full")
    log_message(sprintf("Using chromosome display mode: %s", display_mode))
    
    # Load and prepare karyotype data
    # V3.16: Simplified configuration - direct file path
    genome_regions_path <- config$input_files$genome_regions
    if (is.null(genome_regions_path) || !file.exists(genome_regions_path)) {
      stop(sprintf("Genome regions file not found. Expected path: %s", 
                   genome_regions_path %||% "not specified"))
    }
    
    genome_regions_raw <- readr::read_csv(genome_regions_path, show_col_types = FALSE)
    
    # Generate karyotype dataframe
    karyotype_data <- data.frame(
      Chr = as.character(species_id_mapping),
      Start = rep(0, length(species_list)),
      End = numeric(length(species_list)),
      stringsAsFactors = FALSE
    )
    
    # Fill genome lengths
    for (i in seq_along(species_list)) {
      species_name <- species_list[i]
      original_karyotype <- genome_regions_raw %>% 
        dplyr::filter(species == species_name)
      
      if (nrow(original_karyotype) > 0) {
        is_ir_lacking <- !is.na(original_karyotype$special_handling[1]) && 
                        original_karyotype$special_handling[1] == "IR_lacking_genome"
        
        if (display_mode == "IRA_only" && !is_ir_lacking) {
          karyotype_data$End[i] <- original_karyotype$ssc_end[1]
        } else {
          karyotype_data$End[i] <- original_karyotype$total_length[1]
        }
      } else {
        log_message(sprintf("No karyotype data found for %s", species_name), level = "warning")
        karyotype_data$End[i] <- 150000
      }
    }

    # Merge Chr ID into variant data
    species_mapping_df <- data.frame(
      species = names(species_id_mapping),
      species_id = as.character(species_id_mapping),
      Chr = as.character(species_id_mapping),
      stringsAsFactors = FALSE
    )
    scoped_data <- merge(scoped_data, species_mapping_df, by = "species")
    log_message("M01: Successfully merged Chr ID into main variant data.")

    # Extract window size
    window_size <- get_task_parameter(task_params, config, "window_size", 2000)
    log_message(sprintf("Using window size: %d bp", window_size))

    # =================================================================
    # Part B: Main Processing
    # =================================================================

    log_message("M01 Part B: Processing multiple plot modes")
    
    plot_modes <- get_task_parameter(task_params, config, "plot_modes", c("snp_vs_indel"))
    
    plots_created <- 0
    generated_files <- c()
    
    # Create output directory
    m01_output_dir <- output_dir
    if (!dir.exists(m01_output_dir)) {
      dir.create(m01_output_dir, recursive = TRUE)
    }

    # V20 ENHANCEMENT: Output species-chromosome mapping (supports custom_label)
    species_chr_mapping <- data.frame(
      Chr = as.character(species_id_mapping),
      Species_Name = names(species_id_mapping),
      Karyotype_Length = karyotype_data$End,
      Display_Mode = display_mode,
      stringsAsFactors = FALSE
    )
    
    mapping_file_path <- file.path(m01_output_dir, "M01_species_chromosome_mapping.txt")
    write.table(species_chr_mapping, 
                file = mapping_file_path, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    log_message(sprintf("M01: Saved species-chromosome mapping"))

    # Pre-calculate densities
    density_snp <- scoped_data %>%
      dplyr::filter(var_type == "snp") %>%
      calculate_windowed_density(
        variant_df = .,
        karyotype_df = karyotype_data,
        species_map = species_id_mapping,
        window_size = window_size
      )
    write.csv(density_snp, file.path(m01_output_dir, "M01_snp_density.csv"), row.names = FALSE)
    
    density_indel <- scoped_data %>%
      dplyr::filter(var_type == "INDEL") %>%
      calculate_windowed_density(
        variant_df = .,
        karyotype_df = karyotype_data,
        species_map = species_id_mapping,
        window_size = window_size
      )
    write.csv(density_indel, file.path(m01_output_dir, "M01_indel_density.csv"), row.names = FALSE)
    
    # Create region markers
    marker_data <- NULL
    if (!is.null(genome_regions_data)) {
      marker_data <- create_region_marker_data(
        regions_data = genome_regions_data,
        species_id_mapping = species_mapping_df,
        display_mode = display_mode
      )
      write.csv(marker_data, file.path(m01_output_dir, "M01_region_markers.csv"), row.names = FALSE)
    }
    
    # Main plotting loop
    for (mode in plot_modes) {
      log_message(sprintf("M01: Processing mode - %s", mode))
      
      tryCatch({
        if (mode == "snp_vs_indel") {
          task_params_with_mode <- task_params
          task_params_with_mode$current_plot_mode <- mode
          
          plot_path <- generate_chromosome_distribution_plot(
            karyotype_data = karyotype_data,
            track1_data = density_snp,
            track2_data = density_indel,
            label_type = "heatmap",
            config = config,
            task_params = task_params_with_mode,
            output_path_base = file.path(m01_output_dir, "M01_distribution_snp_vs_indel")
          )
          
          plots_created <- plots_created + 1
          generated_files <- c(generated_files, plot_path)
        }
        else if (mode == "snp_vs_regions" && !is.null(marker_data)) {
          task_params_with_mode <- task_params
          task_params_with_mode$current_plot_mode <- mode
          
          plot_path <- generate_chromosome_distribution_plot(
            karyotype_data = karyotype_data,
            track1_data = density_snp,
            track2_data = marker_data,
            label_type = "marker",
            config = config,
            task_params = task_params_with_mode,
            output_path_base = file.path(m01_output_dir, "M01_distribution_snp_vs_regions")
          )
          
          plots_created <- plots_created + 1
          generated_files <- c(generated_files, plot_path)
        }
        else if (mode == "indel_vs_regions" && !is.null(marker_data)) {
          task_params_with_mode <- task_params
          task_params_with_mode$current_plot_mode <- mode
          
          plot_path <- generate_chromosome_distribution_plot(
            karyotype_data = karyotype_data,
            track1_data = density_indel,
            track2_data = marker_data,
            label_type = "marker",
            config = config,
            task_params_with_mode,
            output_path_base = file.path(m01_output_dir, "M01_distribution_indel_vs_regions")
          )
          
          plots_created <- plots_created + 1
          generated_files <- c(generated_files, plot_path)
        }
      }, error = function(e) {
        log_message(sprintf("M01: Failed to generate plot for mode %s - %s", mode, e$message), level = "error")
      })
    }

    analysis_results$species_processed <- species_list
    analysis_results$success <- TRUE
    analysis_results$plots_generated <- plots_created
    analysis_results$manifest <- generated_files
    
    log_message(sprintf("M01 chromosome analysis completed: %d plots generated", plots_created))
    
  }, error = function(e) {
    analysis_results$success <- FALSE
    analysis_results$error_message <- e$message
    log_message(sprintf("M01 chromosome analysis failed: %s", e$message), level = "error")
  })
  
  return(analysis_results)
}


#' @title Create Region Marker Data for RIdeogram
#' @description Transforms genome region data into point marker format for RIdeogram
#' @param regions_data Data frame from processed_genome_regions.csv
#' @param species_id_mapping Data frame mapping species to chromosome IDs
#' @param display_mode "full" or "IRA_only"
#' @return Data frame formatted for RIdeogram markers
#' @keywords internal
create_region_marker_data <- function(regions_data, species_id_mapping, display_mode = "full") {
  
  merged_data <- merge(regions_data, species_id_mapping, by = "species")
  merged_data <- merged_data[!is.na(merged_data$ir_start), ]

  if (nrow(merged_data) == 0) {
    log_message("M01: No valid IR regions found", level = "warning")
    return(data.frame(Type=character(), Shape=character(), Chr=character(), 
                     Start=numeric(), End=numeric(), color=character()))
  }

  marker_list <- lapply(1:nrow(merged_data), function(i) {
    row <- merged_data[i, ]
    
    if (display_mode == "IRA_only") {
      markers <- data.frame(
        Type  = c("IR", "IR"),
        Shape = c("triangle", "triangle"),
        Chr   = c(row$species_id, row$species_id),
        Start = c(row$ir_start, row$ir_end),
        End   = c(row$ir_start, row$ir_end),
        color = c("#FF7F0EFF", "#FF7F0EFF"),
        stringsAsFactors = FALSE
      )
    } else {
      ir_length <- row$ir_end - row$ir_start
      irb_start <- row$total_length - ir_length
      irb_end <- row$total_length
      
      markers <- data.frame(
        Type  = c("IR", "IR", "IR", "IR"),
        Shape = c("triangle", "triangle", "triangle", "triangle"),
        Chr   = c(row$species_id, row$species_id, row$species_id, row$species_id),
        Start = c(row$ir_start, row$ir_end, irb_start, irb_end),
        End   = c(row$ir_start, row$ir_end, irb_start, irb_end),
        color = c("#FF7F0EFF", "#FF7F0EFF", "#FF7F0EFF", "#FF7F0EFF"),
        stringsAsFactors = FALSE
      )
    }
    
    return(markers)
  })

  final_marker_data <- do.call(rbind, marker_list)
  
  markers_per_species <- if (display_mode == "IRA_only") 2 else 4
  log_message(sprintf("M01: Created %d markers (%d per species) in '%s' mode", 
                      nrow(final_marker_data), markers_per_species, display_mode))
  return(final_marker_data)
}


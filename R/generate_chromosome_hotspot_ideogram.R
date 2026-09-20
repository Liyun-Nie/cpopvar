############################################################
#### generate_chromosome_hotspot_ideogram.R - Chromosome Hotspot Map Visualization ####
############################################################
#
# Generates chromosome ideogram plots with hotspot markers using RIdeogram
# Part of M01-2 Structural & Relational Analysis
#
# This module visualizes CDS and poiGS hotspots on chromosome ideograms
# to reveal spatial clustering and co-localization patterns
#
############################################################

#' Generate Chromosome Hotspot Ideogram
#'
#' @description Main entry function for generating chromosome hotspot map visualization
#' @param m03_hotspots_path Path to M03_core_shared_hotspots.csv
#' @param m04_hotspots_path Path to M04_core_shared_poigs_hotspots.csv
#' @param config Configuration object
#' @param output_dir Output directory for plots and data
#' @param session_id Session identifier
#' @param task_params Task-specific parameters (for custom_label configuration)
#' @return List containing analysis results and file manifest
#' @export
generate_chromosome_hotspot_ideogram <- function(m03_hotspots_path,
                                                m04_hotspots_path,
                                                config,
                                                output_dir,
                                                session_id,
                                                task_params = list()) {
  
  log_message("=== GATE 4.4: CHROMOSOME HOTSPOT MAP VISUALIZATION ===")
  
  # Initialize results structure
  results <- list(
    success = FALSE,
    plots_generated = 0,
    manifest = list(),
    error_message = NULL
  )
  
  tryCatch({
    
    # ================================================================
    # STEP 1: Load region_info_complete.csv
    # ================================================================
    
    log_message("STEP 1: Loading region_info_complete.csv")
    
    session_paths <- get_session_paths(session_id)
    p01_dir <- get_processing_stage_path(session_id, "P01_preprocessed")
    region_info_path <- file.path(p01_dir, "region_info_complete.csv")
    
    if (!file.exists(region_info_path)) {
      stop(sprintf("region_info_complete.csv not found: %s", region_info_path))
    }
    
    # ================================================================
    # STEP 2: Extract hotspot coordinates
    # ================================================================
    
    log_message("STEP 2: Extracting hotspot coordinates")
    
    # Extract CDS hotspots
    cds_hotspots <- extract_cds_hotspot_coordinates(
      m03_hotspots_path = m03_hotspots_path,
      region_info_path = region_info_path,
      session_id = session_id
    )
    
    # Extract poiGS hotspots
    poigs_hotspots <- extract_poigs_hotspot_coordinates(
      m04_hotspots_path = m04_hotspots_path,
      region_info_path = region_info_path,
      session_id = session_id
    )
    
    # ================================================================
    # STEP 3: Load genome regions and create karyotype
    # ================================================================
    
    log_message("STEP 3: Creating karyotype data")
    
    # Load genome regions
    genome_regions_path <- config$input_files$genome_regions
    if (is.null(genome_regions_path) || !file.exists(genome_regions_path)) {
      stop(sprintf("Genome regions file not found: %s", 
                   genome_regions_path %||% "not specified"))
    }
    
    genome_regions_raw <- readr::read_csv(genome_regions_path, show_col_types = FALSE)
    
    # Get species list and apply ordering
    all_species <- unique(c(cds_hotspots$species, poigs_hotspots$species))
    
    ordering_result <- load_and_apply_species_ordering(
      data = data.frame(species = all_species),
      session_id = session_id,
      species_col = "species"
    )
    
    if (ordering_result$success) {
      species_list <- ordering_result$species_order
      log_message(sprintf("Applied user-defined species ordering: %s", 
                         paste(species_list, collapse = ", ")))
    } else {
      species_list <- all_species
      log_message(sprintf("Using default species order: %s", ordering_result$message))
    }
    
    # V20 ENHANCEMENT: Extract custom_label_mapping
    custom_label_mapping <- ordering_result$custom_label_mapping
    
    # V20 FIX: Check if custom labels should be used (task-level parameter only)
    use_custom_labels <- FALSE
    if (!is.null(task_params$use_custom_labels)) {
      use_custom_labels <- task_params$use_custom_labels
      log_message(sprintf("M05: use_custom_labels = %s (from task parameters)", use_custom_labels))
    }
    
    # Create species-to-Chr mapping
    # If custom labels are enabled and available, use them; otherwise use numeric IDs
    if (use_custom_labels && !is.null(custom_label_mapping)) {
      species_id_mapping <- custom_label_mapping[species_list]
      log_message(sprintf("Using custom labels for Chr column: %s", 
                         paste(head(species_id_mapping, 3), collapse = ", ")))
    } else {
      species_id_mapping <- setNames(as.character(seq_along(species_list)), species_list)
      log_message("Using numeric IDs for Chr column")
    }
    log_message(sprintf("Created species mapping for %d species", length(species_id_mapping)))
    
    # Create karyotype data (use IRA_only mode like M01 to avoid IRB blank space)
    karyotype_data <- create_hotspot_karyotype(
      species_list = species_list,
      genome_regions_raw = genome_regions_raw,
      species_id_mapping = species_id_mapping,
      display_mode = "IRA_only"  # consistent with M01 chromosome distribution
    )
    
    # ================================================================
    # STEP 3.5: Prepare hotspot regions heatmap (CDS and poiGS hotspots only)
    # ================================================================
    
    log_message("STEP 3.5: Preparing hotspot regions heatmap")
    
    # V20 FIX: Prepare genome_regions with Chr mapping using species_id_mapping
    # This ensures Chr column matches karyotype (either custom_label or numeric ID)
    genome_regions_for_heatmap <- genome_regions_raw %>%
      dplyr::filter(species %in% species_list) %>%
      dplyr::mutate(Chr = as.character(species_id_mapping[species]))
    
    # Pass hotspot coordinates to prepare_regions_heatmap
    # Only hotspot regions will be shown in heatmap (not all CDS/IGS)
    regions_heatmap <- prepare_hotspot_regions_heatmap(
      cds_hotspots = cds_hotspots,
      poigs_hotspots = poigs_hotspots,
      genome_regions = genome_regions_for_heatmap
    )
    
    # ================================================================
    # STEP 4: Get colors from config
    # ================================================================
    
    log_message("STEP 4: Retrieving hotspot colors")
    
    hotspot_colors <- get_hotspot_colors_from_config(config)
    
    # ================================================================
    # STEP 5: Prepare hotspot markers
    # ================================================================
    
    log_message("STEP 5: Preparing hotspot markers")
    
    hotspot_markers <- prepare_hotspot_markers(
      cds_hotspots = cds_hotspots,
      poigs_hotspots = poigs_hotspots,
      species_id_mapping = species_id_mapping,
      hotspot_colors = hotspot_colors
    )
    
    # ================================================================
    # STEP 5.5: Load IR markers from M01 chromosome distribution
    # ================================================================
    
    log_message("STEP 5.5: Loading IR markers from M01 chromosome distribution")
    
    # Construct path to M01 region markers
    session_paths <- get_session_paths(session_id)
    m01_marker_path <- file.path(
      session_paths$plots,
      "M01_distribution",
      "structural_relational",
      "chromosome_distribution",
      "M01_region_markers.csv"
    )
    
    if (file.exists(m01_marker_path)) {
      log_message(sprintf("Found M01 region markers: %s", m01_marker_path))
      
      # Load M01 markers (already in RIdeogram format with correct Chr labels)
      m01_markers <- readr::read_csv(m01_marker_path, show_col_types = FALSE)
      
      # V20 FIX: M01 markers already use correct Chr labels (custom_label or numeric)
      # CRITICAL: Ensure Chr is character type for bind_rows compatibility
      chr_to_species <- setNames(names(species_id_mapping), species_id_mapping)
      
      # Extract IR markers and add metadata for sorted table
      ir_markers <- m01_markers %>%
        dplyr::mutate(
          Chr = as.character(Chr),  # convert to character (numeric or custom_label)
          species = chr_to_species[Chr],  # Reverse map: Chr -> species
          label = Type,  # "IR"
          hotspot_type = "IR",
          genome_region = Type,
          num_species_shared = NA_integer_,
          position = Start  # For sorted table
        ) %>%
        dplyr::select(Type, Shape, Chr, Start, End, color, species, label, 
                      hotspot_type, genome_region, num_species_shared)
      
      log_message(sprintf("Loaded %d IR markers from M01", nrow(ir_markers)))
      
      # Combine with hotspot markers
      hotspot_markers <- dplyr::bind_rows(hotspot_markers, ir_markers)
      
      log_message(sprintf("Total markers (including IR): %d", nrow(hotspot_markers)))
    } else {
      log_message(sprintf("INFO: M01 region markers not found at: %s", m01_marker_path), 
                 level = "info")
      log_message("IR markers will not be included in the visualization (M01 module not executed)")
      log_message("M05 will proceed with CDS and IGS hotspots only")
    }
    
    # ================================================================
    # STEP 6: Generate visualizations
    # ================================================================
    
    log_message("STEP 6: Generating chromosome hotspot ideogram")
    
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Generate ideogram plot with overlaid heatmap
    plot_path <- plot_hotspot_ideogram(
      karyotype_data = karyotype_data,
      regions_heatmap = regions_heatmap,
      hotspot_markers = hotspot_markers,
      output_dir = output_dir,
      config = config
    )
    
    results$plots_generated <- results$plots_generated + 1
    results$manifest$ideogram_plot <- plot_path
    
    # ================================================================
    # STEP 7: Export data tables
    # ================================================================
    
    log_message("STEP 7: Exporting data tables")
    
    # Export RIdeogram marker data
    markers_path <- file.path(output_dir, "M05_hotspot_markers.csv")
    readr::write_csv(
      hotspot_markers %>% dplyr::select(Type, Shape, Chr, Start, End, color),
      markers_path
    )
    log_message(sprintf("Exported marker data: %s", markers_path))
    results$manifest$markers_data <- markers_path
    
    # Export heatmap data for hotspot regions only
    heatmap_path <- file.path(output_dir, "M05_hotspot_heatmap.csv")
    readr::write_csv(regions_heatmap, heatmap_path)
    log_message(sprintf("Exported heatmap data: %s (%d regions)", 
                       heatmap_path, nrow(regions_heatmap)))
    results$manifest$heatmap_data <- heatmap_path
    
    # Export sorted marker table (for manual labeling)
    # Reorder columns for readability
    sorted_markers_path <- file.path(output_dir, "M05_hotspot_markers_sorted.csv")
    sorted_markers <- hotspot_markers %>%
      dplyr::arrange(Chr, Start) %>%  # Sort by Chr first, then position
      dplyr::select(
        Chr,
        position = Start,
        species,
        hotspot_type,
        genome_region,
        gene_or_region = label,
        num_species_shared,  # track species sharing
        color
        # Removed 'shape' column as requested
      )
    
    readr::write_csv(sorted_markers, sorted_markers_path)
    log_message(sprintf("Exported sorted marker table: %s", sorted_markers_path))
    results$manifest$sorted_markers <- sorted_markers_path
    
    # Generate a text-based figure legend
    legend_path <- file.path(output_dir, "M05_figure_legend.txt")
    generate_figure_legend(
      output_path = legend_path,
      cds_hotspots = cds_hotspots,
      poigs_hotspots = poigs_hotspots,
      hotspot_colors = get_hotspot_colors_from_config(config)
    )
    log_message(sprintf("Generated figure legend: %s", legend_path))
    results$manifest$figure_legend <- legend_path
    
    # Generate summary statistics
    summary_path <- file.path(output_dir, "M05_hotspot_summary.txt")
    write_hotspot_summary(
      cds_hotspots = cds_hotspots,
      poigs_hotspots = poigs_hotspots,
      hotspot_markers = hotspot_markers,
      output_path = summary_path
    )
    log_message(sprintf("Exported summary statistics: %s", summary_path))
    results$manifest$summary <- summary_path
    
    # ================================================================
    # STEP 8: Finalize results
    # ================================================================
    
    results$success <- TRUE
    log_message("=== CHROMOSOME HOTSPOT MAP VISUALIZATION COMPLETED ===")
    log_message(sprintf("Total plots generated: %d", results$plots_generated))
    log_message(sprintf("Output directory: %s", output_dir))
    
    return(results)
    
  }, error = function(e) {
    results$success <- FALSE
    results$error_message <- e$message
    log_message(sprintf("ERROR in chromosome hotspot map: %s", e$message), level = "error")
    return(results)
  })
}


#' Create Hotspot Karyotype Data
#'
#' @description Creates karyotype data frame for RIdeogram
#' @param species_list Vector of species names (ordered)
#' @param genome_regions_raw Raw genome regions data
#' @param species_id_mapping Named vector mapping species to chromosome IDs
#' @param display_mode Display mode: "full" or "IRA_only" (default: "IRA_only")
#' @return Data frame with karyotype information
#' @keywords internal
create_hotspot_karyotype <- function(species_list, genome_regions_raw, species_id_mapping, 
                                     display_mode = "IRA_only") {
  
  log_message("Creating karyotype data for hotspot visualization")
  log_message(sprintf("Display mode: %s", display_mode))
  
  # V20 ENHANCEMENT: species_id_mapping now contains either custom_label or numeric ID
  # Generate karyotype dataframe using the mapping directly
  karyotype_data <- data.frame(
    Chr = as.character(species_id_mapping),  # Uses custom_label if available
    Start = rep(0, length(species_list)),
    End = numeric(length(species_list)),
    stringsAsFactors = FALSE
  )
  
  # Fill genome lengths based on display mode (consistent with M01)
  for (i in seq_along(species_list)) {
    species_name <- species_list[i]
    species_genome <- genome_regions_raw %>% 
      dplyr::filter(species == species_name)
    
    if (nrow(species_genome) > 0) {
      # Check if this species lacks IR
      is_ir_lacking <- !is.na(species_genome$special_handling[1]) && 
                      species_genome$special_handling[1] == "IR_lacking_genome"
      
      # Use IRA_only mode by default (consistent with M01 chromosome distribution)
      if (display_mode == "IRA_only" && !is_ir_lacking) {
        karyotype_data$End[i] <- species_genome$ssc_end[1]
      } else {
        karyotype_data$End[i] <- species_genome$total_length[1]
      }
    } else {
      log_message(sprintf("WARNING: No genome data found for %s", species_name), 
                 level = "warning")
      karyotype_data$End[i] <- 150000  # Fallback
    }
  }
  
  log_message(sprintf("Created karyotype for %d species", nrow(karyotype_data)))
  log_message(sprintf("  - Genome length range: %d - %d bp", 
                     min(karyotype_data$End), 
                     max(karyotype_data$End)))
  
  return(karyotype_data)
}


#' Plot Hotspot Ideogram
#'
#' @description Generates the actual ideogram plot using RIdeogram with overlaid heatmap
#' @param karyotype_data Karyotype data frame
#' @param regions_heatmap CDS and IGS regions heatmap data
#' @param hotspot_markers Hotspot marker data frame
#' @param output_dir Output directory
#' @param config Configuration object
#' @return Path to generated plot
#' @keywords internal
plot_hotspot_ideogram <- function(karyotype_data, regions_heatmap, hotspot_markers, output_dir, config) {
  
  log_message("Generating RIdeogram plot")
  
  # Prepare output paths
  output_base <- file.path(output_dir, "M05_chromosome_hotspot_ideogram")
  svg_path <- paste0(output_base, ".svg")
  pdf_path <- paste0(output_base, ".pdf")
  
  # Call RIdeogram
  tryCatch({
    
    log_message("Calling RIdeogram::ideogram()")
    log_message(sprintf("  - Karyotype: %d chromosomes", nrow(karyotype_data)))
    log_message(sprintf("  - Regions heatmap: %d regions", nrow(regions_heatmap)))
    log_message(sprintf("  - Markers: %d hotspots", nrow(hotspot_markers)))
    
    # Get colors from config (var_type palette)
    color_palette <- get_color_palette(config, "var_type")
    cds_color <- color_palette["snp"]  # Red for CDS
    igs_color <- color_palette["INDEL"]  # Blue for IGS
    
    # Create diverging color palette: Blue (-1) -> White (0) -> Red (1)
    diverging_palette <- c(igs_color, "#FFFFFF", cds_color)
    
    log_message(sprintf("Using diverging color palette: %s (IGS) -> White -> %s (CDS)",
                       igs_color, cds_color))
    
    # RIdeogram call with overlaid heatmap + marker label
    # Hide automatic legends to avoid overlap
    RIdeogram::ideogram(
      karyotype = karyotype_data,
      overlaid = regions_heatmap,
      label = hotspot_markers,
      label_type = "marker",
      colorset1 = diverging_palette,  # Diverging palette for overlaid
      output = svg_path,
      colorset2 = NULL,  # Disable marker legend
      width = 170,       # SVG width
      Lx = 160,          # Legend x position (move far right to hide)
      Ly = 10            # Legend y position
    )
    
    log_message(sprintf("SVG generated: %s", svg_path))
    
    # Convert to PDF with error handling
    log_message("Converting SVG to PDF")
    
    pdf_success <- tryCatch({
      # CRITICAL: Change working directory to ensure PDF is saved in output_dir
      # RIdeogram's convertSVG saves to current working directory by default
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)  # Ensure restoration even on error
      
      setwd(output_dir)
      log_message(sprintf("Changed working directory to: %s", output_dir))
      
      RIdeogram::convertSVG(
        svg = basename(svg_path),  # Use basename only when in output_dir
        device = "pdf",
        width = 10,
        height = max(8, nrow(karyotype_data) * 0.5)  # Dynamic height
      )
      
      # Restore working directory
      setwd(old_wd)
      
      # Verify PDF was created
      if (file.exists(pdf_path)) {
        log_message(sprintf("PDF generated successfully: %s", pdf_path))
        
        # Clean up SVG
        if (file.exists(svg_path)) {
          file.remove(svg_path)
          log_message("Cleaned up intermediate SVG file")
        }
        
        TRUE
      } else {
        log_message("WARNING: convertSVG completed but PDF file not found", level = "warning")
        FALSE
      }
      
    }, error = function(e) {
      log_message(sprintf("WARNING: PDF conversion failed: %s", e$message), level = "warning")
      log_message("SVG file will be kept as fallback", level = "warning")
      FALSE
    })
    
    # Return appropriate path
    if (pdf_success) {
      return(pdf_path)
    } else {
      log_message("Using SVG as final output (PDF conversion failed)")
      return(svg_path)
    }
    
  }, error = function(e) {
    log_message(sprintf("ERROR in RIdeogram plotting: %s", e$message), level = "error")
    stop(sprintf("Failed to generate ideogram plot: %s", e$message))
  })
}


#' Write Hotspot Summary Statistics
#'
#' @description Generates a text summary of hotspot statistics
#' @param cds_hotspots CDS hotspot data
#' @param poigs_hotspots poiGS hotspot data
#' @param hotspot_markers Combined marker data
#' @param output_path Output file path
#' @keywords internal
write_hotspot_summary <- function(cds_hotspots, poigs_hotspots, hotspot_markers, output_path) {
  
  log_message("Generating hotspot summary statistics")
  
  # Calculate statistics
  n_species <- length(unique(c(cds_hotspots$species, poigs_hotspots$species)))
  n_cds_hotspots <- nrow(cds_hotspots)
  n_poigs_hotspots <- nrow(poigs_hotspots)
  n_total_markers <- nrow(hotspot_markers)
  
  # Per-species statistics
  species_stats <- hotspot_markers %>%
    dplyr::group_by(species, Type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Type, values_from = count, values_fill = 0)
  
  # Write summary
  sink(output_path)
  cat("============================================================\n")
  cat("CHROMOSOME HOTSPOT MAP SUMMARY\n")
  cat("============================================================\n\n")
  cat(sprintf("Generated: %s\n", Sys.time()))
  cat(sprintf("Module: M01-2 Structural & Relational Analysis\n\n"))
  
  cat("--- OVERALL STATISTICS ---\n")
  cat(sprintf("Total species: %d\n", n_species))
  cat(sprintf("Total CDS hotspots: %d\n", n_cds_hotspots))
  cat(sprintf("Total poiGS hotspots: %d\n", n_poigs_hotspots))
  # Note: IR markers are structural features, not hotspots
  cat(sprintf("Total markers plotted: %d (CDS + poiGS only)\n\n", n_cds_hotspots + n_poigs_hotspots))
  
  cat("--- PER-SPECIES BREAKDOWN ---\n")
  print(species_stats)
  cat("\n")
  
  cat("--- GENOME REGION DISTRIBUTION ---\n")
  region_stats <- hotspot_markers %>%
    dplyr::group_by(genome_region, Type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Type, values_from = count, values_fill = 0)
  print(region_stats)
  cat("\n")
  
  cat("--- COLOR SCHEME ---\n")
  cat(sprintf("CDS Hotspots: %s (red gradient)\n", 
             unique(hotspot_markers$color[hotspot_markers$Type == "CDS_Hotspot"])[1]))
  cat(sprintf("poiGS Hotspots: %s (blue gradient)\n", 
             unique(hotspot_markers$color[hotspot_markers$Type == "poiGS_Hotspot"])[1]))
  cat("\nNote: Marker color intensity represents species sharing degree.\n")
  cat("      Light colors = Low sharing, Dark colors = High sharing\n")
  cat("      IR markers (orange triangles) are structural features, not hotspots.\n\n")
  
  cat("============================================================\n")
  cat("For manual labeling, refer to: M05_hotspot_markers_sorted.csv\n")
  cat("============================================================\n")
  sink()
  
  log_message("Summary statistics written successfully")
}


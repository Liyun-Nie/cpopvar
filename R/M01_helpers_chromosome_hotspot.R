# ============================================================
# M01 Chromosome Hotspot Map Helper Functions
# ============================================================
#
# This file contains helper functions for M01-2: Structural & Relational Analysis
# Chromosome hotspot map visualization helpers

#' Extract CDS Hotspot Coordinates
#' 
#' @description Extracts CDS hotspot coordinates from M03 output and region_info_complete.csv
#' @param m03_hotspots_path Path to M03_core_shared_hotspots.csv
#' @param region_info_path Path to region_info_complete.csv
#' @param session_id Session identifier
#' @return Data frame with CDS hotspot coordinates
#' @keywords internal
extract_cds_hotspot_coordinates <- function(m03_hotspots_path, region_info_path, session_id) {
  
  log_message("=== EXTRACTING CDS HOTSPOT COORDINATES ===")
  
  # Load M03 hotspots
  if (!file.exists(m03_hotspots_path)) {
    stop(sprintf("M03 hotspots file not found: %s", m03_hotspots_path))
  }
  
  m03_hotspots <- readr::read_csv(m03_hotspots_path, show_col_types = FALSE)
  log_message(sprintf("Loaded M03 hotspots: %d unique genes", nrow(m03_hotspots)))
  log_message(sprintf("  - Total species coverage: %d", sum(m03_hotspots$num_species)))
  
  # Load region_info_complete.csv
  if (!file.exists(region_info_path)) {
    stop(sprintf("region_info_complete.csv not found: %s", region_info_path))
  }
  
  region_info <- readr::read_csv(region_info_path, show_col_types = FALSE)
  log_message(sprintf("Loaded region_info: %d records", nrow(region_info)))
  
  # CRITICAL FIX: Expand hotspots from gene-level to species-gene-level
  # M03 data structure: one row per hotspot gene, with species_list as comma-separated string
  log_message("Expanding hotspots to species-level...")
  
  m03_expanded <- m03_hotspots %>%
    dplyr::mutate(
      # Split species_list string into vector
      species = strsplit(species_list, ", ")
    ) %>%
    tidyr::unnest(species) %>%
    dplyr::mutate(
      # Clean species names (remove whitespace)
      species = trimws(species)
    )
  
  log_message(sprintf("Expanded to species-level: %d records across %d species", 
                     nrow(m03_expanded), 
                     length(unique(m03_expanded$species))))
  
  # Join M03 hotspots with region_info to get coordinates
  # Now we have 'species' column for proper join
  cds_hotspots_with_coords <- m03_expanded %>%
    dplyr::left_join(
      region_info %>% dplyr::filter(region_type == "CDS"),
      by = c("species", "gene"),
      relationship = "many-to-many"  # One gene may have multiple exons
    )
  
  if (nrow(cds_hotspots_with_coords) == 0) {
    stop("No CDS hotspots could be matched with region_info")
  }
  
  log_message(sprintf("Matched M03 hotspots with region_info: %d exon-level records", 
                     nrow(cds_hotspots_with_coords)))
  
  # Keep exon-level coordinates for segment display
  # Each exon will be displayed as a separate segment in the heatmap
  cds_hotspot_coords <- cds_hotspots_with_coords %>%
    dplyr::mutate(
      position = region_start,  # Keep original position for markers
      hotspot_type = "M03_CDS",
      marker_shape = "circle"  # circle is clearer than triangle at this scale
    ) %>%
    dplyr::select(species, gene, num_species, mean_frequency, max_frequency,
                 position, region_start, region_end, genome_region, region_type,
                 hotspot_type, marker_shape)
  
  log_message(sprintf("Prepared CDS hotspot data: %d exon segments", nrow(cds_hotspot_coords)))
  log_message(sprintf("  - Species: %d", length(unique(cds_hotspot_coords$species))))
  log_message(sprintf("  - Unique genes: %d", length(unique(cds_hotspot_coords$gene))))
  log_message(sprintf("  - Species sharing range: %d-%d", 
                     min(cds_hotspot_coords$num_species), 
                     max(cds_hotspot_coords$num_species)))
  
  # Validate coordinates
  invalid_coords <- cds_hotspot_coords %>% 
    dplyr::filter(is.na(position) | position <= 0)
  
  if (nrow(invalid_coords) > 0) {
    log_message(sprintf("WARNING: Found %d CDS hotspots with invalid coordinates", 
                       nrow(invalid_coords)), level = "warning")
  }
  
  return(cds_hotspot_coords)
}


#' Extract poiGS Hotspot Coordinates
#' 
#' @description Extracts poiGS hotspot coordinates from M04 output and region_info_complete.csv
#' @param m04_hotspots_path Path to M04_core_shared_poigs_hotspots.csv
#' @param region_info_path Path to region_info_complete.csv
#' @param session_id Session identifier
#' @return Data frame with poiGS hotspot coordinates
#' @keywords internal
extract_poigs_hotspot_coordinates <- function(m04_hotspots_path, region_info_path, session_id) {
  
  log_message("=== EXTRACTING POIGS HOTSPOT COORDINATES ===")
  
  # Load M04 hotspots
  if (!file.exists(m04_hotspots_path)) {
    stop(sprintf("M04 hotspots file not found: %s", m04_hotspots_path))
  }
  
  m04_hotspots <- readr::read_csv(m04_hotspots_path, show_col_types = FALSE)
  log_message(sprintf("Loaded M04 hotspots: %d unique poiGS regions", nrow(m04_hotspots)))
  log_message(sprintf("  - Total species coverage: %d", sum(m04_hotspots$num_species)))
  
  # Load region_info_complete.csv
  if (!file.exists(region_info_path)) {
    stop(sprintf("region_info_complete.csv not found: %s", region_info_path))
  }
  
  region_info <- readr::read_csv(region_info_path, show_col_types = FALSE)
  
  # CRITICAL FIX: Expand hotspots from poiGS-level to species-poiGS-level
  log_message("Expanding poiGS hotspots to species-level...")
  
  m04_expanded <- m04_hotspots %>%
    dplyr::mutate(
      # Split species_list string into vector
      species = strsplit(species_list, ", ")
    ) %>%
    tidyr::unnest(species) %>%
    dplyr::mutate(
      # Clean species names (remove whitespace)
      species = trimws(species),
      # Extract Region_Name from poiGS_ID (remove "poiGS_" prefix)
      # poiGS_ID format: "poiGS_gene1-gene2" → Region_Name: "IGS_gene1-gene2"
      Region_Name = paste0("IGS_", sub("^poiGS_", "", poiGS_ID))
    )
  
  log_message(sprintf("Expanded to species-level: %d records across %d species", 
                     nrow(m04_expanded), 
                     length(unique(m04_expanded$species))))
  
  # Join M04 hotspots with region_info to get coordinates
  poigs_hotspots_with_coords <- m04_expanded %>%
    dplyr::left_join(
      region_info %>% dplyr::filter(region_type == "IGS"),
      by = c("species", "Region_Name" = "region_name")
    )
  
  if (nrow(poigs_hotspots_with_coords) == 0) {
    stop("No poiGS hotspots could be matched with region_info")
  }
  
  log_message(sprintf("Matched M04 hotspots with region_info: %d IGS records", 
                     nrow(poigs_hotspots_with_coords)))
  
  # Keep region_start/region_end for segment display
  poigs_hotspot_coords <- poigs_hotspots_with_coords %>%
    dplyr::mutate(
      position = region_start,  # Keep for markers
      region_name = Region_Name,
      hotspot_type = "M04_poiGS",
      marker_shape = "circle"  # circle is clearer than triangle at this scale
    ) %>%
    dplyr::select(
      species, 
      region_name,
      poiGS_ID,
      num_species,
      mean_frequency,
      max_frequency,
      genome_region,
      position,
      region_start,
      region_end,
      region_type,
      hotspot_type,
      marker_shape
    )
  
  log_message(sprintf("Extracted poiGS coordinates: %d IGS segments", nrow(poigs_hotspot_coords)))
  log_message(sprintf("  - Species: %d", length(unique(poigs_hotspot_coords$species))))
  log_message(sprintf("  - Unique poiGS regions: %d", length(unique(poigs_hotspot_coords$poiGS_ID))))
  log_message(sprintf("  - Species sharing range: %d-%d", 
                     min(poigs_hotspot_coords$num_species), 
                     max(poigs_hotspot_coords$num_species)))
  
  # Validate coordinates
  invalid_coords <- poigs_hotspot_coords %>% 
    dplyr::filter(is.na(position) | position <= 0)
  
  if (nrow(invalid_coords) > 0) {
    log_message(sprintf("WARNING: Found %d poiGS hotspots with invalid coordinates", 
                       nrow(invalid_coords)), level = "warning")
  }
  
  return(poigs_hotspot_coords)
}


#' Get Hotspot Colors from Configuration
#' 
#' @description Retrieves color palette for CDS and poiGS hotspots from config
#' @param config Configuration object
#' @return Named vector with 'cds' and 'poigs' colors
#' @keywords internal
get_hotspot_colors_from_config <- function(config) {
  
  log_message("=== RETRIEVING HOTSPOT COLORS FROM CONFIG ===")
  
  # Use the var_type color palette
  # CDS hotspots: SNP color (red)
  # poiGS hotspots: INDEL color (blue)
  
  tryCatch({
    # Get var_type color palette
    var_type_palette <- get_color_palette(config, "var_type")
    
    if (is.null(var_type_palette)) {
      stop("var_type color palette not found in config")
    }
    
    # Extract colors (use [[]] to get pure strings without names)
    cds_color <- var_type_palette[["snp"]]
    poigs_color <- var_type_palette[["INDEL"]]
    
    if (is.null(cds_color) || is.null(poigs_color) || is.na(cds_color) || is.na(poigs_color)) {
      stop("Required colors (snp, INDEL) not found in var_type palette")
    }
    
    log_message(sprintf("Retrieved colors from config:"))
    log_message(sprintf("  - CDS hotspots (SNP): %s", cds_color))
    log_message(sprintf("  - poiGS hotspots (INDEL): %s", poigs_color))
    
    return(c(cds = cds_color, poigs = poigs_color))
    
  }, error = function(e) {
    log_message(sprintf("ERROR retrieving colors from config: %s", e$message), level = "error")
    log_message("Using fallback colors", level = "warning")
    
    # Fallback colors
    return(c(cds = "#CE4545", poigs = "#3e9faf"))
  })
}


#' Generate Figure Legend for Chromosome Hotspot Ideogram
#'
#' @description Generates a text-based figure legend explaining the color scheme
#'   and marker meanings in the chromosome hotspot ideogram.
#'   
#'   Marker styling:
#'   - Replaces the ambiguous RIdeogram auto-legend
#'   - Provides a clear biological interpretation
#'   - Explains color gradients for species sharing
#'
#' @param output_path Path to save the legend text file
#' @param cds_hotspots CDS hotspot data frame
#' @param poigs_hotspots poiGS hotspot data frame
#' @param hotspot_colors Named vector of base colors
#' @keywords internal
generate_figure_legend <- function(output_path, cds_hotspots, poigs_hotspots, hotspot_colors) {
  
  # Calculate sharing statistics
  cds_min <- min(cds_hotspots$num_species, na.rm = TRUE)
  cds_max <- max(cds_hotspots$num_species, na.rm = TRUE)
  poigs_min <- min(poigs_hotspots$num_species, na.rm = TRUE)
  poigs_max <- max(poigs_hotspots$num_species, na.rm = TRUE)
  
  # Extract base colors
  cds_base <- hotspot_colors[["cds"]]
  poigs_base <- hotspot_colors[["poigs"]]
  
  # Write legend in journal article format
  sink(output_path)
  cat("FIGURE LEGEND\n\n")
  
  cat("Chromosome-level distribution of variant hotspots across 17 species. ")
  cat("Each vertical bar represents the chloroplast genome of one species (IRA region only, LSC-IRA-SSC). ")
  cat(sprintf("Red regions indicate CDS (coding sequence) hotspots (n=%d), ", 
             nrow(cds_hotspots)))
  cat(sprintf("while blue regions mark IGS (intergenic spacer) hotspots (n=%d). ", 
             nrow(poigs_hotspots)))
  cat("Circular markers denote the precise genomic positions of hotspots, ")
  cat("with color intensity reflecting species sharing degree: ")
  cat(sprintf("light colors represent low sharing (%d-%d species), ", 
             min(cds_min, poigs_min), max(cds_min, poigs_min) + 1))
  cat(sprintf("while dark colors indicate high sharing (%d-%d species). ", 
             min(cds_max, poigs_max) - 1, max(cds_max, poigs_max)))
  cat("Orange triangles mark the boundaries of Inverted Repeat (IR) regions. ")
  cat("Hotspots were identified as genomic regions with variant frequencies exceeding the 75th percentile. ")
  cat(sprintf("CDS hotspots: %d genes shared by %d-%d species. ", 
             length(unique(cds_hotspots$gene)), cds_min, cds_max))
  cat(sprintf("IGS hotspots: %d position-based orthologous intergenic spacers (poiGS) shared by %d-%d species.\n", 
             length(unique(poigs_hotspots$poiGS_ID)), poigs_min, poigs_max))
  sink()
  
  log_message(sprintf("Figure legend generated: %s", output_path))
}


#' Generate Dynamic Color Gradient
#'
#' @description Generates a dynamic color gradient based on actual species sharing range.
#'   Automatically adjusts the number of gradient levels based on the data range.
#'   
#'   Rules:
#'   - Range <= 5 species: Use exact range (e.g., 3-5 species → 3 levels)
#'   - Range 6-10 species: Use 6 levels (balanced)
#'   - Range > 10 species: Use min(10, ceiling(range/2)) levels (avoid over-density)
#'
#' @param base_color Base color (6-digit hex, no prefix) representing high sharing
#' @param min_shared Minimum number of shared species
#' @param max_shared Maximum number of shared species
#' @param n_levels Number of gradient levels (auto-calculated if NULL)
#' @return Vector of color codes (6-digit hex, no prefix)
#' @keywords internal
generate_color_gradient <- function(base_color, min_shared, max_shared, n_levels = NULL) {
  
  # Auto-calculate gradient levels based on sharing range
  if (is.null(n_levels)) {
    sharing_range <- max_shared - min_shared + 1
    
    n_levels <- if (sharing_range <= 5) {
      sharing_range  # Exact match for small ranges
    } else if (sharing_range <= 10) {
      6  # Balanced for medium ranges
    } else {
      min(10, ceiling(sharing_range / 2))  # Avoid over-density for large ranges
    }
  }
  
  log_message(sprintf("Generating %d-level color gradient for sharing range: %d-%d species",
                     n_levels, min_shared, max_shared))
  
  # Convert base_color to RGB
  base_rgb <- col2rgb(paste0("#", base_color))
  
  # Generate light color by interpolating towards white (255, 255, 255)
  # Use 70% interpolation to ensure sufficient contrast
  light_rgb <- base_rgb + (255 - base_rgb) * 0.7
  light_color <- rgb(light_rgb[1], light_rgb[2], light_rgb[3], maxColorValue = 255)
  dark_color <- paste0("#", base_color)
  
  # Generate smooth gradient using colorRampPalette
  gradient_colors <- colorRampPalette(c(light_color, dark_color))(n_levels)
  
  # Remove # prefix for RIdeogram compatibility
  gradient_colors <- substr(gradient_colors, 2, 7)
  
  log_message(sprintf("  Gradient: %s (light, low sharing) -> %s (dark, high sharing)", 
                     gradient_colors[1], gradient_colors[n_levels]))
  
  return(gradient_colors)
}


#' Prepare Hotspot Markers for RIdeogram
#' 
#' @description Combines CDS and poiGS hotspot data into RIdeogram marker format
#'   with color gradients to represent species sharing degree.
#'   
#'   Color and shape:
#'   - Color format: Pure 6-digit hex (no # prefix, no alpha) for RIdeogram compatibility
#'   - Color gradient: Light→Dark based on species sharing degree
#'   - Shape: circle only (simplified from box/circle system)
#'   
#' @param cds_hotspots Data frame with CDS hotspot coordinates
#' @param poigs_hotspots Data frame with poiGS hotspot coordinates
#' @param species_id_mapping Named vector mapping species names to chromosome IDs
#' @param hotspot_colors Named vector with 'cds' and 'poigs' colors
#' @return Data frame in RIdeogram marker format
#' @keywords internal
prepare_hotspot_markers <- function(cds_hotspots, poigs_hotspots, species_id_mapping, hotspot_colors) {
  
  log_message("=== PREPARING HOTSPOT MARKERS FOR RIDEOGRAM ===")
  
  # Extract base colors (remove # prefix for RIdeogram compatibility)
  cds_base <- substr(hotspot_colors[["cds"]], 2, 7)  # CE4545
  poigs_base <- substr(hotspot_colors[["poigs"]], 2, 7)  # 3e9faf
  
  log_message(sprintf("Base colors (6-digit hex, no prefix):"))
  log_message(sprintf("  - CDS base: %s", cds_base))
  log_message(sprintf("  - poiGS base: %s", poigs_base))
  
  # Calculate sharing statistics for CDS
  cds_max_shared <- max(cds_hotspots$num_species, na.rm = TRUE)
  cds_min_shared <- min(cds_hotspots$num_species, na.rm = TRUE)
  
  log_message(sprintf("CDS sharing range: %d-%d species", cds_min_shared, cds_max_shared))
  
  # Generate a dynamic color gradient for CDS
  cds_gradient <- generate_color_gradient(
    base_color = cds_base,
    min_shared = cds_min_shared,
    max_shared = cds_max_shared
  )
  
  # Prepare CDS markers with color gradient
  cds_markers <- cds_hotspots %>%
    dplyr::mutate(
      Chr = as.character(species_id_mapping[species]),
      Start = position,
      End = position,
      Type = "CDS_Hotspot",
      # Shape: circle only (simplified)
      Shape = "circle",
      # Map num_species to color gradient index
      color_index = pmin(pmax(num_species - cds_min_shared + 1, 1), length(cds_gradient)),
      color = cds_gradient[color_index],
      label = gene,
      num_species_shared = num_species,
      hotspot_type = "CDS",
      genome_region = "CDS"
    ) %>%
    dplyr::filter(!is.na(Chr) & !is.na(Start))
  
  log_message(sprintf("Prepared CDS markers: %d records", nrow(cds_markers)))
  log_message(sprintf("  - Color distribution:"))
  color_dist_cds <- table(cds_markers$color)
  for (i in seq_along(color_dist_cds)) {
    log_message(sprintf("    %s: %d markers", names(color_dist_cds)[i], color_dist_cds[i]))
  }
  
  # Calculate sharing statistics for poiGS
  poigs_max_shared <- max(poigs_hotspots$num_species, na.rm = TRUE)
  poigs_min_shared <- min(poigs_hotspots$num_species, na.rm = TRUE)
  
  log_message(sprintf("poiGS sharing range: %d-%d species", poigs_min_shared, poigs_max_shared))
  
  # Generate a dynamic color gradient for poiGS
  poigs_gradient <- generate_color_gradient(
    base_color = poigs_base,
    min_shared = poigs_min_shared,
    max_shared = poigs_max_shared
  )
  
  # Prepare poiGS markers with color gradient
  poigs_markers <- poigs_hotspots %>%
    dplyr::mutate(
      Chr = as.character(species_id_mapping[species]),
      Start = position,
      End = position,
      Type = "poiGS_Hotspot",
      # Shape: circle only
      Shape = "circle",
      # Map num_species to color gradient index
      color_index = pmin(pmax(num_species - poigs_min_shared + 1, 1), length(poigs_gradient)),
      color = poigs_gradient[color_index],
      label = region_name,
      num_species_shared = num_species,
      hotspot_type = "poiGS",
      genome_region = "IGS"
    ) %>%
    dplyr::filter(!is.na(Chr) & !is.na(Start))
  
  log_message(sprintf("Prepared poiGS markers: %d records", nrow(poigs_markers)))
  log_message(sprintf("  - Color distribution:"))
  color_dist_poigs <- table(poigs_markers$color)
  for (i in seq_along(color_dist_poigs)) {
    log_message(sprintf("    %s: %d markers", names(color_dist_poigs)[i], color_dist_poigs[i]))
  }
  
  # Combine markers
  all_markers <- dplyr::bind_rows(
    cds_markers %>% 
      dplyr::select(Type, Shape, Chr, Start, End, color, species, label, hotspot_type, genome_region, num_species_shared),
    poigs_markers %>% 
      dplyr::select(Type, Shape, Chr, Start, End, color, species, label, hotspot_type, genome_region, num_species_shared)
  )
  
  log_message(sprintf("Combined markers: %d total records", nrow(all_markers)))
  log_message(sprintf("  - CDS hotspots: %d", sum(all_markers$Type == "CDS_Hotspot")))
  log_message(sprintf("  - poiGS hotspots: %d", sum(all_markers$Type == "poiGS_Hotspot")))
  log_message(sprintf("  - Species coverage: %d/%d", 
                     length(unique(all_markers$species)), 
                     length(species_id_mapping)))
  log_message(sprintf("  - Species sharing range: %d-%d species per hotspot", 
                     min(all_markers$num_species_shared), 
                     max(all_markers$num_species_shared)))
  
  # Validate marker data
  if (nrow(all_markers) == 0) {
    stop("No valid markers generated")
  }
  
  # Check for coordinate validity
  invalid_markers <- all_markers %>%
    dplyr::filter(is.na(Start) | Start <= 0 | is.na(End) | End <= 0)
  
  if (nrow(invalid_markers) > 0) {
    log_message(sprintf("WARNING: Found %d markers with invalid coordinates", 
                       nrow(invalid_markers)), level = "warning")
  }
  
  return(all_markers)
}


#' Prepare CDS and IGS Regions Heatmap Data
#'
#' @description Extracts all CDS and IGS regions from region_info_complete.csv
#'   and prepares them for RIdeogram overlaid heatmap visualization.
#'   Uses a diverging value scale: CDS = 1 (red), IGS = -1 (blue), dummy = 0 (white)
#' @param region_info_path Path to region_info_complete.csv
#' @param genome_regions Genome regions data with Chr mapping
#' @return Data frame with columns: Chr, Start, End, Value
#' @keywords internal
prepare_regions_heatmap <- function(region_info_path, genome_regions) {
  
  log_message("=== PREPARING CDS AND IGS REGIONS HEATMAP ===")
  
  # Load region_info
  region_info <- readr::read_csv(region_info_path, show_col_types = FALSE)
  
  log_message(sprintf("Loaded region_info: %d total regions", nrow(region_info)))
  
  # Filter for CDS and IGS only
  regions <- region_info %>%
    dplyr::filter(region_type %in% c("CDS", "IGS")) %>%
    dplyr::select(species, region_name, region_type, region_start, region_end)
  
  log_message(sprintf("Filtered regions: %d CDS + %d IGS",
                     sum(regions$region_type == "CDS"),
                     sum(regions$region_type == "IGS")))
  
  # Merge with genome_regions to get Chr
  regions_with_chr <- regions %>%
    dplyr::left_join(
      genome_regions %>% dplyr::select(species, Chr),
      by = "species"
    )
  
  # Assign Value based on region_type
  # CDS = 1 (red end), IGS = -1 (blue end)
  regions_heatmap <- regions_with_chr %>%
    dplyr::mutate(
      Value = dplyr::case_when(
        region_type == "CDS" ~ 1,
        region_type == "IGS" ~ -1,
        TRUE ~ 0
      )
    ) %>%
    dplyr::select(Chr, Start = region_start, End = region_end, Value) %>%
    dplyr::arrange(Chr, Start)
  
  # Add dummy region with Value = 0 to establish color range
  # Use the first chromosome's first position
  first_chr <- regions_heatmap$Chr[1]
  dummy_region <- data.frame(
    Chr = first_chr,
    Start = 1,
    End = 1,
    Value = 0
  )
  
  regions_heatmap <- dplyr::bind_rows(dummy_region, regions_heatmap)
  
  log_message(sprintf("Prepared heatmap data: %d regions (including 1 dummy)", 
                     nrow(regions_heatmap)))
  log_message(sprintf("  - Value range: %.1f to %.1f", 
                     min(regions_heatmap$Value), 
                     max(regions_heatmap$Value)))
  log_message(sprintf("  - CDS regions (Value=1): %d", 
                     sum(regions_heatmap$Value == 1)))
  log_message(sprintf("  - IGS regions (Value=-1): %d", 
                     sum(regions_heatmap$Value == -1)))
  
  return(regions_heatmap)
}


#' Prepare Hotspot Regions Heatmap Data
#'
#' @description Creates heatmap data for ONLY hotspot regions (CDS hotspots and poiGS hotspots)
#'   with simplified value assignments. This differs from prepare_regions_heatmap() which shows
#'   ALL CDS/IGS regions. Only hotspot regions are displayed.
#'   
#'   Value assignments (simplified, no sharing degree encoding):
#'   - CDS hotspots = 1 (red)
#'   - poiGS hotspots = -1 (blue)
#'   - Dummy region = 0 (white, establishes color range)
#'
#' @param cds_hotspots Data frame with CDS hotspot coordinates (from extract_cds_hotspot_coordinates)
#' @param poigs_hotspots Data frame with poiGS hotspot coordinates (from extract_poigs_hotspot_coordinates)
#' @param genome_regions Genome regions data with Chr mapping
#' @return Data frame with columns: Chr, Start, End, Value
#' @keywords internal
prepare_hotspot_regions_heatmap <- function(cds_hotspots, poigs_hotspots, genome_regions) {
  
  log_message("=== PREPARING HOTSPOT REGIONS HEATMAP (CDS + poiGS Hotspots Only) ===")
  
  # Process CDS hotspots: Value = 1 (red)
  # Use region_start/region_end for segment display
  cds_heatmap <- cds_hotspots %>%
    dplyr::left_join(
      genome_regions %>% dplyr::select(species, Chr),
      by = "species"
    ) %>%
    dplyr::mutate(Value = 1) %>%
    dplyr::select(Chr, Start = region_start, End = region_end, Value) %>%
    dplyr::filter(!is.na(Chr), !is.na(Start), !is.na(End))
  
  log_message(sprintf("Processed CDS hotspots: %d exon segments", nrow(cds_heatmap)))
  
  # Process poiGS hotspots: Value = -1 (blue)
  # Use region_start/region_end for segment display
  poigs_heatmap <- poigs_hotspots %>%
    dplyr::left_join(
      genome_regions %>% dplyr::select(species, Chr),
      by = "species"
    ) %>%
    dplyr::mutate(Value = -1) %>%
    dplyr::select(Chr, Start = region_start, End = region_end, Value) %>%
    dplyr::filter(!is.na(Chr), !is.na(Start), !is.na(End))
  
  log_message(sprintf("Processed poiGS hotspots: %d IGS segments", nrow(poigs_heatmap)))
  
  # Combine hotspots
  regions_heatmap <- dplyr::bind_rows(cds_heatmap, poigs_heatmap) %>%
    dplyr::arrange(Chr, Start)
  
  # Add dummy region with Value = 0 to establish color range
  first_chr <- regions_heatmap$Chr[1]
  dummy_region <- data.frame(
    Chr = first_chr,
    Start = 1,
    End = 1,
    Value = 0
  )
  
  regions_heatmap <- dplyr::bind_rows(dummy_region, regions_heatmap)
  
  log_message(sprintf("Prepared hotspot heatmap data: %d hotspots (including 1 dummy)",
                     nrow(regions_heatmap)))
  log_message(sprintf("  - Value range: %.1f to %.1f",
                     min(regions_heatmap$Value),
                     max(regions_heatmap$Value)))
  log_message(sprintf("  - CDS hotspots (Value=1): %d",
                     sum(regions_heatmap$Value == 1)))
  log_message(sprintf("  - poiGS hotspots (Value=-1): %d",
                     sum(regions_heatmap$Value == -1)))
  
  return(regions_heatmap)
}

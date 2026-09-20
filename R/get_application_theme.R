#' Get application theme
#' @param theme_name Name of the theme to apply
#' @param custom_theme_file Path to custom theme file
#' @param base_size Base font size
#' @param base_family Base font family
#' @return ggplot2 theme object
#' @importFrom ggplot2 theme_minimal theme element_text element_rect element_blank margin unit theme_bw theme_classic
#' @export
get_application_theme <- function(theme_name = "professional_light", 
                                 custom_theme_file = NULL, 
                                 base_size = 11, 
                                 base_family = "Arial") {
  
  # Font fallback mechanism for cross-platform compatibility
  if (base_family == "Arial") {
    # Simple fallback strategy: use cross-platform compatible font names
    # This avoids dependency on extrafont package and works more reliably
    if (Sys.info()["sysname"] == "Windows") {
      base_family <- "Arial"  # Arial usually available on Windows
    } else if (Sys.info()["sysname"] == "Darwin") {
      base_family <- "Arial"  # Arial usually available on macOS
    } else {
      # Linux systems may not have Arial, use universally available fonts
      base_family <- "sans"  # Generic sans-serif, always available
      log_message("Using 'sans' font family for better Linux compatibility.", level = "info")
    }
  }
  
  tryCatch({
    
    # If custom theme file is provided, try to load it
    if (!is.null(custom_theme_file) && file.exists(custom_theme_file)) {
      tryCatch({
        source(custom_theme_file)
        if (exists("get_custom_theme", where = parent.frame())) {
          custom_theme <- get("get_custom_theme", mode = "function")()
          if (inherits(custom_theme, "theme")) {
            return(custom_theme)
          }
        }
        warning("Custom theme file loaded but get_custom_theme() function not found or invalid")
      }, error = function(e) {
        warning(sprintf("Failed to load custom theme: %s", e$message))
      })
    }
    
    # Built-in themes
    if (theme_name == "professional_light") {
      # Enhanced professional theme based on mytheme with optimizations
      # - Uses theme_classic() for clean, publication-ready base
      # - Arial font as default (set in function signature)
      # - Adds black border as requested
      # - Maintains mytheme's precision while adding practical enhancements
      return(
  theme_classic(base_size = base_size, base_family = base_family) +
ggplot2::theme(
    # Global text settings with EXPLICIT size calculation
    text =ggplot2::element_text(colour = "black", size = base_size),
    
    # Plot title and subtitle
    plot.title =ggplot2::element_text(size = base_size * 1.2, hjust = 0.5, face = "bold", 
                             margin = margin(b = 20)),
    plot.subtitle =ggplot2::element_text(size = base_size * 1.0, hjust = 0.5, 
                                margin = margin(b = 15)),
    
    # Axis settings with EXPLICIT size calculation
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    axis.title.x =ggplot2::element_text(size = base_size * 1.0, margin = margin(t = 15)),
    axis.title.y =ggplot2::element_text(size = base_size * 1.0, margin = margin(r = 15)),
    axis.text.x =ggplot2::element_text(size = base_size * 0.8, angle = 0),
    axis.text.y =ggplot2::element_text(size = base_size * 0.8),
    
    # Legend with EXPLICIT size calculation
    legend.title =ggplot2::element_text(colour = "black", size = base_size * 1.0, face = "bold"),
    legend.text =ggplot2::element_text(colour = "black", size = base_size * 0.8),
    legend.key.size = unit(4, units = "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 15),
    
    # Panel and border (key user requirement - black border)
    panel.border =ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.background =ggplot2::element_rect(fill = "white", colour = NA),
    
    # Subtle grid for enhanced data readability (compromise between classic and practical)
    panel.grid.major = element_line(color = "grey95", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    
    # Strip (for faceted plots)
    strip.background =ggplot2::element_rect(fill = "grey85", color = "black"),
    strip.text =ggplot2::element_text(size = base_size * 0.9, face = "bold"),
    
    # Plot background
    plot.background =ggplot2::element_rect(fill = "white", color = NA),
    
    # Margins (maintaining good spacing)
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )
      )
      
    } else if (theme_name == "minimal") {
      return(
        theme_minimal(base_size = base_size, base_family = base_family) +
ggplot2::theme(
          plot.title =ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.title =ggplot2::element_text(face = "bold"),
          legend.position = "bottom"
        )
      )
      
    } else if (theme_name == "classic") {
      return(
        theme_classic(base_size = base_size, base_family = base_family) +
ggplot2::theme(
          plot.title =ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.title =ggplot2::element_text(face = "bold")
        )
      )
      
    } else {
      warning(sprintf("Unknown theme '%s', using professional_light", theme_name))
      return(get_application_theme("professional_light", NULL, base_size, base_family))
    }
    
  }, error = function(e) {
    warning(sprintf("Theme generation failed: %s", e$message))
    # Ultimate fallback to basic theme
    return(theme_bw(base_size = base_size, base_family = base_family))
  })
}

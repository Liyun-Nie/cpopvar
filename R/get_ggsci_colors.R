#' Get scientific color palette from ggsci
#' 
#' @title Get scientific color palette from ggsci
#' @description Retrieves color palettes from the ggsci package for scientific publications
#' @param palette Name of the ggsci palette to use
#' @param n_colors Number of colors to retrieve (optional)
#' @return Vector of color codes
#' @export
get_ggsci_colors <- function(palette, n_colors = NULL) {
  if (!requireNamespace("ggsci", quietly = TRUE)) {
    return(NULL)
  }
  
  ggsci_colors <- switch(palette,
    "npg" = ggsci::pal_npg("nrc")(10),
    "aaas" = ggsci::pal_aaas("default")(10),
    "lancet" = ggsci::pal_lancet("lanonc")(9),
    "jco" = ggsci::pal_jco("default")(10),
    "ucscgb" = ggsci::pal_ucscgb("default")(26),
    "d3" = ggsci::pal_d3("category10")(10),
    "startrek" = ggsci::pal_startrek("uniform")(7),
    NULL
  )
  
  if (is.null(ggsci_colors)) {
    return(NULL)
  }
  
  if (!is.null(n_colors)) {
    # Extend colors if needed
    if (n_colors > length(ggsci_colors)) {
      ggsci_colors <- rep(ggsci_colors, ceiling(n_colors / length(ggsci_colors)))
    }
    return(ggsci_colors[1:n_colors])
  }
  
  return(ggsci_colors)
}

# /shiny-workflow/src/utils/generate_flowchart_text.R

#' @title Generate a Text-Based Statistics Workflow Flowchart
#' @description Creates and returns a pre-formatted, plain text string
#' representing the decision logic of the M02 statistical engine.
#' @return A single character string containing the ASCII flowchart.
#' generate_statistics_flowchart_text
#' @export
generate_statistics_flowchart_text <- function() {
  
  flowchart_string <- "
M02 Intelligent Statistical Engine Workflow
===========================================

[ Start: Input Data ]
           |
           v
+--------------------------+
|  1. Diagnostic Tests     |
|   - Homogeneity (Levene) |
|   - Normality (Shapiro)  |
+--------------------------+
           |
           v
+---------------------------------------------------------+
| Decision: Meet assumptions for parametric test?         |
| (Data is Normal AND Variances are Homogeneous)          |
+---------------------------------------------------------+
           |                               |
           |--[ No ]---------------------->|
           |                               |
         (Yes)                           v
           |                  +-------------------------------------+
           |                  | [X] Select Non-Parametric Pathway    |
           |                  +-------------------------------------+
           |                               |
           v                               v
+--------------------------+    [ 2. Kruskal-Wallis Test ]
| [V] Select Parametric      |               |
|    Pathway (ANOVA)       |               v
+--------------------------+    < Is result significant? >
           |                      |        |
           |                      |--[No]-->| [ End: Report No Sig. Difference ]
           |                      |        |
           |                    (Yes)      |
           |                      |        |
           |                      v        |
           |        +--------------------------+
           |        |    3. Post-Hoc Test      |
           |        +--------------------------+
           |                      |
           |                      v
           |             < Are groups > 2? >
           |                |         |
           |                |--[No]-->| [ End: No Post-Hoc Needed ]
           |                |         |
           |              (Yes)       |
           |                |         |
           |                v         |
           |     [ Perform Dunn's Test ]    |
           |                |         |
           `----------------|----.    .----|----------------'
                            |    |    |    |
                            v    v    v    v
                       +-------------------------+
                       |  End: Generate Report   |
                       +-------------------------+
"

  return(trimws(flowchart_string))
}
# ==============================================================================
# HELPER 2: GENERATE PHYLOGENETIC TREE PLOT
# ==============================================================================

#' Generate Phylogenetic Tree Plot
#'
#' Constructs a time-scaled \code{ggtree} plot from numeric tree data and metadata,
#' displaying tip locations and formatted date axes.
#'
#' @param tree_numeric A \code{treedata} object containing tree structure and node dates.
#' @param meta Data frame of specimen metadata with tip labels and locations.
#' @param tree_date_bounds Numeric vector of length 2 specifying plot lower and upper date limits.
#'
#' @return A \code{ggtree} / \code{ggplot} object representing the phylogenetic tree.
#' @import ggplot2
#' @import ggtree
#' @importFrom dplyr left_join
#' @importFrom lubridate date_decimal
#' @export
plot_phylo_tree <- function(tree_numeric, meta, tree_date_bounds) {
  message("[phymapr] Generating Phylogenetic Tree Plot...")

  plot_data <- ggplot2::fortify(tree_numeric)
  plot_data$x <- plot_data$date  

  plot_data_full <- plot_data %>%
    dplyr::left_join(meta, by = c("label" = "SPECIMEN_NUMBER"))

  p_tree <- ggtree::ggtree(plot_data_full) + 
    ggtree::theme_tree2() + 
    ggtree::geom_tippoint(ggplot2::aes(color = location), size = 3) +
    ggtree::geom_tiplab(ggplot2::aes(label = paste(label, location, sep = " - ")), size = 2.5, offset = 0.05) +
    ggplot2::scale_x_continuous(
      limits = tree_date_bounds, 
      labels = function(x) format(lubridate::date_decimal(x), "%Y-%b")
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo", na.value = "grey50") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(x = "Date", color = "Region")

  return(p_tree)
}
# ==============================================================================
# MASTER FUNCTION: RUN FULL PIPELINE
# ==============================================================================
#' Generate Phylogenetic Transmission Pipeline
#' 
#' @import ape
#' @import treeio
#' @import ggtree
#' @import ggplot2
#' @import dplyr
#' @import phytools
#' @import tidygeocoder
#' @import gganimate
#' @import lubridate
#' @import maps
#' @import sf
#' @import gifski
#' @import tidyr
#' @import ggspatial
#' @export
generate_phylo_transmission <- function(
    timetree_file, 
    metadata_file,
    col_specimen = "SPECIMEN_NUMBER", 
    col_date = "date",                
    col_location = "location",        
    n_prune_early = 4, 
    tree_date_bounds = c(NA, NA), 
    map_bounds = c(xmin = NA, xmax = NA, ymin = NA, ymax = NA),
    map_style = "polygon" 
) {
  
  # Step 1: Prep data and establish bounds
  data_objs <- prep_phylo_data(timetree_file, metadata_file, col_specimen, col_date, col_location, tree_date_bounds, map_bounds)
  
  # Step 2: Build tree plot
  p_tree <- plot_phylo_tree(data_objs$tree_numeric, data_objs$meta, data_objs$tree_date_bounds)
  
  # Step 3: Run fastAnc state reconstruction
  mcc_tab <- run_asr(data_objs$tree_numeric, data_objs$meta, data_objs$all_dates)
  
  # Step 4: Map and snap coordinates
  snapped_data <- snap_and_prune_locations(mcc_tab, data_objs$meta, n_prune_early)
  
  # Step 5: Build map using the dynamically updated map_bounds from Step 1
  anim_map <- build_map_animation(snapped_data$mcc_tab_final, snapped_data$location_lookup, data_objs$map_bounds, map_style)
  
  message("Complete! Returning Tree and Map objects.")
  
  # Return final list containing all requested objects
  return(list(
    tree_plot = p_tree,
    map_animation = anim_map,
    mcc_data = snapped_data$mcc_tab_final 
  ))
}
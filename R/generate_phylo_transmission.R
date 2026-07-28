# ==============================================================================
# MASTER FUNCTION: RUN FULL TRANSMISSION PIPELINE
# ==============================================================================

#' Generate End-to-End Phylogenetic Transmission Pipeline
#'
#' Runs the complete end-to-end transmission workflow: pre-processes trees and spatial metadata,
#' generates time-scaled tree plots, performs ancestral state reconstruction (ASR), executes the
#' epidemiologic inference layer to collapse residence episodes and score transition pathways, applies
#' non-linear blended pacing, and renders animated geographic transmission maps.
#'
#' @param timetree_file Path to input tree file (BEAST/TreeTime NEXUS or Newick).
#' @param metadata_file Path to metadata file (CSV or TSV format).
#' @param col_specimen Identifier column name for specimen matching. Default "SPECIMEN_NUMBER".
#' @param col_date Date column name in metadata. Default "date".
#' @param col_location Location column name in metadata. Default "location".
#' @param n_prune_early Number of early lineage tips to prune if requested. Default 4.
#' @param tree_date_bounds Numeric vector for tree plot x-axis date range.
#' @param map_bounds Spatial map bounding box coordinates (xmin, xmax, ymin, ymax).
#' @param map_style Map visual style: "polygon" (fast rendering) or "tile" (dark satellite tiles). Default "polygon".
#' @param threshold_direct_snp Maximum SNP distance threshold for direct transmission classification. Default 2.
#' @param threshold_indirect_snp Maximum SNP distance threshold before categorizing as import/distant. Default 5.
#' @param distant_rendering_style Rendering style for distant/import transitions ("import_arrow", "pulse", "hide", "curve"). Default "import_arrow".
#' @param animation_pace_balance Blended pacing weight between linear time (1.0) and rank order (0.0). Default 0.5.
#'
#' @return A named list containing:
#' \item{tree_plot}{A \code{ggtree} ggplot object}
#' \item{map_animation}{A \code{gganimate} object for spatial map animation}
#' \item{episodes}{Data frame of local circulation residence episodes}
#' \item{pathways}{Data frame of scored transmission pathways with confidence categories}
#'
#' @import ape
#' @import treeio
#' @import ggtree
#' @import ggplot2
#' @import dplyr
#' @import tidygeocoder
#' @import gganimate
#' @import lubridate
#' @import maps
#' @import gifski
#' @import tidyr
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
    map_style = "polygon",
    threshold_direct_snp = 2,
    threshold_indirect_snp = 5,
    distant_rendering_style = "import_arrow",
    animation_pace_balance = 0.5
) {
  # Step 1: Prep data and establish spatial/temporal bounds
  data_objs <- prep_phylo_data(
    timetree_file, metadata_file, col_specimen, col_date, col_location, 
    tree_date_bounds, map_bounds
  )

  # Step 2: Build phylogenetic tree plot
  p_tree <- plot_phylo_tree(data_objs$tree_numeric, data_objs$meta, data_objs$tree_date_bounds)

  # Step 3: Run fast Ancestral State Reconstruction
  mcc_tab <- run_asr(data_objs$tree_numeric, data_objs$meta, data_objs$all_dates)

  # Step 4: Run Epidemiologic Inference Layer (Episode detection + confidence classification)
  inference_res <- run_epidemiologic_inference(
    tree_numeric           = data_objs$tree_numeric,
    meta                   = data_objs$meta,
    threshold_direct_snp   = threshold_direct_snp,
    threshold_indirect_snp = threshold_indirect_snp
  )

  # Step 5: Apply Blended Non-Linear Pacing
  paced_res <- apply_blended_pacing(
    episodes = inference_res$episodes,
    pathways = inference_res$pathways,
    meta     = data_objs$meta,
    lambda   = animation_pace_balance
  )

  episodes_final <- paced_res$episodes %>%
    dplyr::mutate(reveal_time = reveal_start)

  # Step 6: Build map animation using paced results and map_bounds
  anim_map <- build_map_animation(
    episodes                = episodes_final,
    pathways                = paced_res$pathways,
    location_lookup         = inference_res$location_lookup,
    map_bounds              = data_objs$map_bounds,
    map_style               = map_style,
    distant_rendering_style = distant_rendering_style
  )

  message("[phymapr] Complete! Returning Tree, Map, and Epidemiologic Inference objects.")

  return(list(
    tree_plot     = p_tree,
    map_animation = anim_map,
    episodes      = paced_res$episodes,
    pathways      = paced_res$pathways
  ))
}
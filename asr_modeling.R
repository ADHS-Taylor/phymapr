# ==============================================================================
# HELPER 3: ANCESTRAL STATE RECONSTRUCTION
# ==============================================================================
run_asr <- function(tree_numeric, meta, all_dates) {
  message("Running Ancestral State Reconstruction...")
  
  phy_jittered <- tree_numeric@phylo
  phy_jittered$edge.length[phy_jittered$edge.length == 0] <- 1e-06 
  
  tree_labels <- phy_jittered$tip.label
  
  meta_cleaned <- meta %>%
    dplyr::filter(SPECIMEN_NUMBER %in% tree_labels) %>%
    dplyr::distinct(SPECIMEN_NUMBER, .keep_all = TRUE)
  meta_ordered <- meta_cleaned[match(tree_labels, meta_cleaned$SPECIMEN_NUMBER), ]
  
  tip_lat_vec <- setNames(as.numeric(meta_ordered$lat), meta_ordered$SPECIMEN_NUMBER)
  tip_lon_vec <- setNames(as.numeric(meta_ordered$long), meta_ordered$SPECIMEN_NUMBER)
  
  node_lat_reconstructed <- phytools::fastAnc(phy_jittered, tip_lat_vec)
  node_lon_reconstructed <- phytools::fastAnc(phy_jittered, tip_lon_vec)
  
  all_lat <- c(as.numeric(tip_lat_vec), as.numeric(node_lat_reconstructed))
  all_lon <- c(as.numeric(tip_lon_vec), as.numeric(node_lon_reconstructed))
  
  edges <- phy_jittered$edge
  mcc_tab <- data.frame(
    startYear = all_dates[edges[,1]],
    endYear   = all_dates[edges[,2]],
    startLat  = all_lat[edges[,1]],
    startLon  = all_lon[edges[,1]],
    endLat    = all_lat[edges[,2]],
    endLon    = all_lon[edges[,2]]
  ) %>% tidyr::drop_na()
  
  return(mcc_tab)
}
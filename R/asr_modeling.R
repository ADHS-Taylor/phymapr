# ==============================================================================
# HELPER 3: ANCESTRAL STATE RECONSTRUCTION (ASR)
# ==============================================================================

#' Run Brownian Motion Ancestral State Reconstruction
#'
#' Estimates continuous traits (e.g. latitude/longitude coordinates) at internal
#' nodes using a Brownian motion model on the phylogenetic tree.
#'
#' @param phy A \code{phylo} object representing the phylogenetic tree.
#' @param tip_vals A named numeric vector of trait values at tip nodes, where names match tip labels.
#'
#' @return A named numeric vector of estimated trait values for internal nodes (indexed by node ID).
#' @export
run_bm_asr <- function(phy, tip_vals) {
  edges <- phy$edge
  el <- phy$edge.length
  num_nodes <- phy$Nnode + length(phy$tip.label)
  N <- length(phy$tip.label)

  val <- rep(NA, num_nodes)
  val[1:N] <- tip_vals[phy$tip.label]
  var <- rep(0, num_nodes)

  phy_post <- ape::reorder.phylo(phy, "postorder")
  edges_p <- phy_post$edge
  el_p <- phy_post$edge.length

  children <- split(edges_p[, 2], edges_p[, 1])
  edge_idx <- split(seq_len(nrow(edges_p)), edges_p[, 1])
  unique_parents <- unique(edges_p[, 1])

  for (parent in unique_parents) {
    ch <- children[[as.character(parent)]]
    e_idx <- edge_idx[[as.character(parent)]]

    ch_vals <- val[ch]
    ch_vars <- var[ch]
    ch_edges <- el_p[e_idx]

    V <- ch_vars + ch_edges
    w <- 1 / V

    val[parent] <- sum(w * ch_vals) / sum(w)
    var[parent] <- 1 / sum(w)
  }

  res <- val[(N + 1):num_nodes]
  names(res) <- (N + 1):num_nodes
  return(res)
}

#' Run Continuous Ancestral State Reconstruction for Tree Coordinates
#'
#' Extracts continuous geographic coordinates (lat/long) across all nodes and
#' computes per-branch SNP distance metrics.
#'
#' @param tree_numeric A \code{treedata} object containing phylogenetic tree and node metadata.
#' @param meta A data frame containing specimen spatial metadata with \code{SPECIMEN_NUMBER}, \code{lat}, and \code{long}.
#' @param all_dates A numeric vector of node dates corresponding to node IDs.
#'
#' @return A data frame containing MCC table data with edge parent/child node coordinates, dates, and SNP distances.
#' @importFrom dplyr filter distinct
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na
#' @export
run_asr <- function(tree_numeric, meta, all_dates) {
  message("[phymapr] Running Ancestral State Reconstruction...")

  phy_jittered <- tree_numeric@phylo
  phy_jittered$edge.length[phy_jittered$edge.length == 0] <- 1e-06

  tree_labels <- phy_jittered$tip.label

  meta_cleaned <- meta %>%
    dplyr::filter(SPECIMEN_NUMBER %in% tree_labels) %>%
    dplyr::distinct(SPECIMEN_NUMBER, .keep_all = TRUE)
  meta_ordered <- meta_cleaned[match(tree_labels, meta_cleaned$SPECIMEN_NUMBER), ]

  tip_lat_vec <- stats::setNames(as.numeric(meta_ordered$lat), meta_ordered$SPECIMEN_NUMBER)
  tip_lon_vec <- stats::setNames(as.numeric(meta_ordered$long), meta_ordered$SPECIMEN_NUMBER)

  node_lat_reconstructed <- run_bm_asr(phy_jittered, tip_lat_vec)
  node_lon_reconstructed <- run_bm_asr(phy_jittered, tip_lon_vec)

  all_lat <- c(as.numeric(tip_lat_vec), as.numeric(node_lat_reconstructed))
  all_lon <- c(as.numeric(tip_lon_vec), as.numeric(node_lon_reconstructed))

  edges <- phy_jittered$edge

  # Extract SNP distance for each child node
  td <- tibble::as_tibble(tree_numeric)
  snp_distances <- sapply(edges[, 2], function(nd) {
    row <- td[td$node == nd, ]
    if (nrow(row) == 0) return(0)
    if ("mutations" %in% colnames(row)) {
      muts <- row$mutations[[1]]
      if (is.null(muts) || all(is.na(muts)) || (length(muts) == 1 && muts == "")) {
        return(0)
      }
      return(length(muts))
    }
    return(0)
  })

  mcc_tab <- data.frame(
    parent_node  = edges[, 1],
    child_node   = edges[, 2],
    startYear    = all_dates[edges[, 1]],
    endYear      = all_dates[edges[, 2]],
    startLat     = all_lat[edges[, 1]],
    startLon     = all_lon[edges[, 1]],
    endLat       = all_lat[edges[, 2]],
    endLon       = all_lon[edges[, 2]],
    snp_distance = snp_distances
  ) %>% tidyr::drop_na()

  return(mcc_tab)
}
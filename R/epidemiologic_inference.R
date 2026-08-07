# ==============================================================================
# EPIDEMIOLOGIC INFERENCE LAYER
# Collapses redundant local branches into local circulation episodes,
# calculates transition confidence using multivariate features,
# and applies pacing.
# ==============================================================================

#' Run Epidemiologic Inference Layer
#'
#' Evaluates continuous spatial reconstruction, snaps internal nodes to discrete location lookup entries,
#' collapses contiguous regional stays into residence episodes, and categorizes transmission pathways
#' into direct, indirect, or import categories based on SNP distance, duration, and sampling density.
#'
#' @param tree_numeric A \code{treedata} object containing phylogenetic tree and node metadata.
#' @param meta Data frame with specimen spatial and temporal metadata.
#' @param threshold_direct_snp Maximum SNP distance threshold for classifying direct transmission pathways. Default 2.
#' @param threshold_indirect_snp Maximum SNP distance threshold before categorizing as distant/import. Default 5.
#'
#' @return A list containing:
#' \item{episodes}{Data frame of local circulation residence episodes}
#' \item{pathways}{Data frame of categorized transmission pathways with confidence levels}
#' \item{location_lookup}{Distinct location lookup table with spatial coordinates}
#'
#' @importFrom dplyr select drop_na distinct filter mutate group_by summarise n case_when
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble
#' @importFrom lubridate decimal_date as_date
#' @importFrom stats setNames median dist
#' @export
run_epidemiologic_inference <- function(
    tree_numeric, 
    meta, 
    threshold_direct_snp = 2, 
    threshold_indirect_snp = 5
) {
  message("[phymapr] Running Epidemiologic Inference Layer...")

  phy <- tree_numeric@phylo
  td <- tibble::as_tibble(tree_numeric)

  # 1. Establish location lookup table
  location_lookup <- meta %>% 
    dplyr::select(location, lat, long) %>% 
    tidyr::drop_na() %>%
    dplyr::distinct(location, .keep_all = TRUE)

  snap_to_nearest <- function(q_lat, q_lon, lookup) {
    if (is.na(q_lat) || is.na(q_lon)) {
      return(list(location = NA, lat = NA, long = NA, dist = NA))
    }
    dists <- sqrt((lookup$lat - q_lat)^2 + (lookup$long - q_lon)^2)
    idx <- which.min(dists)
    list(
      location = lookup$location[idx],
      lat = lookup$lat[idx],
      long = lookup$long[idx],
      dist = dists[idx]
    )
  }

  # 2. Reconstruct continuous ancestral coordinates
  tree_labels <- phy$tip.label
  meta_cleaned <- meta %>%
    dplyr::filter(SPECIMEN_NUMBER %in% tree_labels) %>%
    dplyr::distinct(SPECIMEN_NUMBER, .keep_all = TRUE)
  meta_ordered <- meta_cleaned[match(tree_labels, meta_cleaned$SPECIMEN_NUMBER), ]

  tip_lat_vec <- stats::setNames(as.numeric(meta_ordered$lat), meta_ordered$SPECIMEN_NUMBER)
  tip_lon_vec <- stats::setNames(as.numeric(meta_ordered$long), meta_ordered$SPECIMEN_NUMBER)

  # Jitter branch lengths to avoid zero lengths for Brownian motion ASR
  phy_jittered <- phy
  phy_jittered$edge.length[phy_jittered$edge.length == 0] <- 1e-06

  node_lat_reconstructed <- run_bm_asr(phy_jittered, tip_lat_vec)
  node_lon_reconstructed <- run_bm_asr(phy_jittered, tip_lon_vec)

  all_lat <- c(as.numeric(tip_lat_vec), as.numeric(node_lat_reconstructed))
  all_lon <- c(as.numeric(tip_lon_vec), as.numeric(node_lon_reconstructed))

  # Extract node dates
  all_dates <- td$date
  if (is.character(all_dates)) {
    all_dates <- as.numeric(all_dates)
  }

  num_nodes <- length(all_dates)

  # 3. Snap coordinates of all nodes to the nearest known location
  node_location <- rep(NA, num_nodes)
  node_lat_snapped <- rep(NA, num_nodes)
  node_lon_snapped <- rep(NA, num_nodes)
  node_snap_dist <- rep(NA, num_nodes)

  for (i in 1:num_nodes) {
    snap <- snap_to_nearest(all_lat[i], all_lon[i], location_lookup)
    node_location[i]     <- snap$location
    node_lat_snapped[i]  <- snap$lat
    node_lon_snapped[i]  <- snap$long
    node_snap_dist[i]    <- snap$dist
  }

  # 4. Collapse nodes into Contiguous Residence Episodes
  parents <- rep(NA, num_nodes)
  parents[phy$edge[, 2]] <- phy$edge[, 1]

  node_order <- order(all_dates)
  episode_ids <- rep(NA, num_nodes)
  next_episode_id <- 1

  for (nd in node_order) {
    p_nd <- parents[nd]
    if (is.na(p_nd)) {
      episode_ids[nd] <- next_episode_id
      next_episode_id <- next_episode_id + 1
    } else {
      # If parent and child share the same location, collapse into the same episode
      if (!is.na(node_location[nd]) && !is.na(node_location[p_nd]) && 
          node_location[nd] == node_location[p_nd]) {
        episode_ids[nd] <- episode_ids[p_nd]
      } else {
        episode_ids[nd] <- next_episode_id
        next_episode_id <- next_episode_id + 1
      }
    }
  }

  # Construct Episode Table
  episode_data <- data.frame(
    node = 1:num_nodes,
    episode_id = episode_ids,
    location = node_location,
    date = all_dates,
    lat = node_lat_snapped,
    long = node_lon_snapped
  ) %>% tidyr::drop_na(location)

  episodes <- episode_data %>%
    dplyr::group_by(episode_id, location, lat, long) %>%
    dplyr::summarise(
      startYear = min(date, na.rm = TRUE),
      endYear   = max(date, na.rm = TRUE),
      num_nodes = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(duration = endYear - startYear)

  # 5. Identify Transition Pathways and Evaluate Confidence
  edges <- phy$edge

  # Calculate SNP distances for each child branch
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

  # Dynamic Auto-Threshold Scaling if SNP distances across tree edges exceed default thresholds
  if (length(snp_distances[snp_distances > 0]) > 0) {
    median_snp <- stats::median(snp_distances[snp_distances > 0], na.rm = TRUE)
    if (!is.na(median_snp) && median_snp > threshold_indirect_snp && threshold_direct_snp == 2 && threshold_indirect_snp == 5) {
      threshold_direct_snp <- max(2, round(median_snp))
      threshold_indirect_snp <- max(5, round(median_snp * 1.5))
      message(sprintf("[phymapr] Auto-scaling SNP confidence thresholds (median branch SNPs = %.1f): direct <= %d, indirect > %d", 
                      median_snp, threshold_direct_snp, threshold_indirect_snp))
    }
  }

  # Identify tip vs internal nodes — tips have reliable locations from metadata,

  # internal nodes have ASR-reconstructed locations that may be unreliable
  n_tips <- length(phy$tip.label)
  is_tip_node <- seq_len(num_nodes) <= n_tips

  # Calculate a snapping distance threshold for internal nodes.
  # Internal nodes with snap distance above this threshold have unreliable

  # location assignments (BM ASR placed them far from any real sample site).
  # Use the median pairwise distance between known locations as a reference.
  if (nrow(location_lookup) >= 2) {
    loc_dists <- as.numeric(stats::dist(location_lookup[, c("lat", "long")]))
    snap_threshold <- stats::median(loc_dists) * 0.4
  } else {
    snap_threshold <- 2.0
  }
  message(sprintf("[phymapr] Internal node snap threshold: %.2f degrees (internal nodes snapped further are treated as unreliable).", snap_threshold))

  # Calculate transition details
  pathways <- data.frame(
    parent_node          = edges[, 1],
    child_node           = edges[, 2],
    startYear            = all_dates[edges[, 1]],
    endYear              = all_dates[edges[, 2]],
    startLat             = node_lat_snapped[edges[, 1]],
    startLon             = node_lon_snapped[edges[, 1]],
    endLat               = node_lat_snapped[edges[, 2]],
    endLon               = node_lon_snapped[edges[, 2]],
    origin_location      = node_location[edges[, 1]],
    destination_location = node_location[edges[, 2]],
    snp_distance         = snp_distances,
    snapping_distance    = node_snap_dist[edges[, 1]] + node_snap_dist[edges[, 2]],
    parent_snap_dist     = node_snap_dist[edges[, 1]],
    child_snap_dist      = node_snap_dist[edges[, 2]],
    parent_is_tip        = is_tip_node[edges[, 1]],
    child_is_tip         = is_tip_node[edges[, 2]]
  ) %>%
    tidyr::drop_na(origin_location, destination_location) %>%
    # Filter for transitions between different locations
    dplyr::filter(origin_location != destination_location) %>%
    # Remove transitions involving internal nodes with unreliable snap locations.
    # An internal node with high snap distance was placed by BM ASR far from any
    # real sample — its "location" is a phantom jump artifact, not a real transition.
    dplyr::filter(
      (parent_is_tip | parent_snap_dist <= snap_threshold) &
      (child_is_tip  | child_snap_dist  <= snap_threshold)
    ) %>%
    dplyr::select(-parent_snap_dist, -child_snap_dist, -parent_is_tip, -child_is_tip) %>%
    dplyr::mutate(time_elapsed = endYear - startYear)

  # Calculate Surveillance Density in metadata
  meta_dates <- as.numeric(lubridate::decimal_date(lubridate::as_date(meta$date)))

  calc_density <- function(t_start, t_end) {
    sum(meta_dates >= (t_start - 0.2) & meta_dates <= (t_end + 0.2), na.rm = TRUE)
  }

  pathways$surveillance_density <- mapply(calc_density, pathways$startYear, pathways$endYear)

  # Classify pathway confidence category
  pathways <- pathways %>%
    dplyr::mutate(
      transmission_type = dplyr::case_when(
        # A. Direct Pathway (Likely)
        snp_distance <= threshold_direct_snp & time_elapsed <= 0.12 & snapping_distance <= 1.5 ~ "Direct Pathway (Likely)",

        # B. Distant / Import
        snp_distance > threshold_indirect_snp | time_elapsed > 0.35 |
          (surveillance_density > 25 & (snp_distance > threshold_direct_snp | time_elapsed > 0.15)) ~ "Distant / Import",

        # C. Indirect Pathway (Missing Intermediates)
        TRUE ~ "Indirect (Missing Intermediates)"
      )
    )

  message(sprintf("[phymapr] Inferred %d Local Residence Episodes and %d Transition Pathways.", 
                  nrow(episodes), nrow(pathways)))

  return(list(
    episodes = episodes,
    pathways = pathways,
    location_lookup = location_lookup
  ))
}

#' Apply Blended Non-Linear Pacing
#'
#' Adjusts temporal reveal coordinates by blending linear timeline position with
#' event rank-order weighting, smoothing out animation gaps in sparse sampling periods.
#'
#' @param episodes Data frame of residence episodes.
#' @param pathways Data frame of transition pathways.
#' @param meta Data frame of metadata.
#' @param lambda Blending weight between linear time (1.0) and rank order (0.0). Default 0.5.
#'
#' @return A list containing updated \code{episodes} and \code{pathways} data frames with \code{reveal_start} and \code{reveal_end} columns.
#' @export
apply_blended_pacing <- function(episodes, pathways, meta, lambda = 0.5) {
  message(sprintf("[phymapr] Applying Blended Non-Linear Pacing (lambda = %.2f)...", lambda))

  # Combine all timepoints
  all_times <- unique(c(episodes$startYear, episodes$endYear, pathways$startYear, pathways$endYear))
  all_times <- sort(all_times)

  t_min <- min(all_times, na.rm = TRUE)
  t_max <- max(all_times, na.rm = TRUE)

  if (length(all_times) <= 1 || t_max == t_min) {
    episodes$reveal_start <- episodes$startYear
    episodes$reveal_end   <- episodes$endYear
    pathways$reveal_start <- pathways$startYear
    pathways$reveal_end   <- pathways$endYear
    return(list(episodes = episodes, pathways = pathways))
  }

  # Map time t using linear-rank blend
  map_time <- function(t) {
    if (is.na(t)) return(NA)
    u <- (t - t_min) / (t_max - t_min)
    rk <- which(all_times == t)
    r <- (rk - 1) / (length(all_times) - 1)
    w <- lambda * u + (1 - lambda) * r
    t_prime <- t_min + w * (t_max - t_min)
    return(t_prime)
  }

  map_time_vec <- function(t_vec) {
    sapply(t_vec, map_time)
  }

  episodes$reveal_start <- map_time_vec(episodes$startYear)
  episodes$reveal_end   <- map_time_vec(episodes$endYear)
  pathways$reveal_start <- map_time_vec(pathways$startYear)
  pathways$reveal_end   <- map_time_vec(pathways$endYear)

  return(list(episodes = episodes, pathways = pathways))
}

# ==============================================================================
# HELPER 1: DATA PREPARATION, GEOCODING & DYNAMIC BOUNDS
# ==============================================================================

#' Prepare Phylogenetic and Spatial Metadata
#'
#' Reads time-scaled phylogenetic trees and metadata files, normalizes dates and specimen identifiers,
#' geocodes locations if coordinates are missing, and computes dynamic map bounding boxes.
#'
#' @param timetree_file Path to the input tree file (BEAST/TreeTime NEXUS or Newick format).
#' @param metadata_file Path to the metadata file (CSV or TSV format).
#' @param col_specimen Column name for specimen identifiers. Default "SPECIMEN_NUMBER".
#' @param col_date Column name for sample collection dates. Default "date".
#' @param col_location Column name for geographic locations. Default "location".
#' @param tree_date_bounds Numeric vector of length 2 giving lower and upper date bounds, or c(NA, NA) to auto-compute.
#' @param map_bounds Named numeric vector with xmin, xmax, ymin, ymax bounds, or NA values to auto-compute.
#'
#' @return A list containing:
#' \item{tree_numeric}{A \code{treedata} object with numeric date metadata}
#' \item{meta}{Cleaned metadata data frame with geocoded lat/long coordinates}
#' \item{all_dates}{Numeric vector of all node dates}
#' \item{tree_date_bounds}{Final expanded date bounds for plotting}
#' \item{map_bounds}{Final spatial map bounding box coordinates}
#'
#' @importFrom dplyr rename all_of mutate
#' @importFrom tibble as_tibble
#' @importFrom treeio read.beast drop.tip as.treedata
#' @importFrom tidygeocoder geocode
#' @export
prep_phylo_data <- function(
    timetree_file, 
    metadata_file, 
    col_specimen = "SPECIMEN_NUMBER", 
    col_date = "date", 
    col_location = "location", 
    tree_date_bounds = c(NA, NA), 
    map_bounds = c(xmin = NA, xmax = NA, ymin = NA, ymax = NA)
) {
  message("[phymapr] Loading tree and metadata...")

  tt_tree <- treeio::read.beast(timetree_file)
  
  # Safely drop tip "Reference" only if present in tree
  if ("Reference" %in% tt_tree@phylo$tip.label) {
    tt_tree <- treeio::drop.tip(tt_tree, "Reference")
  }

  # Support for both TSV and CSV based on file extension
  if (grepl("\\.tsv$|\\.txt$", metadata_file, ignore.case = TRUE)) {
    meta <- utils::read.delim(metadata_file, sep = "\t")
  } else {
    meta <- utils::read.csv(metadata_file)
  }

  meta <- meta %>%
    dplyr::rename(
      SPECIMEN_NUMBER = dplyr::all_of(col_specimen),
      location = dplyr::all_of(col_location)
    )

  if (col_date %in% colnames(meta) && col_date != "date") {
    meta <- meta %>% dplyr::rename(date = dplyr::all_of(col_date))
  }

  td <- tibble::as_tibble(tt_tree) %>%
    dplyr::mutate(date = as.numeric(as.character(date)))
  tt_tree_numeric <- treeio::as.treedata(td)

  all_dates <- td$date
  mrsd <- max(all_dates, na.rm = TRUE)

  if (is.na(tree_date_bounds[2])) tree_date_bounds[2] <- mrsd
  if (is.na(tree_date_bounds[1])) tree_date_bounds[1] <- min(all_dates, na.rm = TRUE)

  # Calculate total date range and add padding
  date_range <- tree_date_bounds[2] - tree_date_bounds[1]
  tree_date_bounds[2] <- tree_date_bounds[2] + (date_range * 0.15)

  if (!all(c("lat", "long") %in% colnames(meta))) {
    message("[phymapr] Geocoding locations via OpenStreetMap...")
    meta <- meta %>% tidygeocoder::geocode(location, method = 'osm', lat = lat, long = long)
  }

  # Dynamic Map Bounds Calculation
  if (any(is.na(map_bounds))) {
    message("[phymapr] Calculating dynamic map bounds with padding...")

    min_lon <- min(as.numeric(meta$long), na.rm = TRUE)
    max_lon <- max(as.numeric(meta$long), na.rm = TRUE)
    min_lat <- min(as.numeric(meta$lat), na.rm = TRUE)
    max_lat <- max(as.numeric(meta$lat), na.rm = TRUE)

    lon_pad <- (max_lon - min_lon) * 0.10
    lat_pad <- (max_lat - min_lat) * 0.10

    if (is.na(lon_pad) || lon_pad == 0) lon_pad <- 2
    if (is.na(lat_pad) || lat_pad == 0) lat_pad <- 2

    if (is.na(map_bounds["xmin"])) map_bounds["xmin"] <- min_lon - lon_pad
    if (is.na(map_bounds["xmax"])) map_bounds["xmax"] <- max_lon + lon_pad
    if (is.na(map_bounds["ymin"])) map_bounds["ymin"] <- min_lat - lat_pad
    if (is.na(map_bounds["ymax"])) map_bounds["ymax"] <- max_lat + lat_pad
  }

  return(list(
    tree_numeric     = tt_tree_numeric,
    meta             = meta,
    all_dates        = all_dates,
    tree_date_bounds = tree_date_bounds,
    map_bounds       = map_bounds
  ))
}
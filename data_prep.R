# ==============================================================================
# HELPER 1: DATA PREPARATION, GEOCODING & DYNAMIC BOUNDS
# ==============================================================================
prep_phylo_data <- function(timetree_file, metadata_file, col_specimen, col_date, col_location, tree_date_bounds, map_bounds) {
  message("Loading tree and metadata...")
  
  tt_tree <- treeio::read.beast(timetree_file) %>% ape::drop.tip("Reference")
  
  # Support for both TSV and CSV based on file extension
  if (grepl("\\.tsv$|\\.txt$", metadata_file, ignore.case = TRUE)) {
    meta <- read.delim(metadata_file, sep = "\t")
  } else {
    meta <- read.csv(metadata_file) 
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
  
  # Calculate total date range
  date_range <- tree_date_bounds[2] - tree_date_bounds[1]
  
  # Add 15% to the upper limit to make room for text
  tree_date_bounds[2] <- tree_date_bounds[2] + (date_range * 0.15)
  
  if (!all(c("lat", "long") %in% colnames(meta))) {
    message("Geocoding locations...")
    meta <- meta %>% tidygeocoder::geocode(location, method = 'osm', lat = lat, long = long)
  }
  
  # Dynamic Map Bounds Calculation
  if (any(is.na(map_bounds))) {
    message("Calculating dynamic map bounds with padding...")
    
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
    tree_numeric = tt_tree_numeric,
    meta = meta,
    all_dates = all_dates,
    tree_date_bounds = tree_date_bounds,
    map_bounds = map_bounds 
  ))
}
# ==============================================================================
# HELPER 4: LOCATION SNAPPING & PRUNING
# ==============================================================================
snap_and_prune_locations <- function(mcc_tab, meta, n_prune_early) {
  message("Mapping and snapping coordinates...")
  
  location_lookup <- meta %>% 
    dplyr::select(location, lat, long) %>% 
    tidyr::drop_na() %>%
    dplyr::distinct(location, .keep_all = TRUE)
  
  find_nearest_location <- function(q_lat, q_lon, lookup) {
    dists <- sqrt((lookup$lat - q_lat)^2 + (lookup$long - q_lon)^2)
    lookup$location[which.min(dists)]
  }
  
  snap_to_nearest <- function(q_lat, q_lon, lookup) {
    dists <- sqrt((lookup$lat - q_lat)^2 + (lookup$long - q_lon)^2)
    lookup[which.min(dists), ]
  }
  
  mcc_tab_final <- mcc_tab %>%
    dplyr::rowwise() %>%
    dplyr::mutate(location = find_nearest_location(startLat, startLon, location_lookup)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(startYear) %>%
    utils::tail(n = -n_prune_early) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(
      start_info = list(snap_to_nearest(startLat, startLon, location_lookup)),
      end_info   = list(snap_to_nearest(endLat, endLon, location_lookup))
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest_wider(start_info, names_sep = "_") %>%
    tidyr::unnest_wider(end_info, names_sep = "_") %>%
    dplyr::mutate(
      startLat = start_info_lat,
      startLon = start_info_long,
      origin_location = start_info_location,
      endLat = end_info_lat,
      endLon = end_info_long
    ) %>%
    dplyr::filter(!(startLat == endLat & startLon == endLon)) 
  
  return(list(mcc_tab_final = mcc_tab_final, location_lookup = location_lookup))
}

# ==============================================================================
# HELPER 5: GENERATE ANIMATED MAP
# ==============================================================================
build_map_animation <- function(mcc_tab_final, location_lookup, map_bounds, map_style) {
  message("Building map animation object...")
  
  anim_map <- ggplot2::ggplot()
  
  if (map_style == "tile") {
    anim_map <- anim_map +
      ggspatial::annotation_map_tile(type = "cartodark", zoom = 3) + 
      ggplot2::geom_curve(data = mcc_tab_final, 
                          ggplot2::aes(x = startLon, y = startLat, xend = endLon, yend = endLat, 
                                       color = origin_location, group = seq_along(startLon)),
                          curvature = 0.2, alpha = 0.9, size = 1.2, 
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +
      ggplot2::geom_point(data = location_lookup, ggplot2::aes(x = long, y = lat), 
                          color = "white", fill = "black", size = 2.5, shape = 21, stroke = 1) +
      ggplot2::coord_sf(xlim = c(map_bounds["xmin"], map_bounds["xmax"]), 
                        ylim = c(map_bounds["ymin"], map_bounds["ymax"]), 
                        crs = 4326, expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "gray10", color = NA),
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(), 
        axis.title = ggplot2::element_blank()
      )
  } else {
    world_map <- ggplot2::map_data("world")
    us_states <- ggplot2::map_data("state") 
    
    anim_map <- anim_map +
      ggplot2::geom_polygon(data = world_map, ggplot2::aes(x = long, y = lat, group = group), 
                            fill = "gray96", color = "gray80", size = 0.2) +
      ggplot2::geom_path(data = us_states, ggplot2::aes(x = long, y = lat, group = group), 
                         color = "white", size = 0.3) +
      ggplot2::geom_curve(data = mcc_tab_final, 
                          ggplot2::aes(x = startLon, y = startLat, xend = endLon, yend = endLat, 
                                       color = origin_location, group = seq_along(startLon)),
                          curvature = 0.2, alpha = 0.8, size = 1.1,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +
      ggplot2::geom_point(data = location_lookup, ggplot2::aes(x = long, y = lat), 
                          color = "black", fill = "white", size = 2, shape = 21) +
      ggplot2::coord_fixed(xlim = c(map_bounds["xmin"], map_bounds["xmax"]), 
                           ylim = c(map_bounds["ymin"], map_bounds["ymax"]), ratio = 1.3) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"))
  }
  
  anim_map <- anim_map +
    gganimate::transition_reveal(endYear) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::labs(subtitle = "Time: {frame_along}", color = "Source Location") +
    ggplot2::theme(
      legend.position = "bottom",                  
      legend.text = ggplot2::element_text(size = 9),        
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.key.size = ggplot2::unit(0.5, "cm")            
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  
  return(anim_map)
}
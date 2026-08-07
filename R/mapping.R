# ==============================================================================
# MAPPING & ARC INTERPOLATION
# Interpolates Bezier curves, imports, and pulses, and builds map animations.
# ==============================================================================

#' Interpolate Bezier Curve Points for Map Arcs
#'
#' @param x1 Origin longitude.
#' @param y1 Origin latitude.
#' @param x2 Destination longitude.
#' @param y2 Destination latitude.
#' @param t_start Animation reveal start time.
#' @param t_end Animation reveal end time.
#' @param curvature Curvature intensity factor. Default 0.2.
#' @param n_points Number of intermediate points to generate along curve. Default 15.
#' @param pathway_id Unique pathway identifier. Default 1.
#'
#' @return Data frame of interpolated curve points with \code{x}, \code{y}, \code{reveal_time}, and \code{pathway_id}.
#' @export
interpolate_bezier <- function(x1, y1, x2, y2, t_start, t_end, curvature = 0.2, n_points = 15, pathway_id = 1) {
  p_seq <- seq(0, 1, length.out = n_points)
  xm <- (x1 + x2) / 2
  ym <- (y1 + y2) / 2
  vx <- -(y2 - y1)
  vy <- (x2 - x1)
  xctrl <- xm - curvature * vx
  yctrl <- ym - curvature * vy
  x <- (1 - p_seq)^2 * x1 + 2 * (1 - p_seq) * p_seq * xctrl + p_seq^2 * x2
  y <- (1 - p_seq)^2 * y1 + 2 * (1 - p_seq) * p_seq * yctrl + p_seq^2 * y2
  t_val <- t_start + p_seq * (t_end - t_start)

  data.frame(
    x = x,
    y = y,
    reveal_time = t_val,
    pathway_id = pathway_id,
    point_index = 1:n_points
  )
}

#' Interpolate Import Arrow Points
#'
#' @param x2 Destination longitude.
#' @param y2 Destination latitude.
#' @param t_start Animation reveal start time.
#' @param t_end Animation reveal end time.
#' @param offset Vertical latitude offset for arrow origin. Default 1.0.
#' @param n_points Number of points to generate. Default 15.
#' @param pathway_id Unique pathway identifier. Default 1.
#'
#' @return Data frame of arrow trajectory points.
#' @export
interpolate_arrow <- function(x2, y2, t_start, t_end, offset = 1.0, n_points = 15, pathway_id = 1) {
  p_seq <- seq(0, 1, length.out = n_points)
  x <- rep(x2, n_points)
  # Downward movement from y2 + offset to y2
  y <- (y2 + offset) - p_seq * offset
  t_val <- t_start + p_seq * (t_end - t_start)

  data.frame(
    x = x,
    y = y,
    reveal_time = t_val,
    pathway_id = pathway_id,
    point_index = 1:n_points
  )
}

#' Interpolate Geographic Pulse Effect Points
#'
#' @param x2 Target longitude.
#' @param y2 Target latitude.
#' @param t_end Reveal trigger time.
#' @param duration Temporal window duration for pulse expanding. Default 0.08.
#' @param n_points Number of points in pulse expansion series. Default 10.
#' @param pathway_id Unique pathway identifier. Default 1.
#'
#' @return Data frame of pulse animation points with \code{pulse_size} and \code{pulse_alpha}.
#' @export
interpolate_pulse <- function(x2, y2, t_end, duration = 0.08, n_points = 10, pathway_id = 1) {
  p_seq <- seq(0, 1, length.out = n_points)
  x <- rep(x2, n_points)
  y <- rep(y2, n_points)
  t_val <- t_end + p_seq * duration
  size <- 2 + p_seq * 8
  alpha <- 1.0 - p_seq * 0.9

  data.frame(
    x = x,
    y = y,
    reveal_time = t_val,
    pathway_id = pathway_id,
    point_index = 1:n_points,
    pulse_size = size,
    pulse_alpha = alpha
  )
}

#' Build Geographic Transmission Map Animation
#'
#' Constructs animated spatial transmission maps with custom styling, pathways,
#' background layers, and gganimate transitions.
#'
#' @param episodes Data frame of residence episodes.
#' @param pathways Data frame of transition pathways.
#' @param location_lookup Lookup table of distinct locations with lat/long.
#' @param map_bounds Named numeric vector of spatial bounds (xmin, xmax, ymin, ymax).
#' @param map_style Map aesthetic style: "polygon" (fast polygon map) or "tile" (dark satellite/tile map). Default "polygon".
#' @param distant_rendering_style How distant/import transitions are rendered: "import_arrow", "pulse", "hide", or "curve". Default "import_arrow".
#'
#' @return A \code{gganimate} object representing the animated transmission map.
#' @import ggplot2
#' @import gganimate
#' @importFrom dplyr group_by filter ungroup
#' @importFrom maps map_data
#' @export
build_map_animation <- function(
    episodes, 
    pathways, 
    location_lookup, 
    map_bounds, 
    map_style = "polygon",
    distant_rendering_style = "import_arrow"
) {
  message("[phymapr] Building map animation object...")

  if (map_style == "tile") {
    if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("ggspatial", quietly = TRUE)) {
      message("[phymapr] Packages 'sf' and/or 'ggspatial' are not available. Falling back to 'polygon' style.")
      map_style <- "polygon"
    }
  }

  # 1. Generate path and pulse coordinates for each transition pathway
  path_list <- list()
  pulse_list <- list()

  if (nrow(pathways) > 0) {
    for (i in seq_len(nrow(pathways))) {
      x1   <- pathways$startLon[i]
      y1   <- pathways$startLat[i]
      x2   <- pathways$endLon[i]
      y2   <- pathways$endLat[i]
      t1   <- pathways$reveal_start[i]
      t2   <- pathways$reveal_end[i]
      type <- pathways$transmission_type[i]
      orig <- pathways$origin_location[i]

      if (type == "Distant / Import") {
        if (distant_rendering_style == "hide") {
          next
        } else if (distant_rendering_style == "pulse") {
          pulse_df <- interpolate_pulse(x2, y2, t2, pathway_id = i)
          pulse_df$transmission_type <- type
          pulse_df$origin_location <- orig
          pulse_list[[length(pulse_list) + 1]] <- pulse_df
        } else if (distant_rendering_style == "import_arrow") {
          arrow_df <- interpolate_arrow(x2, y2, t1, t2, offset = 1.2, pathway_id = i)
          arrow_df$transmission_type <- type
          arrow_df$origin_location <- "External Import"
          path_list[[length(path_list) + 1]] <- arrow_df
        } else {
          bezier_df <- interpolate_bezier(x1, y1, x2, y2, t1, t2, curvature = 0.2, pathway_id = i)
          bezier_df$transmission_type <- type
          bezier_df$origin_location <- orig
          path_list[[length(path_list) + 1]] <- bezier_df
        }
      } else {
        bezier_df <- interpolate_bezier(x1, y1, x2, y2, t1, t2, curvature = 0.2, pathway_id = i)
        bezier_df$transmission_type <- type
        bezier_df$origin_location <- orig
        path_list[[length(path_list) + 1]] <- bezier_df
      }
    }
  }

  path_df <- if (length(path_list) > 0) do.call(rbind, path_list) else data.frame()
  pulse_df <- if (length(pulse_list) > 0) do.call(rbind, pulse_list) else data.frame()

  # 2. Construct ggplot object
  anim_map <- ggplot2::ggplot()

  # Base Map Layers
  if (map_style == "tile") {
    anim_map <- anim_map +
      ggspatial::annotation_map_tile(type = "cartodark", zoom = 3)
  } else {
    world_map <- ggplot2::map_data("world")
    us_states <- ggplot2::map_data("state")

    anim_map <- anim_map +
      ggplot2::geom_polygon(
        data = world_map, 
        ggplot2::aes(x = long, y = lat, group = group), 
        fill = "gray96", color = "gray85", linewidth = 0.2
      ) +
      ggplot2::geom_path(
        data = us_states, 
        ggplot2::aes(x = long, y = lat, group = group), 
        color = "white", linewidth = 0.3
      )
  }

  # Local Circulation Episodes (Static markers only — no animated moving node)
  # Episodes are represented by the static location markers below rather than

  # an animated point that traverses the full timeline.

  # Static Background Markers for all locations
  if (nrow(location_lookup) > 0) {
    anim_map <- anim_map +
      ggplot2::geom_point(
        data = location_lookup, 
        ggplot2::aes(x = long, y = lat), 
        color = "black", fill = "gray80", size = 1.8, shape = 21
      )
  }

  # Inferred Transmission Pathways (Smooth Curve Body + Single Terminal Arrowhead)
  if (nrow(path_df) > 0) {
    path_df$transmission_type <- factor(
      path_df$transmission_type, 
      levels = c("Direct Pathway (Likely)", "Indirect (Missing Intermediates)", "Distant / Import")
    )

    # Split into local pathways (viridis-colored by origin) and distant/import (black)
    local_path_df <- path_df %>% dplyr::filter(transmission_type != "Distant / Import")
    import_path_df <- path_df %>% dplyr::filter(transmission_type == "Distant / Import")

    # A. Draw LOCAL pathway curve bodies (colored by origin location)
    if (nrow(local_path_df) > 0) {
      anim_map <- anim_map +
        ggplot2::geom_path(
          data = local_path_df, 
          ggplot2::aes(
            x = x, y = y, 
            color = origin_location, 
            linetype = transmission_type, 
            alpha = transmission_type, 
            size = transmission_type,
            group = pathway_id
          )
        )
    }

    # B. Draw DISTANT / IMPORT pathway curve bodies in BLACK (fixed color, not mapped)
    if (nrow(import_path_df) > 0) {
      anim_map <- anim_map +
        ggplot2::geom_path(
          data = import_path_df, 
          ggplot2::aes(
            x = x, y = y, 
            linetype = transmission_type, 
            alpha = transmission_type, 
            size = transmission_type,
            group = pathway_id
          ),
          color = "black"
        )
    }

    # C. Extract ONLY terminal segments (last 2 points per pathway) to render EXACTLY ONE arrowhead at destination
    local_end_df <- local_path_df %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::filter(point_index >= (max(point_index) - 1)) %>%
      dplyr::ungroup()

    import_end_df <- import_path_df %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::filter(point_index >= (max(point_index) - 1)) %>%
      dplyr::ungroup()

    if (nrow(local_end_df) > 0) {
      anim_map <- anim_map +
        ggplot2::geom_path(
          data = local_end_df, 
          ggplot2::aes(
            x = x, y = y, 
            color = origin_location, 
            linetype = transmission_type, 
            alpha = transmission_type, 
            size = transmission_type,
            group = pathway_id
          ),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed")
        )
    }

    if (nrow(import_end_df) > 0) {
      anim_map <- anim_map +
        ggplot2::geom_path(
          data = import_end_df, 
          ggplot2::aes(
            x = x, y = y, 
            linetype = transmission_type, 
            alpha = transmission_type, 
            size = transmission_type,
            group = pathway_id
          ),
          color = "black",
          arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed")
        )
    }
  }

  # Import Pulses
  if (nrow(pulse_df) > 0) {
    anim_map <- anim_map +
      ggplot2::geom_point(
        data = pulse_df, 
        ggplot2::aes(
          x = x, y = y, 
          alpha = pulse_alpha,
          group = pathway_id
        ), 
        size = 5, shape = 21, color = "red3", fill = NA, stroke = 1.8
      )
  }

  # Coordinates and Theme Settings
  if (map_style == "tile") {
    anim_map <- anim_map +
      ggplot2::coord_sf(
        xlim = c(map_bounds["xmin"], map_bounds["xmax"]), 
        ylim = c(map_bounds["ymin"], map_bounds["ymax"]), 
        crs = 4326, expand = FALSE
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "gray10", color = NA),
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(), 
        axis.title = ggplot2::element_blank()
      )
  } else {
    anim_map <- anim_map +
      ggplot2::coord_fixed(
        xlim = c(map_bounds["xmin"], map_bounds["xmax"]), 
        ylim = c(map_bounds["ymin"], map_bounds["ymax"]), 
        ratio = 1.3
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"))
  }

  # 3. Apply Scales & Animation parameters
  anim_map <- anim_map +
    gganimate::transition_reveal(reveal_time) +
    gganimate::shadow_mark(past = TRUE, future = FALSE) +
    ggplot2::scale_color_viridis_d(option = "turbo", na.value = "gray50") +
    ggplot2::scale_linetype_manual(
      values = c(
        "Direct Pathway (Likely)" = "solid", 
        "Indirect (Missing Intermediates)" = "longdash", 
        "Distant / Import" = "dotdash"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "Direct Pathway (Likely)" = 0.95, 
        "Indirect (Missing Intermediates)" = 0.75, 
        "Distant / Import" = 0.80
      ),
      drop = FALSE
    ) +
    ggplot2::scale_size_manual(
      values = c(
        "Direct Pathway (Likely)" = 1.6, 
        "Indirect (Missing Intermediates)" = 1.2, 
        "Distant / Import" = 1.0
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      subtitle = "Date: {format(lubridate::date_decimal(frame_along), '%b %d, %Y')}", 
      color = "Source Location",
      linetype = "Pathway Type",
      alpha = "Pathway Type",
      size = "Pathway Type"
    ) +
    ggplot2::theme(
      legend.position = "bottom",                  
      legend.text = ggplot2::element_text(size = 8),        
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.box = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
      linetype = ggplot2::guide_legend(nrow = 1, title.position = "top"),
      alpha = ggplot2::guide_legend(nrow = 1, title.position = "top"),
      size = ggplot2::guide_legend(nrow = 1, title.position = "top")
    )

  return(anim_map)
}


#' Build Static Cumulative Transmission Map
#'
#' Creates a static (non-animated) version of the transmission map showing ALL
#' pathways and episodes simultaneously. This gives a complete overview of
#' transmission activity that is much easier to read than a GIF.
#'
#' @param episodes Data frame of residence episodes.
#' @param pathways Data frame of transition pathways.
#' @param location_lookup Lookup table of distinct locations with lat/long.
#' @param map_bounds Named numeric vector of spatial bounds (xmin, xmax, ymin, ymax).
#' @param map_style Map aesthetic style: "polygon" or "tile". Default "polygon".
#' @param distant_rendering_style How distant/import transitions are rendered. Default "import_arrow".
#'
#' @return A \code{ggplot} object (static, not animated).
#' @import ggplot2
#' @importFrom dplyr group_by filter ungroup mutate n
#' @importFrom maps map_data
#' @export
build_static_map <- function(
    episodes, 
    pathways, 
    location_lookup, 
    map_bounds, 
    map_style = "polygon",
    distant_rendering_style = "import_arrow"
) {
  message("[phymapr] Building static cumulative map...")

  # 1. Build path data for ALL pathways (same logic as animated version)
  path_list <- list()

  if (nrow(pathways) > 0) {
    for (i in seq_len(nrow(pathways))) {
      x1   <- pathways$startLon[i]
      y1   <- pathways$startLat[i]
      x2   <- pathways$endLon[i]
      y2   <- pathways$endLat[i]
      type <- pathways$transmission_type[i]
      orig <- pathways$origin_location[i]

      if (type == "Distant / Import") {
        if (distant_rendering_style == "hide") next
        if (distant_rendering_style == "import_arrow") {
          arrow_df <- interpolate_arrow(x2, y2, 0, 1, offset = 1.2, pathway_id = i)
          arrow_df$transmission_type <- type
          arrow_df$origin_location <- "External Import"
          path_list[[length(path_list) + 1]] <- arrow_df
        } else {
          bezier_df <- interpolate_bezier(x1, y1, x2, y2, 0, 1, curvature = 0.2, pathway_id = i)
          bezier_df$transmission_type <- type
          bezier_df$origin_location <- orig
          path_list[[length(path_list) + 1]] <- bezier_df
        }
      } else {
        bezier_df <- interpolate_bezier(x1, y1, x2, y2, 0, 1, curvature = 0.2, pathway_id = i)
        bezier_df$transmission_type <- type
        bezier_df$origin_location <- orig
        path_list[[length(path_list) + 1]] <- bezier_df
      }
    }
  }

  path_df <- if (length(path_list) > 0) do.call(rbind, path_list) else data.frame()

  # 2. Build static ggplot
  static_map <- ggplot2::ggplot()

  # Base map
  if (map_style == "polygon") {
    world_map <- ggplot2::map_data("world")
    us_states <- ggplot2::map_data("state")

    static_map <- static_map +
      ggplot2::geom_polygon(
        data = world_map, 
        ggplot2::aes(x = long, y = lat, group = group), 
        fill = "gray96", color = "gray85", linewidth = 0.2
      ) +
      ggplot2::geom_path(
        data = us_states, 
        ggplot2::aes(x = long, y = lat, group = group), 
        color = "white", linewidth = 0.3
      )
  }

  # Location markers (sized by number of episodes)
  if (nrow(location_lookup) > 0 && nrow(episodes) > 0) {
    loc_counts <- episodes %>%
      dplyr::group_by(location) %>%
      dplyr::summarise(n_episodes = dplyr::n(), .groups = "drop") %>%
      dplyr::left_join(location_lookup, by = "location")

    static_map <- static_map +
      ggplot2::geom_point(
        data = loc_counts, 
        ggplot2::aes(x = long, y = lat, size = n_episodes), 
        color = "black", fill = "gray60", shape = 21, alpha = 0.7
      ) +
      ggplot2::scale_size_continuous(range = c(2, 8), name = "Episodes")
  }

  # All pathways
  if (nrow(path_df) > 0) {
    path_df$transmission_type <- factor(
      path_df$transmission_type, 
      levels = c("Direct Pathway (Likely)", "Indirect (Missing Intermediates)", "Distant / Import")
    )

    # Split into local and import
    local_path_df <- path_df %>% dplyr::filter(transmission_type != "Distant / Import")
    import_path_df <- path_df %>% dplyr::filter(transmission_type == "Distant / Import")

    if (nrow(local_path_df) > 0) {
      static_map <- static_map +
        ggplot2::geom_path(
          data = local_path_df, 
          ggplot2::aes(
            x = x, y = y, 
            color = origin_location, 
            linetype = transmission_type,
            group = pathway_id
          ),
          linewidth = 0.9, alpha = 0.7
        )

      # Arrowheads for local
      local_end_df <- local_path_df %>%
        dplyr::group_by(pathway_id) %>%
        dplyr::filter(point_index >= (max(point_index) - 1)) %>%
        dplyr::ungroup()

      static_map <- static_map +
        ggplot2::geom_path(
          data = local_end_df, 
          ggplot2::aes(
            x = x, y = y, 
            color = origin_location, 
            linetype = transmission_type,
            group = pathway_id
          ),
          linewidth = 0.9, alpha = 0.7,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")
        )
    }

    if (nrow(import_path_df) > 0) {
      static_map <- static_map +
        ggplot2::geom_path(
          data = import_path_df, 
          ggplot2::aes(
            x = x, y = y, 
            linetype = transmission_type,
            group = pathway_id
          ),
          color = "black", linewidth = 0.8, alpha = 0.7
        )

      import_end_df <- import_path_df %>%
        dplyr::group_by(pathway_id) %>%
        dplyr::filter(point_index >= (max(point_index) - 1)) %>%
        dplyr::ungroup()

      static_map <- static_map +
        ggplot2::geom_path(
          data = import_end_df, 
          ggplot2::aes(x = x, y = y, linetype = transmission_type, group = pathway_id),
          color = "black", linewidth = 0.8, alpha = 0.7,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")
        )
    }
  }

  # Coordinates, scales, theme
  static_map <- static_map +
    ggplot2::coord_fixed(
      xlim = c(map_bounds["xmin"], map_bounds["xmax"]), 
      ylim = c(map_bounds["ymin"], map_bounds["ymax"]), 
      ratio = 1.3
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo", na.value = "gray50") +
    ggplot2::scale_linetype_manual(
      values = c(
        "Direct Pathway (Likely)" = "solid", 
        "Indirect (Missing Intermediates)" = "longdash", 
        "Distant / Import" = "dotdash"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Cumulative Transmission Map",
      subtitle = paste0(nrow(pathways), " pathways across ", 
                        nrow(location_lookup), " locations"),
      color = "Source Location",
      linetype = "Pathway Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "aliceblue"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.box = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
      linetype = ggplot2::guide_legend(nrow = 1, title.position = "top"),
      size = ggplot2::guide_legend(nrow = 1, title.position = "top")
    )

  return(static_map)
}

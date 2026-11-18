# LEAFLET PLOTTING FUNCTION: ----------------------------------------------

library(dplyr)
library(terra)
library(leaflet)
library(sf)   # needed to hand the grid cells to leaflet

# Helper: single-value raster -> one polygon per cell (WGS84 sf) ----------
make_grid_sf <- function(xmin, xmax, ymin, ymax, res, crs_proj, max_cells) {
  n_cells <- ((xmax - xmin) / res) * ((ymax - ymin) / res)
  if (n_cells > max_cells) {
    message(res, " m grid skipped: ", format(n_cells, big.mark = ",", scientific = FALSE),
            " cells exceeds max_grid_cells (", format(max_cells, big.mark = ","),
            "). Filter to a smaller area or raise max_grid_cells.")
    return(NULL)
  }
  r <- rast(extent = ext(xmin, xmax, ymin, ymax), resolution = res,
            crs = crs_proj, vals = 1)               # every cell = 1
  names(r) <- "value"
  # One polygon per cell (aggregate = FALSE), back to WGS84 for leaflet
  st_as_sf(project(as.polygons(r, aggregate = FALSE), "EPSG:4326"))
}

plot_bt_leaflet <- function(data, ..., colour_by = c("Species", "YEAR", "SAMPLE_DESC"),
                            radius = 5, cluster = FALSE, max_legend = 30,
                            jitter = 0,                # degrees, display only (~0.0005 = ~50 m)
                            grids = c(500, 1000),      # cell sizes in metres; NULL = no grids
                            grid_crs = "auto",         # "auto" = UTM zone of the data; or e.g. "EPSG:27700" (British National Grid)
                            grid_buffer = 2000,        # metres of grid beyond the outermost points
                            grid_fill_opacity = 0.04,  # near-transparent fill
                            grid_line_opacity = 0.5,   # cell outlines (what makes the cells visible)
                            max_grid_cells = 50000) {  # safety cap per grid layer
  
  colour_by <- match.arg(colour_by)
  
  # 1. Subset with whatever filter you pass in ----
  sub <- data %>% filter(...)
  if (nrow(sub) == 0) stop("No records match that filter.")
  
  # 2. One row per colour group x original sample coordinate ----
  pts_df <- sub %>%
    summarise(
      n_records = n(),
      n_species = n_distinct(Species),
      n_samples = n_distinct(SAMPLE_DESC),
      years     = paste(sort(unique(YEAR)), collapse = ", "),
      .by = c(assemblageID, STUDY_ID, TAXA, all_of(colour_by), LATITUDE, LONGITUDE)
    ) %>%
    arrange(across(all_of(colour_by)))
  
  # 3. SpatVector (WGS84, true coordinates, never jittered) ----
  pts <- vect(pts_df, geom = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326", keepgeom = TRUE)
  
  # 4. Leaflet map ----
  df <- as.data.frame(pts)
  xy <- crds(pts)
  if (jitter > 0) {                      # nudges the plotted position only
    set.seed(1)
    xy[, 1] <- xy[, 1] + runif(nrow(xy), -jitter, jitter)
    xy[, 2] <- xy[, 2] + runif(nrow(xy), -jitter, jitter)
  }
  df$colour_val <- df[[colour_by]]
  
  if (colour_by == "YEAR") {
    pal <- colorNumeric("viridis", domain = range(df$colour_val))
    n_levels <- 0                        # continuous legend, always shown
  } else {
    lv  <- sort(unique(df$colour_val))
    pal <- colorFactor(grDevices::hcl.colors(length(lv), "Dark 3"), domain = lv)
    n_levels <- length(lv)
  }
  
  m <- leaflet(df) %>%
    addProviderTiles("OpenStreetMap",     group = "Map") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("OpenTopoMap",       group = "Topo")
  
  # 4b. Single-value grid rasters (1 km first so 500 m sits on top) ----
  grid_groups <- character(0)
  if (length(grids) > 0) {
    if (identical(grid_crs, "auto")) {   # UTM zone from the centre of the data
      ctr <- colMeans(crds(pts))
      grid_crs <- paste0("EPSG:", (if (ctr[2] >= 0) 32600 else 32700) + floor((ctr[1] + 180) / 6) + 1)
    }
    e <- ext(project(pts, grid_crs))
    snap <- max(grids)                   # snap to the coarsest size so 500 m cells nest inside 1 km cells
    xmin <- floor((e$xmin - grid_buffer) / snap) * snap
    xmax <- ceiling((e$xmax + grid_buffer) / snap) * snap
    ymin <- floor((e$ymin - grid_buffer) / snap) * snap
    ymax <- ceiling((e$ymax + grid_buffer) / snap) * snap
    
    grid_cols <- c(`500` = "#e6550d", `1000` = "#3182bd")   # orange / blue, so both read when on together
    
    for (res in sort(grids, decreasing = TRUE)) {
      g <- make_grid_sf(xmin, xmax, ymin, ymax, res, grid_crs, max_grid_cells)
      if (is.null(g)) next
      grp <- if (res >= 1000) paste0(res / 1000, " km grid") else paste0(res, " m grid")
      col <- unname(grid_cols[as.character(res)])
      if (is.na(col)) col <- "#444444"
      m <- m %>%
        addPolygons(data = g, group = grp,
                    color = col, weight = if (res >= 1000) 1.5 else 1, opacity = grid_line_opacity,
                    fillColor = col, fillOpacity = grid_fill_opacity,
                    options = pathOptions(interactive = FALSE))   # clicks pass through to the points
      grid_groups <- c(grid_groups, grp)
    }
  }
  
  m <- m %>%
    addCircleMarkers(
      lng = xy[, 1], lat = xy[, 2],
      radius = radius, color = ~pal(colour_val), weight = 1,
      fillColor = ~pal(colour_val), fillOpacity = 0.8,
      clusterOptions = if (cluster) markerClusterOptions() else NULL,
      popup = ~paste0("<b>", colour_by, ": ", colour_val, "</b><br>",
                      "Assemblage: ", assemblageID, "<br>",
                      "Study: ", STUDY_ID, " (", TAXA, ")<br>",
                      "True lat/lon: ", round(LATITUDE, 5), ", ", round(LONGITUDE, 5), "<br>",
                      "Records: ", n_records, " | Species: ", n_species,
                      " | Samples: ", n_samples, "<br>",
                      "Years: ", years)
    ) %>%
    addLayersControl(baseGroups = c("Map", "Satellite", "Topo"),
                     overlayGroups = if (length(grid_groups)) rev(grid_groups) else NULL,   # 500 m listed first
                     options = layersControlOptions(collapsed = FALSE)) %>%
    addScaleBar(position = "bottomleft")
  
  if (n_levels <= max_legend) {
    m <- m %>% addLegend("bottomright", pal = pal, values = ~colour_val,
                         title = colour_by, labFormat = labelFormat(big.mark = ""))
  } else {
    message(n_levels, " levels of ", colour_by, ": legend skipped (raise max_legend to show it).")
  }
  
  print(m)
  invisible(pts)
}

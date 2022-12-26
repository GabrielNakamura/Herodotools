#' Plotting ancestral area character in space
#'
#' @param ada.res An object from ada function
#' @param grid An spatial object containing the cells used to calculate ancestral diversity distribution 
#'     with ada
#' @param patterns Character, a vector containing the names of the metrics to be spatialized or "all" to 
#'     create a map for all metrics. The  characters allowed to be passed are "rich", "Nnodes", "PeakDiv",
#'     "Skewness", "LowDistPeak", "HighDistPeak", "PeakRange". Default is "all"
#' @param coords A two column matrix containing the values of 
#'     Longitude (first column) and Latitude (second column)
#' @param resolution Scalar informing the resolution of grids in grid object
#' @param color_palette Character indicating the palette to be used to plot spatial patterns. it can be one
#'      of the option present in \link[rcartocolor]{scale_fill_carto_c} function. Default is "SunsetDark"
#' @param projection Character indicating the projection to be used in spatial plots.
#'      Default is "+proj=robin"
#'   
#' @return A list containing spatial plots for the metrics calculated in ada.res
plot_ada <- 
  function(ada.res,
           grid,
           coords, 
           resolution, 
           patterns = "all",
           color_palette = "SunsetDark", 
           projection = "+proj=robin")
  {
    
    pkg_req <- c("ggplot2", "raster", "rcartocolor", "sf")
    
    for(pkg in pkg_req) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
          paste0("Package '", pkg, "' must be installed to use this function."),
          call. = FALSE
        )
      }
    }
    
    # preparing data and generating raster ------------------------------------
    
    ada.res <- ada.res$Cell.Metrics
    extend_grid <- grid@bbox
    r <- raster::raster(vals = NA, xmn = extend_grid[1, 1],
                        xmx = extend_grid[1, 2],
                        ymn = extend_grid[2, 1],
                        ymx = extend_grid[2, 2], 
                        resolution = resolution
    )
    cell.r <- raster::cellFromXY(r, coords[rownames(ada.res),])
    values_cell <- rep(NA, raster::ncell(r))
    names(values_cell) <- 1:raster::ncell(r)
    val.cells <- 1:raster::ncell(r) %in% cell.r
    
    values_cell[val.cells] <- ada.res[, 1] # just richness
    r.values <- raster::setValues(r, values = values_cell)
    
    
    # transforming to sf object -----------------------------------------------
    data_sf <- sf::st_as_sf(raster::rasterToPolygons(r.values))
    data_sf$ID <- names(values_cell[val.cells])
    data_sf$Nnodes[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 2]
    data_sf$PeakDiv[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 3]
    data_sf$Skewness[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 4]
    data_sf$LowDistPeak[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 5]
    data_sf$HighDistPeak[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 6]
    data_sf$PeakRange[which(data_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 7]
    
    data_sf_proj <- 
      data_sf %>%
      sf::st_transform(crs = projection)
    
    if(patterns == "all"){
      names_div <- colnames(data_sf)[4:ncol(data_sf)]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            ggplot2::geom_sf(data = data_sf_proj, 
                             ggplot2::aes_string(geometry = "geometry", 
                                        fill = x), 
                             color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    } else{
      names_div <- colnames(data_sf)[4:ncol(data_sf)]
      patterns_plot <- pmatch(patterns, names_div)
      names_div <- names_div[patterns_plot]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            ggplot2::geom_sf(data = data_sf_proj, 
                             ggplot2::aes_string(geometry = "geometry", 
                                                 fill = x), 
                             color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    }
    res_plot
  }
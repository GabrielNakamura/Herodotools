#' Plotting ancestral area character in space
#'
#' @param ada.res An object from ada function
#' @param grid An spatial object containing the cells used to calculate ancestral diversity distribution with ada
#' @param patterns Character, a vector containing the names of the metrics to be spatialized or "all" to create a map for all metrics. The 
#'    characters allowed to be passed are "rich", "Nnodes", "PeakDiv", "Skewness", "LowDistPeak", "HighDistPeak", "PeakRange"
#' @param palette Character. The name of the palette to be used in spatial maps. It can be any of the palette from \link[rcartocolor]{scale_fill_carto_c} function
#'
#' @return A list containing spatial plots for the metrics calculated in ada.res
#' 
#' @export
#'
#' @examples
#' 
plot_ada <- 
  function(ada.res,
           grid,
           coords, 
           resolution, 
           patterns = "all",
           color_palette = "SunsetDark")
  {
    
    
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
    
    
    # transforming to sf object -----------------------------------------------
    test_sf <- sf::st_as_sf(raster::rasterToPolygons(r.n_nodes))
    test_sf$ID <- names(values_cell[val.cells])
    test_sf$rich[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 1]
    test_sf$Nnodes[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 2]
    test_sf$PeakDiv[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 3]
    test_sf$Skewness[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 4]
    test_sf$LowDistPeak[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 5]
    test_sf$HighDistPeak[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 6]
    test_sf$PeakRange[which(test_sf$ID == names(values_cell[val.cells]))] <- ada.res[, 7]
    
    if(patterns == "all"){
      names_div <- colnames(test_sf)[4:ncol(test_sf)]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            geom_sf(data = test_sf_proj, aes_string(geometry = "geometry", 
                                                    fill = x), 
                    color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    } else{
      names_div <- colnames(test_sf)[4:ncol(test_sf)]
      patterns_plot <- pmatch(patterns, names_div)
      names_div <- names_div[patterns_plot]
      res_plot <- 
        lapply(names_div, function(x){
          ggplot2::ggplot() +
            geom_sf(data = test_sf_proj, aes_string(geometry = "geometry", 
                                                    fill = x), 
                    color = "transparent", size = 0.1) +
            rcartocolor::scale_fill_carto_c(palette = color_palette)
        })
      names(res_plot) <- names_div
    }
    res_plot
  }

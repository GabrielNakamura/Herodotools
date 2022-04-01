#' Plotting ancestral area character in space
#'
#' @param ada.res An object from ada function
#' @param grid An spatial object containing the cells used to calculate ancestral diversity distribution with ada
#' @param patterns Character, a vector or a single name containing the name of variable to be plotted
#' @param palette Character. The name of the palette to be used in spatial maps. It can be any of the palettes from rcartocolor package
#'
#' @return An spatial plot
#' @export
#'
#' @examples
#' 
plot_ada <- 
  function(ada.res, grid,  coords, resolution, patterns, color_palette = "SunsetDark"){
    ada.res <- ada.res$Cell.Metrics
    extend_grid <- raster::extend(grid)
    r <- raster::raster(vals = NA, xmn = extend_grid[1],
                        xmx = extend_grid[2],
                        ymn = extend_grid[3],
                        ymx = extend_grid[4]
    )
    cell.r <- raster::cellFromXY(r, coords[rownames(ada.res),])
    values_cell <- rep(NA, raster::ncell(r))
    names(values_cell) <- 1:raster::ncell(r)
    val.cells <- 1:raster::ncell(r) %in% cell.r
    values_cell[val.cells] <- ada.res[, 2]
    r.n_nodes <- raster::setValues(r, values = values_cell)
    projcrs <- "+proj=robin"
    projection(r.n_nodes) <- projcrs
    df_r_nodes <- as.data.frame(r.n_nodes, xy = T)
    spatial_plot <- 
      ggplot() +
      geom_raster(data = na.omit(df_r_nodes), aes(x = x, y = y, fill = layer)) +
      rcartocolor::scale_fill_carto_c(palette = palette
      )
    spatial_plot <- 
      ggplot2::ggplot() +
      ggplot2::geom_raster(data = na.omit(df_r_nodes), aes(x = x, y = y, fill = layer), ) +
      rcartocolor::scale_fill_carto_c(palette = col_palette
      ) +
      labs(fill = "Node Richness")
    spatial_plot
  }

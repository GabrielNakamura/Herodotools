ada.res <- Res
dim(Res$Cell.Metrics)
coords <- coord_tyranidae
dim(coord_tyranidae)
match(rownames(coord_tyranidae), rownames(Res$Cell.Metrics))
resolution <- 1

function(ada.res, coords, resolution, palette = "SunsetDark"){
  ada.res <- ada.res$Cell.Metrics
  # box <- c(xmin = min(coords[, 1]), xmax = max(coords[, 1]), ymin = min(coords[, 2]), ymax = max(coords[, 2]))
  r <- raster::raster(vals = NA, xmn = -170.2166 , xmx = -13.21288, ymn = -55.37714, ymx = 83.6236,
                      resolution = resolution)
  r <- raster::raster(vals = NA, xmn = -170.2166 , xmx = -13.21288, ymn = -55.37714, ymx = 83.6236,
                      resolution = resolution)
  #r <- raster::raster(vals = NA, 
  #                    xmn =  min(coords[, 1]), 
  #                    xmx = max(coords[, 1]), 
  #                    ymn = min(coords[, 2]), 
  #                    ymx = max(coords[, 2]),
  #                    resolution = resolution)
  cell.r <- raster::cellFromXY(r, coords[rownames(ada.res),])
  values_cell <- rep(NA, raster::ncell(r))
  names(values_cell) <- 1:raster::ncell(r)
  val.cells <- 1:raster::ncell(r) %in% cell.r
  values_cell[val.cells] <- ada.res[, 2]
  r.n_nodes <- raster::setValues(r, values = values_cell)
  projcrs <- "+proj=robin"
  projection(r.n_nodes) <- projcrs
  df_r_nodes <- as.data.frame(r.n_nodes, xy = T)
  ggplot() +
    geom_raster(data = na.omit(df_r_nodes), aes(x = x, y = y, fill = layer)) +
    rcartocolor::scale_fill_carto_c(palette = palette
                                    )
}

quartz()
raster::plot(r.n_nodes, xlab = "Longitude", ylab = "Latitude")
library(sf)
quartz()
ggplot2::ggplot() +
  geom_tile(data = df_coord_metrics, aes(x = LON, y = LAT, fill = N_ancestral_nodes))

class(r)
st_as_sf(raster::rasterToPolygons(x = r, na.rm = FALSE))
st_as_sf(r)

length(coords)

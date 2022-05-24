#' Spatial Grid generator
#'
#'
#' @param shape A raster object
#' @param resol Numeric value indicating the resolution of grid cells
#' @param prop Numeric value indicating the minimum proportion to be contained in a cell to build the grids
#'
#' @return A raster object containing the grids with specified values
#' @export
#'
#' @examples
GridFilter <- function(shape,
                       resol = 1,
                       prop = 0)
{
  
  
  
  grid <- raster::raster(raster::extent(shape))
  terra::res(grid) <- resol
  sp::proj4string(grid) <- sp::proj4string(shape)
  gridpolygon <- raster::rasterToPolygons(grid)
  drylandproj <- sp::spTransform(shape, sp::CRS("+proj=laea"))
  gridpolproj <- sp::spTransform(gridpolygon, sp::CRS("+proj=laea"))
  gridpolproj$layer <- c(1:length(gridpolproj$layer))
  areagrid <- rgeos::gArea(gridpolproj, byid=T)
  dry.grid <- raster::intersect(drylandproj, gridpolproj)
  areadrygrid <- rgeos::gArea(dry.grid, byid=T)
  info <- cbind(dry.grid$layer, areagrid[dry.grid$layer], areadrygrid)
  dry.grid$layer <- info[,3]/info[,2]
  dry.grid <- sp::spTransform(dry.grid, sp::CRS(sp::proj4string(shape)))
  dry.grid.filtered <- dry.grid[dry.grid$layer >= prop,]
}


# testing output  ---------------------------------------------------------

test_that("plot_dispersal_from returns a ggplot object (facet = FALSE)", {
  
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")
  
  # ---- fake dispersal data ----
  dispersal_from <- data.frame(
    A = c(0.2, 0.4, 0.6),
    B = c(0.1, 0.3, 0.5)
  )
  
  # ---- coordinates ----
  coords <- data.frame(
    LONG = c(-60, -58, -55),
    LAT  = c(-10, -12, -15)
  )
  
  # ---- minimal sf object ----
  poly <- sf::st_polygon(list(matrix(
    c(-70, -20,
      -40, -20,
      -40,   0,
      -70, -20),
    ncol = 2,
    byrow = TRUE
  )))
  
  shapefile <- sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
  
  p <- plot_dispersal_from(
    dispersal_from = dispersal_from,
    coords = coords,
    shapefile = shapefile,
    area_cols = c(1, 2),
    area_col  = 1,
    facet = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_dispersal_from works with facet = TRUE", {
  
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")
  
  dispersal_from <- data.frame(
    A = c(0.2, 0.4, 0.6),
    B = c(0.1, 0.3, 0.5),
    C = c(0.3, 0.2, 0.1)
  )
  
  coords <- data.frame(
    lon = c(-60, -58, -55),
    lat = c(-10, -12, -15)
  )
  
  poly <- sf::st_polygon(list(matrix(
    c(-70, -20,
      -40, -20,
      -40,   0,
      -70, -20),
    ncol = 2,
    byrow = TRUE
  )))
  
  shapefile <- sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
  
  p <- plot_dispersal_from(
    dispersal_from = dispersal_from,
    coords = coords,
    shapefile = shapefile,
    area_cols = c(1, 3),
    facet = TRUE
  )
  
  expect_s3_class(p, "ggplot")
})

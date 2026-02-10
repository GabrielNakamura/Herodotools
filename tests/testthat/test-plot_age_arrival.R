

# Output validation test --------------------------------------------------


test_that("plot_age_arrival returns a ggplot object", {
  
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")
  
  # ---- fake age_arrival_comm ----
  age_arrival_comm <- list(
    mean_age_per_assemblage = data.frame(mean_age_arrival = c(1.2, 2.4, 3.1))
  )
  
  # ---- coordinates ----
  coords <- data.frame(
    LONG = c(-60, -58, -55),
    LAT  = c(-10, -12, -15)
  )
  
  # ---- minimal sf object (triangle polygon) ----
  poly <- sf::st_polygon(list(matrix(
    c(-70, -20,
      -40, -20,
      -40,   0,
      -70, -20),
    ncol = 2,
    byrow = TRUE
  )))
  
  shapefile <- sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
  
  p <- plot_age_arrival(
    age_arrival_comm = age_arrival_comm,
    coords = coords,
    shapefile = shapefile
  )
  
  expect_s3_class(p, "ggplot")
})


# Input validation test ---------------------------------------------------


test_that("plot_age_arrival errors if coords is not a data.frame", {
  
  age_arrival_comm <- list(
    mean_age_per_assemblage = data.frame(mean_age_arrival = c(1.2, 2.4, 3.1))
  )
  
  expect_error(
    plot_age_arrival(
      age_arrival_comm = age_arrival_comm,
      coords = matrix(1:6, ncol = 2),
      shapefile = NULL
    ),
    "coords must be a data.frame"
  )
})


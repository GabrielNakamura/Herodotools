# Plotting ancestral area character in space

Plotting ancestral area character in space

## Usage

``` r
plot_ada(
  ada.res,
  grid,
  coords,
  resolution,
  patterns = "all",
  color_palette = "SunsetDark",
  projection = "+proj=robin"
)
```

## Arguments

- ada.res:

  An object from ada function

- grid:

  An spatial object containing the cells used to calculate ancestral
  diversity distribution with ada

- coords:

  A two column matrix containing the values of Longitude (first column)
  and Latitude (second column)

- resolution:

  Scalar informing the resolution of grids in grid object

- patterns:

  Character, a vector containing the names of the metrics to be
  spatialized or "all" to create a map for all metrics. The characters
  allowed to be passed are "rich", "Nnodes", "PeakDiv", "Skewness",
  "LowDistPeak", "HighDistPeak", "PeakRange". Default is "all"

- color_palette:

  Character indicating the palette to be used to plot spatial patterns.
  it can be one of the option present in
  [scale_fill_carto_c](https://jakubnowosad.com/rcartocolor/reference/carto_scale.html)
  function. Default is "SunsetDark"

- projection:

  Character indicating the projection to be used in spatial plots.
  Default is "+proj=robin"

## Value

A list containing spatial plots for the metrics calculated in ada.res

# Occupancy probability of each node in each site considering a kernel density funcion

Occupancy probability of each node in each site considering a kernel
density funcion

## Usage

``` r
comp_kernel_gravel(x, phy, coords, w_slope, min_disp_prob)
```

## Arguments

- x:

  Node by communities matrix

- phy:

  phylo object

- coords:

  Data frame with two columns containing spatial coordinate of
  communities

- w_slope:

  A scalar representing the slope of the dispersal kernel function. This
  is the same as passed to `calc_sbears` function

- min_disp_prob:

  A scalar between 0 and 1 representing the minimum probability of
  receiving an immigrant node from a focal cell. Same as passed to
  `calc_sbears` function

## Value

A list with length equal the number of sites. Each element of the list
contains a matrix with nodes in rows and sites in columns and the
probability of occurrence of each ancestral node in each site.

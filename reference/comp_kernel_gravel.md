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

  data frame with spatial coordinate of communities

## Value

A list with length equal the number of sites. Each element of the list
contains a matrix with nodes in rows and sites in columns and the
probability of occurrence of each ancestral node in each site.

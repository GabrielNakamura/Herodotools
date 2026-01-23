# Axiliary function to compute node composition by sites

Axiliary function to compute node composition by sites

## Usage

``` r
comp_ada_nodes_sites(phy, comm, long = FALSE)
```

## Arguments

- phy:

  A phylo object

- comm:

  A community matrix, rows are communities and columns are species

- long:

  Logical, if TRUE dense format of the output matrix is converted to
  long format

## Value

A matrix (or data frame in the long format) with community composition
by nodes

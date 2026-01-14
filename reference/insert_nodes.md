# Insert internal nodes in the phylogenetic tree

It adds non-bifurcating nodes in the branches of the phylogenetic tree
to represent events of ancestral character changes, specifically change
in range states from biogeographic ancestral reconstruction.

## Usage

``` r
insert_nodes(tree, inserts, node_area = NULL)
```

## Arguments

- tree:

  An object of class 'phylo'.

- inserts:

  a data frame create with the function
  [`get_insert_df()`](https://gabrielnakamura.github.io/Herodotools/reference/get_insert_df.md)

- node_area:

  default NULL. One column data frame indicating the range area of each
  node (rows)

## Value

a list with 2 elements:

- `$phylo`: a phylogenetic tree with added non-bifurcating nodes

- `$node_area`: a data frame with the node area for each added node. If
  the argument `node_area` is not NULL, all the node area are included
  in the data frame.

# Auxiliary function to compute information of node path and dispersal path for each species

Auxiliary function to compute information of node path and dispersal
path for each species

## Usage

``` r
get_nodes_info_core(W, tree, ancestral.area, biogeo)
```

## Arguments

- W:

  Numerical matrix, rows are assemblages and columns are species

- tree:

  Phylogenetic tree in newick format

- ancestral.area:

  Single column data frame with nodes in rows and ancestral
  area/Ecoregion of occurrence in columns. If the reconstruction comes
  from BioGeoBears this data frame can be obtained with
  ([`get_node_range_BioGeoBEARS`](https://gabrielnakamura.github.io/Herodotools/reference/get_node_range_BioGeoBEARS.md))

- biogeo:

  Data frame containing the information of the biome/area/ecoregion of
  each assemblage

## Value

A list with auxiliary information on nodes and species

## Examples

``` r
if (FALSE) { # \dontrun{
# hypothetical occurrence matrix with species in columns and assemblages in lines
 W_toy<- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0),
     nrow= 3,ncol= 5,
     dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
     c(paste("s", 1:5, sep=""))))
 
 #toy tree
 data(toy_treeEx)
 
 # hypothetical data indicating the ecoregions of each assemblage
 biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
 
 # hypothetical data indicating the ancestral range of each node
 ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
 
get_nodes_info_core(W=W_toy,tree=toy_treeEx,ancestral.area=ancestral_area_toy,biogeo=biogeo_toy)
} # }
```

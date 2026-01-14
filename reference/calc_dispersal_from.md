# Proportional contribution of an ancestral area to the species in other region

For the species present in a site, the function calculates the
proportional contribution of ancestral areas to dispersal of lineage in
the site's region.

## Usage

``` r
calc_dispersal_from(W, tree, ancestral.area, biogeo)
```

## Arguments

- W:

  Occurrence matrix, rows are assemblages and columns are species

- tree:

  Phylogenetic tree in newick format

- ancestral.area:

  One column data frame with nodes in rows and one column indicating the
  occurrence (biome/ecoregion) area of nodes

- biogeo:

  One column data frame with assemblages in rows and their respective
  biome/ecoregion

## Value

A data frame with assemblages in the rows and regions in the columns.
The values indicates the percentage of contribution of each region to
each assemblage. NA represents no contribution

## See also

[`get_node_range_BioGeoBEARS`](https://gabrielnakamura.github.io/Herodotools/reference/get_node_range_BioGeoBEARS.md)

## Author

Arthur V Rodrigues <rodrigues.arthur.v@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(akodon_sites) # occurrence matrix
akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
data(akodon_newick) # phylogenetic tree
spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
akodon_pa_tree <- akodon_pa[, spp_in_tree]
data(regions) # biogeographic region
data(resDEC) # output from ancestral area reconstruction
node.area <- get_node_range_BioGeoBEARS(resDEC,
                                        phyllip.file = here::here("inst", 
                                        "extdata", "geo_area_akodon.data")
                                        ,akodon_newick,max.range.size = 4) 
calc_dispersal_from(W=akodon_pa_tree,
                    tree=akodon_newick,
                    ancestral.area=node.area,biogeo=regions) # historical dispersal analysis
} # }
```

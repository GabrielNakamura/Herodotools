# Affiliation of assemblages based on phylogenetic turnover

Affiliation of assemblages based on phylogenetic turnover

## Usage

``` r
calc_affiliation_evoreg(phylo.comp.dist, groups)
```

## Arguments

- phylo.comp.dist:

  A distance matrix indicating the phylogenetic (or
  taxonomic/functional) distance composition among assemblages

- groups:

  A character vector indicating the group of each assemblage. This
  object can be obtained with
  [`calc_evoregions`](https://gabrielnakamura.github.io/Herodotools/reference/calc_evoregions.md)

## Value

A list with two matrix, one containing affiliation values and the group
in which each cell is classified and the other containing cell
coordinates

## Details

The function calculates the degree of affiliation of each community to
the region (or evoregion), in which that community was classified. If
used coupled with a analysis of Principal Coordinates of Phylogenetic
Structure (PCPS) to represent the phylogenetic distances the analysis
output of the analysis will represent the degree of membership of
assemblages to each evoregion, as described in Maestri and Duarte 2020.

## References

Maestri, R and Duarte L.d.S. (2020). Evoregions: Mapping shifts in
phylogenetic turnover across biogeographic regions. Methods in Ecology
and Evolution, 11, 1652-1662.

## Examples

``` r
if (FALSE) { # \dontrun{
# First run the classification
#' data(akodon_newick) # phylogenetic tree
data(akodon_sites) # occurrence data 
akodon_pa <- akodon_sites %>% 
    dplyr::select(-LONG, -LAT)
spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
akodon_pa_tree <- akodon_pa[, spp_in_tree]
regions <- calc_evoregions(comm = akodon_pa, phy = akodon_newick)
site_region <- regions$Cluster_Evoregions # classification of each community in regions

axis_sel <- which(regions$PCPS$prop_explainded >= 
    regions$PCPS$tresh_dist) # significant PCPS axis
PCPS_thresh <- regions$PCPS$vectors[, axis_sel] # only significant axis
dist_phylo_PCPS <- vegan::vegdist(PCPS_thresh, method = "euclidean") # distance matrix 
calc_affiliation_evoreg(phylo.comp.dist = dist_phylo_PCPS,
    groups = regions$Cluster_Evoregions) # affiliation
} # }
```

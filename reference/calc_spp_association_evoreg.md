# Quantifying species association with evoregions

Quantifying species association with evoregions

## Usage

``` r
calc_spp_association_evoreg(comm, group.comm, spp.association = 0.7)
```

## Arguments

- comm:

  An assemblage/community matrix with communities in rows and species in
  columns

- group.comm:

  A vector containing the classification of each community in groups

- spp.association:

  A scalar indicating the level of association that will be evaluated to
  assign species to groups. Default is 0.7

## Value

A matrix with species in the rows and one columns indicating the name of
the group in which the species is associated given the level of species
association used in the argument spp.association. If the species fail to
reach the threshold it will be assigned the "fuzzy" status

## References

Maestri, R. and Duarte, L.d.S. (2020). Evoregions: Mapping shifts in
phylogenetic turnover across biogeographic regions. Methods Ecol. Evol.,
11, 1652-1662.

## Author

Leandro Duarte, Renan Maestri and Gabriel Nakamura
<gabriel.nakamura.souza@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data(akodon_sites) # occurrence data 
akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
data(akodon_newick) # phylogenetic tree
spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
akodon_pa_tree <- akodon_pa[, spp_in_tree]
regions <- calc_evoregions(comm = akodon_pa, phy = akodon_newick) # compute evoregions
groups <- regions$cluster_evoregions
calc_spp_association_evoreg(comm = akodon_pa_tree, 
                       group.comm = groups, 
                       spp.association = 0.7) # calculate species association
} # }
```

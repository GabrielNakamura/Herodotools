# Species occurrence in each evoregion

Auxiliary function to produce an occurrence data frame of species in
each region. This object can be with
[`get_tipranges_to_BioGeoBEARS`](https://gabrielnakamura.github.io/Herodotools/reference/get_tipranges_to_BioGeoBEARS.md)
function to produce a phyllip file needed in BioGeoBEARS ancestral area
reconstruction

## Usage

``` r
get_region_occ(comm, site.region)
```

## Arguments

- comm:

  Occurrence data frame with species in colums and rows corresponding to
  assemblages

- site.region:

  a character vector indicating the membership of each assemblage to
  regions

## Value

An occurrence data frame with species as rownames and regions as columns

## See also

[`get_tipranges_to_BioGeoBEARS`](https://gabrielnakamura.github.io/Herodotools/reference/get_tipranges_to_BioGeoBEARS.md)

## Author

Gabriel Nakamura and Arthur Rodrigues

## Examples

``` r
data(akodon_sites)
data(akodon_newick)
akodon_pa <- akodon_sites %>% dplyr::select(-LONG, -LAT)
spp_in_tree <- names(akodon_pa) %in% akodon_newick$tip.label
akodon_pa_tree <- akodon_pa[, spp_in_tree]
data(regions)
site_evoregion <- regions$Cluster_Evoregions
get_region_occ(comm=akodon_pa_tree,site.region=site_evoregion)
#>               A B C D E F G H I
#> A_aerosus     1 1 0 0 0 0 0 0 0
#> A_affinis     0 1 0 0 0 0 0 0 0
#> A_albiventer  1 1 0 0 0 0 0 0 0
#> A_azarae      0 0 1 0 0 1 0 0 1
#> A_boliviensis 1 0 0 0 0 0 0 0 0
#> A_budini      1 0 0 0 0 0 0 0 0
#> A_cursor      0 0 0 0 0 0 0 1 0
#> A_dayi        1 0 0 1 0 0 0 0 0
#> A_dolores     0 0 1 1 0 0 0 0 0
#> A_fumeus      1 0 0 0 0 0 0 0 0
#> A_iniscatus   0 0 0 1 1 0 0 0 0
#> A_juninensis  1 0 0 0 0 0 0 0 1
#> A_kofordi     1 0 0 0 0 0 0 0 0
#> A_lindberghi  1 0 0 0 0 0 0 1 0
#> A_lutescens   1 0 0 0 0 0 0 0 0
#> A_mimus       1 0 0 0 0 0 0 0 0
#> A_mollis      0 1 0 0 0 0 0 0 0
#> A_montensis   0 0 0 0 0 0 0 1 0
#> A_mystax      0 0 0 0 0 0 0 1 0
#> A_orophilus   0 1 0 0 0 0 0 0 0
#> A_paranaensis 0 0 0 0 0 1 0 1 0
#> A_reigi       0 0 0 0 0 1 0 0 0
#> A_siberiae    1 0 0 0 0 0 0 0 0
#> A_simulator   1 0 0 0 0 0 0 0 0
#> A_spegazzinii 1 0 1 0 0 0 0 0 0
#> A_subfuscus   1 0 0 0 0 0 0 0 0
#> A_sylvanus    1 0 0 0 0 0 0 0 0
#> A_toba        0 0 1 0 0 0 1 0 0
#> A_torques     1 0 0 0 0 0 0 0 0
#> A_varius      1 0 0 0 1 0 0 0 0
```

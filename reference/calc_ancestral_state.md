# Ancestral State per Assemblage

Ancestral State per Assemblage

## Usage

``` r
calc_ancestral_state(tree, ancestral.area, prefix = "N")
```

## Arguments

- tree:

  A newick phylogenetic tree object

- ancestral.area:

  One column data frame. Lines are nodes the column are the
  biomes/region of occurrence for each ancestors. Can be obtained by
  using
  ([`get_node_range_BioGeoBEARS`](https://gabrielnakamura.github.io/Herodotools/reference/get_node_range_BioGeoBEARS.md))

- prefix:

  A single character string to be used to name nodes

## Value

A data frame with assemblages in lines and nodes in columns. Each cell
contains the ancestral area/Ecoregion of occurrence for each node and
its respective species.

## Author

Gabriel Nakamura <gabriel.nakamura.souza@gmail.com>

## Examples

``` r
biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
calc_ancestral_state(toy_treeEx, ancestral_area_toy)
#>    s1    s2    s3    s4    s5   
#> N6 "ABC" "ABC" "ABC" "ABC" "ABC"
#> N7 "B"   "B"   "B"   NA    NA   
#> N8 "C"   "C"   NA    NA    NA   
#> N9 NA    NA    NA    "ABC" "ABC"
```

# Get the node area states from the biogeographical stochastic mapping's result

Get the node area states from the biogeographical stochastic mapping's
result

## Usage

``` r
get_bsm_node_area(
  bsm,
  BioGeoBEARS.data,
  phyllip.file,
  tree.path,
  max.range.size
)
```

## Arguments

- bsm:

  results from the function
  [`calc_bsm()`](https://gabrielnakamura.github.io/Herodotools/reference/calc_bsm.md)

- BioGeoBEARS.data:

  a BioGeoBEARS result model object

- phyllip.file:

  path to the phyllip file used in the original model

- tree.path:

  path for the phylogenetic tree used in the original model

- max.range.size:

  maximum range size used in the biogeographical reconstruction

## Value

a list with same length as the number of stochastic mapping. Each list
element is a single column data frame with node area states

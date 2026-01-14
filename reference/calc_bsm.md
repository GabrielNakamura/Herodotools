# calculate Biogeographical Stochastic Mapping (BSM) with BioGeoBEARS

This is a wrapper function to the BioGeoBEARS::runBSM(). It uses the
result of a BioGeoBEARS model as input to generate realized range area
changes across the tree nodes and branches. This mapping is done
stochastically based on the model parameters.

## Usage

``` r
calc_bsm(
  BioGeoBEARS.data,
  phyllip.file,
  tree.path,
  max.maps = 100,
  n.maps.goal = 50,
  seed = NULL,
  save_after_every_try = FALSE,
  ...
)
```

## Arguments

- BioGeoBEARS.data:

  a BioGeoBEARS result model object

- phyllip.file:

  path to the phyllip file used in the original model

- tree.path:

  path for the phylogenetic tree used in the original model

- max.maps:

  maximum number of stochastic maps to try to generate (passed to
  `maxnum_maps_to_try` in `runBSM()` from `{BioGeoBEARS}`)

- n.maps.goal:

  Target number of successfully generated stochastic maps

- seed:

  Optional seed for reproducibility of stochastic simulations.

- save_after_every_try:

  Logical. If `TRUE`, results will be saved after each attempt (default
  is `FALSE`).

- ...:

  Additional arguments passed to `runBSM()` from `{BioGeoBEARS}`

## Value

A list (BSM output) containing the simulated mappings, including
ancestral range estimates and event counts. This output is directly
produced by `runBSM()` from `{BioGeoBEARS}`.

## References

Matzke, N. J. (2013). Probabilistic historical biogeography: new models
for founder-event speciation, imperfect detection, and fossils allow
improved accuracy and model-testing. *Frontiers of Biogeography*, 5(4).

Dupin, J., et al. (2017). Bayesian estimation of the global
biogeographical history of the Solanaceae. *Journal of Biogeography*,
44(4), 887–899.

## Examples

``` r
if (FALSE) { # \dontrun{
data("akodon_newick")
data("resDEC")  # BioGeoBEARS model result

bsm_result <- calc_bsm(
  BioGeoBEARS.data = resDEC,
  phyllip.file = system.file("extdata", "geo_area_akodon.data", package = "Herodotools"),
  tree.path = system.file("extdata", "akodon.new", package = "Herodotools"),
  max.maps = 50,
  n.maps.goal = 10,
  seed = 1234
)
} # }
```

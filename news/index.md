# Changelog

## Herodotools 2.0.1

- New `calc_sbears` function to compute site-based estimation of
  ancestral range of species

- New `calc_bsm` function, a wrapper function used to generate realized
  change areas across the tree nodes and branches extracted from
  BioGeoBEARS model and BioGeoBEARS::runBSM()

- `calc_insitu_diversification` and `calc_age_arrival` now supports
  output from `calc_bsm` function

- `find_max_nclust` now supports parallel computation

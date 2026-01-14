# Tip-based in-situ diversification metrics

Tip-based in-situ diversification metrics

## Usage

``` r
calc_insitu_diversification(
  W,
  tree,
  ancestral.area,
  biogeo,
  type = "equal.splits"
)
```

## Arguments

- W:

  A species occurrence (assemblage) matrix. Rows represent assemblages
  (sites) and columns represent species.

- tree:

  A phylogenetic tree object (class `phylo`).

- ancestral.area:

  A one-column data frame indicating the ecorregion of occurrence of
  each node (rows).

- biogeo:

  A data frame with one column indicating the ecoregion of each
  assemblage (row).

- type:

  Character string indicating the type of evolutionary distinctiveness
  (ED) metric to be used. Options are `"equal.splits"` (default) or
  `"fair.proportion"`.

## Value

A named list containing the following elements (if available):

- `jetz_site_sp`: Site-by-species matrix of DR diversification rates.

- `jetz_comm_mean`: Harmonic mean of DR diversification rates per site.

- `insitu_site_sp`: Site-by-species matrix of \\DR_in-situ\\
  diversification rates.

- `insitu_comm_mean`: Harmonic mean of \\DR_in-situ\\ diversification
  rates per site (ignoring zeros when specified).

- `prop_site_sp`: Site-by-species matrix of proportional diversification
  rates.

- `prop_comm_mean`: Arithmetic mean of proportional diversification
  rates per site.

## Details

This function calculates three related diversification metrics:

- **Jetz diversification rate (DR_jetz):** The inverse of the
  evolutionary distinctiveness (ED) metric for each species across the
  full phylogeny. This is the same as described in [Jetz et al
  (2012)](https://www.nature.com/articles/nature11631)

          \deqn{
                 DR_{i} = (\sum_{j=1}^{N_i}\cdot l_j \cdot 1/2^{j-1})^{-1}
              }

- **In-situ diversification rate (DR_insitu):** In situ measure of
  diversification calculated as the inverse of ED, but restricted to the
  portion of a species' branches that diversified in the region
  (informed in biogeo argument) matching the assemblage.

- **Proportional diversification rate (DR_prop):** The fraction of
  diversification that occurred in a given ecoregion (informed in biogeo
  argument), calculated as the proportion of ED in that region relative
  to the total ED, multiplied by the species’ total diversification rate
  (DR_jetz).

For assemblage-level summaries, site-by-species matrices are generated
for each metric, and site-level means are calculated as:

- **Harmonic mean (Jetz and in-situ):** computed across species in each
  assemblage. Zeros in the in-situ matrix are ignored in the harmonic
  mean calculation to account only for species with diversification in
  situ.

- **Arithmetic mean (proportional rates):** used for proportional
  diversification rates.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with a random tree and toy data
library(ape)
tree <- rcoal(5)
W <- matrix(c(1,0,1,1,1,0,1,1,0,1 ), nrow=2, byrow=TRUE,
            dimnames = list(c("site1","site2"), tree$tip.label))
biogeo <- data.frame(region=c("A","B"))
ancestral.area <- data.frame(region=rep("A", tree$Nnode))

res <- calc_insitu_diversification(W, tree, ancestral.area, biogeo)
str(res)
} # }
```

# Tip-based metrics of trait evolution

Tip-based metrics of trait evolution

## Usage

``` r
calc_tip_based_trait_evo(
  tree,
  trait,
  nsim = 1,
  method = c("transition_rates", "last_transition_time", "stasis_time")
)
```

## Arguments

- tree:

  A phylogenetic tree as an object of class "phylo"

- trait:

  A named vector containing the tip states for a discretely valued
  character. The names must match the tip labels of the tree.

- nsim:

  Number of simulations to stochastic character mapping.

- method:

  Tip-based metric, partial match to "transition_rates",
  "last_transition_time" and "stasis_time".

## Value

A list (length equal to nsim) with tip-based metrics estimated per
species.

## Details

This function calculates three metrics that represent the
macroevolutionary dynamic of traits for each species in a phylogenetic
tree. They are Transition rates, last transition time and stasis time,
all originally proposed in Luza et al (2021). *Transition rates*
indicate how many times the ancestral character has changed over time.
*Stasis time* indicates the maximum branch length (time interval) which
the current tip-character was maintained across the whole phylogeny.
Finally, *last transition time* is the sum of branch lengths from the
tip to the previous node with a reconstructed character equal to the
current tip-character. Originally the function uses a stochastic
character mapping on discrete trait data, but other types of trait
reconstruction can be used.

## Author

André Luza and Vanderlei Debastiani

## Examples

``` r
 
if (FALSE) { # \dontrun{
data(rodent.phylo) # phylogenetic tree
data(trait) # categorical traits on species diet
trans_rates <- calc_tip_based_trait_evo(tree=match_data$phy,trait =trait ,
    nsim = 50,method = c("transition_rates", "last_transition_time", "stasis_time")
    )
} # }
```

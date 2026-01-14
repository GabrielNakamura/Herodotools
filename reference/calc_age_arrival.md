# Computes the arrival ages of each species in the assemblages

Computes the arrival ages of each species in the assemblages

## Usage

``` r
calc_age_arrival(W, tree, ancestral.area, biogeo, age.no.ancestor = "recent")
```

## Arguments

- W:

  Occurrence matrix, rows are assemblages and columns are species

- tree:

  Phylogenetic tree in newick format

- ancestral.area:

  One column data frame, nodes in row and one column containing the
  occurrence of the nodes

- biogeo:

  One column data frame, assemblages in rows and one column containing
  the region in which the assemblage is located

- age.no.ancestor:

  a character string "recent" (default) or "half.edge" indicating how to
  deal with cases in which the most recent aancestor of a tip node are
  not on the same region. "recent" attributes a small age (10e-5), while
  "half.edge" attributes the age as half the length of the branch
  linking the ancestor to the tip.

## Value

A list of length two.

- age_arrival:

  A matrix with the arrival age of each species in each assemblage

- mean_age_arrival:

  mean age values for each assemblage

## Details

This function computes the mean arrival age of species in an assemblage
based in an ancestral area reconstruction. For each assemblage we
calculate the arrival and establishment of each ancestor of present-day
species as showed in Van Dijk et al (2021). We consider an arrival event
as being the arrival of an ancestor and the establishment of this
ancestor until the present day. Arrival events that occurred between the
last speciation event of that lineage and the present time can be
assigned a small value or the age corresponding to the half of the
lenght of the branch linking the ancestor and the present day species

## References

Van Dijk, A.; Nakamura G,; Rodrigues, A.V.; Maestri, R. and Duarte,
L.d.S. (2021). Imprints of tropical niche conservatism and historical
dispersal in the radiation of Tyrannidae (Aves:Passeriformes). Biol.
Journ. Linnean Soc., 134, 57-67.

## Author

Gabriel Nakamura <gabriel.nakamura.souza@gmail.com> and Arthur Rodrigues

## Examples

``` r
 # hypothetical occurrence matrix with species in columns and assemblages in lines
 W_toy<- matrix(c(0, 1, 1, 0, 1,
                  1, 0, 1, 0, 1, 
                  0, 0, 1, 0, 0),nrow= 3,
                  ncol= 5,dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
                  c(paste("s", 1:5, sep=""))))
 
 #toy tree
 data(toy_treeEx)
 
 # hypothetical data indicating the ecoregions of each assemblage
 biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
 
 # hypothetical data indicating the ancestral range of each node
 ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
 
 # caculating age of each assemblage
 age_assemblages <- calc_age_arrival(W_toy, toy_treeEx, ancestral_area_toy, biogeo_toy)
```

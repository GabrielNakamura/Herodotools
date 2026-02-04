# Phylogenetic diversity decomposition

Phylogenetic diversity decomposition

## Usage

``` r
PD_decomposition(comm, sbears.obj, phy, threshold, make.node.label = TRUE)
```

## Arguments

- comm:

  a community matrix with sites in rows and species in columns

- sbears.obj:

  an object calculated from
  [`calc_sbears`](https://gabrielnakamura.github.io/Herodotools/reference/calc_sbears.md)
  containing ancestral community reconstruction

- phy:

  A phylo object containing phylogenetic relationships among species

- threshold:

  Scalar indicating the threshold used to consider the presence of a
  species in a community

- make.node.label:

  Logical. If TRUE (the default) the nodes of the phylogeny will be
  named with the letter N preceding the number of the node

## Value

a list of length two.

- `decomp_potential`: is a list with results of PD decomposition
  considering the range area model reconstruction. This list contains:

- `PD_decomposition`: A data frame with three colums. `partition` is a
  character indicating the component of PD, `value` is a numeric columns
  with respective PD compoent and `community` is a character column with
  the name of the community

       \item \code{decomp_faith}: is a list with the same structure of
           \code{decomp_potential}, but contains the results of PD decomposition
           considering only the current occurrence of lineages in communities

## Examples

``` r
# phylogenetic tree
set.seed(42)
phylo <- ape::rcoal(n = 10)
phylo <- ape::makeNodeLabel(phy = phylo)

# community composition matrix
comm <- 
  matrix(sample(c(0, 1), size = 10*20, replace = TRUE),
         nrow = 20, 
         ncol = 10, 
         dimnames = list(paste("comm", 1:20), phylo$tip.label)
  )

# coordinates - this is necessary for "disperal_assembly" algorithm
xy_coords <- 
  matrix(runif(1:10), 
         nrow = nrow(comm),
         ncol = 2, 
         dimnames = list(rownames(comm), c("lng", "lat")
         )
  )

# running sbears with "single_site" algorithm
out_sbears <- 
  calc_sbears(x = comm,
              phy = phylo, 
              coords = xy_coords, 
              method = "single_site")

# sbears object can be used directly in PD_decomposition function

out_decomp <- 
  PD_decomposition(comm = comm,
                   sbears.obj = out_sbears, 
                   phy = phylo, 
                   threshold = 0.5)
#>   |                                                          |                                                  |   0%  |                                                          |==                                                |   5%  |                                                          |=====                                             |  10%  |                                                          |========                                          |  15%  |                                                          |==========                                        |  20%  |                                                          |============                                      |  25%  |                                                          |===============                                   |  30%  |                                                          |==================                                |  35%  |                                                          |====================                              |  40%  |                                                          |======================                            |  45%  |                                                          |=========================                         |  50%  |                                                          |============================                      |  55%  |                                                          |==============================                    |  60%  |                                                          |================================                  |  65%  |                                                          |===================================               |  70%  |                                                          |======================================            |  75%  |                                                          |========================================          |  80%  |                                                          |==========================================        |  85%  |                                                          |=============================================     |  90%  |                                                          |================================================  |  95%  |                                                          |==================================================| 100%

# Data frame with results from PD decomposition
out_pd_decomposition <- out_decomp$decomp_potential$PD_decomposition
    
    
```

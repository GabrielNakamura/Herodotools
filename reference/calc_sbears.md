# Site-Based Estimation of Ancestral Range of Species

Site-Based Estimation of Ancestral Range of Species

## Usage

``` r
calc_sbears(
  x,
  phy,
  coords,
  method = c("single_site", "disp_assembly"),
  w_slope = 5,
  min_disp_prob = 0.8,
  compute.node.by.sites = TRUE,
  make.node.label = TRUE
)
```

## Arguments

- x:

  Matrix object with community/assemblage composition. Rows represent
  sites and columns represents species.

- phy:

  A phylogenetic tree of `phylo` class object.

- coords:

  A rectangular object (it can be a matrix, data.frame or tibble)
  containing two columns, one with longitude and another with latitude
  coordinates in this order.

- method:

  character indicating how ancestral range probabilities are computed.
  The options are "single_site" or "disp_assembly"

- w_slope:

  A scalar representing the slope of the dispersal kernel function.

- compute.node.by.sites:

  Logical, TRUE (default) computes a matrix of node occurrence by site.

- make.node.label:

  Logical, if TRUE (default) the nodes of the phylogeny will be named as
  the letter "N" preceding node number

## Value

A list with two elements.

- `reconstruction`: a matrix with ancestral nodes in rows and
  assemblages in columns. The numeric values represent the occupation
  probability of each ancestral node in each assemblage.

      \item \code{site_node_composition}: a matrix with assemblage in rows and
          node names in columns. The values represent the presence of a given
          node in a given assemblage based on the occurrence of current species
          in the assemblages.

## Details

SBEARS (Site-Based Estimation of Ancestral Range of Species) is a method
for ancestral state reconstruction of species geographic ranges. It
operates at fine spatial resolution and does not require predefined
discrete biogeographic areas. The method estimates ancestral ranges
under a maximum likelihood framework, using site-level occurrence
information. For large matrix this function supports parallel
computation with package `future`

- the method used to compute ancestral range probabilities can be
  `single_site` or `disp_assembly`. For both methods SBEARS use the
  function [`fastAnc`](https://rdrr.io/pkg/phytools/man/fastAnc.html) to
  reconstruct the the ancestral ranges. In `single_site` each site is
  taken as a binary categorical trait representing the presence
  (state 1) or the absence (state 0) of every species in the site, which
  is reconstructed along the nodes of the phylogeny. In `disp_assembly`
  option the single site reconstruction is carried out, producing a
  matrix describing standardized probabilities of occurrence of each
  node at each site. The difference here is that an additional step is
  performed to compute the probability of a cell being colonized by an
  ancestral species originating from a focal cell. The probability of
  colonization decreases increasing the distance from the focal cell.
  This probability of an ancestral species to disperse from a cell to
  any other can be defined using a dispersal kernel function proposed by
  Gravel et al (2006).

## References

Gravel D., Canhan C.D., Beaudet M. and Messier C. Reconciling niche and
neutrality: the continuum hypothesis. 2006. Ecology Letters
[doi:10.1111/j.1461-0248.2006.00884.x](https://doi.org/10.1111/j.1461-0248.2006.00884.x)

## Examples

``` r
# phylogenetic tree
phylo <- geiger::sim.bdtree(n = 10, seed = 42)
#> Error in loadNamespace(x): there is no package called ‘geiger’
phylo <- ape::makeNodeLabel(phy = phylo)
#> Error: object 'phylo' not found

# community composition matrix
comm <- 
 matrix(sample(c(0, 1), size = 10*20, replace = TRUE),
        nrow = 20, 
        ncol = 10, 
        dimnames = list(paste("comm", 1:20), phylo$tip.label)
 )
#> Error: object 'phylo' not found

# coordinates - this is necessary for "disperal_assembly" algorithm
xy_coords <- 
 matrix(runif(1:10), 
        nrow = nrow(comm),
        ncol = 2, 
        dimnames = list(rownames(comm), c("lng", "lat")
        )
 )
#> Error: object 'comm' not found

# running sbears with "single_site" algorithm
out_sbears <- 
 calc_sbears(x = comm,
             phy = phylo, 
             coords = xy_coords, 
             method = "single_site")
#> Error: object 'phylo' not found

# matrix containing ancestral reconsturction, nodes are rows and columns are sites
anc_reconstruction <- out_sbears$reconstruction
#> Error: object 'out_sbears' not found
```

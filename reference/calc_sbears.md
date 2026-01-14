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
  compute.node.by.sites = FALSE,
  make.node.label = FALSE
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

a list with three elements. reconstruction is the result of ancestral .
area reconstruction; phylogeny is the matrix containing the occurrence
of nodes in sites and joint.phylo.obs is the joint occurrence of nodes
and species in phylogenetic tree.

## Details

SBEARS (Site-Based Estimation of Ancestral Range of Species) is a method
for ancestral state reconstruction of species geographic ranges. It
operates at fine spatial resolution and does not require predefined
discrete biogeographic areas. The method estimates ancestral ranges
under a maximum likelihood framework, using site-level occurrence
information.

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

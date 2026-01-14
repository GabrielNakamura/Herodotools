# Ancestral Diversity Analysis

This function computes different assemblage level metrics using an
ancestral community reconstruction

## Usage

``` r
ada(
  x,
  phy,
  sp.bin = "Sturges",
  marginal = FALSE,
  lik.threshold = FALSE,
  threshold = 0.7,
  compute.fields = F
)
```

## Arguments

- x:

  Community occurrence matrix. Rows are sites and columns are species

- phy:

  Phylogenetic tree

- sp.bin:

  Character indicating the methods to be used to compute the number of
  time slices in which metrics will be computed. Default is "Sturges"

- marginal:

  Logical indicating

- lik.threshold:

  Logical indicating if some threshold value will be used to select the
  likelihood values of ancestor species in a site. Default is TRUE

- threshold:

  Scalar indicating the threshold value used to select the likelihood of
  a given species be present at a site

- compute.fields:

  Logical, indicates if phylogenetic fields will be computed. Default is
  FALSE

## Value

A list of size eight containing the following objects

- Phylogeny:

  Phylogenetic tree

- Root.Age:

  Numeric containing root age of the tree used

- Per.node.ancestral.area:

  Matrix containing the ancestors (nodes) for each species

- Diversity.Through.Time:

  Age of each node in the phylogeny

- Cell.Metric:

  Matrix containing the assemblage level metrics for richness and
  ancestral diversity

## Details

ada function calculates an ancestral assemblage reconstruction for each
assemblage in a matrix using for this reconstruction
[ace](https://rdrr.io/pkg/ape/man/ace.html) function from ape package

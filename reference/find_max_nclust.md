# Estimate the maximum number of groups in DAPC analysis

Estimate the maximum number of groups in DAPC analysis

## Usage

``` r
find_max_nclust(
  x,
  threshold,
  max.nclust,
  nperm = 100,
  method = "kmeans",
  stat = "BIC",
  criterion = "diffNgroup",
  subset = 100,
  confidence.level = c(0.7, 0.8, 0.9, 0.95, 0.99)
)
```

## Arguments

- x:

  A data.frame or matrix object containing eigenvectors by sites.

- threshold:

  Scalar. The number of eigenvectors used to perform classification.

- max.nclust:

  A vector containing values of the maximum number of groups to be
  evaluated.

- nperm:

  Scalar. Number of times classification will be performed.

- method:

  Character, one of c("kmeans","ward"). This will be used in
  `find.clusters` function. See `find.clusters` of adegenet package.
  Default is "kmeans"

- stat:

  Character, one of c("BIC", "AIC", or "WSS"). This will be used in
  `find.clusters` function. See `find.clusters` of adegenet package.
  Default is "BIC".

- criterion:

  Character one of c("diffNgroup", "min","goesup", "smoothNgoesup", or
  "goodfit"). This will be used in `find.clusters` function. Default is
  "diffNgroup". See `find.clusters` of adegenet package.

- subset:

  Scalar. The number of cells used in the analysis. It is particularly
  important whenever the total number of cells is large (\> 1000).

- confidence.level:

  A vector containing values with threshold confidence level used to
  estimate congruence in the classification pattern.

## Value

Matrix containing congruence values ranging between 0-1 for each
max.nclust value (see Arguments) and confidence level.

## Details

This function can be used to find the maximum number of clusters that
maximizes the congruence of the grouping procedure. Basically it
consists in repeat DAPC analysis multiple times using a set of candidate
maximum values (indicated in the argument `max.nclust`)

## Examples

``` r
if (FALSE) { # \dontrun{
data(regions)
evovectors <- regions$PCPS$vectors # eigenvectors by site
find_max_nclust(x = evovectors, threshold = 3, max.nclust = 10)
} # }
```

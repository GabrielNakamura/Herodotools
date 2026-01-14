# Auxiliary function used to define the occurrence threshold used to define node presences

Auxiliary function used to define the occurrence threshold used to
define node presences

## Usage

``` r
find_threshold(
  x,
  threshold.steps = c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
)
```

## Arguments

- x:

  A matri containing the continuous probabilities of occurrences of
  nodes accross the sites

- threshold.steps:

  A numeric vector defining the threshold values that will be compared
  in the procedure

## Value

A list object containing the following information:

    \itemize{
        \item Maximum_site_probability: a vector containing the maximum probability
            of occurrence of nodes for ech site
        \item Removed_sites
        \item Correlation_between_steps
        \item Select_threshold: a scalar indicating the selected threshold
    }

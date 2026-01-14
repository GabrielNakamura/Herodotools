# Calculate the Harmonic Mean

Computes the harmonic mean of a numeric vector. The harmonic mean is
defined as \\n / \sum (1/x_i)\\, where \\n\\ is the number of
observations.

## Usage

``` r
calc_harmonic_mean(x, na.rm = FALSE, ignore.zero = FALSE)
```

## Arguments

- x:

  A numeric vector.

- na.rm:

  Logical; if `TRUE`, remove `NA` values before computation.

- ignore.zero:

  Logical; if `TRUE`, remove zeros before computation.

## Value

A numeric value representing the harmonic mean.

## Details

- If `ignore.zero = FALSE` and the vector contains zeros, the harmonic
  mean is zero by definition.

- If all values are removed (due to `NA`s or zeros), the function
  returns `NA`.

## Examples

``` r
# Simple harmonic mean
calc_harmonic_mean(c(1, 2, 3))
#> [1] 1.636364

# With a zero (harmonic mean collapses to 0)
calc_harmonic_mean(c(1, 2, 0))
#> [1] 0

# Ignoring zeros
calc_harmonic_mean(c(1, 2, 0), ignore.zero = TRUE)
#> [1] 1.333333

# With NA values
calc_harmonic_mean(c(1, 2, NA), na.rm = TRUE)
#> [1] 1.333333

# Combining both options
calc_harmonic_mean(c(1, 2, 0, NA), na.rm = TRUE, ignore.zero = TRUE)
#> [1] 1.333333
```

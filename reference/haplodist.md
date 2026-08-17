# haplodist

Function to extract haplotypes and compute pairwise distances between
haplotypes. This function use the functions `haplotype` of package pegas
and [`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html) of package
ape.

## Usage

``` r
haplodist(x, dist.model = "N", ...)
```

## Arguments

- x:

  A list with the set of DNA sequences (as an object of class "DNAbin"
  or "haplotype") as used by the function `haplotype`.

- dist.model:

  A character string used by the function
  [`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html) to specify the
  evolutionary model to be used to compute pairwise distances from DNA
  sequences (default dist.model = "N").

- ...:

  Additional arguments to the function
  [`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html).

## Value

A list with:

- call:

  Arguments used.

- haplotypes:

  A list with haplotypes indices that identify each observation sharing
  the same haplotype.

- individual.per.haplotype:

  A matrix with individuals per haplotype.

- haplotype.distances:

  A matrix with pairwise distances between haplotypes.

## See also

`HaploVectors`

## Examples

``` r
data(segv)
haplodist(segv$segv.fas)
#> $call
#> haplodist(x = segv$segv.fas)
#> 
#> $haplotypes
#> $haplotypes$haplotype.I
#> [1]  1  3 10 11 14 15
#> 
#> $haplotypes$haplotype.II
#> [1] 2 4 5 6 7 8 9
#> 
#> $haplotypes$haplotype.III
#> [1] 12
#> 
#> $haplotypes$haplotype.IV
#> [1] 13
#> 
#> $haplotypes$haplotype.V
#> [1] 16 17 18
#> 
#> 
#> $individual.per.haplotype
#>          haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> ind01s01           1            0             0            0           0
#> ind02s01           0            1             0            0           0
#> ind03s01           1            0             0            0           0
#> ind01s02           0            1             0            0           0
#> ind02s02           0            1             0            0           0
#> ind03s02           0            1             0            0           0
#> ind01s03           0            1             0            0           0
#> ind02s03           0            1             0            0           0
#> ind03s03           0            1             0            0           0
#> ind01s04           1            0             0            0           0
#> ind02s04           1            0             0            0           0
#> ind03s04           0            0             1            0           0
#> ind01s05           0            0             0            1           0
#> ind02s05           1            0             0            0           0
#> ind03s05           1            0             0            0           0
#> ind01s06           0            0             0            0           1
#> ind02s06           0            0             0            0           1
#> ind03s06           0            0             0            0           1
#> 
#> $haplotype.distances
#>               haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> haplotype.I             0            1             2            3           3
#> haplotype.II            1            0             1            2           2
#> haplotype.III           2            1             0            3           3
#> haplotype.IV            3            2             3            0           2
#> haplotype.V             3            2             3            2           0
#> 
```

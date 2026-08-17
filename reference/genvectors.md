# GenVectors

Function to extract haplotypic/genetic eigenvectors and perform null
model-based tests.

## Usage

``` r
genvectors(
  pop,
  distances,
  checkdata = TRUE,
  dist.model = "N",
  log.frequencies = FALSE,
  method = "euclidean",
  squareroot.dis = TRUE,
  choices = c(1, 2),
  analysis = "none",
  envir,
  formula,
  runs = 999,
  ...
)

# S3 method for class 'genvectors'
print(x, ...)
```

## Arguments

- pop:

  A matrix describing the incidence of each individual (columns) in a
  given locality (rows).

- distances:

  Matrix containing genetic distances between individuals or a list with
  the set of DNA sequences (class "DNAbin" or "haplotype") as used by
  the function `haplotype`.

- checkdata:

  Logical argument (TRUE or FALSE) to check if individual sequences in
  the pop data follow the same order as in the set of DNA sequences
  (Default checkdata = TRUE).

- dist.model:

  A character string used by the function
  [`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html) to specify the
  evolutionary model to be used to computes pairwise distances from DNA
  sequences (default dist.model = "N").

- log.frequencies:

  Logical argument (TRUE or FALSE) to specify if transformation of
  natural logarithms plus one in haplotype per locality data must be
  applied (Default log.frequencies = FALSE).

- method:

  Dissimilarity index to apply in matrix P, which describes localities
  by their haplotypic/genetic composition, as accepted by vegdist
  function in vegan package (Default method = "euclidean").

- squareroot.dis:

  Logical argument (TRUE or FALSE) to specify if use square root of
  dissimilarity index in matrix P (Default squareroot.dis = TRUE).

- choices:

  Axes for re-scaling. Choices must have length equal to two (Default
  choices = c(1, 2)).

- analysis:

  Type of analysis, partial match to "none", "adonis" or "glm" (Default
  analysis = "none").

- envir:

  A matrix with environmental variables for each population, with
  variables as columns and localities as rows. See Details and Examples.

- formula:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html).
  Used in "adonis" or "glm" analysis. See Details and Examples.

- runs:

  Number of permutations for assessing probability of type I error.

- ...:

  Additional arguments to function `matrix.p.sig` and `pcps.sig`.

## Value

A list with:

- call:

  Arguments used.

- haplotypes:

  A list with haplotypes index that identify each observation that share
  the same haplotype.

- genetic.distances:

  A matrix with pairwise genetic/haplotypes distances.

- individual.per.haplotype:

  A matrix with individuals per haplotype.

- genetic.per.locality:

  A matrix with frequency of each genetic/haplotype per locality
  (**\\W\\**).

- vectors:

  Haplotypic/genetic eigenvectors.

- values:

  Eigenvalues, relative eigenvalues and cumulative relative eigenvalues.

- correlations:

  Correlations between haplotypic/genetic eigenvectors and
  haplotypes/alleles.

- P:

  Matrix of haplotypic/genetic composition (**\\P\\**).

- scores:

  Scores for biplots.

- model:

  The observed model.

- fun:

  The funtion used.

- statistic.null.turnover:

  A matrix with null statistic for turnover null model.

- statistic.null.divergence:

  A matrix with null statistic for divergence null model.

- statistic.obs:

  Observed statistic, F value to predefined function.

- p.turnover:

  The p value for the turnover null model.

- p.divergence:

  The p value for the divergence null model.

## Details

genvectors is a function to extract haplotypic/genetic eigenvectors and
perform null model-based tests. Input parameters can be entered in two
ways. The *distances* argument can be supplied as a set of DNA sequences
(class "DNAbin" or "haplotype") as used by the function `haplotype`,
from which pairwise distances are calculated using
[`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html). Alternatively,
it can be provided directly as a genetic distance matrix between
individuals.

The argument *analysis* specifies the type of analysis performed. When
*analysis* is set to "adonis" the analysis is performed on a matrix of
haplotypic/genetic composition (using the `matrix.p.sig` function). The
argument *formula* must be specified, where the left-hand side gives the
resemblance data, right-hand side gives the variables. The resemblance
data is internally named *p.dist*, thus formula is an expression of the
form *p.dist ~ predictors*. If *analysis* is set to "glm" it is
performed with geneticvector (using the `pcps.sig` function). In this
case, the argument *formula* must also be specified, where the left-hand
side gives the vectors used, right-hand side gives the variables. The
vectors are internally named sequentially *geneticvector.1*,
*geneticvector.2*, *geneticvector.3*, and so on. Thus, formula is an
expression of the form *geneticvector.1 ~ predictors*.

## See also

[`haplodist`](https://gabrielnakamura.github.io/Herodotools/reference/haplodist.md),
`matrix.p.sig`, `pcps.sig`

## Examples

``` r
data(segv)

genvectors(segv$segv.pi, segv$segv.fas, envir = segv$segv.envir,
           choices = c(1,2))
#> $call:
#> genvectors(pop = segv$segv.pi, distances = segv$segv.fas, choices = c(1,      2), envir = segv$segv.envir) 
#> 
#> $haplotypes:
#> $haplotype.I
#> [1]  1  3 10 11 14 15
#> 
#> $haplotype.II
#> [1] 2 4 5 6 7 8 9
#> 
#> $haplotype.III
#> [1] 12
#> 
#> $haplotype.IV
#> [1] 13
#> 
#> $haplotype.V
#> [1] 16 17 18
#> 
#> 
#> $individual.per.haplotype:
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
#> $genetic.distances:
#>               haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> haplotype.I             0            1             2            3           3
#> haplotype.II            1            0             1            2           2
#> haplotype.III           2            1             0            3           3
#> haplotype.IV            3            2             3            0           2
#> haplotype.V             3            2             3            2           0
#> $genetic.per.locality:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01           2            1             0            0           0
#> s02           0            3             0            0           0
#> s03           0            3             0            0           0
#> s04           2            0             1            0           0
#> s05           2            0             0            1           0
#> s06           0            0             0            0           3
#> 
#> $vectors:
#>     geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> s01      0.19948623      0.18619387     -0.01180442      0.16492892
#> s02      0.06920902     -0.22038410      0.08709279      0.01032386
#> s03      0.06920902     -0.22038410      0.08709279      0.01032386
#> s04      0.22769103      0.21434391      0.12254724     -0.13229538
#> s05      0.07757363     -0.04102324     -0.30643453     -0.05338445
#> s06     -0.64316893      0.08125368      0.02150611      0.00010318
#> 
#> $values:
#>                 Eigenvalue Relative_eig Cumul_eig
#> geneticvector.1 0.52090168   0.59234074 0.5923407
#> geneticvector.2 0.18603484   0.21154859 0.8038893
#> geneticvector.3 0.12469211   0.14179301 0.9456823
#> geneticvector.4 0.04776669   0.05431765 1.0000000
#> 
#> $correlations:
#>               geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> haplotype.I         0.9286417      0.32658578     -0.17258124      0.03438087
#> haplotype.II        0.9475382     -0.17357486      0.25298497      0.08967594
#> haplotype.III       0.8676974     -0.08828383      0.46398684     -0.15499502
#> haplotype.IV       -0.6826800     -0.43028287     -0.58791392     -0.05623070
#> haplotype.V        -0.9967075      0.01835329      0.06482892      0.04510659
#> 
#> $P:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01   0.4074074    0.3333333     0.1851852   0.03703704  0.03703704
#> s02   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s03   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s04   0.3888889    0.3333333     0.2777778   0.00000000  0.00000000
#> s05   0.3333333    0.2888889     0.1111111   0.20000000  0.06666667
#> s06   0.0000000    0.2000000     0.0000000   0.20000000  0.60000000
#> 
#> $scores:
#>               geneticvector.1 geneticvector.2
#> haplotype.I         0.2114434     0.070001674
#> haplotype.II        0.2157459    -0.038253141
#> haplotype.III       0.1975669    -0.019456352
#> haplotype.IV       -0.4390786    -0.094827505
#> haplotype.V        -0.6410513     0.003933917

genvectors(segv$segv.pi, segv$segv.fas, analysis = "adonis",
           envir = segv$segv.envir, formula = p.dist~R, runs = 99)
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> Set of permutations < 'minperm'. Generating entire set.
#> $call:
#> genvectors(pop = segv$segv.pi, distances = segv$segv.fas, analysis = "adonis",      envir = segv$segv.envir, formula = p.dist ~ R, runs = 99) 
#> 
#> $haplotypes:
#> $haplotype.I
#> [1]  1  3 10 11 14 15
#> 
#> $haplotype.II
#> [1] 2 4 5 6 7 8 9
#> 
#> $haplotype.III
#> [1] 12
#> 
#> $haplotype.IV
#> [1] 13
#> 
#> $haplotype.V
#> [1] 16 17 18
#> 
#> 
#> $individual.per.haplotype:
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
#> $genetic.distances:
#>               haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> haplotype.I             0            1             2            3           3
#> haplotype.II            1            0             1            2           2
#> haplotype.III           2            1             0            3           3
#> haplotype.IV            3            2             3            0           2
#> haplotype.V             3            2             3            2           0
#> $genetic.per.locality:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01           2            1             0            0           0
#> s02           0            3             0            0           0
#> s03           0            3             0            0           0
#> s04           2            0             1            0           0
#> s05           2            0             0            1           0
#> s06           0            0             0            0           3
#> 
#> $vectors:
#>     geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> s01      0.19948623      0.18619387     -0.01180442      0.16492892
#> s02      0.06920902     -0.22038410      0.08709279      0.01032386
#> s03      0.06920902     -0.22038410      0.08709279      0.01032386
#> s04      0.22769103      0.21434391      0.12254724     -0.13229538
#> s05      0.07757363     -0.04102324     -0.30643453     -0.05338445
#> s06     -0.64316893      0.08125368      0.02150611      0.00010318
#> 
#> $values:
#>                 Eigenvalue Relative_eig Cumul_eig
#> geneticvector.1 0.52090168   0.59234074 0.5923407
#> geneticvector.2 0.18603484   0.21154859 0.8038893
#> geneticvector.3 0.12469211   0.14179301 0.9456823
#> geneticvector.4 0.04776669   0.05431765 1.0000000
#> 
#> $correlations:
#>               geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> haplotype.I         0.9286417      0.32658578     -0.17258124      0.03438087
#> haplotype.II        0.9475382     -0.17357486      0.25298497      0.08967594
#> haplotype.III       0.8676974     -0.08828383      0.46398684     -0.15499502
#> haplotype.IV       -0.6826800     -0.43028287     -0.58791392     -0.05623070
#> haplotype.V        -0.9967075      0.01835329      0.06482892      0.04510659
#> 
#> $P:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01   0.4074074    0.3333333     0.1851852   0.03703704  0.03703704
#> s02   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s03   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s04   0.3888889    0.3333333     0.2777778   0.00000000  0.00000000
#> s05   0.3333333    0.2888889     0.1111111   0.20000000  0.06666667
#> s06   0.0000000    0.2000000     0.0000000   0.20000000  0.60000000
#> 
#> $scores:
#>               geneticvector.1 geneticvector.2
#> haplotype.I         0.2114434     0.070001674
#> haplotype.II        0.2157459    -0.038253141
#> haplotype.III       0.1975669    -0.019456352
#> haplotype.IV       -0.4390786    -0.094827505
#> haplotype.V        -0.6410513     0.003933917
#> 
#> $model:
#> Permutation test for adonis under reduced model
#> Marginal effects of terms
#> Permutation: free
#> Number of permutations: 719
#> 
#> vegan::adonis2(formula = formula, data = data.frame(envir), permutations = 2, by = "margin", parallel = NULL)
#>          Df SumOfSqs      R2      F Pr(>F)
#> R         1  0.22208 0.25254 1.3514      1
#> Residual  4  0.65732 0.74746              
#> Total     5  0.87940 1.00000              
#> 
#> $obs.statistic:
#> [1] 1.351422
#> 
#> $p.turnover:
#> [1] 0.32
#> 
#> $p.divergence:
#> [1] 0.9

genvectors(segv$segv.pi, segv$segv.fas, analysis = "glm",
           envir = segv$segv.envir, formula = geneticvector.1~R, runs = 99)
#> $call:
#> genvectors(pop = segv$segv.pi, distances = segv$segv.fas, analysis = "glm",      envir = segv$segv.envir, formula = geneticvector.1 ~ R, runs = 99) 
#> 
#> $haplotypes:
#> $haplotype.I
#> [1]  1  3 10 11 14 15
#> 
#> $haplotype.II
#> [1] 2 4 5 6 7 8 9
#> 
#> $haplotype.III
#> [1] 12
#> 
#> $haplotype.IV
#> [1] 13
#> 
#> $haplotype.V
#> [1] 16 17 18
#> 
#> 
#> $individual.per.haplotype:
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
#> $genetic.distances:
#>               haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> haplotype.I             0            1             2            3           3
#> haplotype.II            1            0             1            2           2
#> haplotype.III           2            1             0            3           3
#> haplotype.IV            3            2             3            0           2
#> haplotype.V             3            2             3            2           0
#> $genetic.per.locality:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01           2            1             0            0           0
#> s02           0            3             0            0           0
#> s03           0            3             0            0           0
#> s04           2            0             1            0           0
#> s05           2            0             0            1           0
#> s06           0            0             0            0           3
#> 
#> $vectors:
#>     geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> s01      0.19948623      0.18619387     -0.01180442      0.16492892
#> s02      0.06920902     -0.22038410      0.08709279      0.01032386
#> s03      0.06920902     -0.22038410      0.08709279      0.01032386
#> s04      0.22769103      0.21434391      0.12254724     -0.13229538
#> s05      0.07757363     -0.04102324     -0.30643453     -0.05338445
#> s06     -0.64316893      0.08125368      0.02150611      0.00010318
#> 
#> $values:
#>                 Eigenvalue Relative_eig Cumul_eig
#> geneticvector.1 0.52090168   0.59234074 0.5923407
#> geneticvector.2 0.18603484   0.21154859 0.8038893
#> geneticvector.3 0.12469211   0.14179301 0.9456823
#> geneticvector.4 0.04776669   0.05431765 1.0000000
#> 
#> $correlations:
#>               geneticvector.1 geneticvector.2 geneticvector.3 geneticvector.4
#> haplotype.I         0.9286417      0.32658578     -0.17258124      0.03438087
#> haplotype.II        0.9475382     -0.17357486      0.25298497      0.08967594
#> haplotype.III       0.8676974     -0.08828383      0.46398684     -0.15499502
#> haplotype.IV       -0.6826800     -0.43028287     -0.58791392     -0.05623070
#> haplotype.V        -0.9967075      0.01835329      0.06482892      0.04510659
#> 
#> $P:
#>     haplotype.I haplotype.II haplotype.III haplotype.IV haplotype.V
#> s01   0.4074074    0.3333333     0.1851852   0.03703704  0.03703704
#> s02   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s03   0.2222222    0.3333333     0.2222222   0.11111111  0.11111111
#> s04   0.3888889    0.3333333     0.2777778   0.00000000  0.00000000
#> s05   0.3333333    0.2888889     0.1111111   0.20000000  0.06666667
#> s06   0.0000000    0.2000000     0.0000000   0.20000000  0.60000000
#> 
#> $scores:
#>               geneticvector.1 geneticvector.2
#> haplotype.I         0.2114434     0.070001674
#> haplotype.II        0.2157459    -0.038253141
#> haplotype.III       0.1975669    -0.019456352
#> haplotype.IV       -0.4390786    -0.094827505
#> haplotype.V        -0.6410513     0.003933917
#> 
#> $model:
#> 
#> Call:  stats::glm(formula = formula, data = data.frame(x, envir))
#> 
#> Coefficients:
#> (Intercept)            R  
#>     -0.2911       0.1455  
#> 
#> Degrees of Freedom: 5 Total (i.e. Null);  4 Residual
#> Null Deviance:       0.5209 
#> Residual Deviance: 0.4362    AIC: 7.299
#> 
#> $obs.statistic:
#> [1] 0.7768179
#> 
#> $p.turnover:
#> [1] 0.55
#> 
#> $p.divergence:
#> [1] 0.91

```

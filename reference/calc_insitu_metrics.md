# Diversification-based PD and PE

Diversification-based PD and PE

## Usage

``` r
calc_insitu_metrics(W, tree, ancestral.area, biogeo, PD = TRUE, PE = TRUE)
```

## Arguments

- W:

  Occurrence matrix, rows are assemblages and columns are species

- tree:

  Phylogenetic hipothesis in newick format

- ancestral.area:

  One column data frame with nodes in rows and one column indicating the
  occurrence area of nodes

- biogeo:

  One columns data frame with rows indicating assemblages and one column
  indicating the biome/ecoregion of each assemblage

- PD:

  Logical, if TRUE (default) PD~in-situ~ will be computed

- PE:

  Logical, if TRUE (default) PE~in-situ~ will be computed

## Value

A data frame containing the values of original PD and PE and also their
in-situ diversification based version

## Details

This function computes two modified versions of Phylogenetic Diversity
(PD) and Phylogenetic endemism by considering in the calculation the
in-situ diversification. the PD~in-situ~ and PE~in-situ~ are calculated
by coupling in the calculation an ancetral area reconstruction that
divide the phylogenetic tree in two parts, one that correspond to
dispersion events and another that correspond to in-situ diversification
events. We use the in-situ diversification part to calculate both Db-PD
and Db-PE

## Author

Gabriel Nakamura <gabriel.nakamura.souza@gmail.com>

## Examples

``` r
W_toy<- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0),
nrow= 3,
ncol= 5,
dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
c(paste("s", 1:5, sep=""))))
 data("toy_treeEx")
 biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
 ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
 assemblage_phylo_metrics <- calc_insitu_metrics(W_toy,
      toy_treeEx,
      ancestral_area_toy, 
      biogeo_toy)
```

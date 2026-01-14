# Species and their respective ancestral nodes

This function computes a matrix containing species and their respective
ancestral nodes

## Usage

``` r
get_spp_nodes(tree, node.prefix = "N")
```

## Arguments

- tree:

  Phylogenetic tree

- node.prefix:

  Character indicating the prefix to be used to name nodes, default is
  the letter "N"

## Value

A matrix with species in lines and nodes in columns. 1 indicates that
the node corresponds to the ancestor of that species and zero indicates
that species does not share the ancestor node

## Details

This is an auxiliary function that, by taking a phylogenetic tree,
allows to represent in a matrix format the relationship between
present-day species (tips) and their ancestors (nodes). The function
provides an easy way to obtain certain objects that can be used in other
functions, like
[`calc_ancestral_state`](https://gabrielnakamura.github.io/Herodotools/reference/calc_ancestral_state.md)

## Examples

``` r
data(akodon_newick)
get_spp_nodes(tree = akodon_newick, node.prefix = "N")
#>     A_mimus A_lindberghi A_subfuscus A_lutescens A_sylvanus A_boliviensis
#> N31       1            1           1           1          1             1
#> N32       0            1           1           1          1             1
#> N33       0            1           1           1          1             1
#> N34       0            0           1           1          1             1
#> N35       0            0           1           1          1             1
#> N36       0            0           1           1          1             1
#> N37       0            0           1           1          1             1
#> N38       0            0           1           1          0             0
#> N39       0            0           0           0          1             1
#> N40       0            0           0           0          0             1
#> N41       0            0           0           0          0             0
#> N42       0            0           0           0          0             0
#> N43       0            0           0           0          0             0
#> N44       0            0           0           0          0             0
#> N45       0            0           0           0          0             0
#> N46       0            0           0           0          0             0
#> N47       0            0           0           0          0             0
#> N48       0            0           0           0          0             0
#> N49       0            0           0           0          0             0
#> N50       0            0           0           0          0             0
#> N51       0            0           0           0          0             0
#> N52       0            0           0           0          0             0
#> N53       0            0           0           0          0             0
#> N54       0            0           0           0          0             0
#> N55       0            0           0           0          0             0
#> N56       0            0           0           0          0             0
#> N57       0            0           0           0          0             0
#> N58       0            0           0           0          0             0
#> N59       0            0           0           0          0             0
#>     A_spegazzinii A_azarae A_fumeus A_kofordi A_juninensis A_cursor A_montensis
#> N31             1        1        1         1            1        1           1
#> N32             1        1        1         1            1        1           1
#> N33             1        1        1         1            1        1           1
#> N34             1        1        1         1            1        1           1
#> N35             1        1        1         1            1        0           0
#> N36             1        1        0         0            0        0           0
#> N37             1        0        0         0            0        0           0
#> N38             0        0        0         0            0        0           0
#> N39             1        0        0         0            0        0           0
#> N40             1        0        0         0            0        0           0
#> N41             0        0        1         1            1        0           0
#> N42             0        0        1         1            0        0           0
#> N43             0        0        0         0            0        1           1
#> N44             0        0        0         0            0        0           1
#> N45             0        0        0         0            0        0           0
#> N46             0        0        0         0            0        0           0
#> N47             0        0        0         0            0        0           0
#> N48             0        0        0         0            0        0           0
#> N49             0        0        0         0            0        0           0
#> N50             0        0        0         0            0        0           0
#> N51             0        0        0         0            0        0           0
#> N52             0        0        0         0            0        0           0
#> N53             0        0        0         0            0        0           0
#> N54             0        0        0         0            0        0           0
#> N55             0        0        0         0            0        0           0
#> N56             0        0        0         0            0        0           0
#> N57             0        0        0         0            0        0           0
#> N58             0        0        0         0            0        0           0
#> N59             0        0        0         0            0        0           0
#>     A_reigi A_paranaensis A_mystax A_aerosus A_affinis A_mollis A_orophilus
#> N31       1             1        1         1         1        1           1
#> N32       1             1        1         1         1        1           1
#> N33       1             1        1         0         0        0           0
#> N34       1             1        1         0         0        0           0
#> N35       0             0        0         0         0        0           0
#> N36       0             0        0         0         0        0           0
#> N37       0             0        0         0         0        0           0
#> N38       0             0        0         0         0        0           0
#> N39       0             0        0         0         0        0           0
#> N40       0             0        0         0         0        0           0
#> N41       0             0        0         0         0        0           0
#> N42       0             0        0         0         0        0           0
#> N43       1             1        1         0         0        0           0
#> N44       1             1        1         0         0        0           0
#> N45       1             1        1         0         0        0           0
#> N46       0             1        1         0         0        0           0
#> N47       0             0        0         1         1        1           1
#> N48       0             0        0         1         1        1           1
#> N49       0             0        0         1         1        1           1
#> N50       0             0        0         1         1        1           1
#> N51       0             0        0         1         1        0           0
#> N52       0             0        0         0         0        1           1
#> N53       0             0        0         0         0        1           1
#> N54       0             0        0         0         0        0           0
#> N55       0             0        0         0         0        0           0
#> N56       0             0        0         0         0        0           0
#> N57       0             0        0         0         0        0           0
#> N58       0             0        0         0         0        0           0
#> N59       0             0        0         0         0        0           0
#>     A_torques A_budini A_siberiae A_simulator A_varius A_albiventer A_iniscatus
#> N31         1        1          1           1        1            1           1
#> N32         1        1          1           1        1            1           1
#> N33         0        0          0           0        0            0           0
#> N34         0        0          0           0        0            0           0
#> N35         0        0          0           0        0            0           0
#> N36         0        0          0           0        0            0           0
#> N37         0        0          0           0        0            0           0
#> N38         0        0          0           0        0            0           0
#> N39         0        0          0           0        0            0           0
#> N40         0        0          0           0        0            0           0
#> N41         0        0          0           0        0            0           0
#> N42         0        0          0           0        0            0           0
#> N43         0        0          0           0        0            0           0
#> N44         0        0          0           0        0            0           0
#> N45         0        0          0           0        0            0           0
#> N46         0        0          0           0        0            0           0
#> N47         1        1          1           1        1            1           1
#> N48         1        1          1           1        1            1           0
#> N49         1        1          1           1        1            0           0
#> N50         1        0          0           0        0            0           0
#> N51         0        0          0           0        0            0           0
#> N52         1        0          0           0        0            0           0
#> N53         0        0          0           0        0            0           0
#> N54         0        1          1           1        1            0           0
#> N55         0        1          1           0        0            0           0
#> N56         0        0          0           1        1            0           0
#> N57         0        0          0           0        0            0           1
#> N58         0        0          0           0        0            0           0
#> N59         0        0          0           0        0            0           0
#>     A_dayi A_toba A_dolores
#> N31      1      1         1
#> N32      1      1         1
#> N33      0      0         0
#> N34      0      0         0
#> N35      0      0         0
#> N36      0      0         0
#> N37      0      0         0
#> N38      0      0         0
#> N39      0      0         0
#> N40      0      0         0
#> N41      0      0         0
#> N42      0      0         0
#> N43      0      0         0
#> N44      0      0         0
#> N45      0      0         0
#> N46      0      0         0
#> N47      1      1         1
#> N48      0      0         0
#> N49      0      0         0
#> N50      0      0         0
#> N51      0      0         0
#> N52      0      0         0
#> N53      0      0         0
#> N54      0      0         0
#> N55      0      0         0
#> N56      0      0         0
#> N57      1      1         1
#> N58      1      1         1
#> N59      0      1         1
```

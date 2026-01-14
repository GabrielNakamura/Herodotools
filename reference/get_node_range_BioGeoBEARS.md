# Get node ranges from BioGeoBEARS biome reconstruction model

This is an auxiliary function that can be with a ancestral area
reconstruction from BioGeoBEARS to obtain a data frame containing the
ancestral area of occurrence (biome range) for all ancestors (nodes)

## Usage

``` r
get_node_range_BioGeoBEARS(
  BioGeoBEARS.data,
  phyllip.file,
  tree,
  max.range.size
)
```

## Arguments

- BioGeoBEARS.data:

  An object containing the result of BioGeoBEARS reconstruction

- phyllip.file:

  Phylip file, the same used in BioGeoBEARS package

- tree:

  A phylogenetic tree

- max.range.size:

  Scalar indicating the maximum number of biomes in which each node can
  occupy

## Value

Data frame with two columns, one indicating the node and other the
ancestral biome range

## Examples

``` r
data(resDEC) # ancestral area reconstruction from biogeobears
data(akodon_newick)
extdata_dir = system.file("extdata", package="Herodotools")
fn = paste(extdata_dir, "/geo_area_akodon.data", sep="")
node_area <- get_node_range_BioGeoBEARS(resDEC,phyllip.file = fn, akodon_newick,
                                        max.range.size = 4) # obtaining area for each node
```

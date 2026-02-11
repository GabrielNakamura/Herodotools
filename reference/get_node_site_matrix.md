# Extract a composition matrix of communities and nodes

Extract a composition matrix of communities and nodes

## Usage

``` r
get_node_site_matrix(simulation_path, save = FALSE, output_file = NULL)
```

## Arguments

- simulation_path:

  A string indicating the path from the root directory to the simulation
  results from gen3sis. Each folder must be separated by a forward slash
  "/"

- save:

  Logical. If FALSE (default), the list with the results won't be saved
  in your machine

- output_file:

  A string indicating the path where the output of this function will be
  saved. If NULL (default), and if save is TRUE the output will be saved
  in the directory root

## Value

A list with four components:

- node_site_matrix:

  A matrix containing the sites in the rows and the nodes of the tree in
  the columns

- tree:

  A phylo object with the species produced in the simulation process.
  The nodes were renamed from the original phylogeny

- nameIndexGuide:

  A data frame with node information

- coordinates:

  A data frame with coordinate information for the sites in
  node_site_matrix

## Examples

``` r
if (FALSE) { # \dontrun{
# Load node-site matrix from gen3sis simulation
result <- get_node_site_matrix(
  simulation_path = "path/to/simulation",
  save = TRUE,
  output_file = "output/node_site_data.rds"
)
} # }
```

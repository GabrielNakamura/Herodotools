# Generates lagrangePHYLIP file for using in BioGeoBEARS

It prepares and save the file to be used as geographical areas in
`BioGeoBEARS` scripts. It was based on
`save_tipranges_to_LagrangePHYLIP` function of `BioGEOBEARS` package.

## Usage

``` r
get_tipranges_to_BioGeoBEARS(
  tip_range,
  filename = "lagrange_area_data_file.data",
  areanames = NULL
)
```

## Arguments

- tip_range:

  data frame. Tip names in rows, geographical areas in columns. Please
  provide the tip names as `row.names`. See example.

- filename:

  filename to store the data. By default it is
  "lagrange_area_data_file.data".

- areanames:

  the names of the areas. By defaut the area names are generated as
  sequencial upper case letters.

## Examples
